//! dbSNP VCF parser for building .osa annotation files.
//!
//! Parses dbSNP's VCF release to extract RS IDs and global MAF.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a dbSNP VCF and produce sorted `AnnotationRecord`s.
pub fn parse_dbsnp_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading dbSNP VCF line")?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let id = fields[2];
        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];

        // Extract RS ID
        let rs_id = if id.starts_with("rs") {
            id.to_string()
        } else {
            // Try INFO field
            let info_map = parse_info(info);
            match info_map.get("RS") {
                Some(rs) => format!("rs{}", rs),
                None => continue, // Skip if no RS ID
            }
        };

        // Parse global MAF if available
        let info_map = parse_info(info);
        let freq = info_map
            .get("CAF")
            .and_then(|caf| {
                // CAF format: ref_freq,alt1_freq,alt2_freq,...
                let parts: Vec<&str> = caf.split(',').collect();
                parts.get(1).and_then(|s| s.parse::<f64>().ok())
            });

        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }

            let mut parts = vec![format!("\"id\":\"{}\"", rs_id)];
            if let Some(f) = freq {
                parts.push(format!("\"globalMaf\":{:.6e}", f));
            }
            let json = format!("{{{}}}", parts.join(","));

            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.to_string(),
                json,
            });
        }
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        if let Some((key, value)) = pair.split_once('=') {
            map.insert(key.to_string(), value.to_string());
        }
    }
    map
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_dbsnp_vcf() {
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;CAF=0.9998,0.0002
1\t10039\trs978760828\tA\tC\t.\t.\tRS=978760828
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 10019);
        assert!(records[0].json.contains("rs775809821"));
        assert!(records[0].json.contains("globalMaf"));

        assert_eq!(records[1].position, 10039);
        assert!(records[1].json.contains("rs978760828"));
        assert!(!records[1].json.contains("globalMaf")); // No CAF
    }
}
