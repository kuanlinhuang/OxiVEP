//! gnomAD VCF parser for building .osa annotation files.
//!
//! Parses gnomAD's sites-only VCF to extract allele frequencies
//! per population, allele counts, and homozygote counts.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Population keys to extract from gnomAD VCF INFO field.
const POPULATIONS: &[&str] = &[
    "afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas",
];

/// Parse a gnomAD sites-only VCF and produce sorted `AnnotationRecord`s.
pub fn parse_gnomad_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading gnomAD VCF line")?;
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

        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];

        let info_map = parse_info(info);

        // Handle multi-allelic: split allele-specific fields by comma
        let alts: Vec<&str> = alt_field.split(',').collect();
        let all_afs = split_info_values(info_map.get("AF").map(|s| s.as_str()));
        let all_ans = split_info_values(info_map.get("AN").map(|s| s.as_str()));
        let all_acs = split_info_values(info_map.get("AC").map(|s| s.as_str()));
        let all_nhomalt = split_info_values(info_map.get("nhomalt").map(|s| s.as_str()));

        for (i, alt) in alts.iter().enumerate() {
            if *alt == "." || *alt == "*" {
                continue;
            }

            let json = build_gnomad_json(
                all_afs.get(i).map(|s| s.as_str()),
                all_ans.first().map(|s| s.as_str()), // AN is typically single-valued
                all_acs.get(i).map(|s| s.as_str()),
                all_nhomalt.get(i).map(|s| s.as_str()),
                &info_map,
                i,
            );

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

fn build_gnomad_json(
    af: Option<&str>,
    an: Option<&str>,
    ac: Option<&str>,
    nhomalt: Option<&str>,
    info_map: &HashMap<String, String>,
    allele_idx: usize,
) -> String {
    let mut parts = Vec::new();

    if let Some(af_str) = af {
        if let Ok(f) = af_str.parse::<f64>() {
            parts.push(format!("\"allAf\":{:.6e}", f));
        }
    }

    if let Some(an_str) = an {
        parts.push(format!("\"allAn\":{}", an_str));
    }

    if let Some(ac_str) = ac {
        parts.push(format!("\"allAc\":{}", ac_str));
    }

    if let Some(nh) = nhomalt {
        parts.push(format!("\"allHc\":{}", nh));
    }

    // Per-population AFs
    for pop in POPULATIONS {
        let key = format!("AF_{}", pop);
        if let Some(val) = info_map.get(&key) {
            let vals = split_info_values(Some(val.as_str()));
            if let Some(af_str) = vals.get(allele_idx) {
                if let Ok(f) = af_str.parse::<f64>() {
                    parts.push(format!("\"{}Af\":{:.6e}", pop, f));
                }
            }
        }
    }

    format!("{{{}}}", parts.join(","))
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

fn split_info_values(value: Option<&str>) -> Vec<String> {
    match value {
        Some(v) => v.split(',').map(|s| s.to_string()).collect(),
        None => Vec::new(),
    }
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
    fn test_parse_gnomad_vcf() {
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF=0.01,0.005;AN=140000;AC=1400,700;nhomalt=10,3;AF_eas=0.02,0.01
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3); // 1 SNV + 2 from multi-allelic

        // First record
        assert_eq!(records[0].position, 10001);
        assert!(records[0].json.contains("\"allAf\":"));
        assert!(records[0].json.contains("\"afrAf\":"));
        assert!(records[0].json.contains("\"nfeAf\":"));

        // Multi-allelic: second alt
        assert_eq!(records[2].position, 20000);
        assert_eq!(records[2].alt_allele, "A");
    }
}
