//! ClinVar VCF parser for building .osa annotation files.
//!
//! Parses ClinVar's VCF release to extract clinical significance,
//! review status, disease names, and accession numbers.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a ClinVar VCF file and produce sorted `AnnotationRecord`s.
///
/// The VCF must be from NCBI's ClinVar release (clinvar.vcf.gz).
/// Records are sorted by (chrom_idx, position) for `SaWriter`.
pub fn parse_clinvar_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading ClinVar VCF line")?;
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
            None => continue, // Skip unrecognized chromosomes
        };

        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];

        // Parse INFO field
        let info = fields[7];
        let info_map = parse_info(info);

        let clnsig = info_map.get("CLNSIG").cloned().unwrap_or_default();
        let clnrevstat = info_map.get("CLNREVSTAT").cloned().unwrap_or_default();
        let clndn = info_map.get("CLNDN").cloned().unwrap_or_default();
        let clnacc = info_map.get("CLNVC").cloned(); // variant class
        let clnid = info_map.get("CLNVCSO").cloned(); // SO accession

        // Handle multi-allelic: each ALT gets its own record
        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }

            let json = build_clinvar_json(&clnsig, &clnrevstat, &clndn, clnacc.as_deref(), clnid.as_deref());

            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.to_string(),
                json,
            });
        }
    }

    // Sort by (chrom_idx, position) as required by SaWriter
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));

    Ok(records)
}

fn build_clinvar_json(
    clnsig: &str,
    clnrevstat: &str,
    clndn: &str,
    clnvc: Option<&str>,
    clnvcso: Option<&str>,
) -> String {
    let mut parts = Vec::new();

    if !clnsig.is_empty() {
        let sigs: Vec<String> = clnsig
            .split('|')
            .map(|s| format!("\"{}\"", escape_json(s)))
            .collect();
        parts.push(format!("\"significance\":[{}]", sigs.join(",")));
    }

    if !clnrevstat.is_empty() {
        parts.push(format!(
            "\"reviewStatus\":\"{}\"",
            escape_json(clnrevstat)
        ));
    }

    if !clndn.is_empty() && clndn != "not_provided" {
        let diseases: Vec<String> = clndn
            .split('|')
            .filter(|s| *s != "not_provided")
            .map(|s| format!("\"{}\"", escape_json(s)))
            .collect();
        if !diseases.is_empty() {
            parts.push(format!("\"phenotypes\":[{}]", diseases.join(",")));
        }
    }

    if let Some(vc) = clnvc {
        parts.push(format!("\"variantClass\":\"{}\"", escape_json(vc)));
    }

    if let Some(vcso) = clnvcso {
        parts.push(format!("\"soAccession\":\"{}\"", escape_json(vcso)));
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

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_clinvar_vcf() {
        let vcf = "\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String>
##INFO=<ID=CLNREVSTAT,Number=.,Type=String>
##INFO=<ID=CLNDN,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t12345\trs123\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=Breast_cancer
1\t67890\trs456\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=not_provided
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_clinvar_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 12345);
        assert!(records[0].json.contains("Pathogenic"));
        assert!(records[0].json.contains("Breast_cancer"));

        assert_eq!(records[1].position, 67890);
        assert!(records[1].json.contains("Benign"));
        // "not_provided" should be filtered out
        assert!(!records[1].json.contains("phenotypes"));
    }
}
