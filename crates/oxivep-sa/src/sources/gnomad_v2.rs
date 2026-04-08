//! gnomAD .osa2 encoder — direct VCF-to-parallel-u32-arrays encoding.
//!
//! Bypasses intermediate JSON for maximum compression and query speed.

use crate::fields::{Field, FieldType};
use crate::writer_v2::{Osa2Metadata, Osa2Record};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Standard gnomAD v2 field configuration.
pub fn gnomad_fields() -> Vec<Field> {
    let pop_fields: Vec<Field> = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
        .iter()
        .map(|pop| Field {
            field: format!("AF_{}", pop),
            alias: format!("{}Af", pop),
            ftype: FieldType::Float,
            multiplier: 2_000_000,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: format!("{} allele frequency", pop.to_uppercase()),
        })
        .collect();

    let mut fields = vec![
        Field {
            field: "AF".into(), alias: "allAf".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Global allele frequency".into(),
        },
        Field {
            field: "AN".into(), alias: "allAn".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Total allele number".into(),
        },
        Field {
            field: "AC".into(), alias: "allAc".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Allele count".into(),
        },
        Field {
            field: "nhomalt".into(), alias: "allHc".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Homozygote count".into(),
        },
    ];
    fields.extend(pop_fields);
    fields
}

/// Standard gnomAD .osa2 metadata.
pub fn gnomad_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("gnomAD population frequencies for {}", assembly),
    }
}

/// Parse gnomAD VCF directly into Osa2Records with parallel u32 values.
pub fn parse_gnomad_to_osa2<R: BufRead>(
    reader: R,
    fields: &[Field],
) -> Result<Vec<Osa2Record>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading gnomAD VCF")?;
        if line.starts_with('#') { continue; }

        let parts: Vec<&str> = line.splitn(9, '\t').collect();
        if parts.len() < 8 { continue; }

        let chrom = normalize_chrom(parts[0]);
        let pos: u32 = match parts[1].parse() { Ok(p) => p, Err(_) => continue };
        let ref_allele = parts[3].as_bytes().to_vec();
        let alt_field = parts[4];
        let info = parts[7];
        let info_map = parse_info(info);

        let alts: Vec<&str> = alt_field.split(',').collect();

        for (ai, alt) in alts.iter().enumerate() {
            if *alt == "." || *alt == "*" { continue; }

            let mut values = Vec::with_capacity(fields.len());
            for field in fields {
                let raw = get_allele_value(&info_map, &field.field, ai);
                let encoded = match field.ftype {
                    FieldType::Float => {
                        match raw.and_then(|v| v.parse::<f64>().ok()) {
                            Some(f) => field.encode_float(f),
                            None => field.missing_value,
                        }
                    }
                    FieldType::Integer => {
                        match raw.and_then(|v| v.parse::<i64>().ok()) {
                            Some(i) => field.encode_int(i),
                            None => field.missing_value,
                        }
                    }
                    _ => field.missing_value,
                };
                values.push(encoded);
            }

            records.push(Osa2Record {
                chrom: chrom.clone(),
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.as_bytes().to_vec(),
                values,
                json_blob: None,
            });
        }
    }

    // Sort by (chrom, position)
    records.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn get_allele_value<'a>(info: &'a HashMap<String, String>, key: &str, allele_idx: usize) -> Option<&'a str> {
    info.get(key).and_then(|v| {
        let parts: Vec<&str> = v.split(',').collect();
        // For per-allele fields, use the allele index; for single-value, use first
        if parts.len() > allele_idx {
            Some(parts[allele_idx])
        } else {
            parts.first().copied()
        }
    })
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for p in info.split(';') {
        if let Some((k, v)) = p.split_once('=') { m.insert(k.into(), v.into()); }
    }
    m
}

fn normalize_chrom(c: &str) -> String {
    if c.starts_with("chr") { c.to_string() } else { format!("chr{}", c) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gnomad_v2_parse() {
        let vcf = "#h\nchr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005\n";
        let fields = gnomad_fields();
        let records = parse_gnomad_to_osa2(vcf.as_bytes(), &fields).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].values.len(), fields.len());

        // Check AF encoding: 0.001 * 2_000_000 = 2000
        assert_eq!(records[0].values[0], 2000);
        // Check AN: 150000
        assert_eq!(records[0].values[1], 150000);
        // Check AC: 150
        assert_eq!(records[0].values[2], 150);
    }
}
