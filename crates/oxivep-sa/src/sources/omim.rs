//! OMIM genemap2.txt parser for building .oga gene annotation files.
//!
//! Extracts gene-phenotype mappings, MIM numbers, and inheritance patterns.

use crate::common::GeneRecord;
use anyhow::{Context, Result};
use std::io::BufRead;

/// Parse OMIM genemap2.txt into GeneRecords.
///
/// Format: tab-separated with columns including:
/// - Column 5: Gene Symbols (may have comma-separated aliases)
/// - Column 8: MIM Number
/// - Column 12: Phenotypes (semicolon-separated, each with MIM and inheritance)
pub fn parse_omim_genemap<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading OMIM line")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 13 {
            continue;
        }

        let gene_symbols = fields[5]; // Approved Gene Symbol
        let mim_number = fields[5 + 3]; // MIM Number (column 8, 0-indexed)

        // Use the primary (first) gene symbol
        let gene_symbol = gene_symbols.split(',').next().unwrap_or("").trim();
        if gene_symbol.is_empty() {
            continue;
        }

        let phenotypes_raw = if fields.len() > 12 { fields[12] } else { "" };

        let mut parts = Vec::new();

        if !mim_number.is_empty() && mim_number != "." {
            parts.push(format!("\"mimNumber\":{}", mim_number));
        }

        if !phenotypes_raw.is_empty() && phenotypes_raw != "." {
            let phenotypes: Vec<String> = phenotypes_raw
                .split(';')
                .filter(|p| !p.trim().is_empty())
                .map(|p| {
                    let p = p.trim();
                    format!("\"{}\"", escape_json(p))
                })
                .collect();
            if !phenotypes.is_empty() {
                parts.push(format!("\"phenotypes\":[{}]", phenotypes.join(",")));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(GeneRecord {
            gene_symbol: gene_symbol.to_string(),
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    Ok(records)
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_omim() {
        // Simplified genemap2 format
        let data = "# Generated\n\
                     # Copyright OMIM\n\
                     1\tp36.33\t1:10001-20000\tGene1\t\tBRCA1\tprotein\t\t113705\t\t\t\tBreast cancer, 114480 (3), Autosomal dominant; Ovarian cancer, 167000 (3)\n";

        let records = parse_omim_genemap(data.as_bytes()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(records[0].json.contains("113705"));
        assert!(records[0].json.contains("Breast cancer"));
    }
}
