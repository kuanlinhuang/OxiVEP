//! gnomAD gene constraint scores parser for building .oga files.
//!
//! Extracts pLI, LOEUF, mis_z, and syn_z per gene.

use crate::common::GeneRecord;
use anyhow::{Context, Result};
use std::io::BufRead;

/// Parse gnomAD constraint metrics TSV into GeneRecords.
///
/// Expected header: gene, transcript, obs_lof, exp_lof, oe_lof, oe_lof_upper, pLI, ...
pub fn parse_gnomad_gene_scores<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    let mut records = Vec::new();
    let mut col_indices: Option<GnomadGeneCols> = None;

    for line in reader.lines() {
        let line = line.context("Reading gnomAD gene scores")?;
        if line.starts_with("gene\t") || line.starts_with("#") {
            col_indices = Some(GnomadGeneCols::from_header(&line));
            continue;
        }
        if line.is_empty() {
            continue;
        }

        let cols = match &col_indices {
            Some(c) => c,
            None => continue,
        };

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= cols.max_idx() {
            continue;
        }

        let gene = fields[cols.gene].trim();
        if gene.is_empty() {
            continue;
        }

        let mut parts = Vec::new();

        if let Some(idx) = cols.pli {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"pLI\":{:.4}", v));
            }
        }
        if let Some(idx) = cols.loeuf {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"loeuf\":{:.4}", v));
            }
        }
        if let Some(idx) = cols.mis_z {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"misZ\":{:.2}", v));
            }
        }
        if let Some(idx) = cols.syn_z {
            if let Ok(v) = fields[idx].parse::<f64>() {
                parts.push(format!("\"synZ\":{:.2}", v));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(GeneRecord {
            gene_symbol: gene.to_string(),
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    Ok(records)
}

struct GnomadGeneCols {
    gene: usize,
    pli: Option<usize>,
    loeuf: Option<usize>,
    mis_z: Option<usize>,
    syn_z: Option<usize>,
}

impl GnomadGeneCols {
    fn from_header(header: &str) -> Self {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |name: &str| fields.iter().position(|f| f.to_lowercase() == name.to_lowercase());

        Self {
            gene: find("gene").unwrap_or(0),
            pli: find("pLI").or_else(|| find("pli")),
            loeuf: find("oe_lof_upper").or_else(|| find("loeuf")),
            mis_z: find("mis_z"),
            syn_z: find("syn_z"),
        }
    }

    fn max_idx(&self) -> usize {
        let mut m = self.gene;
        for opt in [self.pli, self.loeuf, self.mis_z, self.syn_z] {
            if let Some(i) = opt { m = m.max(i); }
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gnomad_gene_scores() {
        let data = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
TP53\tENST00000269305\t0\t25.1\t0.00\t0.05\t0.9999\t5.67\t-0.34
";
        let records = parse_gnomad_gene_scores(data.as_bytes()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(records[0].json.contains("\"pLI\":1.0000"));
        assert!(records[0].json.contains("\"loeuf\":0.0300"));
        assert!(records[0].json.contains("\"misZ\":3.45"));

        assert_eq!(records[1].gene_symbol, "TP53");
        assert!(records[1].json.contains("\"pLI\":0.9999"));
    }
}
