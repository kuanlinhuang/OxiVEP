use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use oxivep_core::Consequence;

/// Format a VCF CSQ INFO field value from a VariationFeature.
///
/// Fields match the standard VEP CSQ format:
/// Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
/// EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
/// Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS
pub fn format_csq(vf: &VariationFeature, fields: &[&str]) -> String {
    let mut result = String::with_capacity(1024);

    let mut first = true;
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            if !first {
                result.push(',');
            }
            first = false;
            let entry = format_csq_entry(vf, tv, aa, fields);
            result.push_str(&entry);
        }
    }

    result
}

fn format_csq_entry(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    aa: &AlleleAnnotation,
    fields: &[&str],
) -> String {
    let mut result = String::with_capacity(512);

    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            result.push('|');
        }
        let value = match *field {
            "Allele" => aa.allele.to_string(),
            "Consequence" => aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join("&"),
            "IMPACT" => format!("{:?}", aa.impact).to_uppercase(),
            "SYMBOL" => tv.gene_symbol.as_deref().unwrap_or_default().to_string(),
            "Gene" => tv.gene_id.to_string(),
            "Feature_type" => "Transcript".to_string(),
            "Feature" => tv.transcript_id.to_string(),
            "BIOTYPE" => tv.biotype.to_string(),
            "EXON" => aa
                .exon
                .map(|(n, t)| format!("{}/{}", n, t))
                .unwrap_or_default(),
            "INTRON" => aa
                .intron
                .map(|(n, t)| format!("{}/{}", n, t))
                .unwrap_or_default(),
            "HGVSc" => aa.hgvsc.clone().unwrap_or_default(),
            "HGVSp" => aa.hgvsp.clone().unwrap_or_default(),
            "cDNA_position" => format_position_range(aa.cdna_position),
            "CDS_position" => format_position_range(aa.cds_position),
            "Protein_position" => format_position_range(aa.protein_position),
            "Amino_acids" => aa
                .amino_acids
                .as_ref()
                .map(|(r, a)| {
                    if r == a {
                        // VEP shows just the amino acid for synonymous variants
                        r.clone()
                    } else {
                        format!("{}/{}", r, a)
                    }
                })
                .unwrap_or_default(),
            "Codons" => aa
                .codons
                .as_ref()
                .map(|(r, a)| format!("{}/{}", r, a))
                .unwrap_or_default(),
            "Existing_variation" => aa.existing_variation.join("&"),
            "REF_ALLELE" => vf.ref_allele.to_string(),
            "UPLOADED_ALLELE" => {
                // Show original VCF alleles if available, otherwise normalized
                if let Some(ref vcf) = vf.vcf_fields {
                    format!("{}/{}", vcf.ref_allele, vcf.alt)
                } else {
                    format!("{}/{}", vf.ref_allele, aa.allele)
                }
            }
            "DISTANCE" => aa
                .distance
                .map(|d| d.to_string())
                .unwrap_or_default(),
            "STRAND" => format!("{}", tv.strand.as_int()),
            "FLAGS" => {
                let flags = &tv.flags;
                if flags.is_empty() {
                    String::new()
                } else {
                    flags.join("&")
                }
            }
            "CANONICAL" => {
                if tv.canonical {
                    "YES".to_string()
                } else {
                    String::new()
                }
            }
            "SYMBOL_SOURCE" => tv.symbol_source.clone().unwrap_or_default(),
            "HGNC_ID" => tv.hgnc_id.clone().unwrap_or_default(),
            "MANE" => {
                // MANE field is the label: "MANE_Select" or "MANE_Plus_Clinical"
                if tv.mane_select.is_some() {
                    "MANE_Select".to_string()
                } else if tv.mane_plus_clinical.is_some() {
                    "MANE_Plus_Clinical".to_string()
                } else {
                    String::new()
                }
            }
            "MANE_SELECT" => tv.mane_select.clone().unwrap_or_default(),
            "MANE_PLUS_CLINICAL" => tv.mane_plus_clinical.clone().unwrap_or_default(),
            "TSL" => tv.tsl.map(|t| t.to_string()).unwrap_or_default(),
            "APPRIS" => tv.appris.clone().unwrap_or_default(),
            "CCDS" => tv.ccds.clone().unwrap_or_default(),
            "ENSP" => tv.protein_id.clone().unwrap_or_default(),
            "SIFT" => aa.sift.clone().unwrap_or_default(),
            "PolyPhen" => aa.polyphen.clone().unwrap_or_default(),
            "AF" => {
                // Use the first matched known variant's gnomAD or minor_allele_freq
                vf.existing_variants.iter().find_map(|kv| {
                    kv.frequencies.get("gnomAD")
                        .or_else(|| kv.frequencies.get("gnomADe"))
                        .or_else(|| kv.frequencies.get("minor_allele_freq"))
                        .map(|f| format!("{}", f))
                }).unwrap_or_default()
            }
            "CLIN_SIG" => {
                vf.existing_variants.iter()
                    .filter_map(|kv| kv.clinical_significance.as_ref())
                    .next()
                    .cloned()
                    .unwrap_or_default()
            }
            "SOMATIC" => {
                let any_somatic = vf.existing_variants.iter().any(|kv| kv.somatic);
                if any_somatic { "1".to_string() } else { String::new() }
            }
            "PHENO" => {
                let any_pheno = vf.existing_variants.iter().any(|kv| kv.phenotype_or_disease);
                if any_pheno { "1".to_string() } else { String::new() }
            }
            "PUBMED" => {
                let pubs: Vec<&str> = vf.existing_variants.iter()
                    .flat_map(|kv| kv.pubmed.iter().map(|s| s.as_str()))
                    .collect();
                pubs.join("&")
            }
            "SOURCE" => tv.source.clone().unwrap_or_default(),
            "HGVS_OFFSET" => aa
                .hgvs_offset
                .map(|o| o.to_string())
                .unwrap_or_default(),
            _ => String::new(),
        };
        escape_csq_into(&value, &mut result);
    }

    result
}

/// Escape special characters in CSQ field values, appending to an existing buffer.
fn escape_csq_into(value: &str, buf: &mut String) {
    for c in value.chars() {
        match c {
            ',' | '|' => buf.push('&'),
            ';' => buf.push_str("%3B"),
            '=' => buf.push_str("%3D"),
            _ => buf.push(c),
        }
    }
}

/// Escape special characters in CSQ field values.
#[cfg(test)]
fn escape_csq_value(value: &str) -> String {
    value
        .replace(',', "&")
        .replace(';', "%3B")
        .replace('=', "%3D")
        .replace('|', "&")
        .replace(' ', "_")
}

fn format_position_range(pos: Option<(u64, u64)>) -> String {
    match pos {
        Some((start, end)) if start == end => start.to_string(),
        Some((start, end)) => format!("{}-{}", start, end),
        None => String::new(),
    }
}

/// Default CSQ fields matching Ensembl VEP's extended output format.
///
/// Includes all standard VEP fields (CANONICAL, CCDS, ENSP, SOURCE,
/// HGVS_OFFSET) plus extended annotations (MANE, SIFT, PolyPhen, etc.).
pub const DEFAULT_CSQ_FIELDS: &[&str] = &[
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "REF_ALLELE",
    "UPLOADED_ALLELE",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SOURCE",
    "HGVS_OFFSET",
    "SIFT",
    "PolyPhen",
    "AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
];

/// Generate the VCF INFO header line for CSQ.
pub fn csq_header_line(fields: &[&str]) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from OxiVEP. Format: {}\">",
        fields.join("|")
    )
}

/// Format a VariationFeature as a tab-delimited VEP output line.
pub fn format_tab_line(vf: &VariationFeature) -> Vec<String> {
    let mut lines = Vec::new();

    let location = if vf.position.start == vf.position.end {
        format!("{}:{}", vf.position.chromosome, vf.position.start)
    } else {
        format!(
            "{}:{}-{}",
            vf.position.chromosome, vf.position.start, vf.position.end
        )
    };

    let uploaded_variation = vf
        .variation_name
        .clone()
        .unwrap_or_else(|| format!("{}_{}", location, vf.allele_string));

    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let consequence_str = aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join(",");

            let impact_str = format!("{:?}", aa.impact).to_uppercase();
            let distance_str = aa.distance.map(|d| d.to_string()).unwrap_or("-".to_string());
            let strand_str = format!("{}", tv.strand.as_int());
            let flags_str = if tv.canonical { "canonical" } else { "-" };

            let line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                uploaded_variation,
                location,
                aa.allele,
                tv.gene_id,
                tv.transcript_id,
                "Transcript",
                consequence_str,
                format_position_range(aa.cdna_position),
                format_position_range(aa.cds_position),
                format_position_range(aa.protein_position),
                aa.amino_acids
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or("-".to_string()),
                aa.codons
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or("-".to_string()),
                if aa.existing_variation.is_empty() { "-".to_string() } else { aa.existing_variation.join(",") },
                impact_str,
                distance_str,
                strand_str,
                flags_str,
            );
            lines.push(line);
        }
    }

    // If no transcript annotations, still output the variant with intergenic
    if vf.transcript_variations.is_empty() {
        for alt in &vf.alt_alleles {
            let line = format!(
                "{}\t{}\t{}\t-\t-\t-\t{}\t-\t-\t-\t-\t-\t-",
                uploaded_variation,
                location,
                alt,
                Consequence::IntergenicVariant.so_term(),
            );
            lines.push(line);
        }
    }

    lines
}

/// Format a VariationFeature as JSON.
pub fn format_json(vf: &VariationFeature) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    obj.insert("id".into(), json_str(&vf.variation_name));
    obj.insert(
        "seq_region_name".into(),
        serde_json::Value::String(vf.position.chromosome.clone()),
    );
    obj.insert("start".into(), serde_json::Value::Number(vf.position.start.into()));
    obj.insert("end".into(), serde_json::Value::Number(vf.position.end.into()));
    obj.insert(
        "allele_string".into(),
        serde_json::Value::String(vf.allele_string.clone()),
    );
    obj.insert("strand".into(), serde_json::Value::Number(vf.position.strand.as_int().into()));

    if let Some(ref msq) = vf.most_severe_consequence {
        obj.insert(
            "most_severe_consequence".into(),
            serde_json::Value::String(msq.so_term().to_string()),
        );
    }

    let transcript_consequences: Vec<serde_json::Value> = vf
        .transcript_variations
        .iter()
        .flat_map(|tv| {
            tv.allele_annotations.iter().map(move |aa| {
                let mut tc = serde_json::Map::new();
                tc.insert(
                    "gene_id".into(),
                    serde_json::Value::String(tv.gene_id.to_string()),
                );
                tc.insert(
                    "transcript_id".into(),
                    serde_json::Value::String(tv.transcript_id.to_string()),
                );
                tc.insert(
                    "biotype".into(),
                    serde_json::Value::String(tv.biotype.to_string()),
                );
                if let Some(ref sym) = tv.gene_symbol {
                    tc.insert(
                        "gene_symbol".into(),
                        serde_json::Value::String(sym.to_string()),
                    );
                }
                tc.insert(
                    "consequence_terms".into(),
                    serde_json::Value::Array(
                        aa.consequences
                            .iter()
                            .map(|c| serde_json::Value::String(c.so_term().to_string()))
                            .collect(),
                    ),
                );
                tc.insert(
                    "impact".into(),
                    serde_json::Value::String(format!("{:?}", aa.impact).to_uppercase()),
                );
                tc.insert(
                    "variant_allele".into(),
                    serde_json::Value::String(aa.allele.to_string()),
                );
                tc.insert(
                    "strand".into(),
                    serde_json::Value::Number(tv.strand.as_int().into()),
                );
                if tv.canonical {
                    tc.insert("canonical".into(), serde_json::Value::Number(1.into()));
                }
                if let Some((s, e)) = aa.cdna_position {
                    tc.insert("cdna_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("cdna_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some((s, e)) = aa.cds_position {
                    tc.insert("cds_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("cds_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some((s, e)) = aa.protein_position {
                    tc.insert("protein_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("protein_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some(ref aas) = aa.amino_acids {
                    tc.insert("amino_acids".into(),
                        serde_json::Value::String(format!("{}/{}", aas.0, aas.1)));
                }
                if let Some(ref cdns) = aa.codons {
                    tc.insert("codons".into(),
                        serde_json::Value::String(format!("{}/{}", cdns.0, cdns.1)));
                }
                if let Some((n, t)) = aa.exon {
                    tc.insert("exon".into(),
                        serde_json::Value::String(format!("{}/{}", n, t)));
                }
                if let Some((n, t)) = aa.intron {
                    tc.insert("intron".into(),
                        serde_json::Value::String(format!("{}/{}", n, t)));
                }
                if let Some(ref h) = aa.hgvsc {
                    tc.insert("hgvsc".into(), serde_json::Value::String(h.clone()));
                }
                if let Some(ref h) = aa.hgvsp {
                    tc.insert("hgvsp".into(), serde_json::Value::String(h.clone()));
                }
                if let Some(d) = aa.distance {
                    tc.insert("distance".into(), serde_json::Value::Number(d.into()));
                }
                serde_json::Value::Object(tc)
            })
        })
        .collect();

    obj.insert(
        "transcript_consequences".into(),
        serde_json::Value::Array(transcript_consequences),
    );

    serde_json::Value::Object(obj)
}

fn json_str(opt: &Option<String>) -> serde_json::Value {
    match opt {
        Some(s) => serde_json::Value::String(s.clone()),
        None => serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_csq_value() {
        assert_eq!(escape_csq_value("hello,world"), "hello&world");
        assert_eq!(escape_csq_value("a;b"), "a%3Bb");
        assert_eq!(escape_csq_value("a|b"), "a&b");
        assert_eq!(escape_csq_value("a b"), "a_b");
        assert_eq!(escape_csq_value("p.Leu153="), "p.Leu153%3D");
    }

    #[test]
    fn test_csq_header() {
        let header = csq_header_line(&["Allele", "Consequence"]);
        assert!(header.contains("Format: Allele|Consequence"));
    }

    #[test]
    fn test_format_position_range() {
        assert_eq!(format_position_range(Some((100, 100))), "100");
        assert_eq!(format_position_range(Some((100, 200))), "100-200");
        assert_eq!(format_position_range(None), "");
    }
}
