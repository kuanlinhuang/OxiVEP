use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_core::Consequence;
use std::fmt::Write as FmtWrite;

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
            format_csq_entry_into(vf, tv, aa, fields, &mut result);
        }
    }

    result
}

/// Write a single CSQ entry directly into the output buffer, avoiding intermediate String allocations.
fn format_csq_entry_into(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    aa: &AlleleAnnotation,
    fields: &[&str],
    buf: &mut String,
) {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            buf.push('|');
        }
        // Write each field value directly into buf via escape_csq_str, avoiding
        // temporary String allocations for most fields.
        match *field {
            "Allele" => escape_csq_str(&aa.allele.to_string(), buf),
            "Consequence" => {
                for (j, c) in aa.consequences.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    buf.push_str(c.so_term());
                }
            }
            "IMPACT" => buf.push_str(aa.impact.as_str()),
            "SYMBOL" => escape_csq_str(tv.gene_symbol.as_deref().unwrap_or_default(), buf),
            "Gene" => escape_csq_str(&tv.gene_id, buf),
            "Feature_type" => buf.push_str("Transcript"),
            "Feature" => escape_csq_str(&tv.transcript_id, buf),
            "BIOTYPE" => escape_csq_str(&tv.biotype, buf),
            "EXON" => {
                if let Some((n, t)) = aa.exon {
                    let _ = write!(buf, "{}/{}", n, t);
                }
            }
            "INTRON" => {
                if let Some((n, t)) = aa.intron {
                    let _ = write!(buf, "{}/{}", n, t);
                }
            }
            "HGVSc" => escape_csq_str(aa.hgvsc.as_deref().unwrap_or_default(), buf),
            "HGVSp" => escape_csq_str(aa.hgvsp.as_deref().unwrap_or_default(), buf),
            "cDNA_position" => write_position_range(aa.cdna_position, buf),
            "CDS_position" => write_position_range(aa.cds_position, buf),
            "Protein_position" => write_position_range(aa.protein_position, buf),
            "Amino_acids" => {
                if let Some((ref r, ref a)) = aa.amino_acids {
                    escape_csq_str(r, buf);
                    if r != a {
                        buf.push('/');
                        escape_csq_str(a, buf);
                    }
                }
            }
            "Codons" => {
                if let Some((ref r, ref a)) = aa.codons {
                    escape_csq_str(r, buf);
                    buf.push('/');
                    escape_csq_str(a, buf);
                }
            }
            "Existing_variation" => {
                for (j, ev) in aa.existing_variation.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    escape_csq_str(ev, buf);
                }
            }
            "REF_ALLELE" => escape_csq_str(&vf.ref_allele.to_string(), buf),
            "UPLOADED_ALLELE" => {
                if let Some(ref vcf) = vf.vcf_fields {
                    escape_csq_str(&vcf.ref_allele, buf);
                    buf.push('/');
                    escape_csq_str(&vcf.alt, buf);
                } else {
                    escape_csq_str(&vf.ref_allele.to_string(), buf);
                    buf.push('/');
                    escape_csq_str(&aa.allele.to_string(), buf);
                }
            }
            "DISTANCE" => {
                if let Some(d) = aa.distance {
                    let _ = write!(buf, "{}", d);
                }
            }
            "STRAND" => { let _ = write!(buf, "{}", tv.strand.as_int()); }
            "FLAGS" => {
                for (j, f) in tv.flags.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    buf.push_str(f);
                }
            }
            "CANONICAL" => { if tv.canonical { buf.push_str("YES"); } }
            "SYMBOL_SOURCE" => escape_csq_str(tv.symbol_source.as_deref().unwrap_or_default(), buf),
            "HGNC_ID" => escape_csq_str(tv.hgnc_id.as_deref().unwrap_or_default(), buf),
            "MANE" => {
                if tv.mane_select.is_some() {
                    buf.push_str("MANE_Select");
                } else if tv.mane_plus_clinical.is_some() {
                    buf.push_str("MANE_Plus_Clinical");
                }
            }
            "MANE_SELECT" => escape_csq_str(tv.mane_select.as_deref().unwrap_or_default(), buf),
            "MANE_PLUS_CLINICAL" => escape_csq_str(tv.mane_plus_clinical.as_deref().unwrap_or_default(), buf),
            "TSL" => { if let Some(t) = tv.tsl { let _ = write!(buf, "{}", t); } }
            "APPRIS" => escape_csq_str(tv.appris.as_deref().unwrap_or_default(), buf),
            "CCDS" => escape_csq_str(tv.ccds.as_deref().unwrap_or_default(), buf),
            "ENSP" => escape_csq_str(tv.protein_id.as_deref().unwrap_or_default(), buf),
            "SIFT" => escape_csq_str(aa.sift.as_deref().unwrap_or_default(), buf),
            "PolyPhen" => escape_csq_str(aa.polyphen.as_deref().unwrap_or_default(), buf),
            "AF" => {
                if let Some(f) = vf.existing_variants.iter().find_map(|kv| {
                    kv.frequencies.get("gnomAD")
                        .or_else(|| kv.frequencies.get("gnomADe"))
                        .or_else(|| kv.frequencies.get("minor_allele_freq"))
                }) {
                    let _ = write!(buf, "{}", f);
                }
            }
            "CLIN_SIG" => {
                if let Some(cs) = vf.existing_variants.iter()
                    .find_map(|kv| kv.clinical_significance.as_deref())
                {
                    escape_csq_str(cs, buf);
                }
            }
            "SOMATIC" => { if vf.existing_variants.iter().any(|kv| kv.somatic) { buf.push('1'); } }
            "PHENO" => { if vf.existing_variants.iter().any(|kv| kv.phenotype_or_disease) { buf.push('1'); } }
            "PUBMED" => {
                let mut first_pub = true;
                for kv in &vf.existing_variants {
                    for p in &kv.pubmed {
                        if !first_pub { buf.push('&'); }
                        first_pub = false;
                        buf.push_str(p);
                    }
                }
            }
            "SOURCE" => escape_csq_str(tv.source.as_deref().unwrap_or_default(), buf),
            "HGVS_OFFSET" => {
                if let Some(o) = aa.hgvs_offset {
                    let _ = write!(buf, "{}", o);
                }
            }
            _ => {}
        }
    }
}

/// Escape special characters in CSQ field values, appending to an existing buffer.
fn escape_csq_str(value: &str, buf: &mut String) {
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

fn write_position_range(pos: Option<(u64, u64)>, buf: &mut String) {
    match pos {
        Some((start, end)) if start == end => { let _ = write!(buf, "{}", start); }
        Some((start, end)) => { let _ = write!(buf, "{}-{}", start, end); }
        None => {}
    }
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
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from fastVEP. Format: {}\">",
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

            let impact_str = aa.impact.as_str().to_string();
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
                    serde_json::Value::String(aa.impact.as_str().to_string()),
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
                // Per-allele supplementary annotations (ClinVar, gnomAD, etc.)
                for (key, json_str) in &aa.supplementary {
                    if let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) {
                        tc.insert(key.clone(), val);
                    }
                }
                serde_json::Value::Object(tc)
            })
        })
        .collect();

    obj.insert(
        "transcript_consequences".into(),
        serde_json::Value::Array(transcript_consequences),
    );

    // Variant type (for SVs)
    if vf.variant_type != fastvep_core::VariantType::Unknown {
        obj.insert(
            "variant_type".into(),
            serde_json::Value::String(format!("{:?}", vf.variant_type)),
        );
    }
    if let Some(sv_end) = vf.sv_end {
        obj.insert("sv_end".into(), serde_json::Value::Number(sv_end.into()));
    }
    if let Some(sv_len) = vf.sv_len {
        obj.insert("sv_len".into(), serde_json::Value::Number(sv_len.into()));
    }

    // Supplementary annotations (pre-serialized JSON from SA providers)
    for sa in &vf.supplementary_annotations {
        if let Ok(val) = serde_json::from_str::<serde_json::Value>(&sa.json_string) {
            obj.insert(sa.json_key.clone(), val);
        }
    }

    // Gene-level annotations
    if !vf.gene_annotations.is_empty() {
        let mut genes_map = serde_json::Map::new();
        for ga in &vf.gene_annotations {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&ga.json_string) {
                genes_map
                    .entry(ga.gene_symbol.clone())
                    .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()))
                    .as_object_mut()
                    .map(|obj| obj.insert(ga.json_key.clone(), val));
            }
        }
        if !genes_map.is_empty() {
            obj.insert("genes".into(), serde_json::Value::Object(genes_map));
        }
    }

    serde_json::Value::Object(obj)
}

/// Format a variant as Nirvana-style structured JSON.
///
/// This is a richer format with nested sections:
/// - `position`: chromosome, position coordinates
/// - `variants`: array of per-allele annotations with supplementary data
/// - `genes`: gene-level annotations keyed by symbol
pub fn format_nirvana_json(vf: &VariationFeature) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    // Position section
    obj.insert("chromosome".into(), serde_json::Value::String(vf.position.chromosome.clone()));
    obj.insert("position".into(), serde_json::Value::Number(vf.position.start.into()));
    obj.insert("end".into(), serde_json::Value::Number(vf.position.end.into()));

    let ref_str = vf.ref_allele.to_string();
    obj.insert("refAllele".into(), serde_json::Value::String(ref_str));

    let alt_strs: Vec<serde_json::Value> = vf
        .alt_alleles
        .iter()
        .map(|a| serde_json::Value::String(a.to_string()))
        .collect();
    obj.insert("altAlleles".into(), serde_json::Value::Array(alt_strs));

    if vf.variant_type != fastvep_core::VariantType::Unknown {
        obj.insert("variantType".into(), serde_json::Value::String(format!("{:?}", vf.variant_type)));
    }

    // Variants section: one per alt allele with all transcript consequences
    let mut variants = Vec::new();
    for alt in &vf.alt_alleles {
        let alt_str = alt.to_string();
        let mut var_obj = serde_json::Map::new();
        var_obj.insert("altAllele".into(), serde_json::Value::String(alt_str.clone()));

        // Collect transcript consequences for this allele
        let mut transcripts = Vec::new();
        for tv in &vf.transcript_variations {
            for aa in &tv.allele_annotations {
                if aa.allele.to_string() == alt_str {
                    let mut tc = serde_json::Map::new();
                    tc.insert("transcriptId".into(), serde_json::Value::String(tv.transcript_id.to_string()));
                    tc.insert("geneId".into(), serde_json::Value::String(tv.gene_id.to_string()));
                    if let Some(ref sym) = tv.gene_symbol {
                        tc.insert("geneSymbol".into(), serde_json::Value::String(sym.to_string()));
                    }
                    tc.insert("biotype".into(), serde_json::Value::String(tv.biotype.to_string()));
                    tc.insert("consequences".into(), serde_json::Value::Array(
                        aa.consequences.iter().map(|c| serde_json::Value::String(c.so_term().to_string())).collect()
                    ));
                    tc.insert("impact".into(), serde_json::Value::String(aa.impact.as_str().to_string()));
                    if tv.canonical { tc.insert("isCanonical".into(), serde_json::Value::Bool(true)); }
                    if let Some(ref h) = aa.hgvsc { tc.insert("hgvsc".into(), serde_json::Value::String(h.clone())); }
                    if let Some(ref h) = aa.hgvsp { tc.insert("hgvsp".into(), serde_json::Value::String(h.clone())); }

                    // Per-allele supplementary annotations
                    for (key, json_str) in &aa.supplementary {
                        if let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) {
                            tc.insert(key.clone(), val);
                        }
                    }

                    transcripts.push(serde_json::Value::Object(tc));
                }
            }
        }
        var_obj.insert("transcripts".into(), serde_json::Value::Array(transcripts));
        variants.push(serde_json::Value::Object(var_obj));
    }
    obj.insert("variants".into(), serde_json::Value::Array(variants));

    // Gene annotations
    if !vf.gene_annotations.is_empty() {
        let mut genes_map = serde_json::Map::new();
        for ga in &vf.gene_annotations {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&ga.json_string) {
                genes_map
                    .entry(ga.gene_symbol.clone())
                    .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()))
                    .as_object_mut()
                    .map(|obj| obj.insert(ga.json_key.clone(), val));
            }
        }
        if !genes_map.is_empty() {
            obj.insert("genes".into(), serde_json::Value::Object(genes_map));
        }
    }

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
