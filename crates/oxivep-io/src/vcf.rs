use anyhow::{Context, Result};
use oxivep_core::{Allele, GenomicPosition, Strand, VariantType};
use std::io::{BufRead, BufReader, Read};

use crate::variant::{VariationFeature, VcfFields};

/// Parse a VCF file and yield VariationFeatures.
pub struct VcfParser<R: Read> {
    reader: BufReader<R>,
    header_lines: Vec<String>,
    line_buf: String,
}

impl<R: Read> VcfParser<R> {
    pub fn new(reader: R) -> Result<Self> {
        let mut buf_reader = BufReader::new(reader);
        let mut header_lines = Vec::new();
        let mut line_buf = String::new();

        // Read header lines
        loop {
            line_buf.clear();
            let bytes = buf_reader.read_line(&mut line_buf)?;
            if bytes == 0 {
                break;
            }
            let trimmed = line_buf.trim_end();
            if trimmed.starts_with('#') {
                header_lines.push(trimmed.to_string());
            } else {
                // This is the first data line; keep it in line_buf
                break;
            }
        }

        Ok(Self {
            reader: buf_reader,
            header_lines,
            line_buf,
        })
    }

    /// Get all VCF header lines (including #CHROM line).
    pub fn header_lines(&self) -> &[String] {
        &self.header_lines
    }

    /// Parse the next variant(s) from the VCF.
    /// Returns None when EOF is reached.
    /// A multi-allelic line may return multiple VariationFeatures if
    /// minimal mode is used in the future.
    pub fn next_variant(&mut self) -> Result<Option<VariationFeature>> {
        if self.line_buf.is_empty() {
            self.line_buf.clear();
            let bytes = self.reader.read_line(&mut self.line_buf)?;
            if bytes == 0 {
                return Ok(None);
            }
        }

        // Trim in-place to avoid allocating a new String
        let trimmed_len = self.line_buf.trim_end().len();
        self.line_buf.truncate(trimmed_len);

        if self.line_buf.is_empty() || self.line_buf.starts_with('#') {
            self.line_buf.clear();
            return self.next_variant();
        }

        let result = parse_vcf_line(&self.line_buf).map(Some);
        self.line_buf.clear(); // consumed
        result
    }

    /// Read all variants into a Vec.
    pub fn read_all(&mut self) -> Result<Vec<VariationFeature>> {
        let mut variants = Vec::new();
        while let Some(vf) = self.next_variant()? {
            variants.push(vf);
        }
        Ok(variants)
    }
}

/// Parse a single VCF data line into a VariationFeature.
pub fn parse_vcf_line(line: &str) -> Result<VariationFeature> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        anyhow::bail!("VCF line has fewer than 8 fields: {}", line);
    }

    let chrom = fields[0];
    let pos: u64 = fields[1]
        .parse()
        .with_context(|| format!("Invalid POS: {}", fields[1]))?;
    let id = fields[2];
    let ref_str = fields[3];
    let alt_str = fields[4];
    let qual = fields[5];
    let filter = fields[6];
    let info = fields[7];

    let rest: Vec<String> = fields[8..].iter().map(|s| s.to_string()).collect();

    // Parse alt alleles (split on comma)
    let raw_alts: Vec<&str> = alt_str.split(',').collect();

    // Determine start/end and normalize alleles
    let mut start = pos;
    let end;
    let mut ref_allele_str = ref_str.to_string();
    let mut alt_allele_strs: Vec<String> = raw_alts.iter().map(|s| s.to_string()).collect();

    // Check if any alt makes this an indel
    let is_indel = alt_allele_strs.iter().any(|alt| {
        alt.starts_with('D')
            || alt.starts_with('I')
            || alt.len() != ref_allele_str.len()
    });

    let is_non_variant = alt_str == "." || alt_str == "<NON_REF>" || alt_str == "<*>";

    // Check for symbolic/structural variant alleles
    let has_symbolic = raw_alts
        .iter()
        .any(|a| a.starts_with('<') && a.ends_with('>') && *a != "<NON_REF>" && *a != "<*>");

    if !is_non_variant && !has_symbolic && is_indel {
        if alt_allele_strs.len() > 1 {
            // Multi-allelic indel: strip shared first base only if ALL non-star alleles share it
            let non_star: Vec<&str> = std::iter::once(ref_allele_str.as_str())
                .chain(alt_allele_strs.iter().filter(|a| !a.contains('*')).map(|s| s.as_str()))
                .collect();

            let all_share_first = non_star.len() > 1
                && non_star
                    .iter()
                    .all(|s| !s.is_empty() && s.as_bytes()[0] == non_star[0].as_bytes()[0]);

            if all_share_first {
                ref_allele_str = if ref_allele_str.len() > 1 {
                    ref_allele_str[1..].to_string()
                } else {
                    "-".to_string()
                };
                start += 1;

                alt_allele_strs = alt_allele_strs
                    .iter()
                    .map(|alt| {
                        if alt.contains('*') {
                            alt.clone()
                        } else if alt.len() > 1 {
                            alt[1..].to_string()
                        } else {
                            "-".to_string()
                        }
                    })
                    .collect();
            }
        } else {
            // Single alt indel: strip shared first base
            let alt = &alt_allele_strs[0];
            if !ref_allele_str.is_empty()
                && !alt.is_empty()
                && ref_allele_str.as_bytes()[0] == alt.as_bytes()[0]
            {
                ref_allele_str = if ref_allele_str.len() > 1 {
                    ref_allele_str[1..].to_string()
                } else {
                    "-".to_string()
                };
                alt_allele_strs[0] = if alt.len() > 1 {
                    alt[1..].to_string()
                } else {
                    "-".to_string()
                };
                start += 1;
            }
        }
    }

    // Calculate end position
    if ref_allele_str == "-" {
        // Insertion: end = start - 1 (zero-length interval in Ensembl convention)
        end = start - 1;
    } else {
        end = start + ref_allele_str.len() as u64 - 1;
    }

    // Build allele string: "REF/ALT1/ALT2"
    let allele_string = if is_non_variant {
        ref_allele_str.clone()
    } else {
        format!(
            "{}/{}",
            ref_allele_str,
            alt_allele_strs.join("/")
        )
    };

    // Convert to Allele enums
    let ref_allele = Allele::from_str(&ref_allele_str);
    let alt_alleles: Vec<Allele> = alt_allele_strs.iter().map(|s| Allele::from_str(s)).collect();

    let variation_name = if id == "." { None } else { Some(id.to_string()) };

    let vcf_fields = VcfFields {
        chrom: chrom.to_string(),
        pos,
        id: id.to_string(),
        ref_allele: ref_str.to_string(),
        alt: alt_str.to_string(),
        qual: qual.to_string(),
        filter: filter.to_string(),
        info: info.to_string(),
        rest,
    };

    // Parse SV-related INFO fields for structural variants
    let (sv_end, sv_len, variant_type) = if has_symbolic {
        let info_map = parse_info_field(info);
        let sv_end = info_map
            .get("END")
            .and_then(|v| v.parse::<u64>().ok());
        let sv_len = info_map
            .get("SVLEN")
            .and_then(|v| v.parse::<i64>().ok());
        let svtype = info_map.get("SVTYPE").map(|s| s.as_str());

        let vtype = classify_sv_type(svtype, &alt_allele_strs);
        (sv_end, sv_len, vtype)
    } else {
        let vtype = classify_small_variant(&ref_allele, &alt_alleles);
        (None, None, vtype)
    };

    // For SVs, use END from INFO to set the genomic end coordinate
    let final_end = if has_symbolic {
        sv_end.unwrap_or(end)
    } else {
        end
    };

    Ok(VariationFeature {
        position: GenomicPosition::new(chrom, start, final_end, Strand::Forward),
        allele_string,
        ref_allele,
        alt_alleles,
        variation_name,
        vcf_fields: Some(vcf_fields),
        transcript_variations: Vec::new(),
        existing_variants: Vec::new(),
        minimised: false,
        most_severe_consequence: None,
        variant_type,
        sv_end,
        sv_len,
        supplementary_annotations: Vec::new(),
        gene_annotations: Vec::new(),
    })
}

/// Parse VCF INFO field into key-value pairs.
fn parse_info_field(info: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for pair in info.split(';') {
        if let Some((key, value)) = pair.split_once('=') {
            map.insert(key.to_string(), value.to_string());
        }
    }
    map
}

/// Classify structural variant type from SVTYPE INFO field or symbolic allele.
fn classify_sv_type(svtype: Option<&str>, alts: &[String]) -> VariantType {
    // Prefer SVTYPE from INFO, fall back to parsing ALT
    let sv = svtype.unwrap_or_else(|| {
        alts.first()
            .map(|a| a.trim_matches(|c| c == '<' || c == '>'))
            .unwrap_or("")
    });

    match sv.to_uppercase().as_str() {
        "DEL" => VariantType::CopyNumberLoss,
        "DUP" | "DUP:TANDEM" => VariantType::TandemDuplication,
        "INV" => VariantType::Inversion,
        "BND" => VariantType::TranslocationBreakend,
        "INS" => VariantType::Insertion,
        "CNV" => VariantType::CopyNumberVariation,
        "STR" => VariantType::ShortTandemRepeatVariation,
        s if s.starts_with("CN") => {
            // <CN0>, <CN1> = loss; <CN3>, <CN4> = gain
            if let Ok(cn) = s.trim_start_matches("CN").parse::<u32>() {
                if cn < 2 { VariantType::CopyNumberLoss }
                else if cn > 2 { VariantType::CopyNumberGain }
                else { VariantType::CopyNumberVariation }
            } else {
                VariantType::CopyNumberVariation
            }
        }
        _ => VariantType::Unknown,
    }
}

/// Classify small variant type from alleles.
fn classify_small_variant(ref_allele: &Allele, alt_alleles: &[Allele]) -> VariantType {
    if alt_alleles.is_empty() {
        return VariantType::Unknown;
    }
    let first_alt = &alt_alleles[0];

    match (ref_allele, first_alt) {
        (Allele::Deletion, Allele::Sequence(_)) => VariantType::Insertion,
        (Allele::Sequence(_), Allele::Deletion) => VariantType::Deletion,
        (Allele::Sequence(r), Allele::Sequence(a)) => {
            if r.len() == 1 && a.len() == 1 {
                VariantType::Snv
            } else if r.len() == a.len() {
                VariantType::Mnv
            } else {
                VariantType::Indel
            }
        }
        (_, Allele::Symbolic(_)) => VariantType::Unknown, // Handled by classify_sv_type
        _ => VariantType::Unknown,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_snv() {
        let line = "1\t100\trs1\tA\tG\t.\tPASS\t.\t";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.chromosome, "1");
        assert_eq!(vf.position.start, 100);
        assert_eq!(vf.position.end, 100);
        assert_eq!(vf.allele_string, "A/G");
        assert_eq!(vf.ref_allele, Allele::Sequence(b"A".to_vec()));
        assert_eq!(vf.alt_alleles, vec![Allele::Sequence(b"G".to_vec())]);
        assert_eq!(vf.variation_name, Some("rs1".to_string()));
    }

    #[test]
    fn test_parse_insertion() {
        // VCF: ref=A, alt=ATCG at pos 100 → Ensembl: ref=-, alt=TCG at pos 101, end=100
        let line = "1\t100\t.\tA\tATCG\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 101);
        assert_eq!(vf.position.end, 100); // insertion
        assert_eq!(vf.allele_string, "-/TCG");
        assert_eq!(vf.ref_allele, Allele::Deletion);
        assert_eq!(vf.alt_alleles, vec![Allele::Sequence(b"TCG".to_vec())]);
        assert!(vf.is_insertion());
    }

    #[test]
    fn test_parse_deletion() {
        // VCF: ref=ATCG, alt=A at pos 100 → Ensembl: ref=TCG, alt=- at pos 101
        let line = "1\t100\t.\tATCG\tA\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 101);
        assert_eq!(vf.position.end, 103);
        assert_eq!(vf.allele_string, "TCG/-");
        assert_eq!(vf.ref_allele, Allele::Sequence(b"TCG".to_vec()));
        assert_eq!(vf.alt_alleles, vec![Allele::Deletion]);
        assert!(vf.is_deletion());
    }

    #[test]
    fn test_parse_multi_allelic_indel() {
        // Multi-allelic: ref=ACG, alt=A,ACGT → strip shared A
        let line = "1\t100\t.\tACG\tA,ACGT\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 101);
        assert_eq!(vf.allele_string, "CG/-/CGT");
    }

    #[test]
    fn test_parse_multi_allelic_snv() {
        let line = "1\t100\t.\tA\tG,T\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 100);
        assert_eq!(vf.position.end, 100);
        assert_eq!(vf.allele_string, "A/G/T");
    }

    #[test]
    fn test_parse_mnv() {
        let line = "1\t100\t.\tAC\tGT\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 100);
        assert_eq!(vf.position.end, 101);
        assert_eq!(vf.allele_string, "AC/GT");
    }

    #[test]
    fn test_parse_non_variant() {
        let line = "1\t100\t.\tA\t.\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.allele_string, "A");
        assert!(vf.alt_alleles.is_empty() || vf.alt_alleles[0] == Allele::from_str("."));
    }

    #[test]
    fn test_vcf_parser_multiple_lines() {
        let vcf = "##fileformat=VCFv4.2\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                    1\t100\trs1\tA\tG\t.\tPASS\t.\n\
                    1\t200\trs2\tC\tT\t.\tPASS\t.\n";
        let mut parser = VcfParser::new(vcf.as_bytes()).unwrap();
        assert_eq!(parser.header_lines().len(), 2);
        let variants = parser.read_all().unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].position.start, 100);
        assert_eq!(variants[1].position.start, 200);
    }

    #[test]
    fn test_star_allele_preserved() {
        // Star allele in multi-allelic: ref=ACG, alt=A,* → strip A from non-star
        let line = "1\t100\t.\tACG\tA,*\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 101);
        assert_eq!(vf.allele_string, "CG/-/*");
    }

    #[test]
    fn test_parse_sv_deletion() {
        let line = "chr1\t10000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20000;SVLEN=-10000";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.position.start, 10000);
        assert_eq!(vf.position.end, 20000); // Uses END from INFO
        assert_eq!(vf.variant_type, VariantType::CopyNumberLoss);
        assert_eq!(vf.sv_end, Some(20000));
        assert_eq!(vf.sv_len, Some(-10000));
        assert!(vf.alt_alleles[0].is_symbolic());
    }

    #[test]
    fn test_parse_sv_duplication() {
        let line = "chr2\t5000\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=15000";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::TandemDuplication);
        assert_eq!(vf.sv_end, Some(15000));
        assert_eq!(vf.position.end, 15000);
    }

    #[test]
    fn test_parse_sv_inversion() {
        let line = "chr3\t1000\t.\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=5000";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::Inversion);
    }

    #[test]
    fn test_parse_sv_breakend() {
        let line = "chr1\t12345\t.\tN\t<BND>\t.\tPASS\tSVTYPE=BND";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::TranslocationBreakend);
    }

    #[test]
    fn test_parse_sv_cnv() {
        let line = "chr1\t100\t.\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=500";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::CopyNumberVariation);
    }

    #[test]
    fn test_small_variant_type_classification() {
        // SNV
        let line = "chr1\t100\t.\tA\tG\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::Snv);

        // Deletion
        let line = "chr1\t100\t.\tAC\tA\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::Deletion);

        // Insertion
        let line = "chr1\t100\t.\tA\tACG\t.\tPASS\t.";
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(vf.variant_type, VariantType::Insertion);
    }
}
