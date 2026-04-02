use anyhow::{Context, Result};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_core::region::Interval;
use noodles_csi::binning_index::BinningIndex;
use noodles_tabix as tabix;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// A parsed record from the VEP variation cache tabix file.
#[derive(Debug, Clone)]
pub struct VariationCacheRecord {
    pub variation_name: String,
    pub failed: bool,
    pub somatic: bool,
    pub start: u64,
    pub end: u64,
    pub allele_string: String,
    pub strand: i8,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    /// Population name → raw frequency string (e.g., "A:0.001,G:0.999")
    pub frequencies: HashMap<String, String>,
}

/// Reader for VEP's tabix-indexed variation cache files.
///
/// Each chromosome has its own `all_vars.gz` + `all_vars.gz.tbi` pair
/// under `{cache_dir}/{chr}/`.
pub struct VariationTabixReader {
    cache_dir: PathBuf,
    /// Column names from info.txt `variation_cols` (includes leading "chr" for tabix format).
    columns: Vec<String>,
    /// Chromosome name mapping: VCF name → cache name (e.g., "chr1" → "1").
    chrom_map: HashMap<String, String>,
    /// Valid chromosome names in the cache.
    _valid_chroms: Vec<String>,
}

impl VariationTabixReader {
    /// Create a new reader for a VEP cache directory.
    ///
    /// `variation_cols` comes from `info.txt` and defines the column order.
    /// `valid_chromosomes` comes from `info.txt` and lists chromosomes present in the cache.
    /// If `valid_chromosomes` is empty, auto-detects from subdirectories containing `all_vars.gz`.
    pub fn new(
        cache_dir: &Path,
        variation_cols: &[String],
        valid_chromosomes: &[String],
    ) -> Result<Self> {
        // If no valid_chromosomes provided, scan for subdirectories with all_vars.gz
        let chroms = if valid_chromosomes.is_empty() {
            let mut detected = Vec::new();
            if let Ok(entries) = std::fs::read_dir(cache_dir) {
                for entry in entries.flatten() {
                    if entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
                        let name = entry.file_name().to_string_lossy().to_string();
                        if entry.path().join("all_vars.gz").exists() {
                            detected.push(name);
                        }
                    }
                }
            }
            detected
        } else {
            valid_chromosomes.to_vec()
        };

        let mut chrom_map = HashMap::new();
        for vc in &chroms {
            // Map both "chr1" → "1" and "1" → "1"
            chrom_map.insert(vc.clone(), vc.clone());
            chrom_map.insert(format!("chr{}", vc), vc.clone());
            // Also map "1" → "1" if cache uses "chr" prefix
            if let Some(stripped) = vc.strip_prefix("chr") {
                chrom_map.insert(stripped.to_string(), vc.clone());
            }
        }

        Ok(Self {
            cache_dir: cache_dir.to_path_buf(),
            columns: variation_cols.to_vec(),
            chrom_map,
            _valid_chroms: chroms,
        })
    }

    /// Normalize a chromosome name from VCF to what the cache uses.
    pub fn normalize_chrom(&self, chrom: &str) -> Option<&str> {
        self.chrom_map.get(chrom).map(|s| s.as_str())
    }

    /// Query the variation cache for a genomic region.
    ///
    /// Returns all variation records that overlap the region.
    pub fn query(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<VariationCacheRecord>> {
        let cache_chrom = match self.normalize_chrom(chrom) {
            Some(c) => c.to_string(),
            None => return Ok(Vec::new()), // Unknown chromosome
        };

        let vars_path = self.cache_dir.join(&cache_chrom).join("all_vars.gz");
        let tbi_path = self.cache_dir.join(&cache_chrom).join("all_vars.gz.tbi");

        if !vars_path.exists() || !tbi_path.exists() {
            log::debug!("No variation cache for chromosome {}", cache_chrom);
            return Ok(Vec::new());
        }

        // Load tabix index
        let index = tabix::fs::read(&tbi_path)
            .with_context(|| format!("Reading tabix index: {}", tbi_path.display()))?;

        // Find the reference sequence ID for this chromosome
        let header = index.header().context("Missing tabix header")?;
        let ref_names = header.reference_sequence_names();
        let ref_id = ref_names
            .iter()
            .position(|name| name == cache_chrom.as_str());
        let ref_id = match ref_id {
            Some(id) => id,
            None => {
                log::debug!("Chromosome {} not found in tabix index", cache_chrom);
                return Ok(Vec::new());
            }
        };

        // Build query interval (noodles Position is 1-based, non-zero)
        // VEP cache uses 1-based coordinates already
        let pos_start = Position::try_from(start as usize)
            .context("Invalid start position")?;
        let pos_end = Position::try_from(end as usize)
            .context("Invalid end position")?;
        let query_interval: Interval = (pos_start..=pos_end).into();

        // Get chunks from index
        let chunks = index.query(ref_id, query_interval)
            .context("Tabix query failed")?;

        if chunks.is_empty() {
            return Ok(Vec::new());
        }

        // Open bgzf reader and read matching lines
        let file = File::open(&vars_path)
            .with_context(|| format!("Opening variation file: {}", vars_path.display()))?;
        let mut reader = bgzf::io::Reader::new(file);
        let mut records = Vec::new();

        for chunk in &chunks {
            reader.seek(chunk.start())?;
            let mut line = String::new();

            loop {
                line.clear();
                let bytes_read = reader.read_line(&mut line)?;
                if bytes_read == 0 {
                    break;
                }

                let trimmed = line.trim_end();
                if trimmed.is_empty() {
                    continue;
                }

                // Parse the line and filter by position
                if let Ok(record) = self.parse_line(trimmed) {
                    // Stop if we've passed the query region
                    if record.start > end {
                        break;
                    }
                    if record.start <= end && record.end >= start {
                        records.push(record);
                    }
                }

                // Stop if we've read past the chunk end
                if reader.virtual_position() >= chunk.end() {
                    break;
                }
            }
        }

        Ok(records)
    }

    /// Parse a single tab-separated line from the variation cache.
    pub fn parse_line(&self, line: &str) -> Result<VariationCacheRecord> {
        let fields: Vec<&str> = line.split('\t').collect();
        let col_map = build_column_map(&self.columns, &fields);

        let variation_name = col_map
            .get("variation_name")
            .and_then(|v| non_dot(v))
            .unwrap_or_default()
            .to_string();

        let failed = col_map
            .get("failed")
            .and_then(|v| non_dot(v))
            .map(|v| v != "0" && !v.is_empty())
            .unwrap_or(false);

        let somatic = col_map
            .get("somatic")
            .and_then(|v| non_dot(v))
            .map(|v| v != "0" && !v.is_empty())
            .unwrap_or(false);

        let start: u64 = col_map
            .get("start")
            .and_then(|v| non_dot(v))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);

        let end: u64 = col_map
            .get("end")
            .and_then(|v| non_dot(v))
            .and_then(|v| v.parse().ok())
            .unwrap_or(start);

        let allele_string = col_map
            .get("allele_string")
            .and_then(|v| non_dot(v))
            .unwrap_or_default()
            .to_string();

        let strand: i8 = col_map
            .get("strand")
            .and_then(|v| non_dot(v))
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);

        let minor_allele = col_map
            .get("minor_allele")
            .and_then(|v| non_dot(v))
            .map(|v| v.to_string());

        let minor_allele_freq = col_map
            .get("minor_allele_freq")
            .and_then(|v| non_dot(v))
            .and_then(|v| v.parse().ok());

        let clin_sig = col_map
            .get("clin_sig")
            .and_then(|v| non_dot(v))
            .map(|v| v.to_string());

        let phenotype_or_disease = col_map
            .get("phenotype_or_disease")
            .and_then(|v| non_dot(v))
            .map(|v| v != "0" && !v.is_empty())
            .unwrap_or(false);

        let pubmed: Vec<String> = col_map
            .get("pubmed")
            .and_then(|v| non_dot(v))
            .map(|v| v.split(',').map(|s| s.to_string()).collect())
            .unwrap_or_default();

        // Collect population frequency columns
        let known_non_freq = [
            "chr", "variation_name", "failed", "somatic", "start", "end",
            "allele_string", "strand", "minor_allele", "minor_allele_freq",
            "clin_sig", "phenotype_or_disease", "pubmed",
            "clin_sig_allele", "clinical_impact", "var_synonyms",
        ];
        let mut frequencies = HashMap::new();
        for (col_name, value) in &col_map {
            if !known_non_freq.contains(&col_name.as_str()) {
                if let Some(v) = non_dot(value) {
                    if !v.is_empty() {
                        frequencies.insert(col_name.clone(), v.to_string());
                    }
                }
            }
        }

        Ok(VariationCacheRecord {
            variation_name,
            failed,
            somatic,
            start,
            end,
            allele_string,
            strand,
            minor_allele,
            minor_allele_freq,
            clin_sig,
            phenotype_or_disease,
            pubmed,
            frequencies,
        })
    }
}

/// Build a mapping from column name → field value for a single line.
fn build_column_map<'a>(columns: &[String], fields: &[&'a str]) -> HashMap<String, &'a str> {
    let mut map = HashMap::new();
    for (i, col) in columns.iter().enumerate() {
        if let Some(&val) = fields.get(i) {
            map.insert(col.clone(), val);
        }
    }
    map
}

/// Return None for "." (VEP's empty marker) or empty strings.
fn non_dot(s: &str) -> Option<&str> {
    if s == "." || s.is_empty() {
        None
    } else {
        Some(s)
    }
}

/// Parse a VEP frequency string like "A:0.001,G:0.999" into allele-frequency pairs.
pub fn parse_freq_string(freq_str: &str) -> Vec<(String, f64)> {
    let mut result = Vec::new();
    for pair in freq_str.split(',') {
        if let Some((allele, freq_s)) = pair.split_once(':') {
            if let Ok(freq) = freq_s.parse::<f64>() {
                result.insert(result.len(), (allele.to_string(), freq));
            }
        }
    }
    result
}

/// Complement a nucleotide base.
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

/// Reverse-complement a DNA string.
fn reverse_complement(seq: &str) -> String {
    seq.bytes().rev().map(|b| complement_base(b) as char).collect()
}

/// Check if a variant's alleles match a known variant record, accounting for
/// strand differences and multi-allelic sites.
///
/// Returns the matched alt allele from the record if found.
pub fn match_alleles(
    input_ref: &str,
    input_alt: &str,
    input_start: u64,
    input_end: u64,
    record: &VariationCacheRecord,
) -> Option<String> {
    // Position must overlap
    if record.start > input_end || record.end < input_start {
        return None;
    }

    // Parse the record's allele string (e.g., "C/T" or "C/T/G")
    let record_alleles: Vec<&str> = record.allele_string.split('/').collect();
    if record_alleles.len() < 2 {
        return None;
    }

    let record_ref = record_alleles[0];
    let record_alts = &record_alleles[1..];

    // Try direct match first
    for &record_alt in record_alts {
        if alleles_equal(input_ref, input_alt, record_ref, record_alt) {
            return Some(record_alt.to_string());
        }
    }

    // Try reverse complement if record is on opposite strand
    if record.strand == -1 {
        let rc_ref = reverse_complement(record_ref);
        for &record_alt in record_alts {
            let rc_alt = reverse_complement(record_alt);
            if alleles_equal(input_ref, input_alt, &rc_ref, &rc_alt) {
                return Some(record_alt.to_string());
            }
        }
    }

    None
}

/// Check if two ref/alt pairs represent the same variant.
fn alleles_equal(ref1: &str, alt1: &str, ref2: &str, alt2: &str) -> bool {
    // Handle deletion markers
    let r1 = if ref1 == "-" { "" } else { ref1 };
    let a1 = if alt1 == "-" { "" } else { alt1 };
    let r2 = if ref2 == "-" { "" } else { ref2 };
    let a2 = if alt2 == "-" { "" } else { alt2 };

    r1.eq_ignore_ascii_case(r2) && a1.eq_ignore_ascii_case(a2)
}

/// Extract the frequency for a specific allele from a population's frequency string.
///
/// Given a freq string like "A:0.001,G:0.999" and allele "A", returns Some(0.001).
pub fn get_allele_freq(freq_str: &str, allele: &str) -> Option<f64> {
    for pair in freq_str.split(',') {
        if let Some((a, f)) = pair.split_once(':') {
            if a.eq_ignore_ascii_case(allele) {
                return f.parse().ok();
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_non_dot() {
        assert_eq!(non_dot("."), None);
        assert_eq!(non_dot(""), None);
        assert_eq!(non_dot("rs123"), Some("rs123"));
    }

    #[test]
    fn test_parse_freq_string() {
        let freqs = parse_freq_string("A:0.001,G:0.999");
        assert_eq!(freqs.len(), 2);
        assert_eq!(freqs[0], ("A".to_string(), 0.001));
        assert_eq!(freqs[1], ("G".to_string(), 0.999));

        // Empty/dot
        let freqs = parse_freq_string("");
        assert!(freqs.is_empty());
    }

    #[test]
    fn test_get_allele_freq() {
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "A"), Some(0.001));
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "G"), Some(0.999));
        assert_eq!(get_allele_freq("A:0.001,G:0.999", "T"), None);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("A"), "T");
        assert_eq!(reverse_complement("ATG"), "CAT");
    }

    #[test]
    fn test_alleles_equal() {
        assert!(alleles_equal("C", "T", "C", "T"));
        assert!(alleles_equal("c", "t", "C", "T"));
        assert!(!alleles_equal("C", "T", "C", "G"));
        // Deletion markers
        assert!(alleles_equal("-", "A", "-", "A"));
    }

    #[test]
    fn test_match_alleles_snv() {
        let record = VariationCacheRecord {
            variation_name: "rs123".into(),
            failed: false,
            somatic: false,
            start: 100,
            end: 100,
            allele_string: "C/T".into(),
            strand: 1,
            minor_allele: Some("T".into()),
            minor_allele_freq: Some(0.01),
            clin_sig: None,
            phenotype_or_disease: false,
            pubmed: vec![],
            frequencies: HashMap::new(),
        };

        // Matching SNV
        assert_eq!(match_alleles("C", "T", 100, 100, &record), Some("T".into()));
        // Non-matching allele
        assert_eq!(match_alleles("C", "G", 100, 100, &record), None);
        // Non-overlapping position
        assert_eq!(match_alleles("C", "T", 200, 200, &record), None);
    }

    #[test]
    fn test_match_alleles_reverse_strand() {
        let record = VariationCacheRecord {
            variation_name: "rs456".into(),
            failed: false,
            somatic: false,
            start: 100,
            end: 100,
            allele_string: "G/A".into(),
            strand: -1,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            phenotype_or_disease: false,
            pubmed: vec![],
            frequencies: HashMap::new(),
        };

        // On minus strand, G/A complements to C/T
        assert_eq!(match_alleles("C", "T", 100, 100, &record), Some("A".into()));
    }

    #[test]
    fn test_match_alleles_multiallelic() {
        let record = VariationCacheRecord {
            variation_name: "rs789".into(),
            failed: false,
            somatic: false,
            start: 100,
            end: 100,
            allele_string: "C/T/G".into(),
            strand: 1,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            phenotype_or_disease: false,
            pubmed: vec![],
            frequencies: HashMap::new(),
        };

        assert_eq!(match_alleles("C", "T", 100, 100, &record), Some("T".into()));
        assert_eq!(match_alleles("C", "G", 100, 100, &record), Some("G".into()));
        assert_eq!(match_alleles("C", "A", 100, 100, &record), None);
    }

    #[test]
    fn test_parse_line() {
        let cols: Vec<String> = "chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,pubmed,AFR,AMR,EAS,EUR,SAS,AA,EA"
            .split(',')
            .map(String::from)
            .collect();

        let reader = VariationTabixReader::new(
            Path::new("/tmp/fake"),
            &cols,
            &["21".to_string()],
        ).unwrap();

        let line = "21\trs559462325\t.\t.\t8522406\t.\tG/A\t.\tA\t0.0002\t.\t.\t.\tA:0\tA:0\tA:0.001\tA:0\tA:0\t.\t.";
        let record = reader.parse_line(line).unwrap();

        assert_eq!(record.variation_name, "rs559462325");
        assert!(!record.failed);
        assert!(!record.somatic);
        assert_eq!(record.start, 8522406);
        assert_eq!(record.end, 8522406); // end defaults to start when "."
        assert_eq!(record.allele_string, "G/A");
        assert_eq!(record.strand, 1); // defaults to 1 when "."
        assert_eq!(record.minor_allele, Some("A".into()));
        assert_eq!(record.minor_allele_freq, Some(0.0002));
        assert!(record.clin_sig.is_none());
        assert!(!record.phenotype_or_disease);

        // Check population frequencies
        assert_eq!(record.frequencies.get("AFR"), Some(&"A:0".to_string()));
        assert_eq!(record.frequencies.get("EAS"), Some(&"A:0.001".to_string()));
    }

    #[test]
    fn test_tabix_query_real_data() {
        // Use real VEP cache test data
        let cache_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent().unwrap()
            .parent().unwrap()
            .join("test_data/vep_cache/homo_sapiens/84_GRCh38");

        if !cache_dir.join("info.txt").exists() {
            eprintln!("Skipping tabix test: no test data at {:?}", cache_dir);
            return;
        }

        let info = crate::info::CacheInfo::from_file(&cache_dir.join("info.txt")).unwrap();

        let reader = VariationTabixReader::new(
            &cache_dir,
            &info.variation_cols,
            &info.valid_chromosomes,
        ).unwrap();

        // Query for a known variant: rs879182429 at chr21:8013350
        let records = reader.query("21", 8013350, 8013350).unwrap();
        assert!(!records.is_empty(), "Expected at least one record at 21:8013350");

        let rs = records.iter().find(|r| r.variation_name == "rs879182429");
        assert!(rs.is_some(), "Expected rs879182429 in results");
        let rs = rs.unwrap();
        assert_eq!(rs.allele_string, "T/G");
        assert_eq!(rs.start, 8013350);

        // Query for a variant with frequency data: rs559462325 at chr21:8522406
        let records = reader.query("21", 8522406, 8522406).unwrap();
        let rs = records.iter().find(|r| r.variation_name == "rs559462325");
        assert!(rs.is_some(), "Expected rs559462325 in results");
        let rs = rs.unwrap();
        assert_eq!(rs.minor_allele, Some("A".to_string()));
        assert_eq!(rs.minor_allele_freq, Some(0.0002));
        assert!(rs.frequencies.contains_key("AFR"));

        // Query with "chr" prefix should also work
        let records = reader.query("chr21", 8013350, 8013350).unwrap();
        assert!(!records.is_empty(), "chr21 prefix should normalize to 21");

        // Query for a non-existent region should return empty
        let records = reader.query("21", 1, 1).unwrap();
        assert!(records.is_empty(), "No variants expected at position 1");
    }
}
