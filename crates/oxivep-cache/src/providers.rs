use anyhow::Result;
use oxivep_genome::Transcript;
use std::collections::HashMap;
use std::sync::Arc;

use crate::info::CacheInfo;
use crate::variation::{self, VariationTabixReader};

/// Trait for providing transcript annotations for a genomic region.
pub trait TranscriptProvider: Send + Sync {
    /// Return all transcripts that overlap the given genomic region.
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>>;

    /// Return all transcripts on a chromosome.
    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>>;
}

/// Trait for providing reference sequences.
pub trait SequenceProvider: Send + Sync {
    /// Fetch reference sequence for a region (1-based, inclusive).
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>>;

    /// Fetch a reference sequence slice. Default delegates to fetch_sequence.
    fn fetch_sequence_slice(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.fetch_sequence(chrom, start, end)
    }
}

/// In-memory transcript provider backed by a Vec<Transcript>.
pub struct MemoryTranscriptProvider {
    transcripts: Vec<Transcript>,
}

impl MemoryTranscriptProvider {
    pub fn new(transcripts: Vec<Transcript>) -> Self {
        Self { transcripts }
    }

    pub fn transcript_count(&self) -> usize {
        self.transcripts.len()
    }
}

impl TranscriptProvider for MemoryTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom && t.start <= end && t.end >= start)
            .collect())
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        Ok(self
            .transcripts
            .iter()
            .filter(|t| &*t.chromosome == chrom)
            .collect())
    }
}

/// High-performance transcript provider using per-chromosome sorted arrays
/// and binary search for O(log n + k) lookups instead of O(n) linear scans.
///
/// Inspired by echtvar's chunked binary search with same-region caching.
pub struct IndexedTranscriptProvider {
    /// Transcripts grouped by chromosome, sorted by start position within each group.
    by_chrom: HashMap<Arc<str>, Vec<Transcript>>,
}

impl IndexedTranscriptProvider {
    pub fn new(mut transcripts: Vec<Transcript>) -> Self {
        let mut by_chrom: HashMap<Arc<str>, Vec<Transcript>> = HashMap::new();
        for tr in transcripts.drain(..) {
            by_chrom
                .entry(Arc::clone(&tr.chromosome))
                .or_default()
                .push(tr);
        }
        // Sort each chromosome's transcripts by start position
        for trs in by_chrom.values_mut() {
            trs.sort_by_key(|t| t.start);
        }
        Self { by_chrom }
    }

    pub fn transcript_count(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }
}

impl TranscriptProvider for IndexedTranscriptProvider {
    fn get_transcripts(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<&Transcript>> {
        let trs = match self.by_chrom.get(chrom) {
            Some(trs) => trs,
            None => return Ok(Vec::new()),
        };

        // Binary search: find the first transcript whose start > end (query end).
        // All transcripts that could overlap must have start <= end, so they're in [0..upper).
        let upper = trs.partition_point(|t| t.start <= end);

        // From [0..upper), filter those whose end >= start (query start).
        // Since transcripts are sorted by start, scan backward from upper until
        // we're sure no more can overlap.
        let mut results = Vec::new();
        for t in trs[..upper].iter().rev() {
            // If this transcript ends before our query start, and since transcripts
            // are sorted by start (not end), we can't early-exit here — a transcript
            // that starts earlier might extend further. So we scan all of [0..upper).
            if t.end >= start {
                results.push(t);
            }
        }
        results.reverse(); // Restore start-position order
        Ok(results)
    }

    fn get_transcripts_by_chrom(&self, chrom: &str) -> Result<Vec<&Transcript>> {
        match self.by_chrom.get(chrom) {
            Some(trs) => Ok(trs.iter().collect()),
            None => Ok(Vec::new()),
        }
    }
}

/// Sequence provider backed by a FASTA reader.
pub struct FastaSequenceProvider {
    reader: crate::fasta::FastaReader,
}

impl FastaSequenceProvider {
    pub fn new(reader: crate::fasta::FastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for FastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.reader.fetch(chrom, start, end)
    }
}

/// Sequence provider backed by a memory-mapped FASTA reader.
/// Uses .fai index for random access without loading the full file into RAM.
pub struct MmapFastaSequenceProvider {
    reader: crate::fasta::MmapFastaReader,
}

impl MmapFastaSequenceProvider {
    pub fn new(reader: crate::fasta::MmapFastaReader) -> Self {
        Self { reader }
    }
}

impl SequenceProvider for MmapFastaSequenceProvider {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        self.reader.fetch(chrom, start, end)
    }
}

/// A matched known variant with its allele-specific frequency data.
#[derive(Debug, Clone)]
pub struct MatchedVariant {
    pub name: String,
    pub matched_allele: String,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub somatic: bool,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    /// Population → allele-specific frequency for the matched allele.
    pub frequencies: HashMap<String, f64>,
}

/// Trait for providing co-located known variant annotations.
pub trait VariationProvider {
    /// Look up known variants overlapping a position that match the given alleles.
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>>;
}

/// Variation provider backed by VEP's tabix-indexed cache files.
pub struct TabixVariationProvider {
    reader: VariationTabixReader,
}

impl TabixVariationProvider {
    /// Create a provider from a VEP cache directory.
    ///
    /// Reads `info.txt` for column definitions and valid chromosomes.
    pub fn new(cache_dir: &std::path::Path, cache_info: &CacheInfo) -> Result<Self> {
        let reader = VariationTabixReader::new(
            cache_dir,
            &cache_info.variation_cols,
            &cache_info.valid_chromosomes,
        )?;
        Ok(Self { reader })
    }
}

impl VariationProvider for TabixVariationProvider {
    fn get_matched_variants(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Vec<MatchedVariant>> {
        let records = self.reader.query(chrom, start, end)?;
        let mut matched = Vec::new();

        for record in &records {
            // Skip failed variants
            if record.failed {
                continue;
            }

            // Check allele match
            if let Some(matched_alt) = variation::match_alleles(
                ref_allele, alt_allele, start, end, record,
            ) {
                // Extract per-allele frequencies for the matched allele
                let mut freqs = HashMap::new();
                for (pop, freq_str) in &record.frequencies {
                    if let Some(f) = variation::get_allele_freq(freq_str, &matched_alt) {
                        freqs.insert(pop.clone(), f);
                    }
                }

                // Also include MAF if minor_allele matches
                if let (Some(ref ma), Some(maf)) = (&record.minor_allele, record.minor_allele_freq) {
                    if ma.eq_ignore_ascii_case(&matched_alt) {
                        freqs.entry("minor_allele_freq".into()).or_insert(maf);
                    }
                }

                matched.push(MatchedVariant {
                    name: record.variation_name.clone(),
                    matched_allele: matched_alt,
                    minor_allele: record.minor_allele.clone(),
                    minor_allele_freq: record.minor_allele_freq,
                    clin_sig: record.clin_sig.clone(),
                    somatic: record.somatic,
                    phenotype_or_disease: record.phenotype_or_disease,
                    pubmed: record.pubmed.clone(),
                    frequencies: freqs,
                });
            }
        }

        Ok(matched)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxivep_core::Strand;
    use oxivep_genome::{Exon, Gene};

    fn make_transcript(chrom: &str, start: u64, end: u64) -> Transcript {
        Transcript {
            stable_id: Arc::from(format!("ENST_{}", start).as_str()),
            version: None,
            gene: Gene {
                stable_id: "ENSG_1".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: chrom.into(),
                start,
                end,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: chrom.into(),
            start,
            end,
            strand: Strand::Forward,
            exons: vec![Exon {
                stable_id: "ENSE_1".into(),
                start,
                end,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 1,
            }],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_memory_transcript_provider() {
        let provider = MemoryTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_indexed_transcript_provider() {
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 2000),
            make_transcript("chr1", 3000, 4000),
            make_transcript("chr2", 1000, 2000),
        ]);

        // Overlapping query
        let results = provider.get_transcripts("chr1", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Non-overlapping
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 0);

        // Different chromosome
        let results = provider.get_transcripts("chr2", 1500, 1600).unwrap();
        assert_eq!(results.len(), 1);

        // By chromosome
        let results = provider.get_transcripts_by_chrom("chr1").unwrap();
        assert_eq!(results.len(), 2);

        // Missing chromosome
        let results = provider.get_transcripts("chr99", 1, 100).unwrap();
        assert_eq!(results.len(), 0);

        // Query spanning two transcripts
        let results = provider.get_transcripts("chr1", 1500, 3500).unwrap();
        assert_eq!(results.len(), 2);

        // Query at exact boundaries
        let results = provider.get_transcripts("chr1", 2000, 2000).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 1000);

        // Transcript count
        assert_eq!(provider.transcript_count(), 3);
    }

    #[test]
    fn test_indexed_provider_overlapping_transcripts() {
        // Test with overlapping transcripts (common in real genomes)
        let provider = IndexedTranscriptProvider::new(vec![
            make_transcript("chr1", 1000, 5000),
            make_transcript("chr1", 2000, 3000),
            make_transcript("chr1", 4000, 6000),
        ]);

        // Query in the overlap region
        let results = provider.get_transcripts("chr1", 2500, 2600).unwrap();
        assert_eq!(results.len(), 2); // Both 1000-5000 and 2000-3000

        // Query spanning all three
        let results = provider.get_transcripts("chr1", 1000, 6000).unwrap();
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_fasta_sequence_provider() {
        let fasta = ">chr1\nACGTACGTAAAACCCC\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);
        let seq = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_fasta_fetch_slice_matches_fetch() {
        let fasta = ">chr1\nacgtACGTaaaa\n>chr2\nTTTTgggg\n";
        let reader = crate::fasta::FastaReader::from_reader(fasta.as_bytes()).unwrap();
        let provider = FastaSequenceProvider::new(reader);

        // fetch_slice returns same data as fetch (both uppercase)
        let slice = provider.fetch_sequence_slice("chr1", 1, 4).unwrap();
        let vec = provider.fetch_sequence("chr1", 1, 4).unwrap();
        assert_eq!(slice, vec.as_slice());
        assert_eq!(slice, b"ACGT");

        // Lowercase input is uppercased at load time
        let slice = provider.fetch_sequence_slice("chr2", 5, 8).unwrap();
        assert_eq!(slice, b"GGGG");
    }
}
