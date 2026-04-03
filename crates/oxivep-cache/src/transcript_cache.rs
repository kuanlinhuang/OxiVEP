//! Binary transcript cache for fast startup.
//!
//! Serializes fully-built `Vec<Transcript>` (including spliced sequences)
//! to a compact binary format using bincode + gzip compression.
//! Subsequent loads skip GFF3 parsing, FASTA loading, and sequence building.

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use oxivep_genome::Transcript;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::time::SystemTime;

const CACHE_MAGIC: &[u8; 8] = b"OXIVEP01";

/// Save transcripts to a binary cache file (bincode + gzip).
pub fn save_cache(transcripts: &[Transcript], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Creating cache file: {}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut gz = GzEncoder::new(writer, Compression::fast());

    // Write magic header
    use std::io::Write;
    gz.write_all(CACHE_MAGIC)?;

    // Serialize with bincode
    bincode::serialize_into(&mut gz, transcripts)
        .with_context(|| "Serializing transcripts to cache")?;

    gz.finish()?;
    Ok(())
}

/// Load transcripts from a binary cache file.
pub fn load_cache(path: &Path) -> Result<Vec<Transcript>> {
    let file = File::open(path)
        .with_context(|| format!("Opening cache file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut gz = GzDecoder::new(reader);

    // Verify magic header
    use std::io::Read;
    let mut magic = [0u8; 8];
    gz.read_exact(&mut magic)
        .with_context(|| "Reading cache header")?;
    if &magic != CACHE_MAGIC {
        anyhow::bail!("Invalid cache file (wrong magic header)");
    }

    // Deserialize with bincode
    let transcripts: Vec<Transcript> = bincode::deserialize_from(&mut gz)
        .with_context(|| "Deserializing transcripts from cache")?;

    Ok(transcripts)
}

/// Check if cache file is newer than source file.
pub fn cache_is_fresh(cache_path: &Path, source_path: &Path) -> bool {
    let cache_mtime = cache_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::UNIX_EPOCH);
    let source_mtime = source_path
        .metadata()
        .and_then(|m| m.modified())
        .unwrap_or(SystemTime::now());
    cache_mtime > source_mtime
}

/// Get the default cache path for a given GFF3 path.
pub fn default_cache_path(gff3_path: &Path) -> std::path::PathBuf {
    let mut cache_path = gff3_path.to_path_buf();
    let name = cache_path
        .file_name()
        .map(|n| format!("{}.oxivep.cache", n.to_string_lossy()))
        .unwrap_or_else(|| "transcripts.oxivep.cache".to_string());
    cache_path.set_file_name(name);
    cache_path
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxivep_core::Strand;
    use oxivep_genome::{Exon, Gene, Transcript};
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn make_test_transcript() -> Transcript {
        Transcript {
            stable_id: Arc::from("ENST00000001"),
            version: Some(1),
            gene: Gene {
                stable_id: Arc::from("ENSG00000001"),
                symbol: Some(Arc::from("TEST")),
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from("protein_coding"),
                chromosome: Arc::from("1"),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: Arc::from("protein_coding"),
            chromosome: Arc::from("1"),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE001".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 1,
                },
            ],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(200),
            coding_region_start: Some(1000),
            coding_region_end: Some(1200),
            spliced_seq: Some("ACGTACGT".into()),
            translateable_seq: Some("ACGT".into()),
            peptide: Some("T".into()),
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: Some("P1".into()),
            ccds: None,
            protein_id: Some("ENSP001".into()),
            protein_version: Some(1),
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
    fn test_cache_roundtrip() {
        let transcripts = vec![make_test_transcript()];
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        save_cache(&transcripts, path).unwrap();
        let loaded = load_cache(path).unwrap();

        assert_eq!(loaded.len(), 1);
        assert_eq!(&*loaded[0].stable_id, "ENST00000001");
        assert_eq!(&**loaded[0].gene.symbol.as_ref().unwrap(), "TEST");
        assert_eq!(loaded[0].spliced_seq.as_deref(), Some("ACGTACGT"));
        assert_eq!(loaded[0].canonical, true);
        assert_eq!(loaded[0].tsl, Some(1));
    }

    #[test]
    fn test_invalid_magic() {
        let tmp = NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"NOTVALID").unwrap();
        assert!(load_cache(tmp.path()).is_err());
    }
}
