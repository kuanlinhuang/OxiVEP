//! Reader for .osa2 format (ZIP-based chunked annotation files).
//!
//! Implements the `AnnotationProvider` trait with O(log n) lookups via
//! Var32 binary search on sorted genomic chunks.

use crate::chunk::{delta_decode, Chunk};
use crate::fields::{Field, FieldType};
use crate::kmer16::LongVariant;
use crate::var32;
use crate::writer_v2::{read_u32_array, Osa2Metadata};
use anyhow::{Context, Result};
use lru::LruCache;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use std::cell::UnsafeCell;
use std::fs::File;
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;

/// Reader for .osa2 annotation files.
///
/// Loads genomic chunks on demand from a ZIP archive, caches recently used
/// chunks in an LRU cache, and performs binary search for variant lookups.
pub struct Osa2Reader {
    /// Path to the .osa2 ZIP file (re-opened for each chunk load).
    zip_path: std::path::PathBuf,
    metadata: Osa2Metadata,
    sa_metadata: SaMetadata,
    fields: Vec<Field>,
    /// Categorical string lookup tables per field.
    string_tables: Vec<Vec<String>>,
    /// LRU cache of loaded chunks, keyed by "chrom/chunk_id".
    /// Written during preload/query (sequential access pattern), read during annotation.
    chunk_cache: UnsafeCell<LruCache<String, Chunk>>,
}

// SAFETY: chunk_cache is written only in single-threaded preload phase
// and individual query paths that are protected by the batch pattern.
unsafe impl Send for Osa2Reader {}
unsafe impl Sync for Osa2Reader {}

impl Osa2Reader {
    /// Open an .osa2 file.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
        let mut archive = zip::ZipArchive::new(file)?;

        // Read metadata
        let metadata: Osa2Metadata = {
            let mut entry = archive.by_name("fastsa/metadata.json")
                .context("Missing fastsa/metadata.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        // Read field config
        let fields: Vec<Field> = {
            let mut entry = archive.by_name("fastsa/config.json")
                .context("Missing fastsa/config.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        // Read string tables
        let mut string_tables: Vec<Vec<String>> = fields.iter().map(|_| Vec::new()).collect();
        for (i, field) in fields.iter().enumerate() {
            if field.ftype == FieldType::Categorical {
                let name = format!("fastsa/strings/{}.txt", field.alias);
                if let Ok(mut entry) = archive.by_name(&name) {
                    let mut buf = String::new();
                    entry.read_to_string(&mut buf)?;
                    string_tables[i] = buf.lines().map(|l| l.to_string()).collect();
                }
            }
        }

        let sa_metadata = SaMetadata {
            name: metadata.name.clone(),
            version: metadata.version.clone(),
            description: metadata.description.clone(),
            assembly: metadata.assembly.clone(),
            json_key: metadata.json_key.clone(),
            match_by_allele: metadata.match_by_allele,
            is_array: metadata.is_array,
            is_positional: metadata.is_positional,
        };

        Ok(Self {
            zip_path: path.to_path_buf(),
            metadata,
            sa_metadata,
            fields,
            string_tables,
            chunk_cache: UnsafeCell::new(LruCache::new(NonZeroUsize::new(8).unwrap())),
        })
    }

    /// Load a chunk from the ZIP archive into the LRU cache.
    fn load_chunk(&self, chrom: &str, chunk_id: u32) -> Result<()> {
        let cache_key = format!("{}/{}", chrom, chunk_id);
        let cache = unsafe { &mut *self.chunk_cache.get() };
        if cache.contains(&cache_key) {
            return Ok(());
        }

        let file = File::open(&self.zip_path)?;
        let mut archive = zip::ZipArchive::new(file)?;
        let prefix = format!("fastsa/{}/{}/", chrom, chunk_id);

        // Read var32 keys
        let var32s = {
            let name = format!("{}var32.bin", prefix);
            match archive.by_name(&name) {
                Ok(mut entry) => {
                    let mut buf = Vec::new();
                    entry.read_to_end(&mut buf)?;
                    let mut keys = read_u32_array(&buf)?;
                    delta_decode(&mut keys); // Reconstruct from deltas
                    keys
                }
                Err(_) => {
                    // Chunk doesn't exist in this archive
                    cache.put(cache_key, Chunk::empty());
                    return Ok(());
                }
            }
        };

        // Read long variants
        let longs: Vec<LongVariant> = {
            let name = format!("{}too-long.enc", prefix);
            match archive.by_name(&name) {
                Ok(mut entry) => {
                    let mut buf = Vec::new();
                    entry.read_to_end(&mut buf)?;
                    bincode::deserialize(&buf).unwrap_or_default()
                }
                Err(_) => Vec::new(),
            }
        };

        // Read parallel value arrays
        let mut values = Vec::new();
        for field in &self.fields {
            if field.ftype == FieldType::JsonBlob {
                continue;
            }
            let name = format!("{}{}.bin", prefix, field.alias);
            match archive.by_name(&name) {
                Ok(mut entry) => {
                    let mut buf = Vec::new();
                    entry.read_to_end(&mut buf)?;
                    values.push(read_u32_array(&buf)?);
                }
                Err(_) => {
                    // Fill with missing values
                    values.push(vec![field.missing_value; var32s.len()]);
                }
            }
        }

        // Read JSON blobs if any
        let json_blobs = {
            let name = format!("{}json_blobs.zst", prefix);
            match archive.by_name(&name) {
                Ok(mut entry) => {
                    let mut buf = Vec::new();
                    entry.read_to_end(&mut buf)?;
                    let decompressed = zstd::decode_all(buf.as_slice())?;
                    let text = String::from_utf8(decompressed)?;
                    Some(text.lines().map(|l| l.to_string()).collect())
                }
                Err(_) => None,
            }
        };

        cache.put(cache_key, Chunk { var32s, longs, values, json_blobs });
        Ok(())
    }

    /// Query a variant in the loaded chunks.
    fn query(&self, chrom: &str, pos: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Result<Option<String>> {
        let chunk_id = pos >> self.metadata.chunk_bits;
        let cache_key = format!("{}/{}", chrom, chunk_id);

        // Ensure chunk is loaded
        self.load_chunk(chrom, chunk_id)?;

        let cache = unsafe { &mut *self.chunk_cache.get() };
        let chunk = match cache.get(&cache_key) {
            Some(c) => c,
            None => return Ok(None),
        };

        if chunk.is_empty() {
            return Ok(None);
        }

        // Try Var32 lookup first
        let chunk_mask = (1u32 << self.metadata.chunk_bits) - 1;
        let within_pos = pos & chunk_mask;

        let idx = if var32::is_long(ref_allele.len(), alt_allele.len()) {
            chunk.find_long(pos, ref_allele, alt_allele)
        } else {
            var32::encode(within_pos, ref_allele, alt_allele)
                .and_then(|key| chunk.find_short(key))
        };

        match idx {
            Some(i) => {
                let json = chunk.reconstruct_json(i, &self.fields, &self.string_tables);
                Ok(Some(json))
            }
            None => Ok(None),
        }
    }
}

impl AnnotationProvider for Osa2Reader {
    fn name(&self) -> &str {
        &self.sa_metadata.name
    }

    fn json_key(&self) -> &str {
        &self.sa_metadata.json_key
    }

    fn metadata(&self) -> &SaMetadata {
        &self.sa_metadata
    }

    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        let ref_bytes = ref_allele.as_bytes();
        let alt_bytes = alt_allele.as_bytes();

        match self.query(chrom, pos as u32, ref_bytes, alt_bytes)? {
            Some(json) => {
                if self.sa_metadata.is_positional {
                    Ok(Some(AnnotationValue::Positional(json)))
                } else {
                    Ok(Some(AnnotationValue::Json(json)))
                }
            }
            None => Ok(None),
        }
    }

    fn preload(&self, chrom: &str, positions: &[u64]) -> Result<()> {
        if positions.is_empty() {
            return Ok(());
        }

        // Determine which chunks need to be loaded
        let mut chunk_ids: Vec<u32> = positions
            .iter()
            .map(|&p| (p as u32) >> self.metadata.chunk_bits)
            .collect();
        chunk_ids.sort_unstable();
        chunk_ids.dedup();

        for cid in chunk_ids {
            self.load_chunk(chrom, cid)?;
        }

        Ok(())
    }
}
