//! Reader for .osa position/allele-level annotation files.
//!
//! Uses memory-mapped I/O for the data file and binary search on the index
//! for O(1) block lookups. Supports preloading for batch annotation.

use crate::block::{BlockEntry, SaBlock};
use crate::common::{ChromMap, OSA_MAGIC};
use crate::index::SaIndex;
use anyhow::Result;
use memmap2::Mmap;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use std::cell::UnsafeCell;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

/// Reader for .osa annotation files.
///
/// Thread-safety: The reader is `Send + Sync`. Preloaded data is stored
/// in an `UnsafeCell<HashMap>` that is written to only during `preload()`
/// (single-threaded phase) and read during `annotate_position()` (parallel phase).
pub struct SaReader {
    mmap: Mmap,
    index: SaIndex,
    metadata: SaMetadata,
    chrom_map: ChromMap,
    /// Cache of decompressed entries keyed by (chrom_idx, block_start_pos).
    /// Written during preload (sequential), read during annotation (parallel).
    preloaded: UnsafeCell<HashMap<(u16, u32), Vec<BlockEntry>>>,
}

// SAFETY: preloaded is written only during preload() which runs sequentially
// before the parallel annotation phase, and read-only during annotate_position().
unsafe impl Send for SaReader {}
unsafe impl Sync for SaReader {}

impl SaReader {
    /// Open an .osa + .osa.idx file pair.
    pub fn open(data_path: &Path) -> Result<Self> {
        let idx_path = data_path.with_extension("osa.idx");

        // Read index
        let mut idx_file = File::open(&idx_path)?;
        let index = SaIndex::read_from(&mut idx_file)?;

        // Memory-map data file
        let data_file = File::open(data_path)?;
        let mmap = unsafe { Mmap::map(&data_file)? };

        // Verify data file magic
        if mmap.len() < 10 || &mmap[..8] != OSA_MAGIC {
            anyhow::bail!("Invalid OSA data file: bad magic");
        }

        let metadata = SaMetadata {
            name: index.header.name.clone(),
            version: index.header.version.clone(),
            description: index.header.description.clone(),
            assembly: index.header.assembly.clone(),
            json_key: index.header.json_key.clone(),
            match_by_allele: index.header.match_by_allele,
            is_array: index.header.is_array,
            is_positional: index.header.is_positional,
        };

        Ok(Self {
            mmap,
            index,
            metadata,
            chrom_map: ChromMap::standard_human(),
            preloaded: UnsafeCell::new(HashMap::new()),
        })
    }

    /// Read and decompress a block at the given file offset.
    fn read_block(&self, file_offset: u64, compressed_len: u32) -> Result<Vec<BlockEntry>> {
        let offset = file_offset as usize;
        // The data file stores: [4-byte compressed_len] [compressed_data]
        let data_start = offset + 4;
        let data_end = data_start + compressed_len as usize;

        if data_end > self.mmap.len() {
            anyhow::bail!("Block extends beyond data file");
        }

        SaBlock::decompress(&self.mmap[data_start..data_end])
    }

    /// Query annotations for a specific position and allele.
    /// First checks preloaded cache, then falls back to direct read.
    fn query(&self, chrom: &str, position: u32, ref_allele: &str, alt_allele: &str) -> Result<Option<String>> {
        // Check preloaded cache first (uses chrom index for fast lookup)
        let preloaded = unsafe { &*self.preloaded.get() };
        if let Some(chrom_idx) = self.chrom_map.get(chrom) {
            if let Some(entries) = preloaded.get(&(chrom_idx, position)) {
                return Ok(self.find_match(entries, position, ref_allele, alt_allele));
            }
        }

        // Fall back to direct block read
        let block_refs = self.index.find_blocks(chrom, position);
        for block_ref in block_refs {
            let entries = self.read_block(block_ref.file_offset, block_ref.compressed_len)?;
            if let Some(json) = self.find_match(&entries, position, ref_allele, alt_allele) {
                return Ok(Some(json));
            }
        }

        Ok(None)
    }

    fn find_match(&self, entries: &[BlockEntry], position: u32, ref_allele: &str, alt_allele: &str) -> Option<String> {
        let allele_ref = if self.metadata.match_by_allele { ref_allele } else { "" };
        let allele_alt = if self.metadata.match_by_allele { alt_allele } else { "" };

        SaBlock::find_by_position(entries, position, allele_ref, allele_alt, self.metadata.is_positional)
            .map(|idx| entries[idx].json.clone())
    }
}

impl AnnotationProvider for SaReader {
    fn name(&self) -> &str {
        &self.metadata.name
    }

    fn json_key(&self) -> &str {
        &self.metadata.json_key
    }

    fn metadata(&self) -> &SaMetadata {
        &self.metadata
    }

    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        match self.query(chrom, pos as u32, ref_allele, alt_allele)? {
            Some(json) => {
                if self.metadata.is_positional {
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

        let chrom_idx = self.chrom_map.get(chrom).unwrap_or(255);
        let cache = unsafe { &mut *self.preloaded.get() };
        cache.clear();

        let min_pos = *positions.iter().min().unwrap() as u32;
        let max_pos = *positions.iter().max().unwrap() as u32;

        let block_refs = self.index.find_blocks_range(chrom, min_pos, max_pos);

        for block_ref in block_refs {
            let entries = self.read_block(block_ref.file_offset, block_ref.compressed_len)?;
            for entry in entries {
                cache
                    .entry((chrom_idx, entry.position))
                    .or_insert_with(Vec::new)
                    .push(entry);
            }
        }

        Ok(())
    }
}
