//! Constants and common types for the fastSA binary annotation format.

/// Magic bytes for position/allele-level annotation files (.osa).
pub const OSA_MAGIC: &[u8; 8] = b"FSTSA_01";

/// Magic bytes for interval-level annotation files (.osi).
pub const OSI_MAGIC: &[u8; 8] = b"FSTSI_01";

/// Magic bytes for gene-level annotation files (.oga).
pub const OGA_MAGIC: &[u8; 8] = b"FSTGA_01";

/// Current schema version. Bump when the binary format changes.
pub const SCHEMA_VERSION: u16 = 1;

/// Default block size for compression (8 MiB).
pub const DEFAULT_BLOCK_SIZE: usize = 8 * 1024 * 1024;

/// Default zstd compression level (3 is a good speed/ratio tradeoff).
pub const ZSTD_LEVEL: i32 = 3;

/// File extension for position/allele-level annotations.
pub const OSA_EXT: &str = "osa";

/// File extension for the index file.
pub const IDX_EXT: &str = "osa.idx";

/// File extension for interval-level annotations.
pub const OSI_EXT: &str = "osi";

/// File extension for gene-level annotations.
pub const OGA_EXT: &str = "oga";

/// A single annotation record ready for writing.
#[derive(Debug, Clone)]
pub struct AnnotationRecord {
    /// Chromosome index (numeric, mapped externally).
    pub chrom_idx: u16,
    /// 1-based genomic position.
    pub position: u32,
    /// Reference allele (empty string for positional annotations).
    pub ref_allele: String,
    /// Alternate allele (empty string for positional annotations).
    pub alt_allele: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// A single interval annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IntervalRecord {
    /// Chromosome name.
    pub chrom: String,
    /// 1-based start position (inclusive).
    pub start: u32,
    /// 1-based end position (inclusive).
    pub end: u32,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// Chromosome name to index mapping for efficient lookups.
#[derive(Debug, Clone)]
pub struct ChromMap {
    name_to_idx: std::collections::HashMap<String, u16>,
}

impl ChromMap {
    /// Create a standard human chromosome mapping (chr1-22, chrX, chrY, chrM).
    pub fn standard_human() -> Self {
        let mut map = std::collections::HashMap::new();
        for i in 1..=22u16 {
            map.insert(format!("chr{}", i), i - 1);
            map.insert(format!("{}", i), i - 1);
        }
        map.insert("chrX".into(), 22);
        map.insert("X".into(), 22);
        map.insert("chrY".into(), 23);
        map.insert("Y".into(), 23);
        map.insert("chrM".into(), 24);
        map.insert("MT".into(), 24);
        map.insert("chrMT".into(), 24);
        Self { name_to_idx: map }
    }

    /// Look up a chromosome index by name.
    #[inline]
    pub fn get(&self, name: &str) -> Option<u16> {
        self.name_to_idx.get(name).copied()
    }
}

/// A single gene annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneRecord {
    /// Gene symbol (e.g., "BRCA1").
    pub gene_symbol: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}
