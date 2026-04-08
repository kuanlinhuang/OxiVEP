//! OxiSA: Supplementary annotation format for OxiVEP.
//!
//! Two format generations:
//! - **v1 (.osa)**: Zstd block-compressed with per-entry JSON strings
//! - **v2 (.osa2)**: ZIP-based chunked format with Var32 encoding, parallel
//!   u32 value arrays, delta encoding, and LRU chunk caching (echtvar-inspired)
//!
//! Additional formats:
//! - **`.osi`** — Interval-level annotations (SV databases, regulatory regions)
//! - **`.oga`** — Gene-level annotations (OMIM, pLI scores, ClinGen)

pub mod block;
pub mod bloom;
pub mod chunk;
pub mod common;
pub mod custom;
pub mod fields;
pub mod gene;
pub mod index;
pub mod interval;
pub mod kmer16;
pub mod reader;
pub mod reader_v2;
pub mod sources;
pub mod var32;
pub mod writer;
pub mod writer_v2;
pub mod zigzag;
