//! Source-specific parsers for building annotation databases.
//!
//! Each submodule implements a parser for a specific data source
//! (ClinVar, gnomAD, dbSNP, etc.) that produces `AnnotationRecord`s
//! for the `SaWriter`.

pub mod clinvar;
pub mod dbsnp;
pub mod gnomad;
