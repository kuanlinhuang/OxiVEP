mod coding;
mod genomic;
mod protein;

pub use coding::{hgvsc, hgvsc_with_seq, hgvsc_intronic, hgvsc_intronic_range, hgvsc_noncoding, hgvsc_noncoding_intronic, hgvsc_noncoding_intronic_range};
pub use genomic::hgvsg;
pub use protein::{hgvsp, hgvsp_frameshift};

/// Full HGVS annotation result.
#[derive(Debug, Clone, Default)]
pub struct HgvsAnnotation {
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsg: Option<String>,
}
