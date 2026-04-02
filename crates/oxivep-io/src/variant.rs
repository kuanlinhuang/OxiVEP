use oxivep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
/// A variant feature ready for annotation.
#[derive(Debug, Clone)]
pub struct VariationFeature {
    pub position: GenomicPosition,
    /// Allele string in Ensembl format: "REF/ALT1/ALT2"
    pub allele_string: String,
    /// The reference allele after normalization.
    pub ref_allele: Allele,
    /// Alternative alleles after normalization.
    pub alt_alleles: Vec<Allele>,
    /// Variant ID (e.g., rs number) from VCF ID column.
    pub variation_name: Option<String>,
    /// Original VCF line for pass-through output.
    pub vcf_line: Option<String>,
    /// Original VCF fields for reconstruction.
    pub vcf_fields: Option<VcfFields>,
    /// Transcript-level annotations (populated during annotation).
    pub transcript_variations: Vec<TranscriptVariation>,
    /// Co-located known variants (populated during annotation).
    pub existing_variants: Vec<KnownVariant>,
    /// Whether the alleles were minimised.
    pub minimised: bool,
    /// Most severe consequence across all transcripts/alleles.
    pub most_severe_consequence: Option<Consequence>,
}

/// Parsed VCF fields for output reconstruction.
#[derive(Debug, Clone)]
pub struct VcfFields {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub rest: Vec<String>,
}

/// Annotation of a variant allele against a specific transcript.
#[derive(Debug, Clone)]
pub struct TranscriptVariation {
    pub transcript_id: String,
    pub gene_id: String,
    pub gene_symbol: Option<String>,
    pub biotype: String,
    pub allele_annotations: Vec<AlleleAnnotation>,
    pub canonical: bool,
    pub strand: Strand,
    pub source: Option<String>,
    pub protein_id: Option<String>,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub symbol_source: Option<String>,
    pub hgnc_id: Option<String>,
    /// Flags like "cds_end_NF", "cds_start_NF"
    pub flags: Vec<String>,
}

/// Annotation for a specific allele against a specific transcript.
#[derive(Debug, Clone)]
pub struct AlleleAnnotation {
    pub allele: Allele,
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub cdna_position: Option<(u64, u64)>,
    pub cds_position: Option<(u64, u64)>,
    pub protein_position: Option<(u64, u64)>,
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    pub exon: Option<(u32, u32)>,
    pub intron: Option<(u32, u32)>,
    pub distance: Option<i64>,
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    pub hgvsg: Option<String>,
    pub existing_variation: Vec<String>,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
}

/// A known/existing variant from the variation cache.
#[derive(Debug, Clone)]
pub struct KnownVariant {
    pub name: String,
    pub allele_string: Option<String>,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clinical_significance: Option<String>,
    pub somatic: bool,
    pub phenotype_or_disease: bool,
    pub pubmed: Vec<String>,
    pub frequencies: std::collections::HashMap<String, f64>,
}

impl VariationFeature {
    /// Compute the most severe consequence across all transcript annotations.
    pub fn compute_most_severe(&mut self) {
        let all_consequences: Vec<Consequence> = self
            .transcript_variations
            .iter()
            .flat_map(|tv| {
                tv.allele_annotations
                    .iter()
                    .flat_map(|aa| aa.consequences.iter().copied())
            })
            .collect();
        self.most_severe_consequence = Consequence::most_severe(&all_consequences);
    }

    /// Check if this is an insertion.
    pub fn is_insertion(&self) -> bool {
        self.ref_allele == Allele::Deletion
    }

    /// Check if this is a deletion.
    pub fn is_deletion(&self) -> bool {
        self.alt_alleles.iter().any(|a| *a == Allele::Deletion)
    }

    /// Check if this is an indel.
    pub fn is_indel(&self) -> bool {
        self.ref_allele == Allele::Deletion
            || self.alt_alleles.iter().any(|a| *a == Allele::Deletion)
            || self.alt_alleles.iter().any(|a| a.len() != self.ref_allele.len())
    }
}
