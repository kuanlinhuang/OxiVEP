//! Structural variant consequence prediction.
//!
//! Predicts consequences for SVs (deletions, duplications, inversions,
//! CNVs, breakends) based on their overlap with transcripts.

use oxivep_core::{Allele, Consequence, Impact, Strand, VariantType};
use oxivep_genome::Transcript;
use std::sync::Arc;

use crate::predictor::{AlleleConsequenceResult, TranscriptConsequence};

/// Predict consequences of a structural variant against a set of transcripts.
///
/// The SV is defined by its genomic coordinates [start, end] and variant type.
/// For each transcript that overlaps, determines whether the SV causes
/// ablation, amplification, or partial disruption.
pub fn predict_sv_consequences(
    chrom: &str,
    start: u64,
    end: u64,
    variant_type: VariantType,
    alt_alleles: &[Allele],
    transcripts: &[&Transcript],
    upstream_distance: u64,
    downstream_distance: u64,
) -> Vec<TranscriptConsequence> {
    let mut results = Vec::new();

    for &transcript in transcripts {
        if *transcript.chromosome != *chrom {
            continue;
        }

        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            .map(|allele| {
                predict_sv_for_transcript(
                    start,
                    end,
                    variant_type,
                    allele,
                    transcript,
                    upstream_distance,
                    downstream_distance,
                )
            })
            .collect();

        if !allele_consequences.is_empty() {
            results.push(TranscriptConsequence {
                transcript_id: Arc::clone(&transcript.stable_id),
                gene_id: transcript.gene.stable_id.clone(),
                gene_symbol: transcript.gene.symbol.clone(),
                biotype: Arc::clone(&transcript.biotype),
                allele_consequences,
                canonical: transcript.canonical,
                strand: transcript.strand,
            });
        }
    }

    results
}

fn predict_sv_for_transcript(
    sv_start: u64,
    sv_end: u64,
    variant_type: VariantType,
    allele: &Allele,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> AlleleConsequenceResult {
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    // Check if SV overlaps the transcript at all
    let overlaps = sv_start <= tx_end && sv_end >= tx_start;

    if !overlaps {
        // Check upstream/downstream
        let (up_dist, down_dist) = if transcript.strand == Strand::Forward {
            (upstream_distance, downstream_distance)
        } else {
            (downstream_distance, upstream_distance)
        };

        let distance = if sv_end < tx_start {
            // SV is before transcript
            if transcript.strand == Strand::Forward {
                Some(tx_start as i64 - sv_end as i64) // upstream
            } else {
                Some(tx_start as i64 - sv_end as i64) // downstream on reverse
            }
        } else if sv_start > tx_end {
            if transcript.strand == Strand::Forward {
                Some(sv_start as i64 - tx_end as i64) // downstream
            } else {
                Some(sv_start as i64 - tx_end as i64) // upstream on reverse
            }
        } else {
            None
        };

        let consequence = if let Some(d) = distance {
            let abs_d = d.unsigned_abs();
            if sv_end < tx_start {
                if (transcript.strand == Strand::Forward && abs_d <= up_dist)
                    || (transcript.strand == Strand::Reverse && abs_d <= down_dist)
                {
                    if transcript.strand == Strand::Forward {
                        Consequence::UpstreamGeneVariant
                    } else {
                        Consequence::DownstreamGeneVariant
                    }
                } else {
                    Consequence::IntergenicVariant
                }
            } else {
                if (transcript.strand == Strand::Forward && abs_d <= down_dist)
                    || (transcript.strand == Strand::Reverse && abs_d <= up_dist)
                {
                    if transcript.strand == Strand::Forward {
                        Consequence::DownstreamGeneVariant
                    } else {
                        Consequence::UpstreamGeneVariant
                    }
                } else {
                    Consequence::IntergenicVariant
                }
            }
        } else {
            Consequence::IntergenicVariant
        };

        return AlleleConsequenceResult {
            allele: allele.clone(),
            consequences: vec![consequence],
            impact: consequence.impact(),
            cdna_start: None,
            cdna_end: None,
            cds_start: None,
            cds_end: None,
            protein_start: None,
            protein_end: None,
            amino_acids: None,
            codons: None,
            exon: None,
            intron: None,
            distance,
        };
    }

    // SV overlaps transcript — determine consequences
    let completely_contains = sv_start <= tx_start && sv_end >= tx_end;
    let consequences = determine_sv_overlap_consequences(
        sv_start,
        sv_end,
        variant_type,
        transcript,
        completely_contains,
    );

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min()
        .unwrap_or(Impact::Modifier);

    AlleleConsequenceResult {
        allele: allele.clone(),
        consequences,
        impact,
        cdna_start: None,
        cdna_end: None,
        cds_start: None,
        cds_end: None,
        protein_start: None,
        protein_end: None,
        amino_acids: None,
        codons: None,
        exon: None,
        intron: None,
        distance: None,
    }
}

fn determine_sv_overlap_consequences(
    sv_start: u64,
    sv_end: u64,
    variant_type: VariantType,
    transcript: &Transcript,
    completely_contains: bool,
) -> Vec<Consequence> {
    let mut consequences = Vec::new();

    match variant_type {
        VariantType::CopyNumberLoss | VariantType::Deletion => {
            if completely_contains {
                consequences.push(Consequence::TranscriptAblation);
            } else {
                // Partial overlap
                consequences.push(Consequence::FeatureTruncation);
                // Check if it hits coding region
                if transcript.is_coding() {
                    if hits_coding_region(sv_start, sv_end, transcript) {
                        consequences.push(Consequence::CodingSequenceVariant);
                    }
                    // Check splice sites
                    if hits_splice_site(sv_start, sv_end, transcript) {
                        consequences.push(Consequence::SpliceAcceptorVariant);
                    }
                }
                consequences.push(Consequence::CopyNumberDecrease);
            }
        }
        VariantType::TandemDuplication | VariantType::CopyNumberGain => {
            if completely_contains {
                consequences.push(Consequence::TranscriptAmplification);
            } else {
                consequences.push(Consequence::FeatureElongation);
                consequences.push(Consequence::CopyNumberIncrease);
            }
        }
        VariantType::CopyNumberVariation => {
            if completely_contains {
                // Could be either ablation or amplification
                consequences.push(Consequence::TranscriptVariant);
            }
            consequences.push(Consequence::CopyNumberChange);
        }
        VariantType::Inversion => {
            if completely_contains {
                consequences.push(Consequence::TranscriptVariant);
            } else {
                consequences.push(Consequence::FeatureTruncation);
                if transcript.is_coding() && hits_coding_region(sv_start, sv_end, transcript) {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
            }
        }
        VariantType::TranslocationBreakend => {
            // Breakend landing in a transcript
            consequences.push(Consequence::TranscriptVariant);
            if transcript.is_coding() && hits_coding_region(sv_start, sv_end, transcript) {
                consequences.push(Consequence::CodingSequenceVariant);
            }
        }
        VariantType::ShortTandemRepeatVariation => {
            consequences.push(Consequence::ShortTandemRepeatChange);
        }
        _ => {
            consequences.push(Consequence::TranscriptVariant);
        }
    }

    if consequences.is_empty() {
        consequences.push(Consequence::TranscriptVariant);
    }

    consequences
}

/// Check if the SV interval overlaps any exon's coding region.
fn hits_coding_region(sv_start: u64, sv_end: u64, transcript: &Transcript) -> bool {
    if let (Some(cds_start), Some(cds_end)) = (transcript.coding_region_start, transcript.coding_region_end) {
        let coding_start = cds_start.min(cds_end);
        let coding_end = cds_start.max(cds_end);
        sv_start <= coding_end && sv_end >= coding_start
    } else {
        false
    }
}

/// Check if the SV overlaps a splice site (2bp at exon boundaries).
fn hits_splice_site(sv_start: u64, sv_end: u64, transcript: &Transcript) -> bool {
    for exon in &transcript.exons {
        // Donor site: 2bp after exon end (exon.end+1, exon.end+2)
        if sv_start <= exon.end + 2 && sv_end >= exon.end + 1 {
            return true;
        }
        // Acceptor site: 2bp before exon start (exon.start-2, exon.start-1)
        if exon.start >= 3 && sv_start <= exon.start - 1 && sv_end >= exon.start - 2 {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxivep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript(start: u64, end: u64) -> Transcript {
        Transcript {
            stable_id: Arc::from("ENST_TEST"),
            version: None,
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start,
                end,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start,
            end,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start,
                    end: start + 200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: start + 500,
                    end: start + 700,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P_TEST".into(),
                genomic_start: start + 50,
                genomic_end: start + 650,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 2,
                end_exon_offset: 150,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(401),
            coding_region_start: Some(start + 50),
            coding_region_end: Some(start + 650),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
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
    fn test_deletion_completely_contains_transcript() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1", 4000, 7000, VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx], 5000, 5000,
        );

        assert_eq!(results.len(), 1);
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::TranscriptAblation));
    }

    #[test]
    fn test_deletion_partial_overlap() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1", 5300, 5800, VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx], 5000, 5000,
        );

        assert_eq!(results.len(), 1);
        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::FeatureTruncation));
        assert!(cons.contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_duplication_completely_contains() {
        let tx = make_coding_transcript(5000, 6000);
        let results = predict_sv_consequences(
            "chr1", 4000, 7000, VariantType::TandemDuplication,
            &[Allele::Symbolic("<DUP>".into())],
            &[&tx], 5000, 5000,
        );

        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_sv_no_overlap_upstream() {
        let tx = make_coding_transcript(10000, 11000);
        let results = predict_sv_consequences(
            "chr1", 5000, 6000, VariantType::CopyNumberLoss,
            &[Allele::Symbolic("<DEL>".into())],
            &[&tx], 5000, 5000,
        );

        let cons = &results[0].allele_consequences[0].consequences;
        assert!(cons.contains(&Consequence::UpstreamGeneVariant));
    }
}
