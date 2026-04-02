use oxivep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
use oxivep_genome::codon::format_codon_change;
use oxivep_genome::{CodonTable, Transcript};

use crate::splice;

/// Result of consequence prediction for a variant against a transcript.
#[derive(Debug, Clone)]
pub struct TranscriptConsequence {
    pub transcript_id: String,
    pub gene_id: String,
    pub gene_symbol: Option<String>,
    pub biotype: String,
    pub allele_consequences: Vec<AlleleConsequenceResult>,
    pub canonical: bool,
    pub strand: Strand,
}

/// Consequence result for a single allele against a transcript.
#[derive(Debug, Clone)]
pub struct AlleleConsequenceResult {
    pub allele: Allele,
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub cdna_start: Option<u64>,
    pub cdna_end: Option<u64>,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub protein_start: Option<u64>,
    pub protein_end: Option<u64>,
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    pub exon: Option<(u32, u32)>,
    pub intron: Option<(u32, u32)>,
    pub distance: Option<i64>,
}

/// Full prediction result for a variant.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub transcript_consequences: Vec<TranscriptConsequence>,
    pub most_severe: Option<Consequence>,
}

/// The consequence prediction engine.
pub struct ConsequencePredictor {
    pub upstream_distance: u64,
    pub downstream_distance: u64,
    codon_table: CodonTable,
}

impl ConsequencePredictor {
    pub fn new(upstream_distance: u64, downstream_distance: u64) -> Self {
        Self {
            upstream_distance,
            downstream_distance,
            codon_table: CodonTable::standard(),
        }
    }

    /// Predict consequences of a variant against a set of transcripts.
    pub fn predict(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcripts: &[&Transcript],
        ref_seq: Option<&[u8]>,
    ) -> PredictionResult {
        let mut transcript_consequences = Vec::new();

        for transcript in transcripts {
            let tc = self.predict_transcript(position, ref_allele, alt_alleles, transcript, ref_seq);
            transcript_consequences.push(tc);
        }

        let all_consequences: Vec<Consequence> = transcript_consequences
            .iter()
            .flat_map(|tc| {
                tc.allele_consequences
                    .iter()
                    .flat_map(|ac| ac.consequences.iter().copied())
            })
            .collect();

        let most_severe = Consequence::most_severe(&all_consequences);

        PredictionResult {
            transcript_consequences,
            most_severe,
        }
    }

    fn predict_transcript(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcript: &Transcript,
        ref_seq: Option<&[u8]>,
    ) -> TranscriptConsequence {
        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            .map(|alt| self.predict_allele(position, ref_allele, alt, transcript, ref_seq))
            .collect();

        TranscriptConsequence {
            transcript_id: transcript.stable_id.clone(),
            gene_id: transcript.gene.stable_id.clone(),
            gene_symbol: transcript.gene.symbol.clone(),
            biotype: transcript.biotype.clone(),
            allele_consequences,
            canonical: transcript.canonical,
            strand: transcript.strand,
        }
    }

    fn predict_allele(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        _ref_seq: Option<&[u8]>,
    ) -> AlleleConsequenceResult {
        let var_start = position.start;
        let var_end = position.end;
        let tr_start = transcript.start;
        let tr_end = transcript.end;

        let mut consequences = Vec::new();
        let mut cds_start = None;
        let mut cds_end = None;
        let mut protein_start = None;
        let mut protein_end = None;
        let mut amino_acids = None;
        let mut codons = None;
        let mut distance = None;

        // 1. Check if variant overlaps the transcript at all
        let overlaps = var_start <= tr_end && var_end >= tr_start;

        if !overlaps {
            // Check upstream/downstream
            let dist = self.distance_to_transcript(var_start, var_end, transcript);
            if let Some(d) = dist {
                distance = Some(d);
                let abs_dist = d.unsigned_abs();
                if abs_dist <= self.upstream_distance {
                    match transcript.strand {
                        Strand::Forward => {
                            if var_end < tr_start {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                        Strand::Reverse => {
                            if var_start > tr_end {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                    }
                }
            }

            if consequences.is_empty() {
                consequences.push(Consequence::IntergenicVariant);
            }

            let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);
            return AlleleConsequenceResult {
                allele: alt_allele.clone(),
                consequences,
                impact,
                cdna_start: None, cdna_end: None,
                cds_start: None, cds_end: None,
                protein_start: None, protein_end: None,
                amino_acids: None, codons: None,
                exon: None, intron: None,
                distance,
            };
        }

        // 2. Map to cDNA coordinates
        let cdna_start = transcript.genomic_to_cdna(var_start);
        let cdna_end = transcript.genomic_to_cdna(var_end);

        // 3. Determine exon/intron location
        let exon_info = transcript.exon_at(var_start).map(|(i, t)| (i as u32 + 1, t as u32));
        let intron_info = transcript.intron_at(var_start).map(|(i, t)| (i as u32 + 1, t as u32));

        let in_exon = exon_info.is_some();
        let in_intron = intron_info.is_some();

        // 4. Check splice sites (always check regardless of coding status)
        if splice::is_splice_donor(transcript, var_start) {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if splice::is_splice_acceptor(transcript, var_start) {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }

        // Only add extended splice consequences if not already a donor/acceptor
        let is_essential_splice = consequences.iter().any(|c| {
            matches!(c, Consequence::SpliceDonorVariant | Consequence::SpliceAcceptorVariant)
        });

        if !is_essential_splice {
            let is_donor_5th = splice::is_splice_donor_5th_base(transcript, var_start);
            let is_donor_region = splice::is_splice_donor_region(transcript, var_start);
            if is_donor_5th {
                consequences.push(Consequence::SpliceDonorFifthBaseVariant);
            }
            if is_donor_region {
                consequences.push(Consequence::SpliceDonorRegionVariant);
            }
            if splice::is_splice_polypyrimidine_tract(transcript, var_start) {
                consequences.push(Consequence::SplicePolypyrimidineTractVariant);
            }
            // VEP excludes splice_region_variant when a more specific splice term is present:
            // splice_donor_region_variant or splice_donor_5th_base_variant
            if !is_donor_5th && !is_donor_region {
                if splice::is_splice_region(transcript, var_start) {
                    consequences.push(Consequence::SpliceRegionVariant);
                }
            }
        }

        // 5. Coding vs non-coding transcript
        if transcript.is_coding() {
            let coding_start = transcript.coding_region_start.unwrap_or(0);
            let coding_end = transcript.coding_region_end.unwrap_or(0);

            // Map to CDS coordinates
            if let Some(cs) = cdna_start {
                cds_start = transcript.cdna_to_cds(cs);
            }
            if let Some(ce) = cdna_end {
                cds_end = transcript.cdna_to_cds(ce);
            }
            if let Some(cds_s) = cds_start {
                protein_start = Some(Transcript::cds_to_protein(cds_s));
            }
            if let Some(cds_e) = cds_end {
                protein_end = Some(Transcript::cds_to_protein(cds_e));
            }

            let in_coding_region = self.is_in_coding_region(var_start, var_end, coding_start, coding_end);
            let in_5_utr = self.is_in_5_utr(var_start, var_end, transcript);
            let in_3_utr = self.is_in_3_utr(var_start, var_end, transcript);

            if in_5_utr && in_exon {
                consequences.push(Consequence::FivePrimeUtrVariant);
            } else if in_3_utr && in_exon {
                consequences.push(Consequence::ThreePrimeUtrVariant);
            } else if in_coding_region && in_exon {
                // Coding exonic variant - determine coding consequence
                let coding_conseq = self.predict_coding_consequence(
                    ref_allele, alt_allele, transcript, cds_start, cds_end,
                );
                if let Some((conseq, aa, cdn)) = coding_conseq {
                    consequences.push(conseq);
                    amino_acids = aa;
                    codons = cdn;
                } else {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
            } else if in_intron && !is_essential_splice {
                // VEP excludes intron_variant for positions at splice donor/acceptor sites
                consequences.push(Consequence::IntronVariant);
            }
        } else {
            // Non-coding transcript
            if in_exon {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else if in_intron {
                if !is_essential_splice {
                    consequences.push(Consequence::IntronVariant);
                }
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }

        // If still no consequences, add catch-all
        if consequences.is_empty() {
            if transcript.is_coding() {
                consequences.push(Consequence::CodingTranscriptVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }

        // Add NMD_transcript_variant modifier for nonsense_mediated_decay transcripts
        if transcript.biotype == "nonsense_mediated_decay" {
            consequences.push(Consequence::NmdTranscriptVariant);
        }

        // Deduplicate
        consequences.sort_by_key(|c| c.rank());
        consequences.dedup();

        let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);

        AlleleConsequenceResult {
            allele: alt_allele.clone(),
            consequences,
            impact,
            cdna_start, cdna_end,
            cds_start, cds_end,
            protein_start, protein_end,
            amino_acids, codons,
            exon: exon_info, intron: intron_info,
            distance,
        }
    }

    /// Predict the coding consequence (missense, synonymous, frameshift, etc.)
    fn predict_coding_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        _cds_end: Option<u64>,
    ) -> Option<(Consequence, Option<(String, String)>, Option<(String, String)>)> {
        let cds_pos = cds_start?;

        let ref_len = ref_allele.len();
        let alt_len = alt_allele.len();

        // Check if this is a frameshift or in-frame indel
        let is_deletion = *ref_allele != Allele::Deletion && *alt_allele == Allele::Deletion;
        let is_insertion = *ref_allele == Allele::Deletion && *alt_allele != Allele::Missing;
        let is_indel = is_deletion || is_insertion || ref_len != alt_len;

        if is_indel {
            let (consequence, is_frameshift) = if is_deletion || is_insertion {
                let indel_len = if is_deletion { ref_len } else { alt_len };
                if indel_len % 3 != 0 {
                    (Consequence::FrameshiftVariant, true)
                } else if is_insertion {
                    (Consequence::InframeInsertion, false)
                } else {
                    (Consequence::InframeDeletion, false)
                }
            } else {
                let len_diff = (ref_len as i64 - alt_len as i64).unsigned_abs() as usize;
                if len_diff % 3 != 0 {
                    (Consequence::FrameshiftVariant, true)
                } else if ref_len > alt_len {
                    (Consequence::InframeDeletion, false)
                } else {
                    (Consequence::InframeInsertion, false)
                }
            };

            // Try to compute amino acids and codons from translateable_seq
            let (aa_pair, codon_pair) = self.compute_indel_amino_acids(
                transcript, cds_pos, ref_allele, alt_allele, is_frameshift,
            );

            return Some((consequence, aa_pair, codon_pair));
        }

        // Same length substitution (SNV or MNV)
        if let Some(ref translateable_seq) = transcript.translateable_seq {
            let seq_bytes = translateable_seq.as_bytes();

            // Get the codon containing this CDS position
            let codon_number = ((cds_pos - 1) / 3) as usize;
            let codon_offset = ((cds_pos - 1) % 3) as usize;
            let codon_start = codon_number * 3;

            if codon_start + 3 <= seq_bytes.len() {
                let ref_codon = [seq_bytes[codon_start], seq_bytes[codon_start + 1], seq_bytes[codon_start + 2]];
                let mut alt_codon = ref_codon;

                // Apply the substitution
                if let Allele::Sequence(alt_bases) = alt_allele {
                    for (i, &base) in alt_bases.iter().enumerate() {
                        let pos = codon_offset + i;
                        if pos < 3 {
                            alt_codon[pos] = match transcript.strand {
                                Strand::Forward => base,
                                Strand::Reverse => complement(base),
                            };
                        }
                    }
                }

                let ref_aa = self.codon_table.translate(&ref_codon);
                let alt_aa = self.codon_table.translate(&alt_codon);

                let (ref_codon_str, alt_codon_str) = format_codon_change(&ref_codon, &alt_codon);

                let ref_aa_str = String::from(ref_aa as char);
                let alt_aa_str = String::from(alt_aa as char);

                let codon_pair = Some((ref_codon_str, alt_codon_str));
                let aa_pair = Some((ref_aa_str.clone(), alt_aa_str.clone()));

                // Check if this affects the start codon
                if codon_number == 0 && !CodonTable::is_start(&alt_codon) {
                    return Some((Consequence::StartLost, aa_pair, codon_pair));
                }

                // Determine consequence type
                if ref_aa == alt_aa {
                    if ref_aa == b'*' {
                        return Some((Consequence::StopRetainedVariant, aa_pair, codon_pair));
                    }
                    return Some((Consequence::SynonymousVariant, aa_pair, codon_pair));
                }

                if alt_aa == b'*' {
                    return Some((Consequence::StopGained, aa_pair, codon_pair));
                }
                if ref_aa == b'*' {
                    return Some((Consequence::StopLost, aa_pair, codon_pair));
                }

                return Some((Consequence::MissenseVariant, aa_pair, codon_pair));
            }
        }

        // Fallback: if we can't determine the exact consequence,
        // classify based on whether it's an in-frame or frameshift change
        if ref_len == alt_len {
            Some((Consequence::MissenseVariant, None, None))
        } else {
            None
        }
    }

    /// Compute amino acids and codons affected by an indel variant.
    /// Returns (amino_acids, codons) tuples.
    /// For frameshifts: ref codon with VEP-style case formatting, truncated alt codon.
    fn compute_indel_amino_acids(
        &self,
        transcript: &Transcript,
        cds_pos: u64,
        ref_allele: &Allele,
        alt_allele: &Allele,
        is_frameshift: bool,
    ) -> (Option<(String, String)>, Option<(String, String)>) {
        let translateable_seq = match transcript.translateable_seq.as_ref() {
            Some(s) => s,
            None => return (None, None),
        };
        let seq_bytes = translateable_seq.as_bytes();
        let cds_idx = (cds_pos - 1) as usize;

        if cds_idx >= seq_bytes.len() {
            return (None, None);
        }

        // Get the codon at the affected position
        let codon_number = cds_idx / 3;
        let codon_offset = cds_idx % 3;
        let codon_start = codon_number * 3;

        if codon_start + 3 > seq_bytes.len() {
            return (None, None);
        }

        let ref_codon = [seq_bytes[codon_start], seq_bytes[codon_start + 1], seq_bytes[codon_start + 2]];
        let ref_aa = self.codon_table.translate(&ref_codon);
        let ref_aa_str = String::from(ref_aa as char);

        if is_frameshift {
            // Build the alt sequence by applying the indel
            let mut alt_seq: Vec<u8> = seq_bytes.to_vec();

            match (ref_allele, alt_allele) {
                (Allele::Sequence(_), Allele::Deletion) => {
                    let del_len = ref_allele.len();
                    let end = (cds_idx + del_len).min(alt_seq.len());
                    alt_seq.drain(cds_idx..end);
                }
                (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                    let mut bases: Vec<u8> = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().map(|&b| complement(b)).collect();
                    }
                    for (i, &b) in bases.iter().enumerate() {
                        alt_seq.insert(cds_idx + i, b);
                    }
                }
                (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) => {
                    let end = (cds_idx + ref_bases.len()).min(alt_seq.len());
                    let mut replacement = alt_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        replacement = replacement.iter().map(|&b| complement(b)).collect();
                    }
                    alt_seq.splice(cds_idx..end, replacement);
                }
                _ => return (None, None),
            }

            // Build codon display: VEP style with deleted base uppercase
            // ref codon: lowercase bases, uppercase at the deleted position(s)
            let mut ref_codon_display = String::with_capacity(3);
            for i in 0..3 {
                if i == codon_offset {
                    ref_codon_display.push((ref_codon[i] as char).to_ascii_uppercase());
                } else {
                    ref_codon_display.push((ref_codon[i] as char).to_ascii_lowercase());
                }
            }

            // alt codon: show only the remaining bases of the original codon after the indel
            // For a deletion at offset 2 in a 3-base codon: show only the 2 remaining bases
            let alt_codon_display: String = {
                let mut original_codon: Vec<u8> = ref_codon.to_vec();
                match (ref_allele, alt_allele) {
                    (Allele::Sequence(_), Allele::Deletion) => {
                        // Remove the deleted base(s) from the codon
                        let del_len = ref_allele.len().min(3 - codon_offset);
                        let end = (codon_offset + del_len).min(original_codon.len());
                        original_codon.drain(codon_offset..end);
                    }
                    (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                        // Insert bases into the codon at the offset
                        let mut bases = ins_bases.clone();
                        if transcript.strand == Strand::Reverse {
                            bases = bases.iter().map(|&b| complement(b)).collect();
                        }
                        for (j, &b) in bases.iter().enumerate() {
                            original_codon.insert(codon_offset + j, b);
                        }
                    }
                    _ => {}
                }
                original_codon.iter()
                    .map(|&b| (b as char).to_ascii_lowercase())
                    .collect()
            };

            // For frameshifts, alt amino acid is always X (unknown/frameshift)
            let aa_pair = Some((ref_aa_str, "X".to_string()));
            let codon_pair = Some((ref_codon_display, alt_codon_display));
            (aa_pair, codon_pair)
        } else {
            // In-frame indel: show affected amino acids
            let aa_pair = match (ref_allele, alt_allele) {
                (Allele::Sequence(_), Allele::Deletion) => {
                    let del_len = ref_allele.len();
                    let end_cds = (cds_idx + del_len).min(seq_bytes.len());
                    let deleted_region = &seq_bytes[codon_start..((end_cds + 2) / 3 * 3).min(seq_bytes.len())];
                    let deleted_aas: String = deleted_region.chunks(3)
                        .filter(|c| c.len() == 3)
                        .map(|c| self.codon_table.translate(&[c[0], c[1], c[2]]) as char)
                        .collect();
                    Some((deleted_aas, "-".to_string()))
                }
                (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                    let mut bases = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().map(|&b| complement(b)).collect();
                    }
                    let inserted_aas: String = bases.chunks(3)
                        .filter(|c| c.len() == 3)
                        .map(|c| self.codon_table.translate(&[c[0], c[1], c[2]]) as char)
                        .collect();
                    Some(("-".to_string(), inserted_aas))
                }
                _ => Some((ref_aa_str, "X".to_string())),
            };
            (aa_pair, None)
        }
    }

    fn distance_to_transcript(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> Option<i64> {
        // For insertions (end < start), use start for distance calculation
        // since start represents the actual insertion position
        let effective_start = var_start.min(var_end);
        let effective_end = var_start.max(var_end);
        if effective_end < transcript.start {
            Some((transcript.start - effective_end) as i64)
        } else if effective_start > transcript.end {
            Some((effective_start - transcript.end) as i64)
        } else {
            Some(0)
        }
    }

    fn is_in_coding_region(&self, var_start: u64, var_end: u64, coding_start: u64, coding_end: u64) -> bool {
        var_start <= coding_end && var_end >= coding_start
    }

    fn is_in_5_utr(&self, var_start: u64, _var_end: u64, transcript: &Transcript) -> bool {
        let coding_start = match transcript.coding_region_start {
            Some(s) => s,
            None => return false,
        };
        let coding_end = match transcript.coding_region_end {
            Some(e) => e,
            None => return false,
        };

        match transcript.strand {
            Strand::Forward => var_start < coding_start && var_start >= transcript.start,
            Strand::Reverse => var_start > coding_end && var_start <= transcript.end,
        }
    }

    fn is_in_3_utr(&self, var_start: u64, _var_end: u64, transcript: &Transcript) -> bool {
        let coding_start = match transcript.coding_region_start {
            Some(s) => s,
            None => return false,
        };
        let coding_end = match transcript.coding_region_end {
            Some(e) => e,
            None => return false,
        };

        match transcript.strand {
            Strand::Forward => var_start > coding_end && var_start <= transcript.end,
            Strand::Reverse => var_start < coding_start && var_start >= transcript.start,
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

impl Default for ConsequencePredictor {
    fn default() -> Self {
        Self::new(5000, 5000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxivep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript() -> Transcript {
        // A simple protein-coding transcript on forward strand:
        // Exon1: 1000-1200 (UTR: 1000-1049, CDS: 1050-1200)
        // Intron: 1201-1999
        // Exon2: 2000-2300 (all CDS)
        // Intron: 2301-3999
        // Exon3: 4000-5000 (CDS: 4000-4500, UTR: 4501-5000)
        //
        // CDS length: 151 + 301 + 501 = 953 bases
        // translateable_seq: from cDNA pos 51 to 953+50=1003
        let translateable = "ATGGCTTCAAAGCCC".to_string() + &"A".repeat(938); // starts with ATG

        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon { stable_id: "E1".into(), start: 1000, end: 1200, strand: Strand::Forward, phase: -1, end_phase: 0, rank: 1 },
                Exon { stable_id: "E2".into(), start: 2000, end: 2300, strand: Strand::Forward, phase: 0, end_phase: 1, rank: 2 },
                Exon { stable_id: "E3".into(), start: 4000, end: 5000, strand: Strand::Forward, phase: 1, end_phase: -1, rank: 3 },
            ],
            translation: Some(Translation {
                stable_id: "ENSP00000001".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(1003),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: Some(translateable),
            peptide: None,
            canonical: true,
            mane_select: None, mane_plus_clinical: None,
            tsl: Some(1), appris: None, ccds: None,
            protein_id: Some("ENSP00000001".into()),
            protein_version: None,
            swissprot: vec![], trembl: vec![], uniparc: vec![],
            refseq_id: None, source: None, gencode_primary: false, flags: vec![], codon_table_start_phase: 0,
        }
    }

    fn make_noncoding_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST_NC".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_NC".into(),
                symbol: Some("NCRNA1".into()),
                symbol_source: None, hgnc_id: None,
                biotype: "lncRNA".into(),
                chromosome: "chr1".into(),
                start: 10000, end: 12000,
                strand: Strand::Forward,
            },
            biotype: "lncRNA".into(),
            chromosome: "chr1".into(),
            start: 10000, end: 12000,
            strand: Strand::Forward,
            exons: vec![
                Exon { stable_id: "E1".into(), start: 10000, end: 10500, strand: Strand::Forward, phase: -1, end_phase: -1, rank: 1 },
                Exon { stable_id: "E2".into(), start: 11500, end: 12000, strand: Strand::Forward, phase: -1, end_phase: -1, rank: 2 },
            ],
            translation: None,
            cdna_coding_start: None, cdna_coding_end: None,
            coding_region_start: None, coding_region_end: None,
            spliced_seq: None, translateable_seq: None, peptide: None,
            canonical: false, mane_select: None, mane_plus_clinical: None,
            tsl: None, appris: None, ccds: None, protein_id: None, protein_version: None,
            swissprot: vec![], trembl: vec![], uniparc: vec![],
            refseq_id: None, source: None, gencode_primary: false, flags: vec![], codon_table_start_phase: 0,
        }
    }

    #[test]
    fn test_upstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 500, 500, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        assert_eq!(result.transcript_consequences.len(), 1);
        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 1);
        assert!(tc.allele_consequences[0].consequences.contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.allele_consequences[0].distance, Some(500));
    }

    #[test]
    fn test_downstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 5500, 5500, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(ac.distance, Some(500));
    }

    #[test]
    fn test_intergenic_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Very far away
        let pos = GenomicPosition::new("chr1", 100000, 100000, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_5_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1020 is in exon1 (1000-1200), before CDS start (1050)
        let pos = GenomicPosition::new("chr1", 1020, 1020, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_3_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 4600 is in exon3 (4000-5000), after CDS end (4500)
        let pos = GenomicPosition::new("chr1", 4600, 4600, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn test_intron_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1500 is in intron1 (1201-1999), away from splice sites
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_splice_donor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1201 is first base of intron1 → splice donor
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("G"), &[Allele::from_str("A")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SpliceDonorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_splice_acceptor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1999 is last base of intron1 → splice acceptor
        let pos = GenomicPosition::new("chr1", 1999, 1999, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("G"), &[Allele::from_str("A")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SpliceAcceptorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_synonymous_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // First CDS position is at genomic 1050, which is cDNA pos 51, CDS pos 1
        // translateable_seq starts with "ATG" (Met)
        // CDS pos 3 (third base of first codon) - change G to A: ATA still codes for... wait
        // ATG -> Met. Let's change position 3 of ATG from G to something that's still Met: not possible
        // Let's use a different codon. CDS pos 4-6 is "GCT" (Ala). GCC also codes for Ala.
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        // Change T at CDS pos 6 to C: GCT -> GCC both = Ala
        let pos = GenomicPosition::new("chr1", 1055, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("C")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::SynonymousVariant),
            "Expected synonymous, got: {:?}", ac.consequences
        );
        assert_eq!(ac.impact, Impact::Low);
    }

    #[test]
    fn test_missense_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT" (Ala). Change first base G to T: TCT = Ser (different!)
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::MissenseVariant),
            "Expected missense, got: {:?}", ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_stop_gained() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT". Change to "TAA" (stop) → need to change pos 4,5,6
        // For simplicity, change CDS pos 4: G->T, pos 5: C->A, pos 6: T->A
        // But our predictor works one SNV at a time. Let's pick a codon that's one base
        // away from a stop. "TCA" (Ser) → change C→A: TAA (stop). But we'd need that codon.
        // Actually, translateable_seq[6..9] = "TCA" (positions 7-9 in 1-based)
        // CDS pos 7 is at genomic 1050+6 = 1056
        // Change T to T (no), we need C at pos 8 to become something.
        // Let's just use translateable[3..6] = "GCT" and change pos 4 (G) to T: "TCT" = Ser
        // That's missense, not stop. Let's try another approach.
        // translateable[9..12] = "AAG" (Lys). Change A at pos 10 to T: TAG = stop!
        // CDS pos 10 is at genomic 1050+9 = 1059
        let pos = GenomicPosition::new("chr1", 1059, 1059, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StopGained),
            "Expected stop_gained, got: {:?}. translateable[9..12]={:?}",
            ac.consequences,
            &tr.translateable_seq.as_ref().unwrap()[9..12]
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_frameshift_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 1 base in CDS → frameshift
        // CDS pos 4 at genomic 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::FrameshiftVariant),
            "Expected frameshift, got: {:?}", ac.consequences
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_inframe_deletion() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 3 bases in CDS → inframe deletion
        let pos = GenomicPosition::new("chr1", 1053, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("GCT"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::InframeDeletion),
            "Expected inframe_deletion, got: {:?}", ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_noncoding_exon_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_noncoding_transcript();
        let pos = GenomicPosition::new("chr1", 10100, 10100, Strand::Forward);
        let result = predictor.predict(&pos, &Allele::from_str("A"), &[Allele::from_str("G")], &[&tr], None);

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn test_start_lost() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 1 (first base of ATG) is at genomic 1050
        // Change A to G: GTG is not a standard start codon
        let pos = GenomicPosition::new("chr1", 1050, 1050, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StartLost),
            "Expected start_lost, got: {:?}", ac.consequences
        );
    }

    #[test]
    fn test_multiple_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr1 = make_coding_transcript();
        let tr2 = make_noncoding_transcript();

        // Position in tr1's intron, not overlapping tr2
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr1, &tr2],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 2);
        // tr1: intron variant
        assert!(result.transcript_consequences[0].allele_consequences[0].consequences.contains(&Consequence::IntronVariant));
        // tr2: 8500bp away (>5000), so intergenic
        assert!(result.transcript_consequences[1].allele_consequences[0].consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_most_severe_across_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Splice donor is more severe than intron variant
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        assert_eq!(result.most_severe, Some(Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_multi_allelic() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Two alt alleles at a coding position
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T"), Allele::from_str("C")],
            &[&tr],
            None,
        );

        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 2);
    }
}
