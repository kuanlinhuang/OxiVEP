use oxivep_core::Allele;

/// Generate HGVSc (coding DNA) notation.
///
/// Uses the transcript ID as the reference and cDNA position numbering.
/// Format: ENST00000001.1:c.123A>G
pub fn hgvsc(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    hgvsc_with_seq(transcript_id, cdna_start, cdna_end, ref_allele, alt_allele, coding_start, coding_end, None, 0)
}

/// Generate HGVSc with optional 3' shifting for deletions/insertions.
///
/// When `spliced_seq` is provided, deletions in repetitive regions are shifted
/// to the most 3' position per HGVS nomenclature standard.
pub fn hgvsc_with_seq(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
    spliced_seq: Option<&str>,
    start_phase: u64,
) -> Option<String> {
    let prefix = format!("{}:c.", transcript_id);

    // Convert cDNA position to CDS position (relative to ATG)
    // For 5'UTR (before ATG): position = cdna - coding_start (negative, no +1)
    // For coding region: position = cdna - coding_start + 1 + start_phase (1-based)
    // start_phase accounts for incomplete start codons (Ensembl convention)
    let cds_pos_start = if (cdna_start as i64) < (coding_start as i64) {
        cdna_start as i64 - coding_start as i64  // 5'UTR: no phase
    } else {
        cdna_start as i64 - coding_start as i64 + 1 + start_phase as i64  // CDS: with phase
    };
    let cds_pos_end = if (cdna_end as i64) < (coding_start as i64) {
        cdna_end as i64 - coding_start as i64
    } else {
        cdna_end as i64 - coding_start as i64 + 1 + start_phase as i64
    };

    let pos_str = if cds_pos_start < 0 {
        // 5' UTR: negative positions
        format!("{}", cds_pos_start)
    } else if coding_end.is_some() && cdna_start > coding_end.unwrap() {
        // 3' UTR: c.*N notation
        let utr_offset_start = cdna_start - coding_end.unwrap();
        if cdna_start == cdna_end {
            format!("*{}", utr_offset_start)
        } else {
            let utr_offset_end = cdna_end - coding_end.unwrap();
            format!("*{}_*{}", utr_offset_start, utr_offset_end)
        }
    } else if cds_pos_start == cds_pos_end {
        format!("{}", cds_pos_start)
    } else {
        format!("{}_{}", cds_pos_start, cds_pos_end)
    };

    let notation = match (ref_allele, alt_allele) {
        // SNV
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix,
                pos_str,
                ref_bases[0] as char,
                alt_bases[0] as char
            )
        }
        // Deletion — apply HGVS 3' shifting in repetitive regions
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            let del_len = ref_bases.len();
            // Try 3' shifting if we have the transcript sequence
            let (shifted_start, shifted_end) = if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let mut s = (cdna_start - 1) as usize; // 0-based start
                let mut e = (cdna_end - 1) as usize;   // 0-based end (inclusive)
                // Shift right while the base at position end+1 matches base at start
                while e + 1 < seq_bytes.len()
                    && seq_bytes[e + 1].to_ascii_uppercase()
                        == seq_bytes[s].to_ascii_uppercase()
                {
                    s += 1;
                    e += 1;
                }
                (s as u64 + 1, e as u64 + 1) // back to 1-based
            } else {
                (cdna_start, cdna_end)
            };
            // Recompute positions for shifted coordinates, with UTR awareness
            let shifted_pos = if coding_end.is_some() && shifted_start > coding_end.unwrap() {
                // 3'UTR deletion
                let utr_s = shifted_start - coding_end.unwrap();
                if del_len == 1 {
                    format!("*{}", utr_s)
                } else {
                    let utr_e = shifted_end - coding_end.unwrap();
                    format!("*{}_*{}", utr_s, utr_e)
                }
            } else if (shifted_start as i64) < (coding_start as i64) {
                // 5'UTR deletion
                let cds_s = shifted_start as i64 - coding_start as i64;
                if del_len == 1 {
                    format!("{}", cds_s)
                } else {
                    let cds_e = shifted_end as i64 - coding_start as i64;
                    if coding_end.is_some() && shifted_end > coding_end.unwrap() {
                        // Spans into 3'UTR
                        let utr_e = shifted_end - coding_end.unwrap();
                        format!("{}_*{}", cds_s, utr_e)
                    } else {
                        format!("{}_{}", cds_s, cds_e)
                    }
                }
            } else {
                let cds_s = shifted_start as i64 - coding_start as i64 + 1;
                if del_len == 1 {
                    format!("{}", cds_s)
                } else {
                    if coding_end.is_some() && shifted_end > coding_end.unwrap() {
                        // Spans from CDS into 3'UTR
                        let utr_e = shifted_end - coding_end.unwrap();
                        format!("{}_*{}", cds_s, utr_e)
                    } else {
                        let cds_e = shifted_end as i64 - coding_start as i64 + 1;
                        format!("{}_{}", cds_s, cds_e)
                    }
                }
            };
            format!("{}{}del", prefix, shifted_pos)
        }
        // Insertion — normalize coordinates: ensure ins_before < ins_after
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_before_cdna = cdna_start.min(cdna_end); // base before insertion
            let ins_after_cdna = cdna_start.max(cdna_end);   // base after insertion
            let ins_before_cds = cds_pos_start.min(cds_pos_end);
            let ins_after_cds = cds_pos_start.max(cds_pos_end);

            // Check for duplication: if inserted bases match the preceding OR following sequence
            let is_dup = if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let ins_len = alt_bases.len();
                let before_pos = (ins_before_cdna - 1) as usize;
                // Check preceding sequence
                let dup_before = if before_pos + 1 >= ins_len && before_pos < seq_bytes.len() {
                    let preceding = &seq_bytes[before_pos + 1 - ins_len..before_pos + 1];
                    preceding.iter().zip(alt_bases.iter())
                        .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                } else {
                    false
                };
                // Check following sequence
                let dup_after = if !dup_before {
                    let after_pos = ins_after_cdna as usize; // 0-based index of base after insertion
                    if after_pos + ins_len <= seq_bytes.len() {
                        let following = &seq_bytes[after_pos..after_pos + ins_len];
                        following.iter().zip(alt_bases.iter())
                            .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                    } else {
                        false
                    }
                } else {
                    false
                };
                dup_before || dup_after
            } else {
                false
            };

            // Position of the base before insertion in CDS or UTR coordinates
            let ins_pos_str = if coding_end.is_some() && ins_before_cdna > coding_end.unwrap() {
                // 3'UTR: both flanking positions are in 3'UTR
                let utr_before = ins_before_cdna - coding_end.unwrap();
                let utr_after = ins_after_cdna - coding_end.unwrap();
                if is_dup {
                    let ins_len = alt_bases.len() as u64;
                    if ins_len == 1 {
                        format!("*{}", utr_before)
                    } else {
                        format!("*{}_*{}", utr_before - ins_len + 1, utr_before)
                    }
                } else {
                    format!("*{}_*{}", utr_before, utr_after)
                }
            } else if ins_before_cds < 0 {
                // 5'UTR insertion
                if is_dup {
                    let ins_len = alt_bases.len() as i64;
                    if ins_len == 1 {
                        format!("{}", ins_before_cds)
                    } else {
                        format!("{}_{}", ins_before_cds - ins_len + 1, ins_before_cds)
                    }
                } else {
                    format!("{}_{}", ins_before_cds, ins_after_cds)
                }
            } else if is_dup {
                let ins_len = alt_bases.len() as i64;
                if ins_len == 1 {
                    format!("{}", ins_before_cds)
                } else {
                    format!("{}_{}", ins_before_cds - ins_len + 1, ins_before_cds)
                }
            } else {
                format!("{}_{}", ins_before_cds, ins_after_cds)
            };

            if is_dup {
                format!("{}{}dup", prefix, ins_pos_str)
            } else {
                format!(
                    "{}{}ins{}",
                    prefix,
                    ins_pos_str,
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            }
        }
        // MNV or complex
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            format!(
                "{}{}delins{}",
                prefix,
                pos_str,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        _ => return None,
    };

    Some(notation)
}

/// Generate HGVSc notation for an intronic variant.
///
/// Uses offset notation from the nearest exon boundary:
/// - Positive offset: c.151+5A>G (donor side)
/// - Negative offset: c.152-3A>G (acceptor side)
pub fn hgvsc_intronic(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    hgvsc_intronic_range(
        transcript_id, nearest_exon_cdna_pos, intron_offset,
        None, // end_cdna_pos
        None, // end_offset
        ref_allele, alt_allele, coding_start, coding_end,
    )
}

/// Generate HGVSc intronic notation with optional end position for multi-base variants.
pub fn hgvsc_intronic_range(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    end_cdna_pos: Option<u64>,
    end_intron_offset: Option<i64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix = format!("{}:c.", transcript_id);

    // Convert cDNA pos to CDS-relative position.
    // For CDS positions: cds_pos = cdna - coding_start + 1 (so position 1 = first CDS base)
    // For 5'UTR positions: cds_pos = cdna - coding_start (no +1, since there's no position 0;
    //   position -1 = last 5'UTR base at cdna = coding_start - 1)
    let raw_cds_pos = nearest_exon_cdna_pos as i64 - coding_start as i64 + 1;
    let cds_pos = if raw_cds_pos <= 0 { raw_cds_pos - 1 } else { raw_cds_pos };

    // Build the position string with offset
    let pos_str = if cds_pos < 0 {
        // 5' UTR
        if intron_offset > 0 {
            format!("{}+{}", cds_pos, intron_offset)
        } else {
            format!("{}{}", cds_pos, intron_offset) // offset is already negative
        }
    } else if coding_end.is_some() && nearest_exon_cdna_pos > coding_end.unwrap() {
        // 3' UTR
        let utr_offset = nearest_exon_cdna_pos - coding_end.unwrap();
        if intron_offset > 0 {
            format!("*{}+{}", utr_offset, intron_offset)
        } else {
            format!("*{}{}", utr_offset, intron_offset)
        }
    } else if intron_offset > 0 {
        format!("{}+{}", cds_pos, intron_offset)
    } else {
        format!("{}{}", cds_pos, intron_offset) // offset is already negative
    };

    // Format the variant
    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!("{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char)
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            if ref_bases.len() == 1 || end_intron_offset.is_none() {
                format!("{}{}del", prefix, pos_str)
            } else {
                // Multi-base intronic deletion: use end position from caller
                let e_cdna = end_cdna_pos.unwrap_or(nearest_exon_cdna_pos);
                let e_offset = end_intron_offset.unwrap();
                let e_raw = e_cdna as i64 - coding_start as i64 + 1;
                let e_cds = if e_raw <= 0 { e_raw - 1 } else { e_raw };
                let end_pos_str = if e_cds < 0 {
                    if e_offset > 0 { format!("{}+{}", e_cds, e_offset) }
                    else { format!("{}{}", e_cds, e_offset) }
                } else if coding_end.is_some() && e_cdna > coding_end.unwrap() {
                    let utr = e_cdna - coding_end.unwrap();
                    if e_offset > 0 { format!("*{}+{}", utr, e_offset) }
                    else { format!("*{}{}", utr, e_offset) }
                } else if e_offset > 0 {
                    format!("{}+{}", e_cds, e_offset)
                } else {
                    format!("{}{}", e_cds, e_offset)
                };

                // For same-sign offsets: smaller absolute value first (closer to exon)
                // For positive: +4035 before +4038
                // For negative: -2625 before -2616
                let (start_str, end_str) = if intron_offset > 0 && e_offset > 0 {
                    if intron_offset < e_offset { (pos_str.clone(), end_pos_str) }
                    else { (end_pos_str, pos_str.clone()) }
                } else if intron_offset < 0 && e_offset < 0 {
                    if intron_offset < e_offset { (pos_str.clone(), end_pos_str) }
                    else { (end_pos_str, pos_str.clone()) }
                } else {
                    (pos_str.clone(), end_pos_str)
                };
                format!("{}{}_{}del", prefix, start_str, end_str)
            }
        }
        // Insertion in intron — show range: c.X+N_X+N+1insB
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_str = std::str::from_utf8(alt_bases).unwrap_or("?");
            let next_offset = if intron_offset < 0 { intron_offset + 1 } else { intron_offset + 1 };
            let build_pos = |cdna: u64, off: i64| -> String {
                let raw = cdna as i64 - coding_start as i64 + 1;
                let cp = if raw <= 0 { raw - 1 } else { raw }; // skip position 0 for 5'UTR
                if cp < 0 {
                    if off > 0 { format!("{}+{}", cp, off) } else { format!("{}{}", cp, off) }
                } else if coding_end.is_some() && cdna > coding_end.unwrap() {
                    let u = cdna - coding_end.unwrap();
                    if off > 0 { format!("*{}+{}", u, off) } else { format!("*{}{}", u, off) }
                } else if off > 0 {
                    format!("{}+{}", cp, off)
                } else {
                    format!("{}{}", cp, off)
                }
            };
            let end_str = build_pos(nearest_exon_cdna_pos, next_offset);
            format!("{}{}_{}ins{}", prefix, pos_str, end_str, ins_str)
        }
        _ => return None,
    };

    Some(notation)
}

/// Generate HGVSc notation for non-coding transcripts using `n.` prefix.
///
/// Non-coding transcripts (lncRNA, retained_intron, etc.) use cDNA position
/// directly with `n.` prefix, e.g. `ENST00000472807.5:n.1234A>G`.
pub fn hgvsc_noncoding(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    let prefix = format!("{}:n.", transcript_id);
    let (pos_min, pos_max) = if cdna_start <= cdna_end { (cdna_start, cdna_end) } else { (cdna_end, cdna_start) };
    let pos_str = if pos_min == pos_max {
        format!("{}", pos_min)
    } else {
        format!("{}_{}", pos_min, pos_max)
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!("{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char)
        }
        (Allele::Sequence(_), Allele::Deletion) => {
            format!("{}{}del", prefix, pos_str)
        }
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_pos = format!("{}_{}", pos_max, pos_max + 1);
            format!("{}{}ins{}",
                prefix, ins_pos, std::str::from_utf8(alt_bases).unwrap_or("?"))
        }
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            format!("{}{}delins{}",
                prefix, pos_str, std::str::from_utf8(alt_bases).unwrap_or("?"))
        }
        _ => return None,
    };
    Some(notation)
}

/// Generate HGVSc intronic notation for non-coding transcripts using `n.` prefix.
pub fn hgvsc_noncoding_intronic(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    hgvsc_noncoding_intronic_range(
        transcript_id, nearest_exon_cdna_pos, intron_offset,
        None, None, ref_allele, alt_allele,
    )
}

/// Generate HGVSc intronic notation for non-coding transcripts with optional end position.
pub fn hgvsc_noncoding_intronic_range(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    end_cdna_pos: Option<u64>,
    end_intron_offset: Option<i64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    let prefix = format!("{}:n.", transcript_id);
    let pos_str = if intron_offset > 0 {
        format!("{}+{}", nearest_exon_cdna_pos, intron_offset)
    } else {
        format!("{}{}", nearest_exon_cdna_pos, intron_offset)
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!("{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char)
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            if ref_bases.len() == 1 || end_intron_offset.is_none() {
                format!("{}{}del", prefix, pos_str)
            } else {
                let e_cdna = end_cdna_pos.unwrap_or(nearest_exon_cdna_pos);
                let e_offset = end_intron_offset.unwrap();
                let end_pos_str = if e_offset > 0 {
                    format!("{}+{}", e_cdna, e_offset)
                } else {
                    format!("{}{}", e_cdna, e_offset)
                };
                // Order: smaller offset first
                let (s, e) = if intron_offset < e_offset {
                    (pos_str.clone(), end_pos_str)
                } else {
                    (end_pos_str, pos_str.clone())
                };
                format!("{}{}_{}del", prefix, s, e)
            }
        }
        // Insertion in intron — show range: n.X+N_X+N+1insB
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_str = std::str::from_utf8(alt_bases).unwrap_or("?");
            let next_offset = intron_offset + 1;
            let end_pos_str = if next_offset > 0 {
                format!("{}+{}", nearest_exon_cdna_pos, next_offset)
            } else {
                format!("{}{}", nearest_exon_cdna_pos, next_offset)
            };
            format!("{}{}_{}ins{}", prefix, pos_str, end_pos_str, ins_str)
        }
        _ => return None,
    };
    Some(notation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsc_snv() {
        let result = hgvsc(
            "ENST00000001",
            151, 151,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51, None,
        );
        assert_eq!(result, Some("ENST00000001:c.101A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion() {
        let result = hgvsc(
            "ENST00000001",
            54, 56,
            &Allele::Sequence(b"ACG".to_vec()),
            &Allele::Deletion,
            51, None,
        );
        assert_eq!(result, Some("ENST00000001:c.4_6del".to_string()));
    }

    #[test]
    fn test_hgvsc_insertion() {
        let result = hgvsc(
            "ENST00000001",
            54, 53,
            &Allele::Deletion,
            &Allele::Sequence(b"TTT".to_vec()),
            51, None,
        );
        assert_eq!(result, Some("ENST00000001:c.3_4insTTT".to_string()));
    }

    #[test]
    fn test_hgvsc_5_utr() {
        let result = hgvsc(
            "ENST00000001",
            10, 10,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51, None,
        );
        // 5'UTR: position = cdna - coding_start = 10 - 51 = -41
        assert_eq!(result, Some("ENST00000001:c.-41A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_3_utr() {
        // coding_end = 1003 in cDNA, variant at cDNA 1021
        let result = hgvsc(
            "ENST00000001",
            1021, 1021,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Sequence(b"A".to_vec()),
            51, Some(1003),
        );
        // utr_offset = 1021 - 1003 = 18
        assert_eq!(result, Some("ENST00000001:c.*18G>A".to_string()));
    }

    #[test]
    fn test_hgvsc_intronic_donor() {
        // Variant 5 bases into intron, near donor (upstream exon)
        // nearest exon cDNA pos = 201 (end of exon 1), offset = +5
        // CDS pos = 201 - 51 + 1 = 151
        let result = hgvsc_intronic(
            "ENST00000001",
            201, 5,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Sequence(b"A".to_vec()),
            51, None,
        );
        assert_eq!(result, Some("ENST00000001:c.151+5G>A".to_string()));
    }

    #[test]
    fn test_hgvsc_intronic_acceptor() {
        // Variant 3 bases before exon 2, near acceptor
        // nearest exon cDNA pos = 202 (start of exon 2), offset = -3
        // CDS pos = 202 - 51 + 1 = 152
        let result = hgvsc_intronic(
            "ENST00000001",
            202, -3,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51, None,
        );
        assert_eq!(result, Some("ENST00000001:c.152-3A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion_3prime_shift() {
        // Sequence: AACTTTTGA at CDS positions 1-9 (coding_start=1, cdna positions 1-9)
        // Deletion of T at cdna position 4 (cds pos 4) should shift to position 7
        // because TTTT is repetitive — shift to the most 3' T
        let seq = "AACTTTTGA";
        let result = hgvsc_with_seq(
            "ENST00000001",
            4, 4, // cdna_start=4, cdna_end=4 (the first T in TTTT)
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Deletion,
            1, None,
            Some(seq), 0,
        );
        // Should shift from pos 4 to pos 7 (last T in the TTTT run)
        assert_eq!(result, Some("ENST00000001:c.7del".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion_no_shift() {
        // Non-repetitive: deletion at position 1 (A) — no shifting possible
        let seq = "ACGTACGT";
        let result = hgvsc_with_seq(
            "ENST00000001",
            1, 1,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Deletion,
            1, None,
            Some(seq), 0,
        );
        assert_eq!(result, Some("ENST00000001:c.1del".to_string()));
    }

    #[test]
    fn test_hgvsc_noncoding_snv() {
        let result = hgvsc_noncoding(
            "ENST00000472807.5",
            100, 100,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
        );
        assert_eq!(result, Some("ENST00000472807.5:n.100A>G".to_string()));
    }
}
