use oxivep_genome::codon::{aa_one_to_three, CodonTable};

/// Generate HGVSp (protein) notation.
///
/// Format: ENSP00000001:p.Arg41Lys (missense)
///         ENSP00000001:p.Arg41Ter (stop gained)
///         ENSP00000001:p.Arg41= (synonymous)
///         ENSP00000001:p.Arg41fs (frameshift)
pub fn hgvsp(
    protein_id: &str,
    protein_pos: u64,
    ref_aa: u8,
    alt_aa: u8,
    is_frameshift: bool,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);
    let ref_aa3 = aa_one_to_three(ref_aa);

    if is_frameshift {
        return Some(format!("{}{}{}fs", prefix, ref_aa3, protein_pos));
    }

    if ref_aa == alt_aa {
        // Synonymous
        return Some(format!("{}{}{}=", prefix, ref_aa3, protein_pos));
    }

    let alt_aa3 = aa_one_to_three(alt_aa);

    if alt_aa == b'*' {
        // Stop gained
        return Some(format!("{}{}{}{}",prefix, ref_aa3, protein_pos, alt_aa3));
    }

    if ref_aa == b'*' {
        // Stop lost - extension
        return Some(format!("{}{}{}ext*?",prefix, alt_aa3, protein_pos));
    }

    // Missense
    Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3))
}

/// Generate HGVSp notation for a frameshift variant.
///
/// Scans the frameshifted sequence to find the first changed amino acid and
/// the position of the new stop codon.
///
/// Format: ENSP00000001:p.Ala498ProfsTer28
///   - Ala498 = first amino acid that changes (ref)
///   - Pro = new amino acid at that position
///   - Ter28 = new stop codon 28 positions downstream
pub fn hgvsp_frameshift(
    protein_id: &str,
    ref_translateable: &[u8],
    alt_translateable: &[u8],
    affected_codon_start: usize, // 0-based codon index where the frameshift starts
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);
    let codon_table = CodonTable::standard();

    // Translate both sequences from the affected codon onwards
    let ref_start = affected_codon_start * 3;
    if ref_start + 3 > ref_translateable.len() {
        return None;
    }

    let ref_peptide: Vec<u8> = ref_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    let alt_peptide: Vec<u8> = alt_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    // Find the first position where amino acids differ
    let mut first_changed_offset = 0;
    for i in 0..ref_peptide.len().min(alt_peptide.len()) {
        if ref_peptide[i] != alt_peptide[i] {
            first_changed_offset = i;
            break;
        }
        // If we reach a stop codon in ref before finding a change,
        // the change starts at this position
        if ref_peptide[i] == b'*' {
            first_changed_offset = i;
            break;
        }
        first_changed_offset = i + 1;
    }

    if first_changed_offset >= ref_peptide.len() && first_changed_offset >= alt_peptide.len() {
        return None;
    }

    let first_changed_pos = affected_codon_start + first_changed_offset + 1; // 1-based
    let ref_aa = if first_changed_offset < ref_peptide.len() {
        ref_peptide[first_changed_offset]
    } else {
        b'X'
    };
    let alt_aa = if first_changed_offset < alt_peptide.len() {
        alt_peptide[first_changed_offset]
    } else {
        b'X'
    };

    let ref_aa3 = aa_one_to_three(ref_aa);
    let alt_aa3 = aa_one_to_three(alt_aa);

    // Find the new stop codon position in the alt sequence.
    // Skip if the sequence is mostly N/X (unresolved) — can't determine real ProfsTer.
    let mut stop_dist = None;
    let unresolved_count = alt_peptide[first_changed_offset..].iter()
        .take(10)
        .filter(|&&aa| aa == b'X')
        .count();
    let mostly_unresolved = unresolved_count > 5;
    if !mostly_unresolved {
        for i in first_changed_offset..alt_peptide.len() {
            if alt_peptide[i] == b'*' {
                stop_dist = Some(i - first_changed_offset + 1);
                break;
            }
        }
    }

    match stop_dist {
        Some(d) => Some(format!("{}{}{}{}fsTer{}", prefix, ref_aa3, first_changed_pos, alt_aa3, d)),
        None => Some(format!("{}{}{}{}fs", prefix, ref_aa3, first_changed_pos, alt_aa3)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsp_missense() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'K', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Lys".to_string()));
    }

    #[test]
    fn test_hgvsp_synonymous() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41=".to_string()));
    }

    #[test]
    fn test_hgvsp_stop_gained() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'*', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Ter".to_string()));
    }

    #[test]
    fn test_hgvsp_frameshift() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'X', true);
        assert_eq!(result, Some("ENSP00000001:p.Arg41fs".to_string()));
    }

    #[test]
    fn test_hgvsp_stop_lost() {
        let result = hgvsp("ENSP00000001", 100, b'*', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg100ext*?".to_string()));
    }
}
