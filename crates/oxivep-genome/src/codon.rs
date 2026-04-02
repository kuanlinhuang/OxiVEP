use std::collections::HashMap;

/// Standard genetic code codon translation table.
pub struct CodonTable {
    table: HashMap<[u8; 3], u8>,
}

impl CodonTable {
    /// Create the standard genetic code (NCBI translation table 1).
    pub fn standard() -> Self {
        let codons = [
            // Phe
            (b"TTT", b'F'), (b"TTC", b'F'),
            // Leu
            (b"TTA", b'L'), (b"TTG", b'L'), (b"CTT", b'L'), (b"CTC", b'L'),
            (b"CTA", b'L'), (b"CTG", b'L'),
            // Ile
            (b"ATT", b'I'), (b"ATC", b'I'), (b"ATA", b'I'),
            // Met (start)
            (b"ATG", b'M'),
            // Val
            (b"GTT", b'V'), (b"GTC", b'V'), (b"GTA", b'V'), (b"GTG", b'V'),
            // Ser
            (b"TCT", b'S'), (b"TCC", b'S'), (b"TCA", b'S'), (b"TCG", b'S'),
            (b"AGT", b'S'), (b"AGC", b'S'),
            // Pro
            (b"CCT", b'P'), (b"CCC", b'P'), (b"CCA", b'P'), (b"CCG", b'P'),
            // Thr
            (b"ACT", b'T'), (b"ACC", b'T'), (b"ACA", b'T'), (b"ACG", b'T'),
            // Ala
            (b"GCT", b'A'), (b"GCC", b'A'), (b"GCA", b'A'), (b"GCG", b'A'),
            // Tyr
            (b"TAT", b'Y'), (b"TAC", b'Y'),
            // Stop
            (b"TAA", b'*'), (b"TAG", b'*'), (b"TGA", b'*'),
            // His
            (b"CAT", b'H'), (b"CAC", b'H'),
            // Gln
            (b"CAA", b'Q'), (b"CAG", b'Q'),
            // Asn
            (b"AAT", b'N'), (b"AAC", b'N'),
            // Lys
            (b"AAA", b'K'), (b"AAG", b'K'),
            // Asp
            (b"GAT", b'D'), (b"GAC", b'D'),
            // Glu
            (b"GAA", b'E'), (b"GAG", b'E'),
            // Cys
            (b"TGT", b'C'), (b"TGC", b'C'),
            // Trp
            (b"TGG", b'W'),
            // Arg
            (b"CGT", b'R'), (b"CGC", b'R'), (b"CGA", b'R'), (b"CGG", b'R'),
            (b"AGA", b'R'), (b"AGG", b'R'),
            // Gly
            (b"GGT", b'G'), (b"GGC", b'G'), (b"GGA", b'G'), (b"GGG", b'G'),
        ];

        let mut table = HashMap::with_capacity(64);
        for (codon, aa) in codons {
            table.insert(*codon, aa);
        }
        Self { table }
    }

    /// Translate a 3-base codon to an amino acid.
    /// Returns 'X' for unknown codons (e.g., containing N).
    pub fn translate(&self, codon: &[u8; 3]) -> u8 {
        // Uppercase the codon
        let upper = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        *self.table.get(&upper).unwrap_or(&b'X')
    }

    /// Translate a DNA sequence to a protein sequence.
    /// The sequence length must be a multiple of 3.
    pub fn translate_seq(&self, dna: &[u8]) -> Vec<u8> {
        dna.chunks_exact(3)
            .map(|chunk| {
                let codon: [u8; 3] = [chunk[0], chunk[1], chunk[2]];
                self.translate(&codon)
            })
            .collect()
    }

    /// Check if a codon is a stop codon.
    pub fn is_stop(&self, codon: &[u8; 3]) -> bool {
        self.translate(codon) == b'*'
    }

    /// Check if a codon is the standard start codon (ATG).
    pub fn is_start(codon: &[u8; 3]) -> bool {
        let upper = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        upper == *b"ATG"
    }
}

impl Default for CodonTable {
    fn default() -> Self {
        Self::standard()
    }
}

/// Convert single-letter amino acid code to three-letter code.
pub fn aa_one_to_three(aa: u8) -> &'static str {
    match aa {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'E' => "Glu",
        b'Q' => "Gln",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'*' => "Ter",
        b'X' => "Xaa",
        _ => "???",
    }
}

/// Format a ref/alt codon pair with changed bases UPPERCASE, unchanged lowercase.
/// Matches VEP convention: e.g., GCA/GAA → "gCa/gAa"
pub fn format_codon_change(ref_codon: &[u8; 3], alt_codon: &[u8; 3]) -> (String, String) {
    let mut ref_display = String::with_capacity(3);
    let mut alt_display = String::with_capacity(3);
    for i in 0..3 {
        if ref_codon[i].to_ascii_uppercase() != alt_codon[i].to_ascii_uppercase() {
            ref_display.push((ref_codon[i] as char).to_ascii_uppercase());
            alt_display.push((alt_codon[i] as char).to_ascii_uppercase());
        } else {
            ref_display.push((ref_codon[i] as char).to_ascii_lowercase());
            alt_display.push((alt_codon[i] as char).to_ascii_lowercase());
        }
    }
    (ref_display, alt_display)
}

/// Compute the reverse complement of a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_codon_table() {
        let table = CodonTable::standard();
        assert_eq!(table.translate(b"ATG"), b'M');
        assert_eq!(table.translate(b"TAA"), b'*');
        assert_eq!(table.translate(b"TAG"), b'*');
        assert_eq!(table.translate(b"TGA"), b'*');
        assert_eq!(table.translate(b"TTT"), b'F');
        assert_eq!(table.translate(b"GGG"), b'G');
    }

    #[test]
    fn test_translate_seq() {
        let table = CodonTable::standard();
        let protein = table.translate_seq(b"ATGTTTTAA");
        assert_eq!(protein, b"MF*");
    }

    #[test]
    fn test_unknown_codon() {
        let table = CodonTable::standard();
        assert_eq!(table.translate(b"NNN"), b'X');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_aa_one_to_three() {
        assert_eq!(aa_one_to_three(b'M'), "Met");
        assert_eq!(aa_one_to_three(b'*'), "Ter");
        assert_eq!(aa_one_to_three(b'X'), "Xaa");
    }

    #[test]
    fn test_format_codon_change() {
        // GCA -> GAA: position 1 changes (C->A)
        let (r, a) = format_codon_change(b"GCA", b"GAA");
        assert_eq!(r, "gCa");
        assert_eq!(a, "gAa");

        // ATG -> GTG: position 0 changes (A->G)
        let (r, a) = format_codon_change(b"ATG", b"GTG");
        assert_eq!(r, "Atg");
        assert_eq!(a, "Gtg");

        // GCT -> GCC: position 2 changes (T->C)
        let (r, a) = format_codon_change(b"GCT", b"GCC");
        assert_eq!(r, "gcT");
        assert_eq!(a, "gcC");

        // All positions change
        let (r, a) = format_codon_change(b"AAA", b"TTT");
        assert_eq!(r, "AAA");
        assert_eq!(a, "TTT");
    }

    #[test]
    fn test_is_stop_and_start() {
        let table = CodonTable::standard();
        assert!(table.is_stop(b"TAA"));
        assert!(table.is_stop(b"TAG"));
        assert!(table.is_stop(b"TGA"));
        assert!(!table.is_stop(b"ATG"));
        assert!(CodonTable::is_start(b"ATG"));
        assert!(!CodonTable::is_start(b"TTT"));
    }
}
