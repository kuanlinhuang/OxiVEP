/// Integration tests validating OxiVEP against Ensembl VEP test patterns.
///
/// These tests verify that OxiVEP's VCF parsing, allele normalization,
/// and variant representation matches the behavior documented in
/// ensembl-vep's Parser_VCF.t test suite.
use oxivep_core::{Allele, Impact, Strand};
use oxivep_io::output;
use oxivep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature, VcfFields};
use oxivep_io::vcf::{parse_vcf_line, VcfParser};

// =============================================================================
// VCF Parsing — matches ensembl-vep Parser_VCF.t assertions
// =============================================================================

#[test]
fn test_vep_snv_basic_parsing() {
    // From ensembl-vep test.vcf line 1: 21 25585733 rs142513484 C T
    // VEP Parser_VCF.t expects: chr=21, start=25585733, end=25585733, allele_string=C/T
    let line = "21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "21");
    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.position.end, 25585733);
    assert_eq!(vf.allele_string, "C/T");
    assert_eq!(vf.variation_name, Some("rs142513484".to_string()));
    assert_eq!(vf.position.strand, Strand::Forward);
    assert!(!vf.is_indel());
}

#[test]
fn test_vep_insertion_parsing() {
    // From test_not_ordered.vcf: 21 30000016 rs202173120 T TCA
    // VEP behavior: strip shared first base T → ref="-", alt="CA", start=30000017, end=30000016
    let line = "21\t30000016\trs202173120\tT\tTCA\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000017);
    assert_eq!(vf.position.end, 30000016); // insertion: end < start
    assert_eq!(vf.allele_string, "-/CA");
    assert_eq!(vf.ref_allele, Allele::Deletion);
    assert_eq!(vf.alt_alleles[0], Allele::Sequence(b"CA".to_vec()));
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_insertion_single_base() {
    // 21 30000105 rs34788396 C CC → insertion of 1 base
    let line = "21\t30000105\trs34788396\tC\tCC\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000106);
    assert_eq!(vf.position.end, 30000105);
    assert_eq!(vf.allele_string, "-/C");
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_large_insertion() {
    // 21 30000094 rs960445955 C CATATTCTCCCCTATT → large insertion
    let line = "21\t30000094\trs960445955\tC\tCATATTCTCCCCTATT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 30000095);
    assert_eq!(vf.position.end, 30000094);
    assert_eq!(vf.allele_string, "-/ATATTCTCCCCTATT");
    assert!(vf.is_insertion());
}

#[test]
fn test_vep_deletion_parsing() {
    // Simulated deletion: ref=GGA, alt=G at chr22:19353532
    // VEP strips shared G: ref=GA, alt=-, start=19353533
    let line = "22\t19353532\trs1162718428\tGGA\tG\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "22");
    assert_eq!(vf.position.start, 19353533);
    assert_eq!(vf.position.end, 19353534);
    assert_eq!(vf.allele_string, "GA/-");
    assert_eq!(vf.ref_allele, Allele::Sequence(b"GA".to_vec()));
    assert_eq!(vf.alt_alleles[0], Allele::Deletion);
    assert!(vf.is_deletion());
}

#[test]
fn test_vep_multi_allelic_snv() {
    // Multi-allelic SNV: two alt alleles
    let line = "21\t25585733\t.\tC\tT,G\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.allele_string, "C/T/G");
    assert_eq!(vf.alt_alleles.len(), 2);
    assert_eq!(vf.alt_alleles[0], Allele::Sequence(b"T".to_vec()));
    assert_eq!(vf.alt_alleles[1], Allele::Sequence(b"G".to_vec()));
}

#[test]
fn test_vep_multi_allelic_indel() {
    // Multi-allelic with indel: ref=ACG, alt=A,ACGT
    // All share first base A, so strip: ref=CG, alt=-,CGT
    let line = "21\t25585733\t.\tACG\tA,ACGT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585734);
    assert_eq!(vf.allele_string, "CG/-/CGT");
}

#[test]
fn test_vep_star_allele_handling() {
    // Star allele should be preserved during normalization
    let line = "21\t25585733\t.\tACG\tA,*\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585734);
    assert!(vf.allele_string.contains("*"));
}

#[test]
fn test_vep_non_variant_site() {
    // REF-only site: alt is "."
    let line = "21\t25585733\t.\tC\t.\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.allele_string, "C");
}

#[test]
fn test_vep_chr_prefix_handling() {
    // Chromosome with "chr" prefix should be preserved
    let line = "chr21\t25585733\trs1\tC\tT\t.\tPASS\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "chr21");
}

#[test]
fn test_vep_mt_chromosome() {
    // Mitochondrial chromosome
    let line = "MT\t4472\t.\tT\tA\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.chromosome, "MT");
    assert_eq!(vf.position.start, 4472);
    assert_eq!(vf.allele_string, "T/A");
}

#[test]
fn test_vep_mnv_parsing() {
    // MNV: multi-nucleotide variant (same ref/alt length >1)
    let line = "21\t25585733\t.\tAC\tGT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.position.end, 25585734);
    assert_eq!(vf.allele_string, "AC/GT");
    assert!(!vf.is_indel());
}

#[test]
fn test_vep_complex_indel() {
    // Complex indel: different ref/alt lengths, not simple ins/del
    let line = "21\t25585733\t.\tACG\tTT\t.\t.\t.";
    let vf = parse_vcf_line(line).unwrap();

    // No shared first base, so no stripping
    assert_eq!(vf.position.start, 25585733);
    assert_eq!(vf.allele_string, "ACG/TT");
}

// =============================================================================
// VCF Parser — multi-line / header handling
// =============================================================================

#[test]
fn test_vep_full_vcf_parser() {
    // Simulates the structure of ensembl-vep test.vcf
    let vcf = "##fileformat=VCFv4.1\n\
##contig=<ID=21,assembly=GCF_000001405.26,length=46709983>\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\n\
21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0\n\
21\t25587701\trs187353664\tT\tC\t.\t.\t.\tGT\t0|0\n\
21\t25587758\trs116645811\tG\tA\t.\t.\t.\tGT\t0|0\n";

    let mut parser = VcfParser::new(vcf.as_bytes()).unwrap();

    // Should have 4 header lines
    assert_eq!(parser.header_lines().len(), 4);
    assert!(parser.header_lines()[0].starts_with("##fileformat"));
    assert!(parser.header_lines()[3].starts_with("#CHROM"));

    let variants = parser.read_all().unwrap();
    assert_eq!(variants.len(), 3);

    assert_eq!(variants[0].variation_name.as_deref(), Some("rs142513484"));
    assert_eq!(variants[0].position.start, 25585733);
    assert_eq!(variants[0].allele_string, "C/T");

    assert_eq!(variants[1].variation_name.as_deref(), Some("rs187353664"));
    assert_eq!(variants[2].variation_name.as_deref(), Some("rs116645811"));
}

#[test]
fn test_vep_unordered_variants() {
    // Variants from test_not_ordered.vcf — mixed chromosomes and positions
    let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t230710034\trs894744940\tC\tT\t.\t.\t.\n\
21\t25587758\trs116645811\tG\tA\t.\t.\t.\n\
21\t25592836\trs1135638\tG\tA\t.\t.\t.\n\
1\t230710045\trs754176245\tC\tT\t.\t.\t.\n";

    let mut parser = VcfParser::new(vcf.as_bytes()).unwrap();
    let variants = parser.read_all().unwrap();

    assert_eq!(variants.len(), 4);
    // Verify all parsed regardless of order
    assert_eq!(variants[0].position.chromosome, "1");
    assert_eq!(variants[1].position.chromosome, "21");
    assert_eq!(variants[2].position.chromosome, "21");
    assert_eq!(variants[3].position.chromosome, "1");
}

#[test]
fn test_vep_indel_variety() {
    // Various indel types from test_not_ordered.vcf
    let test_cases = vec![
        // (line, expected_start, expected_allele_string)
        ("21\t30000016\trs202173120\tT\tTCA\t.\t.\t.", 30000017, "-/CA"),
        ("21\t30000018\trs144102269\tA\tAAT\t.\t.\t.", 30000019, "-/AT"),
        ("21\t30000020\trs35748481\tG\tGAT\t.\t.\t.", 30000021, "-/AT"),
        ("21\t30000105\trs34788396\tC\tCC\t.\t.\t.", 30000106, "-/C"),
    ];

    for (line, expected_start, expected_alleles) in test_cases {
        let vf = parse_vcf_line(line).unwrap();
        assert_eq!(
            vf.position.start, expected_start,
            "Failed for: {} — got start={}", line, vf.position.start
        );
        assert_eq!(
            vf.allele_string, expected_alleles,
            "Failed for: {} — got alleles={}", line, vf.allele_string
        );
    }
}

#[test]
fn test_vep_filter_status_preserved() {
    let line = "21\t25585733\trs1\tC\tT\t30\tPASS\tDP=50";
    let vf = parse_vcf_line(line).unwrap();

    let fields = vf.vcf_fields.unwrap();
    assert_eq!(fields.filter, "PASS");
    assert_eq!(fields.qual, "30");
    assert_eq!(fields.info, "DP=50");
}

#[test]
fn test_vep_genotype_fields_preserved() {
    let line = "21\t25585733\trs1\tC\tT\t.\t.\t.\tGT\t0|0";
    let vf = parse_vcf_line(line).unwrap();

    let fields = vf.vcf_fields.unwrap();
    assert_eq!(fields.rest.len(), 2);
    assert_eq!(fields.rest[0], "GT");
    assert_eq!(fields.rest[1], "0|0");
}

// =============================================================================
// Consequence prediction validation
// =============================================================================

#[test]
fn test_all_consequence_types_have_valid_so_terms() {
    use oxivep_core::Consequence;

    // Verify all SO terms are valid and round-trip
    let all = [
        Consequence::TranscriptAblation,
        Consequence::SpliceAcceptorVariant,
        Consequence::SpliceDonorVariant,
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::StopLost,
        Consequence::StartLost,
        Consequence::TranscriptAmplification,
        Consequence::FeatureElongation,
        Consequence::FeatureTruncation,
        Consequence::InframeInsertion,
        Consequence::InframeDeletion,
        Consequence::MissenseVariant,
        Consequence::ProteinAlteringVariant,
        Consequence::SpliceRegionVariant,
        Consequence::SpliceDonorFifthBaseVariant,
        Consequence::SpliceDonorRegionVariant,
        Consequence::SplicePolypyrimidineTractVariant,
        Consequence::IncompleteTerminalCodonVariant,
        Consequence::StartRetainedVariant,
        Consequence::StopRetainedVariant,
        Consequence::SynonymousVariant,
        Consequence::CodingSequenceVariant,
        Consequence::MatureMirnaVariant,
        Consequence::FivePrimeUtrVariant,
        Consequence::ThreePrimeUtrVariant,
        Consequence::NonCodingTranscriptExonVariant,
        Consequence::IntronVariant,
        Consequence::NmdTranscriptVariant,
        Consequence::NonCodingTranscriptVariant,
        Consequence::CodingTranscriptVariant,
        Consequence::UpstreamGeneVariant,
        Consequence::DownstreamGeneVariant,
        Consequence::TfbsAblation,
        Consequence::TfbsAmplification,
        Consequence::TfBindingSiteVariant,
        Consequence::RegulatoryRegionAblation,
        Consequence::RegulatoryRegionAmplification,
        Consequence::RegulatoryRegionVariant,
        Consequence::IntergenicVariant,
        Consequence::SequenceVariant,
    ];

    for c in &all {
        let term = c.so_term();
        let parsed = Consequence::from_so_term(term);
        assert_eq!(parsed, Some(*c), "SO term round-trip failed for {:?} ({})", c, term);
    }

    // Verify strict ordering by rank
    for i in 0..all.len() - 1 {
        assert!(
            all[i].rank() < all[i + 1].rank(),
            "{:?} (rank {}) should be more severe than {:?} (rank {})",
            all[i], all[i].rank(), all[i + 1], all[i + 1].rank()
        );
    }
}

#[test]
fn test_vep_impact_classification() {
    use oxivep_core::{Consequence, Impact};

    // Verify impact matches VEP's classification
    let high = [
        Consequence::TranscriptAblation,
        Consequence::SpliceAcceptorVariant,
        Consequence::SpliceDonorVariant,
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::StopLost,
        Consequence::StartLost,
    ];
    for c in &high {
        assert_eq!(c.impact(), Impact::High, "{:?} should be HIGH", c);
    }

    let moderate = [
        Consequence::InframeInsertion,
        Consequence::InframeDeletion,
        Consequence::MissenseVariant,
        Consequence::ProteinAlteringVariant,
    ];
    for c in &moderate {
        assert_eq!(c.impact(), Impact::Moderate, "{:?} should be MODERATE", c);
    }

    let low = [
        Consequence::SpliceRegionVariant,
        Consequence::SynonymousVariant,
        Consequence::StopRetainedVariant,
        Consequence::StartRetainedVariant,
    ];
    for c in &low {
        assert_eq!(c.impact(), Impact::Low, "{:?} should be LOW", c);
    }

    let modifier = [
        Consequence::IntronVariant,
        Consequence::UpstreamGeneVariant,
        Consequence::DownstreamGeneVariant,
        Consequence::IntergenicVariant,
        Consequence::FivePrimeUtrVariant,
        Consequence::ThreePrimeUtrVariant,
    ];
    for c in &modifier {
        assert_eq!(c.impact(), Impact::Modifier, "{:?} should be MODIFIER", c);
    }
}

// =============================================================================
// Codon table validation
// =============================================================================

#[test]
fn test_all_64_codons_translate() {
    use oxivep_genome::CodonTable;

    let table = CodonTable::standard();
    let bases = [b'A', b'C', b'G', b'T'];

    let mut count = 0;
    for &b1 in &bases {
        for &b2 in &bases {
            for &b3 in &bases {
                let codon = [b1, b2, b3];
                let aa = table.translate(&codon);
                assert_ne!(aa, b'X', "Codon {:?} should translate", std::str::from_utf8(&codon));
                count += 1;
            }
        }
    }
    assert_eq!(count, 64);
}

#[test]
fn test_codon_table_matches_ncbi() {
    use oxivep_genome::CodonTable;

    let table = CodonTable::standard();

    // NCBI standard genetic code verification
    assert_eq!(table.translate(b"ATG"), b'M'); // Met (start)
    assert_eq!(table.translate(b"TAA"), b'*'); // Stop
    assert_eq!(table.translate(b"TAG"), b'*'); // Stop
    assert_eq!(table.translate(b"TGA"), b'*'); // Stop

    // Verify a selection of other codons
    assert_eq!(table.translate(b"TTT"), b'F'); // Phe
    assert_eq!(table.translate(b"TTC"), b'F'); // Phe
    assert_eq!(table.translate(b"CTG"), b'L'); // Leu (most common)
    assert_eq!(table.translate(b"GAT"), b'D'); // Asp
    assert_eq!(table.translate(b"GAC"), b'D'); // Asp
    assert_eq!(table.translate(b"TGG"), b'W'); // Trp (only codon)
    assert_eq!(table.translate(b"CGT"), b'R'); // Arg
    assert_eq!(table.translate(b"AGA"), b'R'); // Arg
    assert_eq!(table.translate(b"GGG"), b'G'); // Gly
}

// =============================================================================
// CSQ Output Formatting — matches Ensembl VEP output format
// =============================================================================

/// Helper to build a mock VariationFeature for CSQ output testing.
fn mock_vf_missense() -> VariationFeature {
    use oxivep_core::Consequence;

    VariationFeature {
        position: oxivep_core::GenomicPosition::new("1", 65568, 65568, Strand::Forward),
        allele_string: "A/C".into(),
        ref_allele: Allele::from_str("A"),
        alt_alleles: vec![Allele::from_str("C")],
        variation_name: Some("1_65568_A/C".into()),
        vcf_line: None,
        vcf_fields: None,
        transcript_variations: vec![TranscriptVariation {
            transcript_id: "ENST00000641515.2".into(),
            gene_id: "ENSG00000186092".into(),
            gene_symbol: Some("OR4F5".into()),
            biotype: "protein_coding".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::from_str("C"),
                consequences: vec![Consequence::MissenseVariant],
                impact: Impact::Moderate,
                cdna_position: Some((64, 64)),
                cds_position: Some((4, 4)),
                protein_position: Some((2, 2)),
                amino_acids: Some(("K".into(), "Q".into())),
                codons: Some(("Aag".into(), "Cag".into())),
                exon: Some((2, 3)),
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: Some("tolerated_low_confidence(0.06)".into()),
                polyphen: Some("benign(0)".into()),
            }],
            canonical: false,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: Some("NM_001005484.2".into()),
            mane_plus_clinical: None,
            tsl: None,
            appris: Some("P1".into()),
            ccds: None,
            symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:14825".into()),
            flags: vec![],
        }],
        existing_variants: vec![],
        minimised: false,
        most_severe_consequence: Some(Consequence::MissenseVariant),
    }
}

#[test]
fn test_csq_missense_field_values() {
    let vf = mock_vf_missense();
    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);

    // Parse the CSQ string into fields
    let fields: Vec<&str> = csq.split('|').collect();

    // Validate key fields match VEP output format
    assert_eq!(fields[0], "C", "Allele");
    assert_eq!(fields[1], "missense_variant", "Consequence");
    assert_eq!(fields[2], "MODERATE", "IMPACT");
    assert_eq!(fields[3], "OR4F5", "SYMBOL");
    assert_eq!(fields[4], "ENSG00000186092", "Gene");
    assert_eq!(fields[5], "Transcript", "Feature_type");
    assert_eq!(fields[6], "ENST00000641515.2", "Feature");
    assert_eq!(fields[7], "protein_coding", "BIOTYPE");
    assert_eq!(fields[8], "2/3", "EXON");
    assert_eq!(fields[9], "", "INTRON should be empty");
    assert_eq!(fields[12], "64", "cDNA_position");
    assert_eq!(fields[13], "4", "CDS_position");
    assert_eq!(fields[14], "2", "Protein_position");
    assert_eq!(fields[15], "K/Q", "Amino_acids");
    assert_eq!(fields[16], "Aag/Cag", "Codons");
    assert_eq!(fields[18], "A", "REF_ALLELE");
    assert_eq!(fields[19], "A/C", "UPLOADED_ALLELE");
    assert_eq!(fields[21], "1", "STRAND");
    // CANONICAL is now at index 23
    assert_eq!(fields[24], "HGNC", "SYMBOL_SOURCE");
    assert_eq!(fields[25], "HGNC:14825", "HGNC_ID");
    assert_eq!(fields[26], "MANE_Select", "MANE");
    assert_eq!(fields[27], "NM_001005484.2", "MANE_SELECT");
    assert_eq!(fields[30], "P1", "APPRIS");
    // CCDS at 31, ENSP at 32, SOURCE at 33, HGVS_OFFSET at 34
    assert_eq!(fields[35], "tolerated_low_confidence(0.06)", "SIFT");
    assert_eq!(fields[36], "benign(0)", "PolyPhen");
}

#[test]
fn test_csq_frameshift_codon_format() {
    use oxivep_core::Consequence;

    // Simulate VEP output: frameshift with codons "gAg/gg"
    let vf = VariationFeature {
        position: oxivep_core::GenomicPosition::new("3", 319780, 319781, Strand::Forward),
        allele_string: "GA/G".into(),
        ref_allele: Allele::from_str("A"),
        alt_alleles: vec![Allele::Deletion],
        variation_name: Some("3_319781_A/-".into()),
        vcf_line: None,
        vcf_fields: None,
        transcript_variations: vec![TranscriptVariation {
            transcript_id: "ENST00000256509.7".into(),
            gene_id: "ENSG00000134121".into(),
            gene_symbol: Some("CHL1".into()),
            biotype: "protein_coding".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::Deletion,
                consequences: vec![Consequence::FrameshiftVariant],
                impact: Impact::High,
                cdna_position: Some((480, 480)),
                cds_position: Some((5, 5)),
                protein_position: Some((2, 2)),
                amino_acids: Some(("E".into(), "X".into())),
                codons: Some(("gAg".into(), "gg".into())),
                exon: Some((3, 28)),
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: Some("NM_006614.4".into()),
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: Some("P3".into()),
            ccds: None,
            symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:1939".into()),
            flags: vec![],
        }],
        existing_variants: vec![],
        minimised: false,
        most_severe_consequence: Some(Consequence::FrameshiftVariant),
    };

    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let fields: Vec<&str> = csq.split('|').collect();

    assert_eq!(fields[0], "-", "Allele for deletion should be '-'");
    assert_eq!(fields[1], "frameshift_variant", "Consequence");
    assert_eq!(fields[2], "HIGH", "IMPACT");
    assert_eq!(fields[8], "3/28", "EXON");
    assert_eq!(fields[12], "480", "cDNA_position");
    assert_eq!(fields[13], "5", "CDS_position");
    assert_eq!(fields[14], "2", "Protein_position");
    assert_eq!(fields[15], "E/X", "Amino_acids");
    assert_eq!(fields[16], "gAg/gg", "Codons - frameshift format");
    assert_eq!(fields[18], "A", "REF_ALLELE");
    assert_eq!(fields[19], "A/-", "UPLOADED_ALLELE");
    assert_eq!(fields[21], "1", "STRAND");
}

#[test]
fn test_csq_header_includes_new_fields() {
    let header = output::csq_header_line(output::DEFAULT_CSQ_FIELDS);
    assert!(header.contains("REF_ALLELE"), "Header should include REF_ALLELE");
    assert!(header.contains("UPLOADED_ALLELE"), "Header should include UPLOADED_ALLELE");
    assert!(header.contains("FLAGS"), "Header should include FLAGS");
    assert!(header.contains("SYMBOL_SOURCE"), "Header should include SYMBOL_SOURCE");
    assert!(header.contains("HGNC_ID"), "Header should include HGNC_ID");
    assert!(header.contains("MANE|MANE_SELECT"), "Header should include MANE fields");
    assert!(header.contains("TSL"), "Header should include TSL");
    assert!(header.contains("APPRIS"), "Header should include APPRIS");
    assert!(header.contains("TRANSCRIPTION_FACTORS"), "Header should end with TRANSCRIPTION_FACTORS");
    assert!(header.contains("CANONICAL"), "Header should include CANONICAL");
    assert!(header.contains("CCDS"), "Header should include CCDS");
    assert!(header.contains("ENSP"), "Header should include ENSP");
    assert!(header.contains("SOURCE"), "Header should include SOURCE");
    assert!(header.contains("HGVS_OFFSET"), "Header should include HGVS_OFFSET");
}

#[test]
fn test_csq_field_count_matches_vep() {
    // Extended field set includes all VEP fields plus CANONICAL, CCDS, ENSP, SOURCE, HGVS_OFFSET
    assert_eq!(output::DEFAULT_CSQ_FIELDS.len(), 47, "DEFAULT_CSQ_FIELDS should have 47 fields");

    // Verify formatting produces 47 pipe-delimited values
    let vf = mock_vf_missense();
    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let field_count = csq.split('|').count();
    assert_eq!(field_count, 47, "CSQ output should have 47 pipe-delimited fields");
}

#[test]
fn test_csq_missense_full_42_field_match() {
    // Exact expected VEP output for the missense entry from SlDf20upKZNV52SS.vcf:
    // C|missense_variant|MODERATE|OR4F5|ENSG00000186092|Transcript|ENST00000641515.2|
    // protein_coding|2/3||||64|4|2|K/Q|Aag/Cag||A|A/C||1||HGNC|HGNC:14825|
    // MANE_Select|NM_001005484.2|||P1|tolerated_low_confidence(0.06)|benign(0)||||||||||
    let vf = mock_vf_missense();
    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let fields: Vec<&str> = csq.split('|').collect();

    // VEP expected values (all 47 fields)
    let expected: Vec<&str> = vec![
        "C",                                  // 0:  Allele
        "missense_variant",                   // 1:  Consequence
        "MODERATE",                           // 2:  IMPACT
        "OR4F5",                              // 3:  SYMBOL
        "ENSG00000186092",                    // 4:  Gene
        "Transcript",                         // 5:  Feature_type
        "ENST00000641515.2",                  // 6:  Feature
        "protein_coding",                     // 7:  BIOTYPE
        "2/3",                                // 8:  EXON
        "",                                   // 9:  INTRON
        "",                                   // 10: HGVSc
        "",                                   // 11: HGVSp
        "64",                                 // 12: cDNA_position
        "4",                                  // 13: CDS_position
        "2",                                  // 14: Protein_position
        "K/Q",                                // 15: Amino_acids
        "Aag/Cag",                            // 16: Codons
        "",                                   // 17: Existing_variation
        "A",                                  // 18: REF_ALLELE
        "A/C",                                // 19: UPLOADED_ALLELE
        "",                                   // 20: DISTANCE
        "1",                                  // 21: STRAND
        "",                                   // 22: FLAGS
        "",                                   // 23: CANONICAL
        "HGNC",                               // 24: SYMBOL_SOURCE
        "HGNC:14825",                         // 25: HGNC_ID
        "MANE_Select",                        // 26: MANE
        "NM_001005484.2",                     // 27: MANE_SELECT
        "",                                   // 28: MANE_PLUS_CLINICAL
        "",                                   // 29: TSL
        "P1",                                 // 30: APPRIS
        "",                                   // 31: CCDS
        "",                                   // 32: ENSP
        "",                                   // 33: SOURCE
        "",                                   // 34: HGVS_OFFSET
        "tolerated_low_confidence(0.06)",     // 35: SIFT
        "benign(0)",                          // 36: PolyPhen
        "",                                   // 37: AF
        "",                                   // 38: CLIN_SIG
        "",                                   // 39: SOMATIC
        "",                                   // 40: PHENO
        "",                                   // 41: PUBMED
        "",                                   // 42: MOTIF_NAME
        "",                                   // 43: MOTIF_POS
        "",                                   // 44: HIGH_INF_POS
        "",                                   // 45: MOTIF_SCORE_CHANGE
        "",                                   // 46: TRANSCRIPTION_FACTORS
    ];

    assert_eq!(fields.len(), expected.len(),
        "Field count mismatch: got {}, expected {}", fields.len(), expected.len());

    for (i, (got, exp)) in fields.iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, exp,
            "Field {} ({}) mismatch: got {:?}, expected {:?}",
            i, output::DEFAULT_CSQ_FIELDS[i], got, exp);
    }
}

#[test]
fn test_csq_frameshift_full_42_field_match() {
    use oxivep_core::Consequence;

    let vf = VariationFeature {
        position: oxivep_core::GenomicPosition::new("3", 319780, 319781, Strand::Forward),
        allele_string: "GA/G".into(),
        ref_allele: Allele::from_str("A"),
        alt_alleles: vec![Allele::Deletion],
        variation_name: Some("3_319781_A/-".into()),
        vcf_line: None,
        vcf_fields: None,
        transcript_variations: vec![TranscriptVariation {
            transcript_id: "ENST00000421198.5".into(),
            gene_id: "ENSG00000134121".into(),
            gene_symbol: Some("CHL1".into()),
            biotype: "protein_coding".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::Deletion,
                consequences: vec![Consequence::FrameshiftVariant],
                impact: Impact::High,
                cdna_position: Some((258, 258)),
                cds_position: Some((5, 5)),
                protein_position: Some((2, 2)),
                amino_acids: Some(("E".into(), "X".into())),
                codons: Some(("gAg".into(), "gg".into())),
                exon: Some((3, 5)),
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(4),
            appris: None,
            ccds: None,
            symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:1939".into()),
            flags: vec!["cds_end_NF".into()],
        }],
        existing_variants: vec![],
        minimised: false,
        most_severe_consequence: Some(Consequence::FrameshiftVariant),
    };

    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let fields: Vec<&str> = csq.split('|').collect();

    // VEP expected for ENST00000421198.5 frameshift entry:
    // -|frameshift_variant|HIGH|CHL1|ENSG00000134121|Transcript|ENST00000421198.5|
    // protein_coding|3/5||||258|5|2|E/X|gAg/gg||A|A/-||1|cds_end_NF|HGNC|HGNC:1939||||4|||||||||||||
    let expected: Vec<&str> = vec![
        "-",                    // 0:  Allele
        "frameshift_variant",   // 1:  Consequence
        "HIGH",                 // 2:  IMPACT
        "CHL1",                 // 3:  SYMBOL
        "ENSG00000134121",      // 4:  Gene
        "Transcript",           // 5:  Feature_type
        "ENST00000421198.5",    // 6:  Feature
        "protein_coding",       // 7:  BIOTYPE
        "3/5",                  // 8:  EXON
        "",                     // 9:  INTRON
        "",                     // 10: HGVSc
        "",                     // 11: HGVSp
        "258",                  // 12: cDNA_position
        "5",                    // 13: CDS_position
        "2",                    // 14: Protein_position
        "E/X",                  // 15: Amino_acids
        "gAg/gg",               // 16: Codons
        "",                     // 17: Existing_variation
        "A",                    // 18: REF_ALLELE
        "A/-",                  // 19: UPLOADED_ALLELE
        "",                     // 20: DISTANCE
        "1",                    // 21: STRAND
        "cds_end_NF",           // 22: FLAGS
        "",                     // 23: CANONICAL
        "HGNC",                 // 24: SYMBOL_SOURCE
        "HGNC:1939",            // 25: HGNC_ID
        "",                     // 26: MANE
        "",                     // 27: MANE_SELECT
        "",                     // 28: MANE_PLUS_CLINICAL
        "4",                    // 29: TSL
        "",                     // 30: APPRIS
        "",                     // 31: CCDS
        "",                     // 32: ENSP
        "",                     // 33: SOURCE
        "",                     // 34: HGVS_OFFSET
        "",                     // 35: SIFT
        "",                     // 36: PolyPhen
        "",                     // 37: AF
        "",                     // 38: CLIN_SIG
        "",                     // 39: SOMATIC
        "",                     // 40: PHENO
        "",                     // 41: PUBMED
        "",                     // 42: MOTIF_NAME
        "",                     // 43: MOTIF_POS
        "",                     // 44: HIGH_INF_POS
        "",                     // 45: MOTIF_SCORE_CHANGE
        "",                     // 46: TRANSCRIPTION_FACTORS
    ];

    assert_eq!(fields.len(), expected.len(),
        "Field count mismatch: got {}, expected {}", fields.len(), expected.len());

    for (i, (got, exp)) in fields.iter().zip(expected.iter()).enumerate() {
        assert_eq!(got, exp,
            "Field {} ({}) mismatch: got {:?}, expected {:?}",
            i, output::DEFAULT_CSQ_FIELDS[i], got, exp);
    }
}

#[test]
fn test_csq_downstream_variant_match() {
    use oxivep_core::Consequence;

    // From VEP: downstream_gene_variant on ENST00000492842.2
    let vf = VariationFeature {
        position: oxivep_core::GenomicPosition::new("1", 65568, 65568, Strand::Forward),
        allele_string: "A/C".into(),
        ref_allele: Allele::from_str("A"),
        alt_alleles: vec![Allele::from_str("C")],
        variation_name: Some("1_65568_A/C".into()),
        vcf_line: None,
        vcf_fields: None,
        transcript_variations: vec![TranscriptVariation {
            transcript_id: "ENST00000492842.2".into(),
            gene_id: "ENSG00000240361".into(),
            gene_symbol: Some("OR4G11P".into()),
            biotype: "transcribed_unprocessed_pseudogene".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::from_str("C"),
                consequences: vec![Consequence::DownstreamGeneVariant],
                impact: Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: Some(1681),
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:31276".into()),
            flags: vec![],
        }],
        existing_variants: vec![],
        minimised: false,
        most_severe_consequence: Some(Consequence::DownstreamGeneVariant),
    };

    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let fields: Vec<&str> = csq.split('|').collect();

    // VEP expected:
    // C|downstream_gene_variant|MODIFIER|OR4G11P|ENSG00000240361|Transcript|ENST00000492842.2|
    // transcribed_unprocessed_pseudogene|||||||||||A|A/C|1681|1|||HGNC|HGNC:31276|...
    assert_eq!(fields.len(), 47);
    assert_eq!(fields[0], "C");
    assert_eq!(fields[1], "downstream_gene_variant");
    assert_eq!(fields[2], "MODIFIER");
    assert_eq!(fields[3], "OR4G11P");
    assert_eq!(fields[7], "transcribed_unprocessed_pseudogene");
    assert_eq!(fields[18], "A");
    assert_eq!(fields[19], "A/C");
    assert_eq!(fields[20], "1681");
    assert_eq!(fields[21], "1");
    assert_eq!(fields[24], "HGNC");
    assert_eq!(fields[25], "HGNC:31276");
    // All annotation fields should be empty for downstream variant
    assert_eq!(fields[8], "");   // EXON
    assert_eq!(fields[12], "");  // cDNA_position
    assert_eq!(fields[13], "");  // CDS_position
    assert_eq!(fields[15], "");  // Amino_acids
    assert_eq!(fields[16], "");  // Codons
}

#[test]
fn test_csq_intron_variant_match() {
    use oxivep_core::Consequence;

    // From VEP: intron_variant on ENST00000272065.10
    let vf = VariationFeature {
        position: oxivep_core::GenomicPosition::new("2", 265023, 265023, Strand::Forward),
        allele_string: "C/T".into(),
        ref_allele: Allele::from_str("C"),
        alt_alleles: vec![Allele::from_str("T")],
        variation_name: Some("2_265023_C/T".into()),
        vcf_line: None,
        vcf_fields: None,
        transcript_variations: vec![TranscriptVariation {
            transcript_id: "ENST00000272065.10".into(),
            gene_id: "ENSG00000143727".into(),
            gene_symbol: Some("ACP1".into()),
            biotype: "protein_coding".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::from_str("T"),
                consequences: vec![Consequence::IntronVariant],
                impact: Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: Some((1, 5)),
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: Some("NM_004300.4".into()),
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: Some("P3".into()),
            ccds: None,
            symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:122".into()),
            flags: vec![],
        }],
        existing_variants: vec![],
        minimised: false,
        most_severe_consequence: Some(Consequence::IntronVariant),
    };

    let csq = output::format_csq(&vf, output::DEFAULT_CSQ_FIELDS);
    let fields: Vec<&str> = csq.split('|').collect();

    // VEP expected:
    // T|intron_variant|MODIFIER|ACP1|ENSG00000143727|Transcript|ENST00000272065.10|
    // protein_coding||1/5|||||||||C|C/T||1|||HGNC|HGNC:122|MANE_Select|NM_004300.4||1|P3|...
    assert_eq!(fields.len(), 47);
    assert_eq!(fields[0], "T");
    assert_eq!(fields[1], "intron_variant");
    assert_eq!(fields[2], "MODIFIER");
    assert_eq!(fields[3], "ACP1");
    assert_eq!(fields[7], "protein_coding");
    assert_eq!(fields[8], "");     // EXON empty for intron
    assert_eq!(fields[9], "1/5");  // INTRON
    assert_eq!(fields[18], "C");
    assert_eq!(fields[19], "C/T");
    assert_eq!(fields[21], "1");
    assert_eq!(fields[24], "HGNC");     // SYMBOL_SOURCE (shifted +1)
    assert_eq!(fields[25], "HGNC:122");
    assert_eq!(fields[26], "MANE_Select");
    assert_eq!(fields[27], "NM_004300.4");
    assert_eq!(fields[29], "1");   // TSL (shifted +1)
    assert_eq!(fields[30], "P3");  // APPRIS (shifted +1)
}
