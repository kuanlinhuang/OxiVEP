use anyhow::{Context, Result};
use oxivep_cache::fasta::FastaReader;
use oxivep_cache::gff::parse_gff3;
use oxivep_cache::info::CacheInfo;
use oxivep_cache::providers::{
    FastaSequenceProvider, MemoryTranscriptProvider, SequenceProvider, TabixVariationProvider,
    TranscriptProvider, VariationProvider,
};
use oxivep_consequence::ConsequencePredictor;
use oxivep_core::Consequence;
use oxivep_hgvs;
use oxivep_io::output;
use oxivep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use oxivep_io::vcf::VcfParser;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

pub struct AnnotateConfig {
    pub input: String,
    pub output: String,
    pub gff3: Option<String>,
    pub fasta: Option<String>,
    pub output_format: String,
    pub pick: bool,
    pub hgvs: bool,
    pub distance: u64,
    pub cache_dir: Option<String>,
}

pub fn run_annotate(config: AnnotateConfig) -> Result<()> {
    // Load transcript models from GFF3
    let mut transcripts = if let Some(ref gff3_path) = config.gff3 {
        let gff_file = File::open(gff3_path)
            .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
        let trs = parse_gff3(gff_file)?;
        eprintln!("Loaded {} transcripts from {}", trs.len(), gff3_path);
        trs
    } else {
        eprintln!("Warning: No GFF3 file provided. Only intergenic variants will be annotated.");
        Vec::new()
    };

    // Load FASTA reference
    let seq_provider: Option<FastaSequenceProvider> = if let Some(ref fasta_path) = config.fasta {
        let fasta_file = File::open(fasta_path)
            .with_context(|| format!("Opening FASTA file: {}", fasta_path))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        eprintln!("Loaded reference FASTA from {}", fasta_path);
        Some(FastaSequenceProvider::new(reader))
    } else {
        None
    };

    // Build sequences for coding transcripts from FASTA
    if let Some(ref sp) = seq_provider {
        let mut built = 0usize;
        for tr in &mut transcripts {
            if tr.is_coding() {
                if let Err(e) = tr.build_sequences(|chrom, start, end| {
                    sp.fetch_sequence(chrom, start, end)
                        .map_err(|e| e.to_string())
                }) {
                    eprintln!("Warning: could not build sequences for {}: {}", tr.stable_id, e);
                } else {
                    built += 1;
                }
            }
        }
        eprintln!("Built sequences for {} coding transcripts", built);
    }

    let transcript_provider = MemoryTranscriptProvider::new(transcripts);

    // Initialize variation provider from VEP cache if provided
    let var_provider: Option<TabixVariationProvider> = if let Some(ref dir) = config.cache_dir {
        let info_path = Path::new(dir).join("info.txt");
        let cache_info = CacheInfo::from_file(&info_path)
            .with_context(|| format!("Reading cache info: {}", info_path.display()))?;
        eprintln!(
            "Loaded VEP cache info: species={}, assembly={}, {} variation columns",
            cache_info.species, cache_info.assembly, cache_info.variation_cols.len()
        );
        Some(TabixVariationProvider::new(Path::new(dir), &cache_info)?)
    } else {
        None
    };

    // Create consequence predictor
    let predictor = ConsequencePredictor::new(config.distance, config.distance);

    // Open input VCF
    let input_reader: Box<dyn io::Read> = if config.input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(&config.input)
                .with_context(|| format!("Opening input file: {}", config.input))?,
        )
    };
    let mut vcf_parser = VcfParser::new(input_reader)?;

    // Open output
    let output_writer: Box<dyn io::Write> = if config.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(&config.output)
                .with_context(|| format!("Creating output file: {}", config.output))?,
        )
    };
    let mut writer = BufWriter::new(output_writer);

    // Write headers based on output format
    match config.output_format.as_str() {
        "vcf" => {
            // Pass through original VCF headers
            for header_line in vcf_parser.header_lines() {
                if header_line.starts_with("#CHROM") {
                    // Insert CSQ header before #CHROM
                    writeln!(writer, "{}", output::csq_header_line(output::DEFAULT_CSQ_FIELDS))?;
                }
                writeln!(writer, "{}", header_line)?;
            }
        }
        "tab" => {
            writeln!(
                writer,
                "## OxiVEP output"
            )?;
            writeln!(
                writer,
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"
            )?;
        }
        "json" => {
            writeln!(writer, "[")?;
        }
        _ => {}
    }

    // Process variants
    let mut count = 0u64;
    let mut first_json = true;

    while let Some(mut vf) = vcf_parser.next_variant()? {
        // Look up known variants from variation cache
        let matched_by_allele: std::collections::HashMap<String, Vec<oxivep_cache::providers::MatchedVariant>> =
            if let Some(ref vp) = var_provider {
                let mut by_allele = std::collections::HashMap::new();
                for alt in &vf.alt_alleles {
                    let alt_str = alt.to_string();
                    let ref_str = vf.ref_allele.to_string();
                    let matches = vp.get_matched_variants(
                        &vf.position.chromosome,
                        vf.position.start,
                        vf.position.end,
                        &ref_str,
                        &alt_str,
                    ).unwrap_or_default();
                    if !matches.is_empty() {
                        by_allele.insert(alt_str, matches);
                    }
                }
                by_allele
            } else {
                std::collections::HashMap::new()
            };

        // Populate existing_variants on the VF for output access
        for matches in matched_by_allele.values() {
            for m in matches {
                if !vf.existing_variants.iter().any(|kv| kv.name == m.name) {
                    vf.existing_variants.push(oxivep_io::variant::KnownVariant {
                        name: m.name.clone(),
                        allele_string: None,
                        minor_allele: m.minor_allele.clone(),
                        minor_allele_freq: m.minor_allele_freq,
                        clinical_significance: m.clin_sig.clone(),
                        somatic: m.somatic,
                        phenotype_or_disease: m.phenotype_or_disease,
                        pubmed: m.pubmed.clone(),
                        frequencies: m.frequencies.clone(),
                    });
                }
            }
        }

        // Find overlapping transcripts
        let chrom = &vf.position.chromosome;
        let query_start = if vf.position.start > config.distance {
            vf.position.start - config.distance
        } else {
            1
        };
        let query_end = vf.position.end + config.distance;
        let overlapping = transcript_provider.get_transcripts(chrom, query_start, query_end)?;

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(&mut vf);
            // Populate existing_variation on intergenic annotations too
            for tv in &mut vf.transcript_variations {
                for aa in &mut tv.allele_annotations {
                    if let Some(matches) = matched_by_allele.get(&aa.allele.to_string()) {
                        aa.existing_variation = matches.iter().map(|m| m.name.clone()).collect();
                    }
                }
            }
        } else {
            // Get reference sequence if available
            let ref_seq = seq_provider.as_ref().and_then(|sp| {
                sp.fetch_sequence(chrom, query_start, query_end).ok()
            });

            // Run consequence prediction
            let result = predictor.predict(
                &vf.position,
                &vf.ref_allele,
                &vf.alt_alleles,
                &overlapping,
                ref_seq.as_deref(),
            );

            // Convert prediction results to VariationFeature annotations
            for tc in &result.transcript_consequences {
                let transcript = overlapping.iter().find(|t| t.stable_id == tc.transcript_id);

                let allele_annotations: Vec<AlleleAnnotation> = tc
                    .allele_consequences
                    .iter()
                    .map(|ac| {
                        let mut ann = AlleleAnnotation {
                            allele: ac.allele.clone(),
                            consequences: ac.consequences.clone(),
                            impact: ac.impact,
                            cdna_position: zip_positions(ac.cdna_start, ac.cdna_end),
                            cds_position: zip_positions(ac.cds_start, ac.cds_end),
                            protein_position: zip_positions(ac.protein_start, ac.protein_end),
                            amino_acids: ac.amino_acids.clone(),
                            codons: ac.codons.clone(),
                            exon: ac.exon,
                            intron: ac.intron,
                            distance: ac.distance,
                            hgvsc: None,
                            hgvsp: None,
                            hgvsg: None,
                            existing_variation: matched_by_allele
                                .get(&ac.allele.to_string())
                                .map(|matches| matches.iter().map(|m| m.name.clone()).collect())
                                .unwrap_or_default(),
                            sift: None,
                            polyphen: None,
                        };

                        // Generate HGVS if requested
                        if config.hgvs {
                            ann.hgvsg = Some(oxivep_hgvs::hgvsg(
                                chrom,
                                vf.position.start,
                                vf.position.end,
                                &vf.ref_allele,
                                &ac.allele,
                            ));

                            if let Some(tr) = transcript {
                                // Build versioned IDs for HGVS notation
                                let versioned_tid = match tr.version {
                                    Some(v) => format!("{}.{}", tc.transcript_id, v),
                                    None => tc.transcript_id.clone(),
                                };

                                // Determine alleles for HGVS - complement for minus strand
                                let (hgvs_ref, hgvs_alt) = if tr.strand == oxivep_core::Strand::Reverse {
                                    (complement_allele(&vf.ref_allele), complement_allele(&ac.allele))
                                } else {
                                    (vf.ref_allele.clone(), ac.allele.clone())
                                };

                                if let Some(coding_start) = tr.cdna_coding_start {
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        // Normalize cDNA positions (minus-strand can reverse order)
                                        let (cs, ce) = (cs.min(ce), cs.max(ce));
                                        // Exonic variant: standard HGVSc with 3' shifting
                                        ann.hgvsc = oxivep_hgvs::hgvsc_with_seq(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                            coding_start,
                                            tr.cdna_coding_end,
                                            tr.spliced_seq.as_deref(),
                                            tr.codon_table_start_phase,
                                        );
                                    } else if ac.intron.is_some() {
                                        // Intronic variant: offset notation
                                        // Note: intronic HGVS uses original coding_start (no phase adjustment)
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(vf.position.start) {
                                            // For multi-base variants, compute end position too
                                            let (end_cdna, end_offset) = if vf.position.start != vf.position.end {
                                                tr.genomic_to_intronic_cdna(vf.position.end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = oxivep_hgvs::hgvsc_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                            );
                                            // For intronic insertions, check if it's a dup.
                                            // A dup occurs when the inserted bases match the
                                            // immediately preceding OR following reference sequence.
                                            // We use the original (non-complemented) alt allele for
                                            // comparison since we're checking against genomic ref.
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    // Get the original (genomic-strand) insertion bases
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            // Check if insertion duplicates the preceding bases
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                                                            } else { false };
                                                            // Check if insertion duplicates the following bases
                                                            let dup_after = if !dup_before {
                                                                let check_start = vf.position.start;
                                                                let check_end = check_start + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence(chrom, check_start, check_end) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                hgvsc = convert_ins_to_dup(h, offset, ins_len, cdna_pos, coding_start, tr.cdna_coding_end);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                } else {
                                    // Non-coding transcript: use n. notation
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        ann.hgvsc = oxivep_hgvs::hgvsc_noncoding(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                        );
                                    } else if ac.intron.is_some() {
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(vf.position.start) {
                                            // For multi-base variants, compute end position
                                            let (end_cdna, end_offset) = if vf.position.start != vf.position.end {
                                                tr.genomic_to_intronic_cdna(vf.position.end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = oxivep_hgvs::hgvsc_noncoding_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &hgvs_alt,
                                            );
                                            // Dup detection for non-coding intronic insertions
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                                                            } else { false };
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.to_ascii_uppercase() == b.to_ascii_uppercase())
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                hgvsc = convert_ins_to_dup_noncoding(h, offset, ins_len, cdna_pos);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                }
                            }

                            if let (Some(ref aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
                                if let Some(tr) = transcript {
                                    if let Some(ref pid) = tr.protein_id {
                                        let versioned_pid = match tr.protein_version {
                                            Some(v) => {
                                                let suffix = format!(".{}", v);
                                                if pid.ends_with(&suffix) {
                                                    pid.clone()
                                                } else {
                                                    format!("{}.{}", pid, v)
                                                }
                                            }
                                            None => pid.clone(),
                                        };
                                        let is_fs = ac.consequences.contains(&Consequence::FrameshiftVariant);

                                        if is_fs {
                                            // Frameshift: build alt sequence and scan for first changed AA + new stop
                                            // Use spliced_seq from CDS start onwards (includes 3'UTR for stop codon search)
                                            if let (Some(ref spliced), Some(coding_start), Some(cds_s)) =
                                                (&tr.spliced_seq, tr.cdna_coding_start, ac.cds_start)
                                            {
                                                // Extract from CDS start to end of spliced seq (includes 3'UTR)
                                                let coding_start_idx = (coding_start - 1) as usize;
                                                let ref_from_cds = &spliced.as_bytes()[coding_start_idx..];
                                                let cds_idx = (cds_s - 1) as usize;
                                                let mut alt_from_cds = ref_from_cds.to_vec();

                                                // Apply the indel to build the frameshifted sequence
                                                if ac.allele == Allele::Deletion {
                                                    let del_len = vf.ref_allele.len();
                                                    let end = (cds_idx + del_len).min(alt_from_cds.len());
                                                    alt_from_cds.drain(cds_idx..end);
                                                } else if let Allele::Sequence(ins_bases) = &ac.allele {
                                                    let mut bases = ins_bases.clone();
                                                    if tr.strand == oxivep_core::Strand::Reverse {
                                                        bases = bases.iter().map(|&b| match b {
                                                            b'A' => b'T', b'T' => b'A',
                                                            b'C' => b'G', b'G' => b'C',
                                                            o => o,
                                                        }).collect();
                                                    }
                                                    for (j, &b) in bases.iter().enumerate() {
                                                        if cds_idx + j <= alt_from_cds.len() {
                                                            alt_from_cds.insert(cds_idx + j, b);
                                                        }
                                                    }
                                                }

                                                let codon_start = cds_idx / 3;
                                                ann.hgvsp = oxivep_hgvs::hgvsp_frameshift(
                                                    &versioned_pid,
                                                    ref_from_cds,
                                                    &alt_from_cds,
                                                    codon_start,
                                                );
                                            }
                                        } else {
                                            let ref_aa_byte = aa.0.as_bytes().first().copied().unwrap_or(b'X');
                                            let alt_aa_byte = aa.1.as_bytes().first().copied().unwrap_or(b'X');
                                            ann.hgvsp = oxivep_hgvs::hgvsp(
                                                &versioned_pid, ps, ref_aa_byte, alt_aa_byte, false,
                                            );
                                        }
                                    }
                                }
                            }
                        }

                        ann
                    })
                    .collect();

                // Apply pick filter if needed
                let should_include = !config.pick || tc.canonical || vf.transcript_variations.is_empty();

                if should_include {
                    vf.transcript_variations.push(TranscriptVariation {
                        transcript_id: tc.transcript_id.clone(),
                        gene_id: tc.gene_id.clone(),
                        gene_symbol: tc.gene_symbol.clone(),
                        biotype: tc.biotype.clone(),
                        allele_annotations,
                        canonical: tc.canonical,
                        strand: tc.strand,
                        source: None,
                        protein_id: transcript.and_then(|t| t.protein_id.clone()),
                        mane_select: transcript.and_then(|t| t.mane_select.clone()),
                        mane_plus_clinical: transcript.and_then(|t| t.mane_plus_clinical.clone()),
                        tsl: transcript.and_then(|t| t.tsl),
                        appris: transcript.and_then(|t| t.appris.clone()),
                        ccds: transcript.and_then(|t| t.ccds.clone()),
                        symbol_source: transcript.and_then(|t| t.gene.symbol_source.clone()),
                        hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                        flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                    });
                }
            }

            vf.compute_most_severe();
        }

        // Write output
        match config.output_format.as_str() {
            "vcf" => write_vcf_line(&mut writer, &vf)?,
            "tab" => {
                for line in output::format_tab_line(&vf) {
                    writeln!(writer, "{}", line)?;
                }
            }
            "json" => {
                if !first_json {
                    writeln!(writer, ",")?;
                }
                first_json = false;
                let json = output::format_json(&vf);
                write!(writer, "{}", serde_json::to_string_pretty(&json)?)?;
            }
            _ => {}
        }

        count += 1;
    }

    // Close JSON array
    if config.output_format == "json" {
        writeln!(writer, "\n]")?;
    }

    writer.flush()?;
    eprintln!("Annotated {} variants", count);

    Ok(())
}

/// Convert an intronic insertion HGVSc notation to duplication notation for coding transcripts.
///
/// For single-base dups: `c.X+N_X+N+1insB` → `c.X+Ndup`
/// For multi-base dups: `c.X+N_X+N+1insBBB` → `c.X+(N-len+1)_X+Ndup`
fn convert_ins_to_dup(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    // Find the prefix (e.g., "ENST00000334363.14:c." or "ENST...:n.")
    // The prefix ends at "c." or "n." — everything up to and including that
    let prefix_end = hgvsc.find(":c.").map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

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

    let end_pos = build_pos(nearest_exon_cdna_pos, intron_offset);
    if ins_len == 1 {
        Some(format!("{}{}dup", prefix, end_pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert an intronic insertion HGVSc notation to duplication notation for non-coding transcripts.
fn convert_ins_to_dup_noncoding(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
) -> Option<String> {
    let prefix_end = hgvsc.find(":n.").map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |off: i64| -> String {
        if off > 0 {
            format!("{}+{}", nearest_exon_cdna_pos, off)
        } else {
            format!("{}{}", nearest_exon_cdna_pos, off)
        }
    };

    let end_pos = build_pos(intron_offset);
    if ins_len == 1 {
        Some(format!("{}{}dup", prefix, end_pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(start_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

fn annotate_intergenic(vf: &mut VariationFeature) {
    for alt in &vf.alt_alleles {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![Consequence::IntergenicVariant],
                impact: oxivep_core::Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
            }],
            canonical: false,
            strand: oxivep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        });
    }
    vf.most_severe_consequence = Some(Consequence::IntergenicVariant);
}

fn write_vcf_line(writer: &mut impl Write, vf: &VariationFeature) -> Result<()> {
    if let Some(ref fields) = vf.vcf_fields {
        let csq = output::format_csq(vf, output::DEFAULT_CSQ_FIELDS);
        let info = if fields.info == "." && !csq.is_empty() {
            format!("CSQ={}", csq)
        } else if !csq.is_empty() {
            format!("{};CSQ={}", fields.info, csq)
        } else {
            fields.info.clone()
        };

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fields.chrom, fields.pos, fields.id, fields.ref_allele, fields.alt,
            fields.qual, fields.filter, info
        )?;

        for rest_field in &fields.rest {
            write!(writer, "\t{}", rest_field)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn zip_positions(start: Option<u64>, end: Option<u64>) -> Option<(u64, u64)> {
    match (start, end) {
        (Some(s), Some(e)) => {
            // Normalize order: smaller position first (minus-strand can reverse start/end)
            Some((s.min(e), s.max(e)))
        }
        (Some(s), None) => Some((s, s)),
        (None, Some(e)) => Some((e, e)),
        _ => None,
    }
}

use oxivep_core::Allele;

/// Complement an allele for minus strand HGVS notation.
fn complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => {
            let comp: Vec<u8> = bases.iter().map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                other => other,
            }).collect();
            Allele::Sequence(comp)
        }
        other => other.clone(),
    }
}

use serde_json;
