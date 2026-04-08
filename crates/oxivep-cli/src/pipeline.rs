use anyhow::{Context, Result};
use oxivep_cache::annotation::{AnnotationProvider, AnnotationValue};
use oxivep_cache::fasta::FastaReader;
use oxivep_cache::gff::parse_gff3;
use oxivep_cache::info::CacheInfo;
use oxivep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MatchedVariant, SequenceProvider,
    TabixVariationProvider, TranscriptProvider, VariationProvider,
};
use oxivep_consequence::ConsequencePredictor;
use oxivep_core::Consequence;
use oxivep_hgvs;
use oxivep_io::output;
use oxivep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use oxivep_io::vcf::VcfParser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

const BATCH_SIZE: usize = 1024;

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
    pub transcript_cache: Option<String>,
    /// Directory containing supplementary annotation files (.osa, .osi, .oga).
    pub sa_dir: Option<String>,
}

pub fn run_annotate(config: AnnotateConfig) -> Result<()> {
    // Extract the GFF3 source name (filename) for the SOURCE field
    let gff3_source: Option<String> = config.gff3.as_ref().map(|p| {
        Path::new(p).file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| p.clone())
    });

    // Load transcript models: try binary cache first, fall back to GFF3
    let cache_path = config.transcript_cache.as_ref().map(|p| Path::new(p).to_path_buf())
        .or_else(|| config.gff3.as_ref().map(|p| oxivep_cache::transcript_cache::default_cache_path(Path::new(p))));

    let mut transcripts = 'load: {
        // Try loading from binary cache
        if let Some(ref cp) = cache_path {
            if cp.exists() {
                let is_fresh = config.gff3.as_ref()
                    .map(|gff| oxivep_cache::transcript_cache::cache_is_fresh(cp, Path::new(gff)))
                    .unwrap_or(true);
                if is_fresh {
                    match oxivep_cache::transcript_cache::load_cache(cp) {
                        Ok(trs) => {
                            eprintln!("Loaded {} transcripts from cache {}", trs.len(), cp.display());
                            break 'load trs;
                        }
                        Err(e) => {
                            eprintln!("Warning: cache load failed ({}), falling back to GFF3", e);
                        }
                    }
                } else {
                    eprintln!("Cache is stale, rebuilding from GFF3");
                }
            }
        }

        // Fall back to GFF3 parsing
        if let Some(ref gff3_path) = config.gff3 {
            let gff_path = Path::new(gff3_path);
            let tbi_path = format!("{}.tbi", gff3_path);

            // Use indexed loading if .gff3.gz + .tbi available
            if gff3_path.ends_with(".gz") && Path::new(&tbi_path).exists() {
                // Pre-scan VCF to collect regions
                let regions = prescan_vcf_regions(&config.input, config.distance)?;
                eprintln!("Pre-scanned {} variant regions from {}", regions.len(), config.input);
                let trs = oxivep_cache::gff::parse_gff3_indexed(gff_path, &regions)?;
                eprintln!("Loaded {} transcripts from indexed {}", trs.len(), gff3_path);
                trs
            } else {
                let gff_file = File::open(gff3_path)
                    .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                let trs = parse_gff3(gff_file)?;
                eprintln!("Loaded {} transcripts from {}", trs.len(), gff3_path);
                trs
            }
        } else {
            eprintln!("Warning: No GFF3 file provided. Only intergenic variants will be annotated.");
            Vec::new()
        }
    };

    // Load FASTA reference (prefer mmap with .fai index, fall back to in-memory)
    let seq_provider: Option<Box<dyn SequenceProvider>> = if let Some(ref fasta_path) = config.fasta {
        let fai_path = format!("{}.fai", fasta_path);
        if Path::new(&fai_path).exists() {
            let reader = oxivep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
            eprintln!("Memory-mapped reference FASTA from {} (using .fai index)", fasta_path);
            Some(Box::new(oxivep_cache::providers::MmapFastaSequenceProvider::new(reader)))
        } else {
            let fasta_file = File::open(fasta_path)
                .with_context(|| format!("Opening FASTA file: {}", fasta_path))?;
            let reader = FastaReader::from_reader(fasta_file)?;
            eprintln!("Loaded reference FASTA from {}", fasta_path);
            Some(Box::new(FastaSequenceProvider::new(reader)))
        }
    } else {
        None
    };

    // Build sequences for coding transcripts from FASTA (skip if loaded from cache with sequences)
    let needs_seq_build = transcripts.iter().any(|t| t.is_coding() && t.spliced_seq.is_none());
    if needs_seq_build {
        if let Some(ref sp) = seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.is_coding() && tr.spliced_seq.is_none() {
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
    }

    // Save cache if we parsed from GFF3 and have a cache path
    if let Some(ref cp) = cache_path {
        if !cp.exists() || !config.transcript_cache.is_some() {
            // Auto-save cache after GFF3 parse + sequence build
            if config.gff3.is_some() {
                if let Err(e) = oxivep_cache::transcript_cache::save_cache(&transcripts, cp) {
                    eprintln!("Warning: could not save cache: {}", e);
                } else {
                    eprintln!("Saved transcript cache to {}", cp.display());
                }
            }
        }
    }

    let transcript_provider = IndexedTranscriptProvider::new(transcripts);

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

    // Load supplementary annotation providers from --sa-dir
    let sa_providers: Vec<Box<dyn AnnotationProvider>> = if let Some(ref dir) = config.sa_dir {
        load_sa_providers(Path::new(dir))?
    } else {
        Vec::new()
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

    // Process variants in batches for parallel annotation
    let mut count = 0u64;
    let mut first_json = true;

    loop {
        // Phase 1: Read a batch of variants (sequential - VCF parser is not Sync)
        let mut batch: Vec<(VariationFeature, HashMap<String, Vec<MatchedVariant>>)> = Vec::with_capacity(BATCH_SIZE);
        for _ in 0..BATCH_SIZE {
            match vcf_parser.next_variant()? {
                Some(mut vf) => {
                    // Variation lookup (sequential - TabixVariationProvider is not Sync)
                    let matched_by_allele: HashMap<String, Vec<MatchedVariant>> =
                        if let Some(ref vp) = var_provider {
                            let mut by_allele = HashMap::new();
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
                            HashMap::new()
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

                    batch.push((vf, matched_by_allele));
                }
                None => break,
            }
        }

        if batch.is_empty() {
            break;
        }

        // Phase 1.5: Preload SA providers for this batch (sequential)
        if !sa_providers.is_empty() && !batch.is_empty() {
            // Collect all positions in this batch, grouped by chromosome
            let mut chrom_positions: HashMap<&str, Vec<u64>> = HashMap::new();
            for (vf, _) in &batch {
                chrom_positions
                    .entry(&vf.position.chromosome)
                    .or_default()
                    .push(vf.position.start);
            }
            for sa in &sa_providers {
                for (chrom, positions) in &chrom_positions {
                    let _ = sa.preload(chrom, positions);
                }
            }
        }

        // Phase 2: Annotate batch in parallel (transcript lookup + consequence prediction + HGVS)
        batch.par_iter_mut().for_each(|(vf, matched_by_allele)| {
            let chrom = &vf.position.chromosome;
            let query_start = if vf.position.start > config.distance {
                vf.position.start - config.distance
            } else {
                1
            };
            let query_end = vf.position.end + config.distance;
            let overlapping = transcript_provider.get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(vf);
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

            // Run consequence prediction — dispatch SVs to SV predictor
            let transcript_consequences = if vf.variant_type.is_structural() {
                oxivep_consequence::sv_predictor::predict_sv_consequences(
                    chrom,
                    vf.position.start,
                    vf.position.end,
                    vf.variant_type,
                    &vf.alt_alleles,
                    &overlapping,
                    config.distance,
                    config.distance,
                )
            } else {
                let result = predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );
                result.transcript_consequences
            };

            // Convert prediction results to VariationFeature annotations
            for tc in &transcript_consequences {
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
                            hgvs_offset: None,
                            existing_variation: matched_by_allele
                                .get(&ac.allele.to_string())
                                .map(|matches| matches.iter().map(|m| m.name.clone()).collect())
                                .unwrap_or_default(),
                            sift: None,
                            polyphen: None,
                            supplementary: Vec::new(),
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
                                    None => tc.transcript_id.to_string(),
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
                                        // Apply HGVS 3' normalization for intronic indels
                                        let (shifted_start, shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    // Use genomic-strand alleles for ref comparison
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };
                                        // For insertions, build the rotated insertion bases
                                        // after 3' shifting (bases rotate as position shifts)
                                        let shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                // Calculate how many positions we shifted
                                                let shift_amount = if tr.strand == oxivep_core::Strand::Forward {
                                                    (shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - shifted_start as i64) as usize
                                                };
                                                // Rotate: for forward strand, each shift moves first base to end
                                                // For reverse strand, each shift moves last base to front
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        oxivep_core::Strand::Forward => {
                                                            rotated.rotate_left(effective_shift);
                                                        }
                                                        oxivep_core::Strand::Reverse => {
                                                            rotated.rotate_right(effective_shift);
                                                        }
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        // For insertions, use position before insertion
                                        // for the primary HGVS coordinate (ins is BETWEEN two bases).
                                        // On reverse strand, the insertion is between P and P+1 in
                                        // genomic coords, but P+1 is 5' in transcript order, so we
                                        // use P+1 as the HGVS start coordinate.
                                        let is_insertion = matches!((&hgvs_ref, &shifted_hgvs_alt), (Allele::Deletion, Allele::Sequence(_)));
                                        let hgvs_pos = if is_insertion {
                                            if tr.strand == oxivep_core::Strand::Reverse {
                                                shifted_end + 1
                                            } else {
                                                shifted_end // base before insertion
                                            }
                                        } else {
                                            shifted_start
                                        };
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(hgvs_pos) {
                                            // For multi-base variants, compute end position too
                                            let (end_cdna, end_offset) = if shifted_start != shifted_end && hgvs_pos == shifted_start {
                                                tr.genomic_to_intronic_cdna(shifted_end)
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
                                                &shifted_hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                            );
                                            // For intronic insertions, check if it's a dup.
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
                                                            // Check dup_before: base(s) before insertion match
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            // Check dup_after: base(s) after insertion match
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                // For dups, determine the dup base position and 3' shift it
                                                                let dup_base_pos = if dup_before {
                                                                    // Dup base is before insertion: position.end
                                                                    vf.position.end
                                                                } else {
                                                                    // Dup base is after insertion: position.start
                                                                    vf.position.start
                                                                };
                                                                // 3' shift the dup position within the intron
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                // Use shifted_dup (start of dup region) for offset computation
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup(h, dup_offset, ins_len, dup_cdna, coding_start, tr.cdna_coding_end);
                                                                }
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
                                        // Apply 3' normalization for non-coding intronic indels
                                        let (nc_shifted_start, nc_shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };

                                        // Rotate insertion bases for non-coding
                                        let nc_shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if nc_shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                let shift_amount = if tr.strand == oxivep_core::Strand::Forward {
                                                    (nc_shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - nc_shifted_start as i64) as usize
                                                };
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        oxivep_core::Strand::Forward => rotated.rotate_left(effective_shift),
                                                        oxivep_core::Strand::Reverse => rotated.rotate_right(effective_shift),
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(nc_shifted_start) {
                                            let (end_cdna, end_offset) = if nc_shifted_start != nc_shifted_end {
                                                tr.genomic_to_intronic_cdna(nc_shifted_end)
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
                                                &nc_shifted_hgvs_alt,
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
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                let dup_base_pos = if dup_before { vf.position.end } else { vf.position.start };
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup_noncoding(h, dup_offset, ins_len, dup_cdna);
                                                                }
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
                        source: gff3_source.clone(),
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

            // Supplementary annotation: query SA providers for each allele
            if !sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        let alt_str = aa.allele.to_string();
                        let ref_str = vf.ref_allele.to_string();
                        for sa in &sa_providers {
                            if let Ok(Some(ann)) =
                                sa.annotate_position(chrom, vf.position.start, &ref_str, &alt_str)
                            {
                                let json_str = match ann {
                                    AnnotationValue::Json(j) => j,
                                    AnnotationValue::Positional(j) => j,
                                    AnnotationValue::Interval(v) => {
                                        format!("[{}]", v.join(","))
                                    }
                                };
                                aa.supplementary.push((
                                    sa.json_key().to_string(),
                                    json_str,
                                ));
                            }
                        }
                    }
                }
            }

            vf.compute_most_severe();
        }
        }); // end par_iter_mut

        // Phase 3: Write output sequentially (preserves VCF order)
        for (vf, _) in &batch {
            match config.output_format.as_str() {
                "vcf" => write_vcf_line(&mut writer, vf)?,
                "tab" => {
                    for line in output::format_tab_line(vf) {
                        writeln!(writer, "{}", line)?;
                    }
                }
                "json" => {
                    if !first_json {
                        writeln!(writer, ",")?;
                    }
                    first_json = false;
                    let json = output::format_json(vf);
                    write!(writer, "{}", serde_json::to_string_pretty(&json)?)?;
                }
                _ => {}
            }
        }

        count += batch.len() as u64;
    } // end batch loop

    // Close JSON array
    if config.output_format == "json" {
        writeln!(writer, "\n]")?;
    }

    writer.flush()?;
    eprintln!("Annotated {} variants", count);

    Ok(())
}

/// Pre-loaded annotation context for the web server.
/// Holds transcript models and sequence provider so they're loaded once at startup.
pub struct AnnotationContext {
    pub transcript_provider: IndexedTranscriptProvider,
    pub seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>>,
    pub predictor: ConsequencePredictor,
    pub gff3_source: Option<String>,
    pub distance: u64,
    pub hgvs: bool,
}

impl AnnotationContext {
    /// Build a context from GFF3 (and optional FASTA) paths — called once at server startup.
    pub fn new(
        gff3: Option<&str>,
        fasta: Option<&str>,
        distance: u64,
    ) -> Result<Self> {
        let gff3_source: Option<String> = gff3.map(|p| {
            Path::new(p).file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| p.to_string())
        });

        // Load transcripts (same logic as run_annotate)
        let mut transcripts = if let Some(gff3_path) = gff3 {
            let cache_path = oxivep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
            let from_cache = if cache_path.exists() {
                let is_fresh = oxivep_cache::transcript_cache::cache_is_fresh(&cache_path, Path::new(gff3_path));
                if is_fresh {
                    oxivep_cache::transcript_cache::load_cache(&cache_path).ok()
                } else { None }
            } else { None };
            if let Some(trs) = from_cache {
                eprintln!("[web] Loaded {} transcripts from cache", trs.len());
                trs
            } else {
                let tbi_path = format!("{}.tbi", gff3_path);
                let trs = if gff3_path.ends_with(".gz") && Path::new(&tbi_path).exists() {
                    // No region pre-scan for web mode — load all
                    let gff_file = File::open(gff3_path)
                        .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                    parse_gff3(gff_file)?
                } else {
                    let gff_file = File::open(gff3_path)
                        .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                    parse_gff3(gff_file)?
                };
                eprintln!("[web] Loaded {} transcripts from {}", trs.len(), gff3_path);
                // Save cache for next time
                if let Err(e) = oxivep_cache::transcript_cache::save_cache(&trs, &cache_path) {
                    eprintln!("[web] Warning: could not save cache: {}", e);
                }
                trs
            }
        } else {
            Vec::new()
        };

        // Load FASTA
        let seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>> = if let Some(fasta_path) = fasta {
            let fai_path = format!("{}.fai", fasta_path);
            if Path::new(&fai_path).exists() {
                let reader = oxivep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
                eprintln!("[web] Memory-mapped FASTA from {}", fasta_path);
                Some(Box::new(oxivep_cache::providers::MmapFastaSequenceProvider::new(reader)))
            } else {
                let fasta_file = File::open(fasta_path)
                    .with_context(|| format!("Opening FASTA: {}", fasta_path))?;
                let reader = FastaReader::from_reader(fasta_file)?;
                eprintln!("[web] Loaded FASTA from {}", fasta_path);
                Some(Box::new(FastaSequenceProvider::new(reader)))
            }
        } else {
            None
        };

        // Build sequences for coding transcripts
        if let Some(ref sp) = seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if tr.build_sequences(|chrom, start, end| {
                        sp.fetch_sequence(chrom, start, end).map_err(|e| e.to_string())
                    }).is_ok() {
                        built += 1;
                    }
                }
            }
            eprintln!("[web] Built sequences for {} coding transcripts", built);
        }

        let transcript_provider = IndexedTranscriptProvider::new(transcripts);
        let predictor = ConsequencePredictor::new(distance, distance);

        Ok(Self {
            transcript_provider,
            seq_provider,
            predictor,
            gff3_source,
            distance,
            hgvs: true,
        })
    }

    /// Replace the transcript models by parsing GFF3 text uploaded from the browser.
    /// Returns (gene_count, transcript_count).
    pub fn update_gff3_text(&mut self, gff3_text: &str) -> Result<(usize, usize)> {
        let mut transcripts = parse_gff3(gff3_text.as_bytes())?;
        let gene_count = {
            let mut genes = std::collections::HashSet::new();
            for t in &transcripts {
                genes.insert(t.gene.stable_id.clone());
            }
            genes.len()
        };

        // Build sequences if FASTA is available
        if let Some(ref sp) = self.seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if tr.build_sequences(|chrom, start, end| {
                        sp.fetch_sequence(chrom, start, end).map_err(|e| e.to_string())
                    }).is_ok() {
                        built += 1;
                    }
                }
            }
            if built > 0 {
                eprintln!("[web] Built sequences for {} coding transcripts", built);
            }
        }

        let tr_count = transcripts.len();
        self.transcript_provider = IndexedTranscriptProvider::new(transcripts);
        self.gff3_source = Some("user-upload".to_string());
        eprintln!("[web] Updated GFF3: {} genes, {} transcripts", gene_count, tr_count);
        Ok((gene_count, tr_count))
    }

    /// Annotate VCF text and return JSON results.
    pub fn annotate_vcf_text(&self, vcf_text: &str, pick: bool) -> Result<Vec<serde_json::Value>> {
        let mut vcf_parser = VcfParser::new(vcf_text.as_bytes())?;
        let mut variants = vcf_parser.read_all()?;

        // Annotate all variants (single-threaded for simplicity in web context;
        // requests are typically small)
        for vf in &mut variants {
            let chrom = &vf.position.chromosome;
            let query_start = vf.position.start.saturating_sub(self.distance).max(1);
            let query_end = vf.position.end + self.distance;
            let overlapping = self.transcript_provider.get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

            if overlapping.is_empty() {
                annotate_intergenic(vf);
            } else {
                let ref_seq = self.seq_provider.as_ref().and_then(|sp| {
                    sp.fetch_sequence(chrom, query_start, query_end).ok()
                });

                let result = self.predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );

                for tc in &result.transcript_consequences {
                    let transcript = overlapping.iter().find(|t| t.stable_id == tc.transcript_id);

                    let allele_annotations: Vec<AlleleAnnotation> = tc.allele_consequences.iter().map(|ac| {
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
                            hgvs_offset: None,
                            existing_variation: vec![],
                            sift: None,
                            polyphen: None,
                            supplementary: Vec::new(),
                        };

                        if self.hgvs {
                            ann.hgvsg = Some(oxivep_hgvs::hgvsg(
                                chrom, vf.position.start, vf.position.end,
                                &vf.ref_allele, &ac.allele,
                            ));
                            if let Some(tr) = transcript {
                                let versioned_tid = match tr.version {
                                    Some(v) => format!("{}.{}", tc.transcript_id, v),
                                    None => tc.transcript_id.to_string(),
                                };
                                let (hgvs_ref, hgvs_alt) = if tr.strand == oxivep_core::Strand::Reverse {
                                    (complement_allele(&vf.ref_allele), complement_allele(&ac.allele))
                                } else {
                                    (vf.ref_allele.clone(), ac.allele.clone())
                                };
                                if let Some(coding_start) = tr.cdna_coding_start {
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        let (cs, ce) = (cs.min(ce), cs.max(ce));
                                        ann.hgvsc = oxivep_hgvs::hgvsc_with_seq(
                                            &versioned_tid, cs, ce,
                                            &hgvs_ref, &hgvs_alt,
                                            coding_start, tr.cdna_coding_end,
                                            tr.spliced_seq.as_deref(), tr.codon_table_start_phase,
                                        );
                                    } else if ac.intron.is_some() {
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(vf.position.start) {
                                            let (end_cdna, end_offset) = if vf.position.start != vf.position.end {
                                                tr.genomic_to_intronic_cdna(vf.position.end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else { (None, None) };
                                            ann.hgvsc = oxivep_hgvs::hgvsc_intronic_range(
                                                &versioned_tid, cdna_pos, offset,
                                                end_cdna, end_offset,
                                                &hgvs_ref, &hgvs_alt,
                                                coding_start, tr.cdna_coding_end,
                                            );
                                        }
                                    }
                                } else if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                    ann.hgvsc = oxivep_hgvs::hgvsc_noncoding(
                                        &versioned_tid, cs, ce, &hgvs_ref, &hgvs_alt,
                                    );
                                } else if ac.intron.is_some() {
                                    if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(vf.position.start) {
                                        let (end_cdna, end_offset) = if vf.position.start != vf.position.end {
                                            tr.genomic_to_intronic_cdna(vf.position.end)
                                                .map(|(c, o)| (Some(c), Some(o)))
                                                .unwrap_or((None, None))
                                        } else { (None, None) };
                                        ann.hgvsc = oxivep_hgvs::hgvsc_noncoding_intronic_range(
                                            &versioned_tid, cdna_pos, offset,
                                            end_cdna, end_offset,
                                            &hgvs_ref, &hgvs_alt,
                                        );
                                    }
                                }

                                // HGVSp
                                if let (Some(ref aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
                                    if let Some(ref pid) = tr.protein_id {
                                        let versioned_pid = match tr.protein_version {
                                            Some(v) => {
                                                let suffix = format!(".{}", v);
                                                if pid.ends_with(&suffix) { pid.clone() } else { format!("{}.{}", pid, v) }
                                            }
                                            None => pid.clone(),
                                        };
                                        let is_fs = ac.consequences.contains(&Consequence::FrameshiftVariant);
                                        if is_fs {
                                            if let (Some(ref spliced), Some(coding_start), Some(cds_s)) =
                                                (&tr.spliced_seq, tr.cdna_coding_start, ac.cds_start)
                                            {
                                                let coding_start_idx = (coding_start - 1) as usize;
                                                let ref_from_cds = &spliced.as_bytes()[coding_start_idx..];
                                                let cds_idx = (cds_s - 1) as usize;
                                                let mut alt_from_cds = ref_from_cds.to_vec();
                                                if ac.allele == oxivep_core::Allele::Deletion {
                                                    let del_len = vf.ref_allele.len();
                                                    let end = (cds_idx + del_len).min(alt_from_cds.len());
                                                    alt_from_cds.drain(cds_idx..end);
                                                } else if let oxivep_core::Allele::Sequence(ins_bases) = &ac.allele {
                                                    let mut bases = ins_bases.clone();
                                                    if tr.strand == oxivep_core::Strand::Reverse {
                                                        bases = bases.iter().map(|&b| match b {
                                                            b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', o => o,
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
                                                    &versioned_pid, ref_from_cds, &alt_from_cds, codon_start,
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
                    }).collect();

                    let should_include = !pick || tc.canonical || vf.transcript_variations.is_empty();
                    if should_include {
                        vf.transcript_variations.push(TranscriptVariation {
                            transcript_id: tc.transcript_id.clone(),
                            gene_id: tc.gene_id.clone(),
                            gene_symbol: tc.gene_symbol.clone(),
                            biotype: tc.biotype.clone(),
                            allele_annotations,
                            canonical: tc.canonical,
                            strand: tc.strand,
                            source: self.gff3_source.clone(),
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
            }
            vf.compute_most_severe();
        }

        Ok(variants.iter().map(|vf| output::format_json(vf)).collect())
    }
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

    if ins_len == 1 {
        let pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        // The dup range is ins_len bases ending at intron_offset
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        let end_pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (coding).
#[allow(dead_code)]
fn convert_ins_to_dup_range(
    hgvsc: &str,
    start_offset: i64,
    end_offset: i64,
    nearest_exon_cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc.find(":c.").map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
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

    if start_offset == end_offset {
        let pos = build_pos(nearest_exon_cdna_pos, start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        let end_pos = build_pos(nearest_exon_cdna_pos, end_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (non-coding).
#[allow(dead_code)]
fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start_offset: i64,
    end_offset: i64,
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

    if start_offset == end_offset {
        let pos = build_pos(start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(start_offset);
        let end_pos = build_pos(end_offset);
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

    if ins_len == 1 {
        let pos = build_pos(intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(start_offset);
        let end_pos = build_pos(intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// 3' shift an intronic indel along the transcript direction.
///
/// HGVS requires variants to be described at the most 3' position.
/// For intronic deletions and insertions/dups in repetitive regions,
/// the position must be shifted toward the 3' end of the transcript.
///
/// Returns the shifted genomic start and end positions.
fn three_prime_shift_intronic(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &oxivep_core::Allele,
    alt_allele: &oxivep_core::Allele,
    strand: oxivep_core::Strand,
    intron_genomic_start: u64,
    intron_genomic_end: u64,
) -> (u64, u64) {
    use oxivep_core::Allele;

    match (ref_allele, alt_allele) {
        // Deletion: shift the deleted bases toward 3' end
        (Allele::Sequence(ref_bases), Allele::Deletion) if !ref_bases.is_empty() => {
            let _del_len = ref_bases.len() as u64;
            let mut s = start;
            let mut e = end;

            match strand {
                oxivep_core::Strand::Forward => {
                    // 3' is toward higher genomic coords
                    // Shift while base after deletion matches first deleted base
                    loop {
                        let next_pos = e + 1;
                        if next_pos > intron_genomic_end { break; }
                        // Get the base after the current deletion end
                        let next_base = match seq_provider.fetch_sequence(chrom, next_pos, next_pos) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        // Get the first deleted base
                        let first_base = match seq_provider.fetch_sequence(chrom, s, s) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        if next_base == first_base {
                            s += 1;
                            e += 1;
                        } else {
                            break;
                        }
                    }
                }
                oxivep_core::Strand::Reverse => {
                    // 3' is toward lower genomic coords
                    loop {
                        if s == 0 || s - 1 < intron_genomic_start { break; }
                        let prev_pos = s - 1;
                        let prev_base = match seq_provider.fetch_sequence(chrom, prev_pos, prev_pos) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        let last_base = match seq_provider.fetch_sequence(chrom, e, e) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        if prev_base == last_base {
                            s -= 1;
                            e -= 1;
                        } else {
                            break;
                        }
                    }
                }
            }
            (s, e)
        }
        // Insertion/dup: shift toward 3' end using the actual inserted bases
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let mut pos = start;
            // Use the actual inserted bases (genomic strand) for comparison.
            // ins_bases here are on the genomic strand since they come from VCF normalization.
            let genomic_ins: Vec<u8> = ins_bases.iter().map(|b| b.to_ascii_uppercase()).collect();

            match strand {
                oxivep_core::Strand::Forward => {
                    // 3' is toward higher genomic coords
                    // Compare base at position with cycling insertion bases
                    let mut shift_count = 0u64;
                    loop {
                        if pos > intron_genomic_end { break; }
                        let check_base = match seq_provider.fetch_sequence(chrom, pos, pos) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        let idx = (shift_count as usize) % ins_len;
                        if check_base == genomic_ins[idx] {
                            pos += 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
                oxivep_core::Strand::Reverse => {
                    // 3' is toward lower genomic coords
                    // Compare base before position with cycling insertion bases (from end)
                    let mut shift_count = 0u64;
                    loop {
                        if pos == 0 || pos - 1 < intron_genomic_start { break; }
                        let check_pos = pos - 1;
                        let check_base = match seq_provider.fetch_sequence(chrom, check_pos, check_pos) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        // Compare from end of insertion, cycling backwards
                        let idx = ins_len - 1 - (shift_count as usize % ins_len);
                        if check_base == genomic_ins[idx] {
                            pos -= 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
            }
            (pos, pos.saturating_sub(1))
        }
        _ => (start, end),
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
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
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

/// Build a binary transcript cache from GFF3 + optional FASTA.
pub fn run_cache_build(gff3_path: &str, fasta_path: Option<&str>, output_path: &str) -> Result<()> {
    let gff_file = File::open(gff3_path)
        .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
    let mut transcripts = parse_gff3(gff_file)?;
    eprintln!("Loaded {} transcripts from {}", transcripts.len(), gff3_path);

    if let Some(fasta) = fasta_path {
        let fasta_file = File::open(fasta)
            .with_context(|| format!("Opening FASTA file: {}", fasta))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        let sp = FastaSequenceProvider::new(reader);
        eprintln!("Loaded reference FASTA from {}", fasta);

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

    oxivep_cache::transcript_cache::save_cache(&transcripts, Path::new(output_path))?;
    eprintln!("Saved transcript cache to {}", output_path);
    Ok(())
}

/// Quick VCF pre-scan to collect variant regions for indexed GFF3 loading.
/// Returns merged (chrom, start, end) regions expanded by the given distance.
fn prescan_vcf_regions(vcf_path: &str, distance: u64) -> Result<Vec<(String, u64, u64)>> {
    let file = File::open(vcf_path)
        .with_context(|| format!("Pre-scanning VCF: {}", vcf_path))?;
    let reader = io::BufReader::new(file);

    let mut regions: HashMap<String, (u64, u64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() { continue; }
        let mut fields = line.split('\t');
        let chrom = match fields.next() { Some(c) => c.to_string(), None => continue };
        let pos: u64 = match fields.next().and_then(|p| p.parse().ok()) { Some(p) => p, None => continue };
        let start = pos.saturating_sub(distance);
        let end = pos + distance;

        let entry = regions.entry(chrom).or_insert((start, end));
        entry.0 = entry.0.min(start);
        entry.1 = entry.1.max(end);
    }

    Ok(regions.into_iter().map(|(chrom, (s, e))| (chrom, s, e)).collect())
}

// =============================================================================
// SA Build: Build supplementary annotation databases from source VCFs
// =============================================================================

/// Standard chromosome ordering for SA builds.
fn standard_chrom_map() -> (Vec<String>, std::collections::HashMap<String, u16>) {
    let chroms: Vec<String> = (1..=22)
        .map(|i| format!("chr{}", i))
        .chain(["chrX", "chrY", "chrM"].iter().map(|s| s.to_string()))
        .collect();
    let map: std::collections::HashMap<String, u16> = chroms
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i as u16))
        .collect();
    (chroms, map)
}

/// Build a supplementary annotation .osa file from a source VCF.
pub fn run_sa_build(source: &str, input: &str, output: &str, assembly: &str) -> Result<()> {
    use oxivep_sa::index::IndexHeader;
    use oxivep_sa::writer::SaWriter;

    let (chrom_list, chrom_map) = standard_chrom_map();

    let header = match source {
        "clinvar" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "clinvar".into(),
            name: "ClinVar".into(),
            version: "latest".into(),
            description: format!("ClinVar annotations for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: true,
            is_positional: false,
        },
        "gnomad" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "gnomad".into(),
            name: "gnomAD".into(),
            version: "latest".into(),
            description: format!("gnomAD population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbsnp" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "dbsnp".into(),
            name: "dbSNP".into(),
            version: "latest".into(),
            description: format!("dbSNP RS IDs for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "phylop" | "gerp" | "dann" => {
            let json_key = match source {
                "phylop" => "phylopScore",
                "gerp" => "gerpScore",
                "dann" => "dannScore",
                _ => unreachable!(),
            };
            IndexHeader {
                schema_version: oxivep_sa::common::SCHEMA_VERSION,
                json_key: json_key.into(),
                name: source.to_uppercase().into(),
                version: "latest".into(),
                description: format!("{} conservation/prediction scores for {}", source, assembly),
                assembly: assembly.into(),
                match_by_allele: false,
                is_array: false,
                is_positional: true,
            }
        },
        "revel" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "revel".into(),
            name: "REVEL".into(),
            version: "latest".into(),
            description: format!("REVEL missense pathogenicity scores for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "spliceai" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "spliceAI".into(),
            name: "SpliceAI".into(),
            version: "latest".into(),
            description: format!("SpliceAI splice site predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "primateai" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "primateAI".into(),
            name: "PrimateAI".into(),
            version: "latest".into(),
            description: format!("PrimateAI pathogenicity predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbnsfp" => IndexHeader {
            schema_version: oxivep_sa::common::SCHEMA_VERSION,
            json_key: "dbnsfp".into(),
            name: "dbNSFP".into(),
            version: "latest".into(),
            description: format!("dbNSFP SIFT/PolyPhen predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        _ => anyhow::bail!(
            "Unknown source: {}. Supported: clinvar, gnomad, dbsnp, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp",
            source
        ),
    };

    eprintln!("Building {} .osa from: {}", source, input);

    let file = File::open(input)
        .with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records = match source {
        "clinvar" => oxivep_sa::sources::clinvar::parse_clinvar_vcf(buf_reader, &chrom_map)?,
        "gnomad" => oxivep_sa::sources::gnomad::parse_gnomad_vcf(buf_reader, &chrom_map)?,
        "dbsnp" => oxivep_sa::sources::dbsnp::parse_dbsnp_vcf(buf_reader, &chrom_map)?,
        "phylop" => oxivep_sa::sources::scores::parse_wigfix(buf_reader, &chrom_map)?,
        "gerp" | "dann" => oxivep_sa::sources::scores::parse_score_tsv(buf_reader, &chrom_map, false)?,
        "revel" => oxivep_sa::sources::revel::parse_revel(buf_reader, &chrom_map, 2)?,
        "spliceai" => oxivep_sa::sources::spliceai::parse_spliceai_vcf(buf_reader, &chrom_map)?,
        "primateai" => oxivep_sa::sources::primateai::parse_primateai(buf_reader, &chrom_map)?,
        "dbnsfp" => oxivep_sa::sources::dbnsfp::parse_dbnsfp(buf_reader, &chrom_map)?,
        _ => unreachable!(),
    };

    eprintln!("Parsed {} records from {}", records.len(), source);

    let output_path = Path::new(output);
    let mut writer = SaWriter::new(header);
    writer.write_to_files(output_path, records.into_iter(), &chrom_list)?;

    eprintln!(
        "Wrote: {} and {}",
        output_path.with_extension("osa").display(),
        output_path.with_extension("osa.idx").display()
    );

    Ok(())
}

// =============================================================================
// SA Provider Loading: discover and load .osa files from --sa-dir
// =============================================================================

/// Load all .osa annotation providers from a directory.
pub fn load_sa_providers(
    sa_dir: &Path,
) -> Result<Vec<Box<dyn oxivep_cache::annotation::AnnotationProvider>>> {
    use oxivep_sa::reader::SaReader;

    let mut providers: Vec<Box<dyn oxivep_cache::annotation::AnnotationProvider>> = Vec::new();

    if !sa_dir.is_dir() {
        anyhow::bail!("SA directory does not exist: {}", sa_dir.display());
    }

    for entry in std::fs::read_dir(sa_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("osa") {
            match SaReader::open(&path) {
                Ok(reader) => {
                    eprintln!(
                        "Loaded SA provider: {} ({})",
                        reader.name(),
                        path.display()
                    );
                    providers.push(Box::new(reader));
                }
                Err(e) => {
                    eprintln!(
                        "Warning: could not load SA file {}: {}",
                        path.display(),
                        e
                    );
                }
            }
        }
    }

    Ok(providers)
}
