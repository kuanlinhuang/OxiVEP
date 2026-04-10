use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, SequenceProvider, TranscriptProvider,
};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence};
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use std::fs::File;
use std::path::Path;
use std::sync::Mutex;

/// Pre-loaded annotation context for the web server.
/// Holds transcript models, sequence provider, and supplementary annotation providers.
pub struct AnnotationContext {
    pub transcript_provider: IndexedTranscriptProvider,
    pub seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>>,
    pub predictor: ConsequencePredictor,
    pub gff3_source: Option<String>,
    pub distance: u64,
    pub hgvs: bool,
    /// Supplementary annotation providers (ClinVar, gnomAD, etc.)
    /// Wrapped in Mutex because SA readers use internal caches that need &mut.
    pub sa_providers: Vec<Mutex<Box<dyn AnnotationProvider>>>,
}

impl AnnotationContext {
    /// Build a context from GFF3, optional FASTA, and optional SA directory.
    pub fn new(
        gff3: Option<&str>,
        fasta: Option<&str>,
        sa_dir: Option<&str>,
        distance: u64,
    ) -> Result<Self> {
        let gff3_source: Option<String> = gff3.map(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| p.to_string())
        });

        let mut transcripts = if let Some(gff3_path) = gff3 {
            let cache_path =
                fastvep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
            let from_cache = if cache_path.exists() {
                let is_fresh = fastvep_cache::transcript_cache::cache_is_fresh(
                    &cache_path,
                    Path::new(gff3_path),
                );
                if is_fresh {
                    fastvep_cache::transcript_cache::load_cache(&cache_path).ok()
                } else {
                    None
                }
            } else {
                None
            };
            if let Some(trs) = from_cache {
                tracing::info!("Loaded {} transcripts from cache", trs.len());
                trs
            } else {
                let gff_file = File::open(gff3_path)
                    .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                let trs = parse_gff3(gff_file)?;
                tracing::info!("Loaded {} transcripts from {}", trs.len(), gff3_path);
                if let Err(e) = fastvep_cache::transcript_cache::save_cache(&trs, &cache_path) {
                    tracing::warn!("Could not save cache: {}", e);
                }
                trs
            }
        } else {
            Vec::new()
        };

        let seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>> =
            if let Some(fasta_path) = fasta {
                let fai_path = format!("{}.fai", fasta_path);
                if Path::new(&fai_path).exists() {
                    let reader =
                        fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
                    tracing::info!("Memory-mapped FASTA from {}", fasta_path);
                    Some(Box::new(
                        fastvep_cache::providers::MmapFastaSequenceProvider::new(reader),
                    ))
                } else {
                    let fasta_file = File::open(fasta_path)
                        .with_context(|| format!("Opening FASTA: {}", fasta_path))?;
                    let reader = FastaReader::from_reader(fasta_file)?;
                    tracing::info!("Loaded FASTA from {}", fasta_path);
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
                    if tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                    {
                        built += 1;
                    }
                }
            }
            tracing::info!("Built sequences for {} coding transcripts", built);
        }

        let transcript_provider = IndexedTranscriptProvider::new(transcripts);
        let predictor = ConsequencePredictor::new(distance, distance);

        // Load supplementary annotation providers (.osa, .osa2 files)
        let sa_providers = if let Some(dir) = sa_dir {
            load_sa_providers(Path::new(dir))?
        } else {
            Vec::new()
        };

        Ok(Self {
            transcript_provider,
            seq_provider,
            predictor,
            gff3_source,
            distance,
            hgvs: true,
            sa_providers,
        })
    }

    pub fn transcript_count(&self) -> usize {
        self.transcript_provider.transcript_count()
    }

    /// Names of loaded supplementary annotation sources.
    pub fn sa_source_names(&self) -> Vec<String> {
        self.sa_providers
            .iter()
            .filter_map(|m| {
                let guard = m.lock().ok()?;
                Some(guard.name().to_string())
            })
            .collect()
    }

    /// Load a genome from GFF3 path (+ optional FASTA + optional SA directory).
    /// Replaces transcripts, sequence provider, and SA providers.
    pub fn load_genome(
        &mut self,
        gff3_path: &str,
        fasta_path: Option<&str>,
        sa_dir: Option<&str>,
    ) -> Result<usize> {
        let new_ctx = Self::new(Some(gff3_path), fasta_path, sa_dir, self.distance)?;
        let tr_count = new_ctx.transcript_provider.transcript_count();
        self.transcript_provider = new_ctx.transcript_provider;
        self.seq_provider = new_ctx.seq_provider;
        self.predictor = new_ctx.predictor;
        self.gff3_source = new_ctx.gff3_source;
        self.sa_providers = new_ctx.sa_providers;
        Ok(tr_count)
    }

    /// Replace the transcript models by parsing GFF3 text uploaded from the browser.
    pub fn update_gff3_text(&mut self, gff3_text: &str) -> Result<(usize, usize)> {
        let mut transcripts = parse_gff3(gff3_text.as_bytes())?;
        let gene_count = {
            let mut genes = std::collections::HashSet::new();
            for t in &transcripts {
                genes.insert(t.gene.stable_id.clone());
            }
            genes.len()
        };

        if let Some(ref sp) = self.seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                    {
                        built += 1;
                    }
                }
            }
            if built > 0 {
                tracing::info!("Built sequences for {} coding transcripts", built);
            }
        }

        let tr_count = transcripts.len();
        self.transcript_provider = IndexedTranscriptProvider::new(transcripts);
        self.gff3_source = Some("user-upload".to_string());
        tracing::info!(
            "Updated GFF3: {} genes, {} transcripts",
            gene_count,
            tr_count
        );
        Ok((gene_count, tr_count))
    }

    /// Annotate VCF text and return JSON results.
    pub fn annotate_vcf_text(&self, vcf_text: &str, pick: bool) -> Result<Vec<serde_json::Value>> {
        let mut vcf_parser = VcfParser::new(vcf_text.as_bytes())?;
        let mut variants = vcf_parser.read_all()?;

        for vf in &mut variants {
            let chrom = &vf.position.chromosome;
            let query_start = vf.position.start.saturating_sub(self.distance).max(1);
            let query_end = vf.position.end + self.distance;
            let overlapping = self
                .transcript_provider
                .get_transcripts(chrom, query_start, query_end)
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
                                protein_position: zip_positions(
                                    ac.protein_start,
                                    ac.protein_end,
                                ),
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
                                ann.hgvsg = Some(fastvep_hgvs::hgvsg(
                                    chrom,
                                    vf.position.start,
                                    vf.position.end,
                                    &vf.ref_allele,
                                    &ac.allele,
                                ));
                                if let Some(tr) = transcript {
                                    let versioned_tid = match tr.version {
                                        Some(v) => {
                                            format!("{}.{}", tc.transcript_id, v)
                                        }
                                        None => tc.transcript_id.to_string(),
                                    };
                                    let (hgvs_ref, hgvs_alt) =
                                        if tr.strand == fastvep_core::Strand::Reverse {
                                            (
                                                complement_allele(&vf.ref_allele),
                                                complement_allele(&ac.allele),
                                            )
                                        } else {
                                            (vf.ref_allele.clone(), ac.allele.clone())
                                        };
                                    if let Some(coding_start) = tr.cdna_coding_start {
                                        if let (Some(cs), Some(ce)) =
                                            (ac.cdna_start, ac.cdna_end)
                                        {
                                            let (cs, ce) = (cs.min(ce), cs.max(ce));
                                            ann.hgvsc = fastvep_hgvs::hgvsc_with_seq(
                                                &versioned_tid,
                                                cs,
                                                ce,
                                                &hgvs_ref,
                                                &hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                                tr.spliced_seq.as_deref(),
                                                tr.codon_table_start_phase,
                                            );
                                        } else if ac.intron.is_some() {
                                            if let Some((cdna_pos, offset)) =
                                                tr.genomic_to_intronic_cdna(vf.position.start)
                                            {
                                                let (end_cdna, end_offset) =
                                                    if vf.position.start != vf.position.end {
                                                        tr.genomic_to_intronic_cdna(
                                                            vf.position.end,
                                                        )
                                                        .map(|(c, o)| (Some(c), Some(o)))
                                                        .unwrap_or((None, None))
                                                    } else {
                                                        (None, None)
                                                    };
                                                ann.hgvsc = fastvep_hgvs::hgvsc_intronic_range(
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
                                            }
                                        }
                                    } else if let (Some(cs), Some(ce)) =
                                        (ac.cdna_start, ac.cdna_end)
                                    {
                                        ann.hgvsc = fastvep_hgvs::hgvsc_noncoding(
                                            &versioned_tid,
                                            cs,
                                            ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                        );
                                    } else if ac.intron.is_some() {
                                        if let Some((cdna_pos, offset)) =
                                            tr.genomic_to_intronic_cdna(vf.position.start)
                                        {
                                            let (end_cdna, end_offset) =
                                                if vf.position.start != vf.position.end {
                                                    tr.genomic_to_intronic_cdna(vf.position.end)
                                                        .map(|(c, o)| (Some(c), Some(o)))
                                                        .unwrap_or((None, None))
                                                } else {
                                                    (None, None)
                                                };
                                            ann.hgvsc =
                                                fastvep_hgvs::hgvsc_noncoding_intronic_range(
                                                    &versioned_tid,
                                                    cdna_pos,
                                                    offset,
                                                    end_cdna,
                                                    end_offset,
                                                    &hgvs_ref,
                                                    &hgvs_alt,
                                                );
                                        }
                                    }

                                    // HGVSp
                                    if let (Some(ref aa), Some(ps)) =
                                        (&ac.amino_acids, ac.protein_start)
                                    {
                                        if let Some(ref pid) = tr.protein_id {
                                            let versioned_pid: String = match tr.protein_version {
                                                Some(v) => {
                                                    let suffix = format!(".{}", v);
                                                    if pid.ends_with(suffix.as_str()) {
                                                        pid.clone()
                                                    } else {
                                                        format!("{}.{}", pid, v)
                                                    }
                                                }
                                                None => pid.clone(),
                                            };
                                            let is_fs = ac
                                                .consequences
                                                .contains(&Consequence::FrameshiftVariant);
                                            if is_fs {
                                                if let (
                                                    Some(ref spliced),
                                                    Some(coding_start),
                                                    Some(cds_s),
                                                ) = (
                                                    &tr.spliced_seq,
                                                    tr.cdna_coding_start,
                                                    ac.cds_start,
                                                ) {
                                                    let coding_start_idx =
                                                        (coding_start - 1) as usize;
                                                    let spliced_bytes: &[u8] = spliced.as_bytes();
                                                    let ref_from_cds =
                                                        &spliced_bytes[coding_start_idx..];
                                                    let cds_idx = (cds_s - 1) as usize;
                                                    let mut alt_from_cds =
                                                        ref_from_cds.to_vec();
                                                    if ac.allele == Allele::Deletion {
                                                        let del_len = vf.ref_allele.len();
                                                        let end = (cds_idx + del_len)
                                                            .min(alt_from_cds.len());
                                                        alt_from_cds.drain(cds_idx..end);
                                                    } else if let Allele::Sequence(ins_bases) =
                                                        &ac.allele
                                                    {
                                                        let mut bases = ins_bases.clone();
                                                        if tr.strand
                                                            == fastvep_core::Strand::Reverse
                                                        {
                                                            bases = bases
                                                                .iter()
                                                                .map(|&b| match b {
                                                                    b'A' => b'T',
                                                                    b'T' => b'A',
                                                                    b'C' => b'G',
                                                                    b'G' => b'C',
                                                                    o => o,
                                                                })
                                                                .collect();
                                                        }
                                                        for (j, &b) in
                                                            bases.iter().enumerate()
                                                        {
                                                            if cds_idx + j
                                                                <= alt_from_cds.len()
                                                            {
                                                                alt_from_cds
                                                                    .insert(cds_idx + j, b);
                                                            }
                                                        }
                                                    }
                                                    let codon_start = cds_idx / 3;
                                                    ann.hgvsp =
                                                        fastvep_hgvs::hgvsp_frameshift(
                                                            &versioned_pid,
                                                            ref_from_cds,
                                                            &alt_from_cds,
                                                            codon_start,
                                                        );
                                                }
                                            } else {
                                                let ref_aa_byte = aa
                                                    .0
                                                    .as_bytes()
                                                    .first()
                                                    .copied()
                                                    .unwrap_or(b'X');
                                                let alt_aa_byte = aa
                                                    .1
                                                    .as_bytes()
                                                    .first()
                                                    .copied()
                                                    .unwrap_or(b'X');
                                                ann.hgvsp = fastvep_hgvs::hgvsp(
                                                    &versioned_pid,
                                                    ps,
                                                    ref_aa_byte,
                                                    alt_aa_byte,
                                                    false,
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                            ann
                        })
                        .collect();

                    let should_include =
                        !pick || tc.canonical || vf.transcript_variations.is_empty();
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
                            mane_plus_clinical: transcript
                                .and_then(|t| t.mane_plus_clinical.clone()),
                            tsl: transcript.and_then(|t| t.tsl),
                            appris: transcript.and_then(|t| t.appris.clone()),
                            ccds: transcript.and_then(|t| t.ccds.clone()),
                            symbol_source: transcript
                                .and_then(|t| t.gene.symbol_source.clone()),
                            hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                            flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                        });
                    }
                }
            }

            // Supplementary annotation: query SA providers for each allele
            if !self.sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        let alt_str = aa.allele.to_string();
                        let ref_str = vf.ref_allele.to_string();
                        for sa in &self.sa_providers {
                            let sa_guard = sa.lock().unwrap();
                            if let Ok(Some(ann)) = sa_guard.annotate_position(
                                chrom,
                                vf.position.start,
                                &ref_str,
                                &alt_str,
                            ) {
                                let json_str = match ann {
                                    AnnotationValue::Json(j) => j,
                                    AnnotationValue::Positional(j) => j,
                                    AnnotationValue::Interval(v) => {
                                        format!("[{}]", v.join(","))
                                    }
                                };
                                aa.supplementary
                                    .push((sa_guard.json_key().to_string(), json_str));
                            }
                        }
                    }
                }
            }

            vf.compute_most_severe();
        }

        Ok(variants.iter().map(|vf| output::format_json(vf)).collect())
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
                impact: fastvep_core::Impact::Modifier,
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
            strand: fastvep_core::Strand::Forward,
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

fn zip_positions(start: Option<u64>, end: Option<u64>) -> Option<(u64, u64)> {
    match (start, end) {
        (Some(s), Some(e)) => Some((s.min(e), s.max(e))),
        (Some(s), None) => Some((s, s)),
        (None, Some(e)) => Some((e, e)),
        _ => None,
    }
}

fn complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => {
            let comp: Vec<u8> = bases
                .iter()
                .map(|&b| match b {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    other => other,
                })
                .collect();
            Allele::Sequence(comp)
        }
        other => other.clone(),
    }
}

/// Load supplementary annotation providers (.osa, .osa2 files) from a directory.
fn load_sa_providers(
    sa_dir: &Path,
) -> Result<Vec<Mutex<Box<dyn AnnotationProvider>>>> {
    use fastvep_sa::reader::SaReader;
    use fastvep_sa::reader_v2::Osa2Reader;

    let mut providers: Vec<Mutex<Box<dyn AnnotationProvider>>> = Vec::new();

    if !sa_dir.is_dir() {
        anyhow::bail!("SA directory does not exist: {}", sa_dir.display());
    }

    for entry in std::fs::read_dir(sa_dir)? {
        let entry = entry?;
        let path = entry.path();
        let ext = path.extension().and_then(|e| e.to_str());

        match ext {
            Some("osa2") => match Osa2Reader::open(&path) {
                Ok(reader) => {
                    tracing::info!("Loaded SA v2: {} ({})", reader.name(), path.display());
                    providers.push(Mutex::new(Box::new(reader)));
                }
                Err(e) => {
                    tracing::warn!("Could not load {}: {}", path.display(), e);
                }
            },
            Some("osa") => match SaReader::open(&path) {
                Ok(reader) => {
                    tracing::info!("Loaded SA: {} ({})", reader.name(), path.display());
                    providers.push(Mutex::new(Box::new(reader)));
                }
                Err(e) => {
                    tracing::warn!("Could not load {}: {}", path.display(), e);
                }
            },
            _ => {}
        }
    }

    Ok(providers)
}
