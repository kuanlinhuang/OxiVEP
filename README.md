# OxiVEP

A high-performance Variant Effect Predictor written in Rust. OxiVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions, structural variants) on genes, transcripts, and protein sequences, with direct integration of clinical and population databases.

OxiVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and [Illumina Nirvana](https://github.com/Illumina/Nirvana), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

## Features

- **Variant Consequence Prediction** — Classifies variants using 49 [Sequence Ontology](http://www.sequenceontology.org/) terms (missense, frameshift, splice donor, copy_number_change, transcript_ablation, etc.)
- **Structural Variant Support** — Full SV pipeline: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<BND>`, `<INS>`, `<STR>` with SV-specific consequence prediction
- **Supplementary Annotations** — Direct integration with ClinVar, gnomAD, dbSNP, COSMIC, 1000 Genomes, TOPMed, MitoMap via the native OxiSA format (v1: zstd block compression; v2: echtvar-inspired chunked ZIP with Var32 encoding, parallel u32 value arrays, delta encoding, and LRU caching)
- **Prediction Scores** — PhyloP, GERP, REVEL, SpliceAI, PrimateAI, DANN conservation and pathogenicity scores; SIFT/PolyPhen via dbNSFP
- **Gene-Level Annotations** — OMIM phenotypes, gnomAD gene constraint (pLI, LOEUF), ClinGen gene-disease validity
- **Filter Engine** — Expression-based filtering compatible with VEP's filter_vep syntax
- **HGVS Nomenclature** — Generates HGVSg, HGVSc, and HGVSp notations with 3' normalization
- **Multiple Output Formats** — VCF (with 47-field CSQ), tab-delimited, JSON (including Nirvana-style structured output)
- **Multi-Sample Support** — Parse FORMAT/GT/DP/GQ/AD fields per sample with genotype classification
- **Regulatory Region Detection** — Promoters, enhancers, CTCF binding sites, TF binding sites from Ensembl regulatory build
- **Mitochondrial Support** — Circular coordinate handling, vertebrate mitochondrial codon table (NCBI table 2)
- **Custom Annotations** — User-provided VCF and BED annotation files
- **Web Interface** — Built-in web GUI for interactive variant annotation
- **GFF3 Annotation Support** — Load gene models from standard GFF3 files (any organism)

## Quick Start

### 1. Install Rust (if you don't have it)

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

### 2. Build and install OxiVEP

```bash
git clone https://github.com/kuanlinhuang/OxiVEP.git
cd OxiVEP

# Build and install to your PATH (~/.cargo/bin/oxivep)
cargo install --path crates/oxivep-cli

# Verify it works
oxivep --version
```

> **Note:** `cargo install` places the binary in `~/.cargo/bin/`. If `oxivep` is not found after install, run `source "$HOME/.cargo/env"` or add this line to your `~/.zshrc` (or `~/.bashrc`):
> ```bash
> source "$HOME/.cargo/env"
> ```

### 3. Try it — annotate the included test data

OxiVEP ships with a small test VCF and GFF3 so you can try it immediately:

```bash
# Annotate 12 test variants covering SNVs, indels, splice sites, UTRs, and intergenic regions
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format tab
```

### 4. Build supplementary annotation databases

```bash
# Build ClinVar annotation database
oxivep sa-build --source clinvar --input clinvar.vcf.gz --output clinvar

# Build gnomAD population frequency database
oxivep sa-build --source gnomad --input gnomad.genomes.v4.vcf.bgz --output gnomad

# Build PhyloP conservation scores
oxivep sa-build --source phylop --input hg38.phyloP100way.wigFix.gz --output phylop

# Build SpliceAI predictions
oxivep sa-build --source spliceai --input spliceai_scores.vcf.gz --output spliceai
```

### 5. Annotate with supplementary databases

```bash
# Annotate with all databases in a directory
oxivep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 Homo_sapiens.GRCh38.112.gff3 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir /path/to/annotation_databases/ \
  --hgvs
```

### 6. Filter annotated variants

```bash
# Filter for high-impact or rare missense variants
oxivep filter \
  -i annotated.vcf \
  --filter "IMPACT is HIGH or (Consequence in missense_variant and AF < 0.001)"
```

### 7. Launch the web interface

```bash
oxivep web --port 8080
```

Open http://localhost:8080 in your browser.

## Using Your Own Data

### Annotate a real VCF against Ensembl gene models

```bash
# Download Ensembl GFF3 (human GRCh38)
wget https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz
gunzip Homo_sapiens.GRCh38.112.gff3.gz

# Annotate your VCF
oxivep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 Homo_sapiens.GRCh38.112.gff3 \
  --hgvs

# Or pipe from bcftools
bcftools view sample.vcf.gz chr21 | oxivep annotate -i - --gff3 Homo_sapiens.GRCh38.112.gff3 --output-format tab
```

### Mouse, zebrafish, or other organisms

```bash
# Mouse
wget https://ftp.ensembl.org/pub/release-112/gff3/mus_musculus/Mus_musculus.GRCm39.112.gff3.gz

# Zebrafish
wget https://ftp.ensembl.org/pub/release-112/gff3/danio_rerio/Danio_rerio.GRCz11.112.gff3.gz
```

OxiVEP works with any organism — just provide the matching GFF3.

## Supplementary Annotation Sources

OxiVEP supports direct integration with clinical and population databases through its native OxiSA binary format. Build once with `oxivep sa-build`, then use `--sa-dir` to annotate:

| Source | Type | Description | Build Command |
|--------|------|-------------|---------------|
| **ClinVar** | Allele-specific | Clinical significance, review status, phenotypes | `--source clinvar` |
| **gnomAD** | Allele-specific | Population frequencies (8 populations), allele counts | `--source gnomad` |
| **dbSNP** | Allele-specific | RS IDs, global minor allele frequency | `--source dbsnp` |
| **COSMIC** | Allele-specific | Somatic mutations, gene, sample counts | `--source cosmic` |
| **1000 Genomes** | Allele-specific | Population frequencies (AFR, AMR, EAS, EUR, SAS) | `--source onekg` |
| **TOPMed** | Allele-specific | Population frequencies, allele counts | `--source topmed` |
| **MitoMap** | Allele-specific | Mitochondrial disease associations | `--source mitomap` |
| **PhyloP** | Positional | Phylogenetic conservation scores | `--source phylop` |
| **GERP** | Positional | Evolutionary rate profiling | `--source gerp` |
| **DANN** | Positional | Deleterious annotation scores | `--source dann` |
| **REVEL** | Allele-specific | Missense pathogenicity predictions | `--source revel` |
| **SpliceAI** | Allele-specific | Splice site effect predictions (delta scores) | `--source spliceai` |
| **PrimateAI** | Allele-specific | Primate-based pathogenicity | `--source primateai` |
| **dbNSFP** | Allele-specific | SIFT/PolyPhen predictions | `--source dbnsfp` |

## Command Reference

### `oxivep annotate`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VCF file (`-` for stdin) | *required* |
| `-o, --output` | Output file (`-` for stdout) | `-` |
| `--gff3` | GFF3 gene annotation file | -- |
| `--fasta` | Reference FASTA file | -- |
| `--output-format` | `vcf`, `tab`, or `json` | `vcf` |
| `--hgvs` | Include HGVS notations | off |
| `--pick` | Report only the most severe consequence per variant | off |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--sa-dir` | Directory containing .osa supplementary annotation files | -- |
| `--cache-dir` | Path to VEP cache directory for known variant annotation | -- |
| `--transcript-cache` | Path to binary transcript cache file | -- |

### `oxivep sa-build`

| Flag | Description | Default |
|------|-------------|---------|
| `--source` | Source type (clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp) | *required* |
| `-i, --input` | Input file (VCF/TSV/wigFix, supports .gz) | *required* |
| `-o, --output` | Output base path (creates .osa and .osa.idx) | *required* |
| `--assembly` | Genome assembly | `GRCh38` |

### `oxivep filter`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VEP-annotated VCF | *required* |
| `-o, --output` | Output file | `-` |
| `--filter` | Filter expression (filter_vep-compatible syntax) | *required* |

Filter syntax examples:
```
IMPACT is HIGH
Consequence in missense_variant,stop_gained,frameshift_variant
AF < 0.001
IMPACT is HIGH and AF < 0.01
(IMPACT is HIGH or IMPACT is MODERATE) and not Consequence is synonymous_variant
```

### `oxivep web`

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 gene annotation file | -- |
| `--fasta` | Reference FASTA file | -- |
| `--port` | HTTP port | `8080` |

### `oxivep cache`

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 annotation file | *required* |
| `--fasta` | Reference FASTA (for pre-building sequences) | -- |
| `-o, --output` | Output cache file path | *required* |

## Output Formats

### VCF Output

Annotations are added as a `CSQ` field in the INFO column with 47 pipe-delimited fields matching Ensembl VEP's extended format.

### Tab Output

One line per variant-transcript-allele combination with 17 columns.

### JSON Output

Structured JSON with `transcript_consequences` array per variant, including supplementary annotations from SA providers (ClinVar, gnomAD, etc.) and gene-level annotations.

## Consequence Types

OxiVEP predicts 49 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification, TFBS_ablation, regulatory_region_ablation |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant, regulatory_region_amplification, TFBS_amplification |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, copy_number_change, copy_number_increase, copy_number_decrease, short_tandem_repeat_change, transcript_variant, and others |

## Architecture

```
crates/
  oxivep-core/         # Core types: Consequence (49 SO terms), VariantType, Allele, Impact
  oxivep-genome/       # Transcript, Exon, Gene, CodonTable, mitochondrial codon table
  oxivep-cache/        # GFF3 parser, FASTA reader, annotation providers, regulatory regions
  oxivep-consequence/  # Consequence prediction: small variants + SV predictor
  oxivep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  oxivep-io/           # VCF parser (incl. SVs), output formatters, multi-sample parsing
  oxivep-filter/       # Filter engine: lexer, parser, evaluator (filter_vep-compatible)
  oxivep-sa/           # Supplementary annotation format (OxiSA):
                       #   v1 (.osa): zstd block compression, binary search
                       #   v2 (.osa2): echtvar-inspired chunked ZIP with Var32 encoding,
                       #     parallel u32 value arrays, delta encoding, LRU chunk cache,
                       #     Bloom filters for negative lookups
                       # Source parsers: ClinVar, gnomAD, dbSNP, COSMIC, 1000G, TOPMed,
                       # MitoMap, PhyloP, GERP, DANN, REVEL, SpliceAI, PrimateAI, dbNSFP
                       # Custom VCF/BED annotation providers
  oxivep-cli/          # CLI binary, annotation pipeline, sa-build, web server
web/                   # Web GUI (HTML/CSS/JS)
tests/                 # Test VCF and GFF3 files
```

## Running Tests

```bash
cargo test --workspace          # 233 tests
cargo test -p oxivep-consequence  # Consequence prediction tests (incl. SV)
cargo test -p oxivep-filter       # Filter engine tests
cargo test -p oxivep-sa           # Supplementary annotation format tests
```

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), release build with LTO. Median of 3 runs, full Ensembl annotations with FASTA and HGVS.

### Multi-Organism Throughput

| Organism | Transcripts | Variants | Time | Throughput |
|----------|-------------|----------|------|------------|
| Yeast (R64, full genome) | 7,036 | 260,526 | 1.47s | **176,796 v/s** |
| C. elegans (WBcel235, full) | 44,365 | 100,000 | 2.70s | **37,023 v/s** |
| Drosophila (BDGP6, full) | 35,442 | 100,000 | 2.76s | **36,286 v/s** |
| Arabidopsis (TAIR10, full) | 54,013 | 500,000 | 5.95s | **83,967 v/s** |
| Human chr22 (GRCh38) | 11,605 | 500,000 | 4.76s | **105,117 v/s** |
| Mouse (GRCm39, full genome) | 142,626 | 500,000 | 12.32s | **40,577 v/s** |
| Human full WGS (GRCh38) | 508,530 | 3,893,341 | 29.69s | **131,155 v/s** |

### vs. Ensembl VEP

| Metric | Ensembl VEP (Perl) | OxiVEP (Rust) |
|--------|-------------------|---------------|
| Full WGS (3.9M variants) | est. ~109 min | **29.7s** |
| 500K variants (chr22) | 674s | **4.76s (142x)** |
| Peak memory (100K variants) | ~500 MB | **2.8 MB** |
| Binary size | ~200 MB installed | **2.3 MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None** |

## License

Apache License 2.0

## Related Tools

- **[RastQC](https://doi.org/10.64898/2026.03.31.715630)** — High-performance sequencing quality control written in Rust (Huang, 2026)

## Acknowledgements

OxiVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI and [Illumina Nirvana](https://github.com/Illumina/Nirvana). The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework.
