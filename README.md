# fastVEP

A high-performance Variant Effect Predictor written in Rust. fastVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions, structural variants) on genes, transcripts, and protein sequences, with direct integration of clinical and population databases.

fastVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and [Illumina Nirvana](https://github.com/Illumina/Nirvana), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

**Try it now:** A hosted web server is available at [fastVEP.org](https://fastVEP.org) — paste VCF data and get annotated results instantly, no installation required.

## Features

- **Variant Consequence Prediction** — Classifies variants using 49 [Sequence Ontology](http://www.sequenceontology.org/) terms (missense, frameshift, splice donor, copy_number_change, transcript_ablation, etc.)
- **Structural Variant Support** — Full SV pipeline: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<BND>`, `<INS>`, `<STR>` with SV-specific consequence prediction
- **Supplementary Annotations** — Direct integration with ClinVar, gnomAD, dbSNP, COSMIC, 1000 Genomes, TOPMed, MitoMap via the native fastSA format (v1: zstd block compression; v2: echtvar-inspired chunked ZIP with Var32 encoding, parallel u32 value arrays, delta encoding, and LRU caching)
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

### 2. Build and install fastVEP

```bash
git clone https://github.com/Huang-lab/fastVEP.git
cd fastVEP

# Build and install both binaries to ~/.cargo/bin/
cargo install --path crates/fastvep-cli   # fastvep (CLI annotator)
cargo install --path crates/fastvep-web   # fastvep-web (production web server)

# Verify it works
fastvep --version
```

> **Note:** `cargo install` places the binary in `~/.cargo/bin/`. If `fastvep` is not found after install, run `source "$HOME/.cargo/env"` or add this line to your `~/.zshrc` (or `~/.bashrc`):
> ```bash
> source "$HOME/.cargo/env"
> ```

### 3. Try it — annotate the included test data

fastVEP ships with a small test VCF and GFF3 so you can try it immediately:

```bash
# Annotate 12 test variants covering SNVs, indels, splice sites, UTRs, and intergenic regions
fastvep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format tab
```

### 4. Build supplementary annotation databases

```bash
# Build ClinVar annotation database
fastvep sa-build --source clinvar --input clinvar.vcf.gz --output clinvar

# Build gnomAD population frequency database
fastvep sa-build --source gnomad --input gnomad.genomes.v4.vcf.bgz --output gnomad

# Build PhyloP conservation scores
fastvep sa-build --source phylop --input hg38.phyloP100way.wigFix.gz --output phylop

# Build SpliceAI predictions
fastvep sa-build --source spliceai --input spliceai_scores.vcf.gz --output spliceai
```

### 5. Annotate with supplementary databases

```bash
# Annotate with all databases in a directory
fastvep annotate \
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
fastvep filter \
  -i annotated.vcf \
  --filter "IMPACT is HIGH or (Consequence in missense_variant and AF < 0.001)"
```

### 7. Launch the web interface

```bash
# Quick start — uses a built-in example gene model (OR4F5, chr1)
fastvep-web

# With your own data
fastvep-web --gff3 Homo_sapiens.GRCh38.115.gff3 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa

# With supplementary annotations (ClinVar, gnomAD, etc.)
fastvep-web --gff3 genes.gff3 --fasta ref.fa --sa-dir /path/to/sa_databases/
```

Open http://localhost:8080 in your browser. The web interface lets you paste VCF data, switch gene models, and view results in an interactive table.

> **Note:** `fastvep-web` is a separate production-quality binary (axum/tokio, async, multi-connection). The legacy `fastvep web` command still works but is single-threaded.

## Local Setup Guide

This section walks through setting up fastVEP with full annotation capabilities (gene models, reference sequence, and supplementary databases like ClinVar and gnomAD).

### Step 1: Download reference data

```bash
mkdir -p data && cd data

# Gene models (GFF3) — pick your organism
# Human GRCh38
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA (needed for HGVS and sequence context)
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Create FASTA index (enables memory-mapped access — important for large genomes)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Step 2: Build supplementary annotation databases

Download source files and build fastSA databases. Each build is a one-time operation.

```bash
mkdir -p sa_databases

# ClinVar — clinical variant significance (~30MB download)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
fastvep sa-build --source clinvar -i clinvar.vcf.gz -o sa_databases/clinvar --assembly GRCh38

# gnomAD v4 — population allele frequencies (~25GB download per chromosome)
# Download from https://gnomad.broadinstitute.org/downloads
fastvep sa-build --source gnomad -i gnomad.genomes.v4.0.sites.vcf.bgz -o sa_databases/gnomad --assembly GRCh38

# dbSNP — variant identifiers
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
fastvep sa-build --source dbsnp -i GCF_000001405.40.gz -o sa_databases/dbsnp --assembly GRCh38

# COSMIC — somatic mutations (requires license, see https://cancer.sanger.ac.uk/cosmic/download)
fastvep sa-build --source cosmic -i CosmicCodingMuts.vcf.gz -o sa_databases/cosmic --assembly GRCh38
```

### Step 3: Run the CLI annotator

```bash
fastvep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --hgvs
```

### Step 4: Run the web server

```bash
# Install the web server binary
cargo install --path crates/fastvep-web

# Run with all annotation sources
fastvep-web \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --port 8080
```

All flags also accept environment variables (`FASTVEP_GFF3`, `FASTVEP_FASTA`, `FASTVEP_SA_DIR`, `FASTVEP_PORT`) for container deployments.

### Multi-organism setup

To serve multiple genomes from the web interface, organize data into subdirectories and use `--data-dir`. Each subdirectory is one genome — the server auto-detects GFF3, FASTA, and SA files inside.

**Directory layout:**

```
genomes/
  human_grch38/
    Homo_sapiens.GRCh38.115.gff3       # gene models (required)
    Homo_sapiens.GRCh38.dna.primary_assembly.fa   # reference (optional, for HGVS)
    Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
    sa/                                 # supplementary annotations (optional)
      clinvar.osa2
      gnomad.osa2
      dbsnp.osa2
  mouse_grcm39/
    Mus_musculus.GRCm39.115.gff3
    mouse.fa
    mouse.fa.fai
  zebrafish/
    Danio_rerio.GRCz11.115.gff3
```

**Setup:**

```bash
mkdir -p genomes/human_grch38/sa genomes/mouse_grcm39 genomes/zebrafish

# Human: GFF3 + FASTA + SA databases
cp data/Homo_sapiens.GRCh38.115.gff3 genomes/human_grch38/
cp data/Homo_sapiens.GRCh38.dna.primary_assembly.fa* genomes/human_grch38/
cp sa_databases/*.osa2 genomes/human_grch38/sa/   # ClinVar, gnomAD, etc.

# Mouse
wget -O- https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz | gunzip > genomes/mouse_grcm39/Mus_musculus.GRCm39.115.gff3

# Zebrafish
wget -O- https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz | gunzip > genomes/zebrafish/Danio_rerio.GRCz11.115.gff3
```

**Run:**

```bash
fastvep-web --data-dir genomes/
```

Users can switch between organisms from the dropdown in the web UI. When a genome has a `sa/` subdirectory, its SA databases are automatically loaded. The dropdown shows "(FASTA + SA)" labels for genomes that have these resources.

`--sa-dir` is optional — if provided, it serves as a fallback for genomes that don't have their own `sa/` folder. If the directory doesn't exist, the server starts without SA (no error).

fastVEP works with any organism — just provide the matching GFF3 (and optionally FASTA for HGVS).

## Supplementary Annotation Sources

fastVEP supports direct integration with clinical and population databases through its native fastSA binary format. Build once with `fastvep sa-build`, then use `--sa-dir` to annotate:

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

### `fastvep annotate`

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

### `fastvep sa-build`

| Flag | Description | Default |
|------|-------------|---------|
| `--source` | Source type (clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp) | *required* |
| `-i, --input` | Input file (VCF/TSV/wigFix, supports .gz) | *required* |
| `-o, --output` | Output base path (creates .osa and .osa.idx) | *required* |
| `--assembly` | Genome assembly | `GRCh38` |

### `fastvep filter`

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

### `fastvep-web` (production web server)

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 gene annotation file | -- |
| `--fasta` | Reference FASTA file | -- |
| `--sa-dir` | Directory containing .osa/.osa2 supplementary annotation files | -- |
| `--data-dir` | Directory of genome subdirectories (for multi-organism switching) | -- |
| `--port` | HTTP port (also `FASTVEP_PORT` env) | `8080` |
| `--bind` | Bind address (also `FASTVEP_BIND` env) | `0.0.0.0` |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--max-body-size` | Max request body in bytes | `10485760` |
| `--max-concurrent` | Max concurrent annotation requests | `64` |

### `fastvep cache`

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

fastVEP predicts 49 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification, TFBS_ablation, regulatory_region_ablation |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant, regulatory_region_amplification, TFBS_amplification |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, copy_number_change, copy_number_increase, copy_number_decrease, short_tandem_repeat_change, transcript_variant, and others |

## Architecture

```
crates/
  fastvep-core/         # Core types: Consequence (49 SO terms), VariantType, Allele, Impact
  fastvep-genome/       # Transcript, Exon, Gene, CodonTable, mitochondrial codon table
  fastvep-cache/        # GFF3 parser, FASTA reader, annotation providers, regulatory regions
  fastvep-consequence/  # Consequence prediction: small variants + SV predictor
  fastvep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  fastvep-io/           # VCF parser (incl. SVs), output formatters, multi-sample parsing
  fastvep-filter/       # Filter engine: lexer, parser, evaluator (filter_vep-compatible)
  fastvep-sa/           # Supplementary annotation format (fastSA):
                       #   v1 (.osa): zstd block compression, binary search
                       #   v2 (.osa2): echtvar-inspired chunked ZIP with Var32 encoding,
                       #     parallel u32 value arrays, delta encoding, LRU chunk cache,
                       #     Bloom filters for negative lookups
                       # Source parsers: ClinVar, gnomAD, dbSNP, COSMIC, 1000G, TOPMed,
                       # MitoMap, PhyloP, GERP, DANN, REVEL, SpliceAI, PrimateAI, dbNSFP
                       # Custom VCF/BED annotation providers
  fastvep-cli/          # CLI binary: annotation pipeline, sa-build, legacy web server
  fastvep-web/          # Production web server (axum/tokio): async, multi-connection,
                       #   genome switching, SA integration, rate limiting
web/                   # Web GUI (HTML/CSS/JS, embedded in both server binaries)
tests/                 # Test data: chr1 (OR4F5) and chr17 (BRCA1) VCF + GFF3
```

## Running Tests

```bash
cargo test --workspace          # 233 tests
cargo test -p fastvep-consequence  # Consequence prediction tests (incl. SV)
cargo test -p fastvep-filter       # Filter engine tests
cargo test -p fastvep-sa           # Supplementary annotation format tests
```

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), release build with LTO. Median of 3 runs, full Ensembl annotations with FASTA and HGVS.

### Multi-Organism Throughput (Gold-Standard Datasets)

| Organism | Transcripts | Variants | Source | Time | Throughput |
|----------|-------------|----------|--------|------|------------|
| Yeast (R64, full genome) | 7,036 | 260,526 | Ensembl/SGD | 3.0s | **85,934 v/s** |
| Drosophila (BDGP6, full) | 35,442 | 4,438,427 | DGRP2 | 57.3s | **77,486 v/s** |
| Arabidopsis (TAIR10, full) | 54,013 | 12,883,854 | 1001 Genomes | 168.7s | **76,378 v/s** |
| Mouse (GRCm39, full genome) | 142,626 | 26,062,054 | MGP CAST/EiJ | 338.0s | **77,113 v/s** |
| Human full WGS (GRCh38) | 508,530 | 4,048,342 | GIAB HG002 | 86.3s | **46,917 v/s** |

### vs. Ensembl VEP v115.1 (head-to-head, GIAB HG002 chr22)

| Variants | fastVEP | VEP | Speedup |
|----------|---------|-----|---------|
| 1,000 | 0.40s | 1.06s | **2.6x** |
| 5,000 | 0.47s | 13.9s | **29x** |
| 10,000 | 0.67s | 30.3s | **45x** |
| 50,000 | 1.59s | 206.1s | **130x** |
| 4,048,342 (full WGS) | 86.3s | cannot complete | -- |
| Peak memory (100K variants) | ~500 MB | **2.8 MB** |
| Binary size | ~200 MB installed | **3.3 MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None** |

## License

Apache License 2.0

## Related Tools

- **[RastQC](https://doi.org/10.64898/2026.03.31.715630)** — High-performance sequencing quality control written in Rust (Huang, 2026)

## Acknowledgements

fastVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI and [Illumina Nirvana](https://github.com/Illumina/Nirvana). The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework.
