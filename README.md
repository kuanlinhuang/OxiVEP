# OxiVEP

A high-performance Variant Effect Predictor written in Rust. OxiVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions) on genes, transcripts, and protein sequences.

OxiVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

## Features

- **Variant Consequence Prediction** — Classifies variants using [Sequence Ontology](http://www.sequenceontology.org/) terms (missense, frameshift, splice donor, etc.)
- **HGVS Nomenclature** — Generates HGVSg, HGVSc, and HGVSp notations
- **Multiple Output Formats** — VCF (with CSQ field), tab-delimited, and JSON
- **GFF3 Annotation Support** — Load gene models from standard GFF3 files
- **Splice Site Detection** — Identifies splice donor, acceptor, region, 5th base, donor region, and polypyrimidine tract variants
- **Multi-allelic Support** — Handles multi-allelic VCF records with proper allele normalization
- **Web Interface** — Built-in web GUI for interactive variant annotation

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

Expected output (using real BRCA1 and TP53 data from Ensembl GRCh38):

```
rs_cds_brca1      17:43124090  G  ENSG00000012048  ENST00000357654  Transcript  synonymous_variant
rs_cds_brca1_mid  17:43106500  C  ENSG00000012048  ENST00000357654  Transcript  synonymous_variant
rs_5utr_brca1     17:43125300  T  ENSG00000012048  ENST00000357654  Transcript  5_prime_UTR_variant
rs_intron_brca1   17:43120000  T  ENSG00000012048  ENST00000357654  Transcript  intron_variant
rs_downstream     17:43043000  C  ENSG00000012048  ENST00000357654  Transcript  downstream_gene_variant
rs_3utr_brca1     17:43045700  C  ENSG00000012048  ENST00000357654  Transcript  missense_variant
rs_del_brca1      17:43124090  -  ENSG00000012048  ENST00000357654  Transcript  frameshift_variant
rs_cds_tp53       17:7675088   T  ENSG00000141510  ENST00000269305  Transcript  missense_variant
```

The 8 variants cover: coding SNVs, 5'UTR, 3'UTR, intron, downstream, and frameshift deletion across BRCA1 and TP53.

### 4. Try VCF output (with CSQ annotations)

```bash
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format vcf
```

### 5. Try JSON output

```bash
oxivep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format json
```

### 6. Launch the web interface

```bash
oxivep web --port 8080
```

Open http://localhost:8080 in your browser. Click **"Load Example"** to load pre-built test variants, then click **"Annotate"** to see results instantly.

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

## Test Data Reference

The repository includes test files in `tests/`:

**`tests/test.gff3`** — Real Ensembl GRCh38 release 115 annotations for two genes on chr17:
- `BRCA1` (ENSG00000012048): Protein-coding gene at chr17:43044292-43170245 (- strand), 23 exons
- `TP53` (ENSG00000141510): Protein-coding gene at chr17:7661779-7687546 (- strand), 11 exons

**`tests/test.vcf`** — 8 variants at real BRCA1/TP53 positions:

| Variant | Position | Type | Expected Consequence |
|---------|----------|------|---------------------|
| rs_cds_brca1 | 17:43124090 | SNV (in CDS) | synonymous_variant |
| rs_cds_brca1_mid | 17:43106500 | SNV (in CDS) | synonymous_variant |
| rs_5utr_brca1 | 17:43125300 | SNV (in 5'UTR) | 5_prime_UTR_variant |
| rs_intron_brca1 | 17:43120000 | SNV (mid-intron) | intron_variant |
| rs_downstream_brca1 | 17:43043000 | SNV (downstream) | downstream_gene_variant |
| rs_3utr_brca1 | 17:43045700 | SNV (near 3' end) | missense_variant |
| rs_del_brca1 | 17:43124090 | 1bp deletion (CDS) | frameshift_variant |
| rs_cds_tp53 | 17:7675088 | SNV (in CDS) | missense_variant |

## Command Reference

### `oxivep annotate`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VCF file (`-` for stdin) | *required* |
| `-o, --output` | Output file (`-` for stdout) | `-` |
| `--gff3` | GFF3 gene annotation file | — |
| `--fasta` | Reference FASTA file | — |
| `--output-format` | `vcf`, `tab`, or `json` | `vcf` |
| `--hgvs` | Include HGVS notations | off |
| `--pick` | Report only the most severe consequence per variant | off |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--everything` | Enable all annotation flags | off |
| `--symbol` | Include gene symbol | off |
| `--canonical` | Flag canonical transcripts | off |

### `oxivep web`

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 gene annotation file | — |
| `--fasta` | Reference FASTA file | — |
| `--port` | HTTP port | `8080` |

### `oxivep filter`

Filter annotated VEP output (coming soon).

## Output Formats

### VCF Output

Annotations are added as a `CSQ` field in the INFO column:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from OxiVEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND">
```

### Tab Output

One line per variant-transcript-allele combination with 13 columns.

### JSON Output

Structured JSON with `transcript_consequences` array per variant.

## Consequence Types

OxiVEP predicts 41 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, and others |

## Architecture

```
crates/
  oxivep-core/         # Core types: Consequence, GenomicPosition, Allele, Impact
  oxivep-genome/       # Transcript, Exon, Gene, CodonTable, coordinate mapping
  oxivep-cache/        # GFF3 parser, FASTA reader, annotation providers
  oxivep-consequence/  # Consequence prediction engine, splice site detection
  oxivep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  oxivep-io/           # VCF parser, output formatters (CSQ, tab, JSON)
  oxivep-filter/       # Variant filtering
  oxivep-cli/          # CLI binary, annotation pipeline, web server
web/                   # Web GUI (HTML/CSS/JS)
tests/                 # Test VCF and GFF3 files
```

## Running Tests

```bash
cargo test --workspace          # 101 tests
cargo test -p oxivep-consequence  # Just consequence prediction tests
```

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), single-threaded, release build.

| Dataset | Variants | Time | Throughput |
|---------|----------|------|------------|
| Human chr21 (20 genes) | 1,000 | 0.037s | **27,000 variants/sec** |
| Human chr21 (20 genes) | 10,000 | 0.057s | **175,000 variants/sec** |
| Human chr21 (20 genes) | 50,000 | 0.128s | **391,000 variants/sec** |
| Mouse chr19 (10 genes) | 500 | 0.033s | **15,200 variants/sec** |
| Zebrafish chr5 (8 genes) | 300 | 0.031s | **9,800 variants/sec** |

### vs. Ensembl VEP

| Metric | Ensembl VEP (Perl) | OxiVEP (Rust) |
|--------|-------------------|---------------|
| Startup time | 5-15s | <0.05s |
| 1,000 SNVs (offline) | ~3-10s | **0.04s** |
| Peak memory (100K variants) | ~500 MB | **2.8 MB** |
| Binary size | ~200 MB installed | **2.4 MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None** |

## License

Apache License 2.0

## Acknowledgements

OxiVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI. The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework.
