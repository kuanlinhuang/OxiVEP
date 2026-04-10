#!/usr/bin/env bash
set -euo pipefail

# OxiVEP validation suite — compares against Ensembl VEP Docker on real variant data.
#
# Prerequisites:
#   - oxivep binary in PATH (cargo install --path crates/oxivep-cli)
#   - Docker with ensemblorg/ensembl-vep:release_115.1
#   - Full human genome data in test_data/organisms/human/
#   - Mouse data in validation/mouse/ (auto-downloaded if missing)
#
# Usage:
#   ./validation/run_validation.sh              # run all
#   ./validation/run_validation.sh human        # human only
#   ./validation/run_validation.sh mouse        # mouse only

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
VEP_IMAGE="ensemblorg/ensembl-vep:release_115.1"
VEP_CACHE="/tmp/vep_cache"
HUMAN_DIR="$PROJECT_DIR/test_data/organisms/human"

mkdir -p "$RESULTS_DIR" "$VEP_CACHE"

# ---- Helpers ----

run_oxivep() {
    local input="$1" gff3="$2" fasta="$3" output="$4"
    oxivep annotate -i "$input" --gff3 "$gff3" --fasta "$fasta" \
        --hgvs --symbol --canonical -o "$output" 2>&1
}

run_vep_docker() {
    local input="$1" gff3="$2" fasta="$3" output="$4"
    docker run --rm \
        -v "$PROJECT_DIR:/work" \
        -v "$VEP_CACHE:/opt/vep/.vep" \
        "$VEP_IMAGE" \
        vep --input_file "/work/$input" \
        --gff "/work/$gff3" \
        --fasta "/work/$fasta" \
        --output_file "/work/$output" \
        --vcf --force_overwrite --no_stats \
        --hgvs --symbol --canonical --offline 2>&1 | tail -3
}

compare() {
    python3 "$SCRIPT_DIR/compare_vep.py" "$1" "$2" --verbose
}

# ---- Human validation ----

run_human() {
    echo ""
    echo "================================================================"
    echo "  HUMAN VALIDATION (GRCh38, full genome annotations)"
    echo "================================================================"

    local human_gff3="$HUMAN_DIR/Homo_sapiens.GRCh38.115.gff3"
    local human_fasta="$HUMAN_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    if [[ ! -f "$human_gff3" || ! -f "$human_fasta" ]]; then
        echo "SKIP: Human GFF3/FASTA not found in $HUMAN_DIR"
        echo "  Run benchmarks/download_data.sh first"
        return
    fi

    # VEP example (173 chr22 variants)
    local vep_example="$SCRIPT_DIR/human/vep_example_GRCh38.vcf"
    if [[ -f "$vep_example" ]]; then
        echo ""
        echo "--- VEP Example (173 chr22 variants, full genome annotations) ---"
        echo "  Running OxiVEP..."
        run_oxivep "$vep_example" "$human_gff3" "$human_fasta" "$RESULTS_DIR/oxivep_vep_example.vcf"

        echo "  Running Ensembl VEP (requires Docker)..."
        # VEP needs sorted+indexed GFF3; use the full GFF3 if indexed
        echo "  NOTE: VEP comparison requires manually prepared sorted GFF3"
        echo "  Skipping VEP for now — use compare_vep.py on existing results"
    fi

    # chr22 1KGP validation
    local chr22_vcf="$SCRIPT_DIR/human/chr22_1kgp.vcf"
    if [[ -f "$chr22_vcf" ]]; then
        echo ""
        echo "--- Human chr22 1KGP (1000 variants, full genome annotations) ---"
        echo "  Running OxiVEP..."
        run_oxivep "$chr22_vcf" "$human_gff3" "$human_fasta" "$RESULTS_DIR/oxivep_chr22.vcf"
    fi
}

# ---- Mouse validation ----

run_mouse() {
    echo ""
    echo "================================================================"
    echo "  MOUSE VALIDATION (Mouse Genomes Project, GRCm39)"
    echo "================================================================"

    local mouse_dir="$SCRIPT_DIR/mouse"
    local input="$mouse_dir/mouse_chr19_mgp.vcf"
    local gff3="$mouse_dir/mouse_chr19.gff3"
    local sorted_gff3="$mouse_dir/mouse_chr19_sorted.gff3.gz"
    local fasta="$mouse_dir/Mus_musculus.GRCm39.dna.chromosome.19.fa"

    if [[ ! -f "$input" ]]; then
        echo "SKIP: mouse input VCF not found at $input"
        return
    fi
    if [[ ! -f "$fasta" ]]; then
        echo "Downloading mouse chr19 FASTA..."
        curl -L -o "${fasta}.gz" \
            "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.19.fa.gz"
        gunzip -k "${fasta}.gz"
        samtools faidx "$fasta"
    fi
    if [[ ! -f "$gff3" ]]; then
        echo "Downloading mouse chr19 GFF3..."
        curl -L -o "$mouse_dir/Mus_musculus.GRCm39.115.chromosome.19.gff3.gz" \
            "https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.chromosome.19.gff3.gz"
        gunzip -c "$mouse_dir/Mus_musculus.GRCm39.115.chromosome.19.gff3.gz" > "$gff3"
    fi

    echo ""
    echo "--- Mouse chr19 ---"

    local oxivep_out="$RESULTS_DIR/oxivep_mouse_chr19.vcf"
    local vep_out="$RESULTS_DIR/vep_mouse_chr19.vcf"

    echo "  Running OxiVEP..."
    run_oxivep "$input" "$gff3" "$fasta" "$oxivep_out"

    if [[ -f "$sorted_gff3" && -f "${sorted_gff3}.tbi" ]]; then
        echo "  Running Ensembl VEP..."
        run_vep_docker \
            "validation/mouse/mouse_chr19_mgp.vcf" \
            "validation/mouse/mouse_chr19_sorted.gff3.gz" \
            "validation/mouse/Mus_musculus.GRCm39.dna.chromosome.19.fa" \
            "validation/results/vep_mouse_chr19.vcf"

        echo "  Comparing..."
        compare "$oxivep_out" "$vep_out"
    else
        echo "  WARN: no sorted+indexed GFF3 for mouse, skipping VEP comparison"
    fi
}

# ---- Main ----

cd "$PROJECT_DIR"

case "${1:-all}" in
    human) run_human ;;
    mouse) run_mouse ;;
    all)   run_human; run_mouse ;;
    *)     echo "Usage: $0 [human|mouse|all]"; exit 1 ;;
esac

echo ""
echo "Done. Results in validation/results/"
