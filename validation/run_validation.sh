#!/usr/bin/env bash
set -euo pipefail

# OxiVEP validation suite — compares against Ensembl VEP Docker on real variant data.
#
# Prerequisites:
#   - oxivep binary in PATH (cargo install --path crates/oxivep-cli)
#   - Docker with ensemblorg/ensembl-vep:release_115.1
#   - Human GFF3/FASTA in test_data/ (chr7, chr11, chr17, chr22)
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

mkdir -p "$RESULTS_DIR" "$VEP_CACHE/homo_sapiens/115_GRCh38" "$VEP_CACHE/mus_musculus/115_GRCm39"

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
    echo "  HUMAN VALIDATION (1000 Genomes Project, GRCh38)"
    echo "================================================================"

    # Also run VEP's own example (173 chr22 variants)
    local vep_example="validation/human/vep_example_GRCh38.vcf"
    if [[ -f "$PROJECT_DIR/$vep_example" ]]; then
        echo ""
        echo "--- VEP Example (173 chr22 variants) ---"
        echo "  Running OxiVEP..."
        run_oxivep "$PROJECT_DIR/$vep_example" "$PROJECT_DIR/test_data/chr22.gff3" \
            "$PROJECT_DIR/test_data/chr22.fa" "$PROJECT_DIR/validation/results/oxivep_vep_example.vcf"
        local sorted_gff22=""
        for c in "test_data/chr22_sorted.gff3.gz"; do
            if [[ -f "$PROJECT_DIR/$c" && -f "$PROJECT_DIR/${c}.tbi" ]]; then sorted_gff22="$c"; break; fi
        done
        if [[ -n "$sorted_gff22" ]]; then
            echo "  Running Ensembl VEP..."
            run_vep_docker "$vep_example" "$sorted_gff22" "test_data/chr22.fa" "validation/results/vep_vep_example.vcf"
            echo "  Comparing..."
            compare "$PROJECT_DIR/validation/results/oxivep_vep_example.vcf" "$PROJECT_DIR/validation/results/vep_vep_example.vcf"
        fi
    fi

    for chr_num in 7 11 17 22; do
        local input="validation/human/chr${chr_num}_1kgp.vcf"
        local oxivep_out="validation/results/oxivep_chr${chr_num}.vcf"
        local vep_out="validation/results/vep_chr${chr_num}.vcf"

        if [[ ! -f "$PROJECT_DIR/$input" ]]; then
            echo "SKIP chr${chr_num}: input $input not found"
            continue
        fi

        # Find the GFF3/FASTA
        local gff3="test_data/chr${chr_num}.gff3"
        local fasta="test_data/chr${chr_num}.fa"
        if [[ ! -f "$PROJECT_DIR/$gff3" ]]; then
            echo "SKIP chr${chr_num}: $gff3 not found"
            continue
        fi

        echo ""
        echo "--- Human chr${chr_num} ---"

        # Find sorted GFF3 for VEP (needs tabix index)
        local sorted_gff3=""
        for candidate in "test_data/chr${chr_num}_sorted.gff3.gz" "validation/chr${chr_num}_sorted.gff3.gz"; do
            if [[ -f "$PROJECT_DIR/$candidate" && -f "$PROJECT_DIR/${candidate}.tbi" ]]; then
                sorted_gff3="$candidate"
                break
            fi
        done
        if [[ -z "$sorted_gff3" ]]; then
            echo "WARN: no sorted+indexed GFF3 for chr${chr_num}, skipping VEP"
            continue
        fi

        echo "  Running OxiVEP..."
        run_oxivep "$PROJECT_DIR/$input" "$PROJECT_DIR/$gff3" "$PROJECT_DIR/$fasta" "$PROJECT_DIR/$oxivep_out"

        echo "  Running Ensembl VEP..."
        run_vep_docker "$input" "$sorted_gff3" "$fasta" "$vep_out"

        echo "  Comparing..."
        compare "$PROJECT_DIR/$oxivep_out" "$PROJECT_DIR/$vep_out"
    done
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
        echo "  Download with: bcftools view -G -r 19:3000000-3500000,... mgp_REL2021_snps.vcf.gz"
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
