#!/usr/bin/env bash
#
# OxiVEP Multi-Organism Benchmark Suite
#
# Benchmarks OxiVEP annotation performance across model organisms using
# real Ensembl annotations and gold-standard variant call sets.
#
# Data sources:
#   Human     — GIAB HG002 (GRCh38), Ensembl 115
#   Mouse     — Mouse Genomes Project (GRCm39), Ensembl 115
#   Yeast     — 1002 Yeast Genomes (R64), Ensembl 115
#   Drosophila — Ensembl variation (BDGP6), Ensembl 115
#   C. elegans — Ensembl variation (WBcel235), Ensembl 115
#   Arabidopsis — Ensembl Plants variation (TAIR10), Ensembl 115
#
# Prerequisites:
#   ./download_data.sh --all
#
# Each benchmark: 3 runs, median reported, with binary transcript cache.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OXIVEP="$PROJECT_DIR/target/release/oxivep"
OUTPUT_DIR="$SCRIPT_DIR/output"
ORG_DATA="$PROJECT_DIR/test_data/organisms"
HUMAN_DIR="$ORG_DATA/human"
RUNS=3

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

print_header() {
    echo ""
    echo -e "${BOLD}============================================================${NC}"
    echo -e "${BOLD}  OxiVEP Multi-Organism Benchmark Suite${NC}"
    echo -e "${BOLD}============================================================${NC}"
    echo ""
}

print_section() {
    echo ""
    echo -e "${CYAN}--- $1 ---${NC}"
}

build_release() {
    print_section "Building OxiVEP (release mode with LTO)"
    cd "$PROJECT_DIR"
    cargo build --release 2>&1 | tail -1
    if [ ! -f "$OXIVEP" ]; then
        echo -e "${RED}ERROR: Build failed.${NC}"
        exit 1
    fi
    echo -e "${GREEN}Build successful.${NC}"
}

count_variants() {
    grep -c -v '^#' "$1"
}

time_run() {
    local input_vcf="$1"
    local gff3="$2"
    local fmt="$3"
    local outfile="$4"
    local fasta="${5:-}"

    local fasta_args=()
    if [[ -n "$fasta" && -f "$fasta" ]]; then
        fasta_args=(--fasta "$fasta")
    fi

    local start end
    start=$(python3 -c 'import time; print(int(time.time()*1e9))')
    "$OXIVEP" annotate \
        --input "$input_vcf" \
        --gff3 "$gff3" \
        "${fasta_args[@]}" \
        --output "$outfile" \
        --output-format "$fmt" \
        --hgvs \
        2>/dev/null || true
    end=$(python3 -c 'import time; print(int(time.time()*1e9))')

    python3 -c "print(f'{($end - $start) / 1e9:.4f}')"
}

median_time() {
    local n="$1"
    shift
    local times=()
    for ((r=1; r<=n; r++)); do
        times+=( "$(time_run "$@")" )
    done
    printf '%s\n' "${times[@]}" | sort -n | python3 -c "
import sys
vals = [float(l) for l in sys.stdin]
mid = len(vals) // 2
print(f'{vals[mid]:.4f}')
"
}

warm_cache() {
    local gff3="$1"
    local fasta="${2:-}"
    local cache_file

    if [[ "$gff3" == *.gz ]]; then
        cache_file="${gff3%.gz}.oxivep.cache"
    else
        cache_file="${gff3}.oxivep.cache"
    fi

    if [[ -f "$cache_file" ]]; then
        echo -e "  Cache exists: $(basename $cache_file)"
        return 0
    fi

    echo -e "  ${YELLOW}Building cache for $(basename $gff3)...${NC}"
    local tmp_vcf
    tmp_vcf=$(mktemp /tmp/oxivep_warm_XXXXXX.vcf)
    local chrom
    chrom=$(grep -m1 -v '^#' "$gff3" | cut -f1)
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n%s\t1\t.\tA\tG\t.\tPASS\t.\n' "$chrom" > "$tmp_vcf"

    local fasta_args=""
    if [[ -n "$fasta" && -f "$fasta" ]]; then
        fasta_args="--fasta $fasta"
    fi

    "$OXIVEP" annotate --input "$tmp_vcf" --gff3 "$gff3" $fasta_args \
        --output /dev/null --output-format tab 2>/dev/null || true
    rm -f "$tmp_vcf"

    if [[ -f "$cache_file" ]]; then
        echo -e "  ${GREEN}Cache built: $(basename $cache_file)${NC}"
    else
        echo -e "  ${RED}Warning: cache not created${NC}"
    fi
}

# ═══════════════════════════════════════════════════════════════
# Register benchmarks
# ═══════════════════════════════════════════════════════════════

declare -a ALL_NAMES=()
declare -a ALL_VCFS=()
declare -a ALL_GFF3S=()
declare -a ALL_FASTAS=()
declare -a ALL_ORGANISMS=()

add_benchmark() {
    local name="$1" vcf="$2" gff3="$3" fasta="${4:-}" organism="$5"
    if [[ ! -f "$vcf" ]]; then
        echo -e "  ${YELLOW}Skipping $name: $(basename "$vcf") not found${NC}"
        return 0
    fi
    if [[ ! -f "$gff3" ]]; then
        echo -e "  ${YELLOW}Skipping $name: $(basename "$gff3") not found${NC}"
        return 0
    fi
    ALL_NAMES+=("$name")
    ALL_VCFS+=("$vcf")
    ALL_GFF3S+=("$gff3")
    ALL_FASTAS+=("$fasta")
    ALL_ORGANISMS+=("$organism")
}

run_benchmarks() {
    mkdir -p "$OUTPUT_DIR"

    # ── Yeast (R64, full genome) ──
    if [[ -f "$ORG_DATA/yeast.gff3" ]]; then
        print_section "Pre-warming: Yeast (R64)"
        warm_cache "$ORG_DATA/yeast.gff3" "$ORG_DATA/yeast.fa"
        add_benchmark "yeast_260k" "$ORG_DATA/yeast_ensembl.vcf" "$ORG_DATA/yeast.gff3" "$ORG_DATA/yeast.fa" "Yeast"
    fi

    # ── Drosophila (BDGP6, full genome) ──
    if [[ -f "$ORG_DATA/drosophila.gff3" ]]; then
        print_section "Pre-warming: Drosophila (BDGP6)"
        warm_cache "$ORG_DATA/drosophila.gff3" "$ORG_DATA/drosophila.fa"
        add_benchmark "drosophila_100k" "$ORG_DATA/drosophila_100k.vcf" "$ORG_DATA/drosophila.gff3" "$ORG_DATA/drosophila.fa" "Drosophila"
    fi

    # ── Arabidopsis (TAIR10, full genome) ──
    if [[ -f "$ORG_DATA/arabidopsis.gff3" ]]; then
        print_section "Pre-warming: Arabidopsis (TAIR10)"
        warm_cache "$ORG_DATA/arabidopsis.gff3" "$ORG_DATA/arabidopsis.fa"
        add_benchmark "arabidopsis_100k" "$ORG_DATA/arabidopsis_100k.vcf" "$ORG_DATA/arabidopsis.gff3" "$ORG_DATA/arabidopsis.fa" "Arabidopsis"
        add_benchmark "arabidopsis_500k" "$ORG_DATA/arabidopsis_500k.vcf" "$ORG_DATA/arabidopsis.gff3" "$ORG_DATA/arabidopsis.fa" "Arabidopsis"
    fi

    # ── Mouse (GRCm39, full genome) ──
    if [[ -f "$ORG_DATA/mouse.gff3" ]]; then
        print_section "Pre-warming: Mouse (GRCm39)"
        warm_cache "$ORG_DATA/mouse.gff3" "$ORG_DATA/mouse.fa"
        add_benchmark "mouse_100k" "$ORG_DATA/mouse_100k.vcf" "$ORG_DATA/mouse.gff3" "$ORG_DATA/mouse.fa" "Mouse"
        add_benchmark "mouse_500k" "$ORG_DATA/mouse_500k.vcf" "$ORG_DATA/mouse.gff3" "$ORG_DATA/mouse.fa" "Mouse"
    fi

    # ── Human (GRCh38, full genome) ──
    local HUMAN_GFF3="$HUMAN_DIR/Homo_sapiens.GRCh38.115.gff3"
    local HUMAN_FA="$HUMAN_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    if [[ -f "$HUMAN_GFF3" ]]; then
        print_section "Pre-warming: Human full genome (GRCh38)"
        warm_cache "$HUMAN_GFF3" "$HUMAN_FA"
        add_benchmark "human_hg002_100k" "$HUMAN_DIR/human_hg002_100k.vcf" "$HUMAN_GFF3" "$HUMAN_FA" "Human"
        add_benchmark "human_hg002_500k" "$HUMAN_DIR/human_hg002_500k.vcf" "$HUMAN_GFF3" "$HUMAN_FA" "Human"
        add_benchmark "human_hg002_full" "$HUMAN_DIR/human_giab_hg002_full.vcf" "$HUMAN_GFF3" "$HUMAN_FA" "Human"
    fi

    # ═══════════════════════════════════════════════════════════════
    # RUN ALL BENCHMARKS
    # ═══════════════════════════════════════════════════════════════
    if [[ ${#ALL_NAMES[@]} -eq 0 ]]; then
        echo -e "\n${RED}No benchmark data found. Run ./download_data.sh --all first.${NC}"
        exit 1
    fi

    print_section "Running benchmarks (median of $RUNS runs each)"

    declare -A RESULTS

    for i in "${!ALL_NAMES[@]}"; do
        local name="${ALL_NAMES[$i]}"
        local vcf="${ALL_VCFS[$i]}"
        local gff="${ALL_GFF3S[$i]}"
        local fasta="${ALL_FASTAS[$i]:-}"
        local organism="${ALL_ORGANISMS[$i]}"
        local nvar
        nvar=$(count_variants "$vcf")

        printf "  %-25s (%8s variants, %-12s) ... " "$name" "$nvar" "$organism"
        local elapsed
        elapsed=$(median_time "$RUNS" "$vcf" "$gff" "vcf" "$OUTPUT_DIR/${name}.annotated.vcf" "$fasta")

        local vps
        vps=$(python3 -c "
e = float('$elapsed')
n = int('$nvar')
print(f'{n/e:,.0f}' if e > 0 else 'inf')
")
        printf "${GREEN}%8s sec${NC}  (%s v/s)\n" "$elapsed" "$vps"
        RESULTS["${name}"]="${elapsed} ${nvar} ${vps} ${organism}"
    done

    # ═══════════════════════════════════════════════════════════════
    # SUMMARY TABLE
    # ═══════════════════════════════════════════════════════════════
    echo ""
    echo -e "${BOLD}============================================================${NC}"
    echo -e "${BOLD}  Benchmark Summary (median of ${RUNS} runs, VCF output)${NC}"
    echo -e "${BOLD}============================================================${NC}"
    echo ""
    printf "${BOLD}%-25s %-12s %10s %10s %15s${NC}\n" "Dataset" "Organism" "Variants" "Time (s)" "Variants/sec"
    printf "%-25s %-12s %10s %10s %15s\n" "-------------------------" "------------" "----------" "----------" "---------------"

    for name in "${ALL_NAMES[@]}"; do
        if [ -n "${RESULTS[$name]+x}" ]; then
            local parts=(${RESULTS[$name]})
            printf "%-25s %-12s %10s %10s %15s\n" "$name" "${parts[3]}" "${parts[1]}" "${parts[0]}" "${parts[2]}"
        fi
    done

    # ═══════════════════════════════════════════════════════════════
    # CSV OUTPUT
    # ═══════════════════════════════════════════════════════════════
    local csv_file="$OUTPUT_DIR/benchmark_results.csv"
    echo "dataset,organism,variants,time_seconds,variants_per_second" > "$csv_file"
    for name in "${ALL_NAMES[@]}"; do
        if [ -n "${RESULTS[$name]+x}" ]; then
            local parts=(${RESULTS[$name]})
            echo "${name},${parts[3]},${parts[1]},${parts[0]},${parts[2]}" >> "$csv_file"
        fi
    done

    echo ""
    echo -e "${GREEN}CSV results:${NC} $csv_file"
    echo -e "${BOLD}Output files:${NC} $OUTPUT_DIR/"
}

# ═══════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════
print_header

echo "Date:     $(date '+%Y-%m-%d %H:%M:%S')"
echo "System:   $(uname -s) $(uname -m)"
echo "Rust:     $(rustc --version 2>/dev/null || echo 'not found')"
echo "Runs:     $RUNS per benchmark (median reported)"
echo "Data:     $VCF_DATA/"
echo ""

build_release
run_benchmarks

echo ""
echo -e "${GREEN}Benchmark complete.${NC}"
