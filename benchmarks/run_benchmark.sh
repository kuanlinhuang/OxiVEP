#!/usr/bin/env bash
#
# OxiVEP Multi-Organism Benchmark Suite
#
# Benchmarks OxiVEP annotation performance across 7 model organisms using
# real Ensembl annotations and gold-standard/Ensembl variation VCF data.
#
# Organisms:
#   Human (GRCh38 chr22)    - GIAB NA12878 + scaled Ensembl variation
#   Mouse (GRCm39 full)     - Ensembl variation database
#   Zebrafish (GRCz11 chr5) - Generated benchmark variants
#   Drosophila (BDGP6 full) - Generated benchmark variants
#   C. elegans (WBcel235)   - Generated benchmark variants
#   Arabidopsis (TAIR10)    - Generated benchmark variants
#   Yeast (R64 full)        - Ensembl variation database
#
# Each benchmark: 3 runs, median reported, with binary transcript cache.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OXIVEP="$PROJECT_DIR/target/release/oxivep"
OUTPUT_DIR="$SCRIPT_DIR/output"
TEST_DATA="$PROJECT_DIR/test_data"
ORG_DATA="$TEST_DATA/organisms"
RUNS=3

# Colors
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
        --symbol --hgvs --canonical \
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

    # Determine cache path
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
    # Get first chromosome from the GFF3
    local chrom
    chrom=$(grep -m1 -v '^#' "$gff3" | cut -f1)
    echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n${chrom}\t1\t.\tA\tG\t.\tPASS\t." > "$tmp_vcf"

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

# Main benchmark logic
run_benchmarks() {
    mkdir -p "$OUTPUT_DIR"

    # Arrays for collecting all benchmarks
    declare -a ALL_NAMES=()
    declare -a ALL_VCFS=()
    declare -a ALL_GFF3S=()
    declare -a ALL_FASTAS=()
    declare -a ALL_ORGANISMS=()

    # ═══════════════════════════════════════════════════════════════
    # HUMAN (GRCh38 chr22) - GIAB + scaled Ensembl variation
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$TEST_DATA/chr22.gff3" && -f "$TEST_DATA/chr22.fa" ]]; then
        print_section "Pre-warming: Human chr22 (GRCh38)"
        warm_cache "$TEST_DATA/chr22.gff3" "$TEST_DATA/chr22.fa"

        # Scaled real chr22 datasets
        for size_vcf in "$SCRIPT_DIR"/real_chr22_*.vcf; do
            [[ -f "$size_vcf" ]] || continue
            local bn=$(basename "$size_vcf" .vcf)
            ALL_NAMES+=("human_$bn")
            ALL_VCFS+=("$size_vcf")
            ALL_GFF3S+=("$TEST_DATA/chr22.gff3")
            ALL_FASTAS+=("$TEST_DATA/chr22.fa")
            ALL_ORGANISMS+=("Human")
        done

        # 1KGP chr22
        if [[ -f "$SCRIPT_DIR/human_chr22_1kgp.vcf" ]]; then
            ALL_NAMES+=("human_1kgp_1k")
            ALL_VCFS+=("$SCRIPT_DIR/human_chr22_1kgp.vcf")
            ALL_GFF3S+=("$TEST_DATA/chr22.gff3")
            ALL_FASTAS+=("$TEST_DATA/chr22.fa")
            ALL_ORGANISMS+=("Human")
        fi

        # GIAB NA12878 chr22
        if [[ -f "$SCRIPT_DIR/giab/NA12878_GRCh38_chr22.vcf" ]]; then
            ALL_NAMES+=("human_giab_chr22")
            ALL_VCFS+=("$SCRIPT_DIR/giab/NA12878_GRCh38_chr22.vcf")
            ALL_GFF3S+=("$TEST_DATA/chr22.gff3")
            ALL_FASTAS+=("$TEST_DATA/chr22.fa")
            ALL_ORGANISMS+=("Human")
        fi
    fi

    # ═══════════════════════════════════════════════════════════════
    # MOUSE (GRCm39 full genome)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$ORG_DATA/mouse.gff3" && -f "$ORG_DATA/mouse.fa" ]]; then
        print_section "Pre-warming: Mouse full genome (GRCm39)"
        warm_cache "$ORG_DATA/mouse.gff3" "$ORG_DATA/mouse.fa"

        for sz in 10k 100k 500k; do
            if [[ -f "$ORG_DATA/mouse_real_${sz}.vcf" ]]; then
                ALL_NAMES+=("mouse_${sz}")
                ALL_VCFS+=("$ORG_DATA/mouse_real_${sz}.vcf")
                ALL_GFF3S+=("$ORG_DATA/mouse.gff3")
                ALL_FASTAS+=("$ORG_DATA/mouse.fa")
                ALL_ORGANISMS+=("Mouse")
            fi
        done
    elif [[ -f "$TEST_DATA/mouse_chr19.gff3" && -f "$TEST_DATA/mouse_chr19.fa" ]]; then
        # Fallback: mouse chr19 only
        print_section "Pre-warming: Mouse chr19 (GRCm39)"
        warm_cache "$TEST_DATA/mouse_chr19.gff3" "$TEST_DATA/mouse_chr19.fa"

        for sz in 10k 100k 500k; do
            if [[ -f "$ORG_DATA/mouse_real_${sz}.vcf" ]]; then
                ALL_NAMES+=("mouse_chr19_${sz}")
                ALL_VCFS+=("$ORG_DATA/mouse_real_${sz}.vcf")
                ALL_GFF3S+=("$TEST_DATA/mouse_chr19.gff3")
                ALL_FASTAS+=("$TEST_DATA/mouse_chr19.fa")
                ALL_ORGANISMS+=("Mouse")
            fi
        done
    fi

    # ═══════════════════════════════════════════════════════════════
    # ZEBRAFISH (GRCz11 chr5)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$TEST_DATA/zebrafish_chr5.gff3" && -f "$TEST_DATA/zebrafish_chr5.fa" ]]; then
        print_section "Pre-warming: Zebrafish chr5 (GRCz11)"
        warm_cache "$TEST_DATA/zebrafish_chr5.gff3" "$TEST_DATA/zebrafish_chr5.fa"
    fi

    # ═══════════════════════════════════════════════════════════════
    # YEAST (S. cerevisiae R64 full genome)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$ORG_DATA/yeast.gff3" && -f "$ORG_DATA/yeast.fa" ]]; then
        print_section "Pre-warming: Yeast full genome (R64)"
        warm_cache "$ORG_DATA/yeast.gff3" "$ORG_DATA/yeast.fa"

        for sz in 10k 50k; do
            if [[ -f "$ORG_DATA/yeast_real_${sz}.vcf" ]]; then
                ALL_NAMES+=("yeast_real_${sz}")
                ALL_VCFS+=("$ORG_DATA/yeast_real_${sz}.vcf")
                ALL_GFF3S+=("$ORG_DATA/yeast.gff3")
                ALL_FASTAS+=("$ORG_DATA/yeast.fa")
                ALL_ORGANISMS+=("Yeast")
            fi
        done
        if [[ -f "$ORG_DATA/yeast_real_all.vcf" ]]; then
            ALL_NAMES+=("yeast_real_all")
            ALL_VCFS+=("$ORG_DATA/yeast_real_all.vcf")
            ALL_GFF3S+=("$ORG_DATA/yeast.gff3")
            ALL_FASTAS+=("$ORG_DATA/yeast.fa")
            ALL_ORGANISMS+=("Yeast")
        fi
    fi

    # ═══════════════════════════════════════════════════════════════
    # C. ELEGANS (WBcel235 full genome)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$ORG_DATA/elegans.gff3" && -f "$ORG_DATA/elegans.fa" ]]; then
        print_section "Pre-warming: C. elegans full genome (WBcel235)"
        warm_cache "$ORG_DATA/elegans.gff3" "$ORG_DATA/elegans.fa"

        for sz in 10k 50k 100k; do
            if [[ -f "$ORG_DATA/elegans_${sz}.vcf" ]]; then
                ALL_NAMES+=("elegans_${sz}")
                ALL_VCFS+=("$ORG_DATA/elegans_${sz}.vcf")
                ALL_GFF3S+=("$ORG_DATA/elegans.gff3")
                ALL_FASTAS+=("$ORG_DATA/elegans.fa")
                ALL_ORGANISMS+=("C.elegans")
            fi
        done
    fi

    # ═══════════════════════════════════════════════════════════════
    # DROSOPHILA (BDGP6 full genome)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$ORG_DATA/drosophila.gff3" && -f "$ORG_DATA/drosophila.fa" ]]; then
        print_section "Pre-warming: Drosophila full genome (BDGP6)"
        warm_cache "$ORG_DATA/drosophila.gff3" "$ORG_DATA/drosophila.fa"

        for sz in 10k 50k 100k; do
            if [[ -f "$ORG_DATA/drosophila_${sz}.vcf" ]]; then
                ALL_NAMES+=("drosophila_${sz}")
                ALL_VCFS+=("$ORG_DATA/drosophila_${sz}.vcf")
                ALL_GFF3S+=("$ORG_DATA/drosophila.gff3")
                ALL_FASTAS+=("$ORG_DATA/drosophila.fa")
                ALL_ORGANISMS+=("Drosophila")
            fi
        done
    fi

    # ═══════════════════════════════════════════════════════════════
    # ARABIDOPSIS (TAIR10 full genome)
    # ═══════════════════════════════════════════════════════════════
    if [[ -f "$ORG_DATA/arabidopsis.gff3" && -f "$ORG_DATA/arabidopsis.fa" ]]; then
        print_section "Pre-warming: Arabidopsis full genome (TAIR10)"
        warm_cache "$ORG_DATA/arabidopsis.gff3" "$ORG_DATA/arabidopsis.fa"

        for sz in 10k 100k 500k; do
            if [[ -f "$ORG_DATA/arabidopsis_real_${sz}.vcf" ]]; then
                ALL_NAMES+=("arabidopsis_${sz}")
                ALL_VCFS+=("$ORG_DATA/arabidopsis_real_${sz}.vcf")
                ALL_GFF3S+=("$ORG_DATA/arabidopsis.gff3")
                ALL_FASTAS+=("$ORG_DATA/arabidopsis.fa")
                ALL_ORGANISMS+=("Arabidopsis")
            fi
        done
    fi

    # ═══════════════════════════════════════════════════════════════
    # RUN ALL BENCHMARKS (VCF format only for speed)
    # ═══════════════════════════════════════════════════════════════
    declare -A RESULTS
    declare -A VARIANT_COUNTS

    for i in "${!ALL_NAMES[@]}"; do
        local name="${ALL_NAMES[$i]}"
        local vcf="${ALL_VCFS[$i]}"
        local gff="${ALL_GFF3S[$i]}"
        local fasta="${ALL_FASTAS[$i]:-}"
        local organism="${ALL_ORGANISMS[$i]}"
        local nvar
        nvar=$(count_variants "$vcf")
        VARIANT_COUNTS["$name"]="$nvar"

        printf "  %-30s (%8s variants, %s) ... " "$name" "$nvar" "$organism"
        local elapsed
        elapsed=$(median_time "$RUNS" "$vcf" "$gff" "vcf" "$OUTPUT_DIR/${name}.annotated.vcf" "$fasta")

        local vps
        vps=$(python3 -c "
e = float('$elapsed')
n = int('$nvar')
print(f'{n/e:.1f}' if e > 0 else 'inf')
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
    printf "${BOLD}%-30s %-12s %10s %10s %15s${NC}\n" "Dataset" "Organism" "Variants" "Time (s)" "Variants/sec"
    printf "%-30s %-12s %10s %10s %15s\n"   "------------------------------" "------------" "----------" "----------" "---------------"

    for name in "${ALL_NAMES[@]}"; do
        if [ -n "${RESULTS[$name]+x}" ]; then
            local parts=(${RESULTS[$name]})
            printf "%-30s %-12s %10s %10s %15s\n" "$name" "${parts[3]}" "${parts[1]}" "${parts[0]}" "${parts[2]}"
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
echo ""

build_release
run_benchmarks

echo ""
echo -e "${GREEN}Benchmark complete.${NC}"
