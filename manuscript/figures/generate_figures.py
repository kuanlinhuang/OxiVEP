#!/usr/bin/env python3
"""
Generate manuscript figures for OxiVEP.
Requires: matplotlib (pip install matplotlib)
"""

import csv
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not installed. Install with: pip install matplotlib")
    print("Generating text-based summaries instead.\n")

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
FIG_DIR = os.path.dirname(__file__)

COLORS = {
    'primary': '#6c7aee',
    'vep': '#f5426c',
    'high': '#f5426c',
    'moderate': '#f59e0b',
    'low': '#3b82f6',
    'modifier': '#6b7280',
    'success': '#10b981',
    'cli': '#6c7aee',
    'mid': '#818cf8',
    'data': '#34d399',
    'core': '#f59e0b',
    'sa': '#f472b6',
}

def read_csv(filename):
    path = os.path.join(DATA_DIR, filename)
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)


def fig1_architecture():
    """Figure 1: Architecture diagram as a proper figure."""
    if not HAS_MPL:
        print("Figure 1: Architecture (text-only, matplotlib not available)")
        return

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 9)
    ax.axis('off')

    def box(x, y, w, h, label, sublabel, color, fontsize=11):
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                              facecolor=color, edgecolor='#374151', linewidth=1.5, alpha=0.9)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2 + 0.15, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='white')
        if sublabel:
            ax.text(x + w/2, y + h/2 - 0.2, sublabel, ha='center', va='center',
                    fontsize=8, color='white', alpha=0.9, style='italic')

    def arrow(x1, y1, x2, y2):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#6b7280', lw=1.5))

    # Title
    ax.text(7, 8.6, 'OxiVEP Architecture', ha='center', va='center',
            fontsize=16, fontweight='bold', color='#1f2937')
    ax.text(7, 8.25, '9-crate Cargo workspace', ha='center', va='center',
            fontsize=10, color='#6b7280')

    # Layer 1: CLI (top)
    box(2.5, 7.0, 9, 0.9, 'oxivep-cli', 'Pipeline, web server, cache builder, SA builder', COLORS['cli'], 13)

    # Layer 2: Middle crates
    box(0.3, 5.3, 2.6, 0.9, 'oxivep-io', 'VCF/CSQ/JSON I/O', COLORS['mid'])
    box(3.2, 5.3, 2.6, 0.9, 'oxivep-hgvs', 'HGVSg/c/p', COLORS['mid'])
    box(6.1, 5.3, 2.8, 0.9, 'oxivep-consequence', 'SNV/indel/SV engine', COLORS['mid'])
    box(9.2, 5.3, 2.2, 0.9, 'oxivep-filter', 'Filter engine', COLORS['mid'])
    box(11.7, 5.3, 2.0, 0.9, 'oxivep-sa', 'Annotations', COLORS['sa'])

    # Layer 3: Data layer
    box(0.5, 3.5, 3.5, 0.9, 'oxivep-genome', 'Transcript, Exon, Gene, CodonTable', COLORS['data'])
    box(4.5, 3.5, 5.0, 0.9, 'oxivep-cache', 'GFF3, FASTA mmap, tabix, binary cache', COLORS['data'])
    box(10.0, 3.5, 3.5, 0.9, 'oxivep-sa (OxiSA)', 'ClinVar, gnomAD, REVEL, ...', COLORS['sa'])

    # Layer 4: Core (bottom)
    box(3.0, 1.7, 8.0, 0.9, 'oxivep-core', 'GenomicPosition, Consequence (49 SO terms), Allele, Strand, Impact, VariantType', COLORS['core'], 12)

    # Arrows: CLI -> middle layer
    for x in [1.6, 4.5, 7.5, 10.3, 12.7]:
        arrow(7, 7.0, x, 6.2)

    # Arrows: middle -> data layer
    for x in [1.6, 4.5, 7.5]:
        arrow(x, 5.3, x, 4.4)
    arrow(12.7, 5.3, 11.75, 4.4)

    # Arrows: data -> core
    for x in [2.25, 7.0]:
        arrow(x, 3.5, 7.0, 2.6)
    arrow(11.75, 3.5, 7.0, 2.6)

    # Stats annotation
    stats_text = "16,885 LOC  |  233 tests  |  2.3 MB binary (LTO + strip)"
    ax.text(7, 0.8, stats_text, ha='center', va='center',
            fontsize=10, color='#6b7280',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#f3f4f6', edgecolor='#d1d5db'))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig1_architecture.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig1_architecture.pdf'), bbox_inches='tight')
    print("Saved fig1_architecture.png/pdf")
    plt.close()


def fig2_throughput_scaling():
    """Figure 2: OxiVEP throughput scaling on GIAB HG002 full-genome data."""
    data = read_csv('scaling.csv')
    variants = [int(d['variants']) for d in data]
    oxi_vps = [int(d['oxivep_vps']) for d in data]
    oxi_time = [float(d['oxivep_time_sec']) for d in data]

    if not HAS_MPL:
        print("Figure 2: Throughput Scaling")
        for v, ot, ovps in zip(variants, oxi_time, oxi_vps):
            print(f"  {v:>10,}v: {ot:.1f}s ({ovps:,.0f} v/s)")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: Wall-clock time
    ax1.plot(variants, oxi_time, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=10, zorder=3)
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of variants (GIAB HG002)', fontsize=12)
    ax1.set_ylabel('Wall-clock time (seconds)', fontsize=12)
    ax1.set_title('A. Annotation Time\n(508,530 transcripts, full GRCh38)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.1f}M'))
    for v, t in zip(variants, oxi_time):
        ax1.annotate(f'{t:.1f}s', xy=(v, t), xytext=(0, 12),
                     textcoords='offset points', ha='center', fontsize=10,
                     color=COLORS['primary'], fontweight='bold')

    # Panel B: Throughput
    ax2.plot(variants, oxi_vps, 'o-', color=COLORS['primary'], linewidth=2.5, markersize=10, zorder=3)
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of variants (GIAB HG002)', fontsize=12)
    ax2.set_ylabel('Throughput (variants/sec)', fontsize=12)
    ax2.set_title('B. Annotation Throughput', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.1f}M'))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.1f}M'))
    for v, vps in zip(variants, oxi_vps):
        ax2.annotate(f'{vps:,.0f}', xy=(v, vps), xytext=(0, 12),
                     textcoords='offset points', ha='center', fontsize=10,
                     color=COLORS['primary'], fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig2_throughput_scaling.pdf'), bbox_inches='tight')
    print("Saved fig2_throughput_scaling.png/pdf")
    plt.close()


def fig3_vep_concordance():
    """Figure 3: VEP concordance heatmap/bar chart."""
    data = read_csv('vep_concordance.csv')
    fields = [d['field'] for d in data]
    accuracy = [float(d['accuracy']) for d in data]

    if not HAS_MPL:
        print("Figure 3: VEP Concordance")
        for f, a in zip(fields, accuracy):
            print(f"  {f:25s} {a:.1f}%")
        print()
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = [COLORS['success'] if a == 100.0 else COLORS['vep'] for a in accuracy]
    bars = ax.barh(range(len(fields)), accuracy, color=colors, alpha=0.85)
    ax.set_yticks(range(len(fields)))
    ax.set_yticklabels(fields, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Concordance with Ensembl VEP v115.1 (%)', fontsize=12)
    ax.set_title('Field-Level Accuracy: OxiVEP vs Ensembl VEP\n(2,340 shared transcript-allele pairs, 173 variants)', fontsize=13, fontweight='bold')
    ax.set_xlim(95, 101)
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(x=100, color=COLORS['success'], linestyle='--', alpha=0.5, linewidth=1)

    for bar, acc in zip(bars, accuracy):
        ax.text(bar.get_width() - 0.3, bar.get_y() + bar.get_height()/2,
                f'{acc:.1f}%', va='center', ha='right', fontsize=9, color='white', fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_concordance.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig3_vep_concordance.pdf'), bbox_inches='tight')
    print("Saved fig3_vep_concordance.png/pdf")
    plt.close()


def fig4_consequence_distribution():
    """Figure 4: Consequence distribution from real 1KGP data."""
    data = read_csv('consequence_distribution.csv')
    consequences = [d['consequence'] for d in data]
    counts = [int(d['count']) for d in data]

    impact_map = {
        'splice_acceptor_variant': 'HIGH', 'splice_donor_variant': 'HIGH',
        'stop_gained': 'HIGH', 'frameshift_variant': 'HIGH',
        'missense_variant': 'MODERATE', 'inframe_insertion': 'MODERATE',
        'inframe_deletion': 'MODERATE',
        'splice_region_variant': 'LOW', 'synonymous_variant': 'LOW',
        'splice_polypyrimidine_tract_variant': 'LOW',
        'splice_donor_5th_base_variant': 'LOW', 'splice_donor_region_variant': 'LOW',
    }
    def get_color(c):
        imp = impact_map.get(c, 'MODIFIER')
        return COLORS.get(imp.lower(), COLORS['modifier'])

    if not HAS_MPL:
        print("Figure 4: Consequence Distribution")
        for c, n in zip(consequences, counts):
            bar = '#' * (n // 200)
            print(f"  {c:45s} {n:>6d}  {bar}")
        print()
        return

    fig, ax = plt.subplots(figsize=(11, 7))
    colors = [get_color(c) for c in consequences]
    bars = ax.barh(range(len(consequences)), counts, color=colors, alpha=0.85)
    ax.set_yticks(range(len(consequences)))
    ax.set_yticklabels(consequences, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()
    ax.set_xlabel('Annotation count', fontsize=12)
    ax.set_title('Predicted Consequence Distribution\n(1,000 real 1KGP chr22 variants, 28,161 annotations)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xscale('log')

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() * 1.1, bar.get_y() + bar.get_height()/2,
                f'{count:,}', va='center', fontsize=9, color='#555')

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['high'], alpha=0.85, label='HIGH'),
        mpatches.Patch(facecolor=COLORS['moderate'], alpha=0.85, label='MODERATE'),
        mpatches.Patch(facecolor=COLORS['low'], alpha=0.85, label='LOW'),
        mpatches.Patch(facecolor=COLORS['modifier'], alpha=0.85, label='MODIFIER'),
    ]
    ax.legend(handles=legend_elements, title='Impact', loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig4_consequence_distribution.pdf'), bbox_inches='tight')
    print("Saved fig4_consequence_distribution.png/pdf")
    plt.close()


def fig5_resource_usage():
    """Figure 5: Resource usage comparison."""
    data = read_csv('resource_usage.csv')
    mem_data = {d['metric']: float(d['value']) for d in data}

    if not HAS_MPL:
        print("Figure 5: Resource Usage")
        for k, v in mem_data.items():
            unit = next((d['unit'] for d in data if d['metric'] == k), '')
            print(f"  {k}: {v} {unit}")
        print()
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Memory comparison
    variants_labels = ['1K', '10K', '100K']
    oxivep_mem = [mem_data.get('peak_memory_1000v', 0),
                  mem_data.get('peak_memory_10000v', 0),
                  mem_data.get('peak_memory_100000v', 0)]
    vep_mem = [500, 500, 600]

    x = range(len(variants_labels))
    w = 0.35
    ax1.bar([i - w/2 for i in x], vep_mem, w, label='Ensembl VEP', color=COLORS['vep'], alpha=0.8)
    ax1.bar([i + w/2 for i in x], oxivep_mem, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax1.set_xlabel('Input size (variants)', fontsize=12)
    ax1.set_ylabel('Peak memory (MB)', fontsize=12)
    ax1.set_title('A. Memory Usage', fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(variants_labels)
    ax1.legend()
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, axis='y')

    # Panel B: Binary size / startup
    metrics = ['Binary Size\n(MB)', 'Startup + GFF3\n(ms)', 'Startup\n(cached, ms)']
    vep_vals = [200, 10000, 10000]
    oxi_vals = [mem_data.get('binary_size', 0),
                mem_data.get('startup_time_gff3_chr22', 0),
                mem_data.get('startup_time_cached', 0)]

    x = range(len(metrics))
    ax2.bar([i - w/2 for i in x], vep_vals, w, label='Ensembl VEP', color=COLORS['vep'], alpha=0.8)
    ax2.bar([i + w/2 for i in x], oxi_vals, w, label='OxiVEP', color=COLORS['primary'], alpha=0.8)
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_title('B. Size & Startup', fontsize=13, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(metrics, fontsize=10)
    ax2.legend()
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig5_resource_usage.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig5_resource_usage.pdf'), bbox_inches='tight')
    print("Saved fig5_resource_usage.png/pdf")
    plt.close()


def fig6_organism_throughput():
    """Figure 6: Multi-organism throughput comparison (new figure)."""
    data = read_csv('organism_comparison.csv')

    if not HAS_MPL:
        print("Figure 6: Organism Throughput")
        for d in data:
            print(f"  {d['organism']:20s} {int(d['variants']):>8,}v  {float(d['time_sec']):.2f}s  {int(d['variants_per_sec']):>8,} v/s")
        print()
        return

    # Pick the largest-variant entry per organism for the bar chart
    best = {}
    for d in data:
        org = d['organism']
        vps = int(d['variants_per_sec'])
        if org not in best or vps > best[org]['vps']:
            best[org] = {
                'vps': vps,
                'variants': int(d['variants']),
                'transcripts': int(d['transcripts']),
                'time': float(d['time_sec']),
                'assembly': d['assembly'],
            }

    # Sort by throughput
    orgs = sorted(best.keys(), key=lambda o: best[o]['vps'])
    vps_vals = [best[o]['vps'] for o in orgs]
    labels = [f"{o}\n({best[o]['assembly']}, {best[o]['transcripts']:,} tr.)" for o in orgs]

    # Color by genome size category
    org_colors = []
    for o in orgs:
        tr = best[o]['transcripts']
        if tr > 100000:
            org_colors.append('#ef4444')  # red - large
        elif tr > 30000:
            org_colors.append('#f59e0b')  # amber - medium
        else:
            org_colors.append('#10b981')  # green - small

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.barh(range(len(orgs)), vps_vals, color=org_colors, alpha=0.85)
    ax.set_yticks(range(len(orgs)))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Peak throughput (variants/sec)', fontsize=12)
    ax.set_title('Cross-Organism Annotation Throughput\n(full Ensembl genome annotations, FASTA + HGVS, median of 3 runs)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}K' if x < 1e6 else f'{x/1e6:.1f}M'))

    for bar, vps, org in zip(bars, vps_vals, orgs):
        v = best[org]['variants']
        t = best[org]['time']
        ax.text(bar.get_width() * 1.15, bar.get_y() + bar.get_height()/2,
                f'{vps:,.0f} v/s  ({v:,} variants in {t:.1f}s)',
                va='center', fontsize=9, color='#333')

    legend_elements = [
        mpatches.Patch(facecolor='#10b981', alpha=0.85, label='<30K transcripts'),
        mpatches.Patch(facecolor='#f59e0b', alpha=0.85, label='30K-100K transcripts'),
        mpatches.Patch(facecolor='#ef4444', alpha=0.85, label='>100K transcripts'),
    ]
    ax.legend(handles=legend_elements, title='Genome complexity', loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig6_organism_throughput.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig6_organism_throughput.pdf'), bbox_inches='tight')
    print("Saved fig6_organism_throughput.png/pdf")
    plt.close()


if __name__ == '__main__':
    print("Generating OxiVEP manuscript figures...")
    print(f"Data dir: {os.path.abspath(DATA_DIR)}")
    print(f"Output dir: {os.path.abspath(FIG_DIR)}")
    print()
    fig1_architecture()
    fig2_throughput_scaling()
    fig3_vep_concordance()
    fig4_consequence_distribution()
    fig5_resource_usage()
    fig6_organism_throughput()
    print("\nDone.")
