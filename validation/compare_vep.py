#!/usr/bin/env python3
"""Compare OxiVEP and Ensembl VEP annotated VCF outputs.

Usage:
    python compare_vep.py <oxivep.vcf> <vep.vcf> [--coding] [--verbose]

Compares CSQ annotations field by field on shared (allele, transcript) pairs.
Reports per-field accuracy, consequence type distribution, and example mismatches.
"""
import sys
import re
import argparse
from collections import defaultdict, Counter


def parse_csq_header(line):
    m = re.search(r'Format: ([^"]+)', line)
    return m.group(1).split('|') if m else []


def parse_vcf(filename):
    """Parse VCF with CSQ INFO field. Returns {(chrom,pos,ref,alt): [csq_dicts]}."""
    csq_fields = []
    variants = {}
    with open(filename) as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                csq_fields = parse_csq_header(line)
                continue
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            csq_str = ''
            for field in parts[7].split(';'):
                if field.startswith('CSQ='):
                    csq_str = field[4:]
                    break
            key = (chrom, pos, ref, alt)
            entries = []
            if csq_str:
                for entry in csq_str.split(','):
                    values = entry.split('|')
                    d = {csq_fields[i]: values[i] if i < len(values) else ''
                         for i in range(len(csq_fields))}
                    entries.append(d)
            variants[key] = entries
    return variants


COMPARE_FIELDS = [
    'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature', 'BIOTYPE',
    'EXON', 'INTRON', 'HGVSc', 'HGVSp',
    'cDNA_position', 'CDS_position', 'Protein_position',
    'Amino_acids', 'Codons', 'DISTANCE', 'STRAND',
]

CODING_CONSEQS = {
    'missense_variant', 'synonymous_variant', 'frameshift_variant',
    'stop_gained', 'stop_lost', 'start_lost', 'inframe_insertion',
    'inframe_deletion', 'splice_donor_variant', 'splice_acceptor_variant',
    '5_prime_UTR_variant', '3_prime_UTR_variant', 'splice_region_variant',
    'stop_retained_variant', 'protein_altering_variant',
    'incomplete_terminal_codon_variant', 'coding_sequence_variant',
}


def compare(oxi, vep, coding_only=False, verbose=False):
    field_match = defaultdict(int)
    field_mismatch = defaultdict(int)
    field_examples = defaultdict(list)
    csq_types_oxi = Counter()
    csq_types_vep = Counter()
    extra_oxi_bt = Counter()
    missing_oxi_bt = Counter()
    total_shared = 0
    total_variants = 0

    for key in sorted(set(oxi.keys()) | set(vep.keys())):
        oxi_csq = oxi.get(key, [])
        vep_csq = vep.get(key, [])
        if not oxi_csq or not vep_csq:
            continue
        total_variants += 1

        oxi_map = {(e.get('Allele', ''), e.get('Feature', '')): e for e in oxi_csq}
        vep_map = {(e.get('Allele', ''), e.get('Feature', '')): e for e in vep_csq}

        # Track extra/missing transcripts
        for atx in set(oxi_map) - set(vep_map):
            bt = oxi_map[atx].get('BIOTYPE', '?')
            extra_oxi_bt[bt] += 1
        for atx in set(vep_map) - set(oxi_map):
            bt = vep_map[atx].get('BIOTYPE', '?')
            missing_oxi_bt[bt] += 1

        chrom, pos, ref, alt = key
        loc = f"{chrom}:{pos} {ref}>{alt}"

        for atx in set(oxi_map) & set(vep_map):
            o = oxi_map[atx]
            v = vep_map[atx]

            # Filter to coding if requested
            if coding_only:
                oc = set(o.get('Consequence', '').split('&'))
                vc = set(v.get('Consequence', '').split('&'))
                if not (oc & CODING_CONSEQS or vc & CODING_CONSEQS):
                    continue

            total_shared += 1
            allele, tx = atx

            for f in COMPARE_FIELDS:
                ov = o.get(f, '')
                vv = v.get(f, '')
                if ov == vv:
                    field_match[f] += 1
                else:
                    field_mismatch[f] += 1
                    if len(field_examples[f]) < 5:
                        field_examples[f].append(
                            f"    {loc} tx={tx}: oxi=\"{ov}\" vep=\"{vv}\""
                        )

            # Count consequence types
            for c in o.get('Consequence', '').split('&'):
                if c:
                    csq_types_oxi[c] += 1
            for c in v.get('Consequence', '').split('&'):
                if c:
                    csq_types_vep[c] += 1

    # ---- Report ----
    print(f"\n{'=' * 70}")
    print(f"VEP CONCORDANCE REPORT")
    print(f"{'=' * 70}")
    print(f"  Variants compared:         {total_variants}")
    print(f"  Shared (allele, tx) pairs: {total_shared}")
    extra = sum(extra_oxi_bt.values())
    missing = sum(missing_oxi_bt.values())
    print(f"  Extra transcripts (OxiVEP only):  {extra}")
    print(f"  Missing transcripts (VEP only):   {missing}")

    print(f"\n--- FIELD-LEVEL ACCURACY ---")
    print(f"{'Field':<20} {'Match':>8} {'Mismatch':>10} {'Accuracy':>10}")
    print(f"{'-' * 52}")
    for f in COMPARE_FIELDS:
        m = field_match[f]
        mm = field_mismatch[f]
        total = m + mm
        pct = 100 * m / total if total > 0 else 0
        flag = '' if mm == 0 else '  <--'
        print(f"{f:<20} {m:>8} {mm:>10} {pct:>9.2f}%{flag}")

    if any(field_mismatch[f] > 0 for f in COMPARE_FIELDS):
        print(f"\n--- MISMATCH EXAMPLES ---")
        for f in COMPARE_FIELDS:
            if field_mismatch[f] > 0:
                print(f"\n  {f} ({field_mismatch[f]} mismatches):")
                for ex in field_examples[f]:
                    print(ex)

    if verbose:
        all_types = sorted(set(csq_types_oxi) | set(csq_types_vep))
        print(f"\n--- CONSEQUENCE TYPE DISTRIBUTION ---")
        print(f"{'Type':<45} {'OxiVEP':>8} {'VEP':>8} {'Diff':>6}")
        for t in all_types:
            o = csq_types_oxi.get(t, 0)
            v = csq_types_vep.get(t, 0)
            d = o - v
            flag = ' ***' if d != 0 else ''
            print(f"{t:<45} {o:>8} {v:>8} {d:>+6}{flag}")

        if extra_oxi_bt:
            print(f"\n--- EXTRA TRANSCRIPTS BY BIOTYPE ---")
            for bt, cnt in extra_oxi_bt.most_common(10):
                print(f"  {bt}: {cnt}")
        if missing_oxi_bt:
            print(f"\n--- MISSING TRANSCRIPTS BY BIOTYPE ---")
            for bt, cnt in missing_oxi_bt.most_common(10):
                print(f"  {bt}: {cnt}")

    print()


def main():
    parser = argparse.ArgumentParser(description='Compare OxiVEP vs Ensembl VEP outputs')
    parser.add_argument('oxivep_vcf', help='OxiVEP annotated VCF')
    parser.add_argument('vep_vcf', help='Ensembl VEP annotated VCF')
    parser.add_argument('--coding', action='store_true',
                        help='Only compare coding/functional variant annotations')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show consequence type distribution and biotype breakdown')
    args = parser.parse_args()

    oxi = parse_vcf(args.oxivep_vcf)
    vep = parse_vcf(args.vep_vcf)
    compare(oxi, vep, coding_only=args.coding, verbose=args.verbose)


if __name__ == '__main__':
    main()
