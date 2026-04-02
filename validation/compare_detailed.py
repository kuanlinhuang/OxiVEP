#!/usr/bin/env python3
"""Detailed comparison of OxiVEP vs Ensembl VEP outputs."""
import sys
import re
from collections import defaultdict

def parse_csq_header(line):
    m = re.search(r'Format: ([^"]+)', line)
    if m:
        return m.group(1).split('|')
    return []

def parse_vcf(filename):
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
            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7] if len(parts) > 7 else ''
            
            csq_str = ''
            for field in info.split(';'):
                if field.startswith('CSQ='):
                    csq_str = field[4:]
                    break
            
            key = (chrom, pos, ref, alt)
            csq_entries = []
            if csq_str:
                for entry in csq_str.split(','):
                    values = entry.split('|')
                    d = {}
                    for i, name in enumerate(csq_fields):
                        d[name] = values[i] if i < len(values) else ''
                    csq_entries.append(d)
            variants[key] = csq_entries
    return variants, csq_fields

def main():
    oxi, oxi_fields = parse_vcf(sys.argv[1])
    vep, vep_fields = parse_vcf(sys.argv[2])
    
    compare_fields = ['Consequence', 'IMPACT', 'SYMBOL', 'BIOTYPE', 'EXON', 'INTRON',
                       'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 
                       'Protein_position', 'Amino_acids', 'Codons', 'DISTANCE', 'STRAND']
    
    all_keys = sorted(set(oxi.keys()) | set(vep.keys()))
    
    stats = {
        'total_variants': len(all_keys),
        'perfect_match': 0,
        'transcript_count_mismatch': 0,
        'consequence_mismatch': 0,
    }
    field_mismatches = defaultdict(int)
    field_examples = defaultdict(list)
    extra_oxi_transcripts = 0
    missing_oxi_transcripts = 0
    consequence_mismatches = []
    
    for key in all_keys:
        oxi_csq = oxi.get(key, [])
        vep_csq = vep.get(key, [])
        
        # Build lookup by (allele, transcript)
        def by_allele_tx(csq_list):
            d = {}
            for e in csq_list:
                k = (e.get('Allele',''), e.get('Feature',''))
                d[k] = e
            return d
        
        oxi_map = by_allele_tx(oxi_csq)
        vep_map = by_allele_tx(vep_csq)
        
        all_tx = set(oxi_map.keys()) | set(vep_map.keys())
        variant_ok = True
        
        for atx in all_tx:
            if atx not in vep_map:
                extra_oxi_transcripts += 1
                variant_ok = False
                continue
            if atx not in oxi_map:
                missing_oxi_transcripts += 1
                variant_ok = False
                continue
            
            o = oxi_map[atx]
            v = vep_map[atx]
            
            for field in compare_fields:
                oval = o.get(field, '')
                vval = v.get(field, '')
                if oval != vval:
                    variant_ok = False
                    field_mismatches[field] += 1
                    if len(field_examples[field]) < 5:
                        chrom, pos, ref, alt = key
                        allele, tx = atx
                        field_examples[field].append(
                            f"    {chrom}:{pos} {ref}>{alt} allele={allele} tx={tx}: oxi='{oval}' vep='{vval}'"
                        )
                    if field == 'Consequence' and len(consequence_mismatches) < 20:
                        chrom, pos, ref, alt = key
                        allele, tx = atx
                        consequence_mismatches.append(
                            f"  {chrom}:{pos} {ref}>{alt} tx={tx}:\n    OxiVEP: {oval}\n    VEP:    {vval}"
                        )
        
        if variant_ok:
            stats['perfect_match'] += 1
    
    print(f"\n{'='*70}")
    print(f"DETAILED COMPARISON: {stats['total_variants']} variants")
    print(f"{'='*70}")
    print(f"  Perfect match (all transcripts identical):  {stats['perfect_match']}")
    print(f"  Extra transcripts only in OxiVEP:           {extra_oxi_transcripts}")
    print(f"  Transcripts missing from OxiVEP:            {missing_oxi_transcripts}")
    
    print(f"\nFIELD-LEVEL MISMATCHES (on shared transcripts):")
    for field, count in sorted(field_mismatches.items(), key=lambda x: -x[1]):
        print(f"  {field}: {count}")
        for ex in field_examples[field]:
            print(ex)
    
    if consequence_mismatches:
        print(f"\nCONSEQUENCE MISMATCHES (first 20):")
        for m in consequence_mismatches:
            print(m)

if __name__ == '__main__':
    main()
