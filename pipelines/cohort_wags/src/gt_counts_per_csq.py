#!/usr/bin/env python3

import argparse
import csv
from cyvcf2 import VCF
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Add simplified per-site genotype class counts (hom_ref, het, hom_var, missing) to each CSQ line.")
    parser.add_argument("--vcf", required=True, help="Input VCF (bgzipped or not)")
    parser.add_argument("--vep_split", required=True, help="Input split-vep TSV file from bcftools +split-vep")
    parser.add_argument("--output", required=True, help="Output TSV with per-site genotype counts for each CSQ line")
    return parser.parse_args()

def normalize_alt_for_csq(alt, ref):
    return alt[len(ref):] if alt.startswith(ref) and len(ref) > 0 else alt

def count_genotypes_per_site(record):
    genotypes = record.genotypes
    counts = {'hom_ref': 0, 'het': 0, 'hom_var': 0, 'missing': 0}
    for gt in genotypes:
        a1, a2 = gt[0], gt[1]
        if a1 is None or a2 is None or a1 == -1 or a2 == -1:
            counts['missing'] += 1
        elif a1 == 0 and a2 == 0:
            counts['hom_ref'] += 1
        elif a1 == a2 and a1 > 0:
            counts['hom_var'] += 1
        else:
            counts['het'] += 1
    return counts

def main():
    args = parse_args()
    vcf_by_pos = defaultdict(list)
    vcf = VCF(args.vcf)
    for rec in vcf:
        vcf_by_pos[(rec.CHROM, rec.POS)].append(rec)

    with open(args.vep_split) as fin, open(args.output, 'w') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        for row in reader:
            if len(row) < 6:
                writer.writerow(row + ['.', '.', '.', '.', '.', '.'])
                continue

            chrom = row[0]
            pos = int(row[1])
            ref = row[2]
            alt_field = row[3]
            ac_field = row[4]
            csq = row[5]

            alts = alt_field.split(',')
            acs = ac_field.split(',') if ac_field != '.' else ['.'] * len(alts)
            csq_allele = csq.split('|')[0]
            normalized_alts = [normalize_alt_for_csq(a, ref) for a in alts]

            try:
                alt_index = normalized_alts.index(csq_allele)
            except ValueError:
                alt_index = None

            vcf_record_list = vcf_by_pos.get((chrom, pos), [])
            if not vcf_record_list:
                writer.writerow(row + [csq_allele, '.', '.', '.', '.', '.'])
                continue

            record = vcf_record_list[0]
            site_counts = count_genotypes_per_site(record)

            writer.writerow(row + [
                csq_allele,
                acs[alt_index] if alt_index is not None and alt_index < len(acs) else '.',
                site_counts['hom_ref'],
                site_counts['het'],
                site_counts['hom_var'],
                site_counts['missing']
            ])

if __name__ == "__main__":
    main()

