#!/usr/bin/env python 

import argparse
import vcf

def filter_vcf_for_overlaps(in_vcf):
    vcf_reader = vcf.Reader(open(in_vcf, 'r'))
    out_vcf = vcf.Writer(open(in_vcf.replace(".vcf",".filtered.vcf"),'w'), vcf_reader)

    last_record = None
    for record in vcf_reader:
        if last_record == None or record.POS > last_record.POS + len(last_record.REF):
            out_vcf.write_record(record)
            last_record = record
    out_vcf.close()

parser = argparse.ArgumentParser(description='Filter out VCF entries which overlap with previous')
parser.add_argument('--vcf', type=str,
                    help='VCF of calls to filter')
args = parser.parse_args()

filter_vcf_for_overlaps(args.vcf)
