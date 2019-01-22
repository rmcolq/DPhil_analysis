#!/usr/bin/env python 

import argparse
import os
import shutil
import glob
import subprocess
import sys
import operator
from Bio import SeqIO
import vcf
import gzip


class Error (Exception): pass

def syscall(command, allow_fail=False):
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if (not allow_fail) and completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('\nOutput from stdout:', completed_process.stdout, sep='\n', file=sys.stderr)
        print('\nOutput from stderr:', completed_process.stderr, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')
    print(completed_process.stdout)
    return completed_process

def filter_vcf_by_ref_pos(in_vcf, ref_fasta, flank_size):
    vcf_reader = vcf.Reader(open(in_vcf, 'r'))
    out_vcf_good = vcf.Writer(open(in_vcf.replace(".vcf",".good.vcf"),'w'), vcf_reader)
    out_vcf_bad = vcf.Writer(open(in_vcf.replace(".vcf",".bad.vcf"),'w'), vcf_reader)

    if ref_fasta[:-3] == ".gz":
        with gzip.open(ref_fasta, 'rt') as input_handle:
            vcf_ref = SeqIO.to_dict(SeqIO.parse(input_handle, "fasta"))

            for record in vcf_reader:
                try:
                    if (record.POS < flank_size or len(vcf_ref[record.CHROM]) - (record.POS + len(record.REF)) < flank_size):
                        out_vcf_bad.write_record(record)
                    else:
                        out_vcf_good.write_record(record)
                except:
                    print(record)
                    print(len(vcf_ref[record.CHROM]))
    else:
        with open(ref_fasta, 'r') as input_handle:
            vcf_ref = SeqIO.to_dict(SeqIO.parse(input_handle, "fasta"))

            for record in vcf_reader:
                try:
                    if (record.POS < flank_size or len(vcf_ref[record.CHROM]) - (record.POS + len(record.REF)) < flank_size):
                        out_vcf_bad.write_record(record)
                    else:
                        out_vcf_good.write_record(record)
                except:
                    print(record)
                    print(len(vcf_ref[record.CHROM]))
    out_vcf_bad.close()
    out_vcf_good.close()

def restrict_to_snps(in_vcf):
    vcf_reader = vcf.Reader(open(in_vcf, 'r'))
    out_vcf = vcf.Writer(open(in_vcf.replace(".vcf",".snps.vcf"),'w'), vcf_reader)

    for record in vcf_reader:
        if len(record.REF) > 1:
            continue
        for alt in record.ALT: 
            if len(alt) > 1:
                continue
        out_vcf.write_record(record)

def run_filter(vcf, vcf_ref, flank, snps):
    '''Filters VCF by ref pos, and optionally restricts to SNPs'''
    filter_vcf_by_ref_pos(vcf, vcf_ref, flank)
    filtered_vcf = vcf.replace(".vcf",".good.vcf")
    if snps:
        restrict_to_snps(filtered_vcf)
        filtered_vcf = filtered_vcf.replace(".vcf",".snps.vcf")
    os.rename(filtered_vcf, vcf.replace(".vcf",".filtered.vcf"))

parser = argparse.ArgumentParser(description='Identifies pairs of VCF files called against the same reference in 2 sample directories and runs a minos system call to evaluate how many dnadiff snp differences are captued between each pair of VCFs.')
parser.add_argument('--vcf', type=str,
                    help='VCF of calls to filter')
parser.add_argument('--vcf_ref', type=str,
                    help='Reference FASTA for VCF')
parser.add_argument('--flank', '-f', type=int, default=5,
                    help='Size of flank sequence to use when comparing true alleles to vcf alleles')
parser.add_argument('--snps', action='store_true',
                    help='Only evaluates SNP VCF records')
args = parser.parse_args()

vcf_ref = args.vcf_ref
if vcf_ref.endswith(".gz"):
    command = ' '.join(['zcat', vcf_ref, '> truth.fa'])
    syscall(command)
    vcf_ref = "truth.fa"

run_filter(args.vcf, vcf_ref, args.flank, args.snps)
