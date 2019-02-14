#!/usr/bin/env python 

import argparse
import os
import shutil
import glob
import subprocess
import pandas as pd
import numpy as np
import _pickle as pk
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import operator
from Bio import SeqIO
import vcf
import gzip
import itertools
import seaborn as sns


class Error (Exception): pass

def find_binary(program, allow_fail=False):
    which_output = shutil.which(program)
    if which_output is None and not allow_fail:
        raise Error('Error finding ' + program + ' in $PATH.')
    return which_output

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

def index_ref(ref):
    bwa_binary = find_binary('bwa')
    command = ' '.join([bwa_binary, 'index', ref])
    syscall(command)

def filter_vcf_by_ref_pos(in_vcf, ref_fasta, flank_size):
    vcf_reader = vcf.Reader(open(in_vcf, 'r'))
    out_vcf_good = vcf.Writer(open(in_vcf.replace(".vcf",".good.vcf"),'w'), vcf_reader)
    out_vcf_bad = vcf.Writer(open(in_vcf.replace(".vcf",".bad.vcf"),'w'), vcf_reader)

    if ref_fasta[-3:] == ".gz":
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

def run_dnadiff(ref, query):
    '''Runs dnadiff to map ref fasta to query fasta and create a SNPs file of differences.'''
    dnadiff_binary = find_binary('dnadiff')
    command = ' '.join([dnadiff_binary, ref, query, '-p', 'tmp.dnadiff'])
    syscall(command)
    files = glob.iglob('tmp.dnadiff*')
    for f in files:
        if not f.endswith('.snps'):
            os.unlink(f)

def run_compare(snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, name, flank):
    '''Runs minos command to compare pair of VCFs to dnadiff snps file'''
    minos_binary = find_binary('minos')
    command = ' '.join([minos_binary, 'check_snps', snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, "tmp.compare." + name, "--flank_length", str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--allow_flank_mismatches", "--max_soft_clipped", str(5)])
    syscall(command)

def run_individual(truth, vcf, vcf_ref, name, flank):
    '''Runs minos command to check each VCF with truth'''
    minos_binary = find_binary('minos')
    index_ref(truth)
    filter_vcf_by_ref_pos(vcf, vcf_ref, flank)
    filtered_vcf = vcf.replace(".vcf",".good.vcf")
    command = ' '.join([minos_binary, 'check_with_ref', filtered_vcf, vcf_ref, truth, "tmp.individual." + name, "--allow_flank_mismatches", "--flank_length",  str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--max_soft_clipped", str(20)])
    syscall(command)

def minos_stat_recall(name):
    stats = pd.read_table("tmp.compare." + name + ".stats.tsv")
    if sum(stats['total'].values) > 0:
        return float(sum(stats['found_vars'].values))/sum(stats['total'].values)
    else:
        return float(0)

def minos_stat_indiv(name1, name2):
    total = 0
    total_right = 0
    for name in [name1, name2]:
        stats = pd.read_table("tmp.individual." + name + ".stats.tsv")
        if sum(stats['total'].values) > 0:
            total += sum(stats['total'].values) - sum(stats['gt_excluded'].values)
            total_right += sum(stats['gt_correct'].values)
    if total > 0:
        return float(total_right)/total
    else:
        return float(0)

def find_diffs(vcf_file1, vcf_file2):
    vcf_reader1 = vcf.Reader(open(vcf_file1, 'r'))
    vcf_reader2 = vcf.Reader(open(vcf_file2, 'r'))
    out_vcf1 = vcf.Writer(open(vcf_file1.replace(".vcf",".filtered.vcf"),'w'), vcf_reader1)
    out_vcf2 = vcf.Writer(open(vcf_file2.replace(".vcf",".filtered.vcf"),'w'), vcf_reader2)
    for (record1, record2) in zip(vcf_reader1, vcf_reader2):
        rec1 = record1.samples[0]
        gt1 = rec1['GT'][0]
        rec2 = record2.samples[0]
        gt2 = rec2['GT'][0]
        if gt1 != "." and gt2 != ".": #gt1 != gt2
            out_vcf1.write_record(record1)
            out_vcf2.write_record(record2)

def compare_results(result_dir, recall_flank, precision_flank):
    refs = glob.glob("%s/genomes*.fa" %result_dir)
    vcfs = glob.glob("%s/pandora_illumina*.vcf" %result_dir)
    vcf_refs = glob.glob("%s/pandora_illumina*.vcf_ref.fa" %result_dir)
    df = pd.DataFrame()
    df['id'] = ['genome1', 'genome2', 'genome3', 'genome4']
    df['truth'] = refs
    df['vcf'] = vcfs
    df['vcf_ref'] = vcf_refs
    
    results_df = pd.DataFrame(index=df['id'].values, columns=df['id'].values)
    
    for pair in itertools.combinations(df.itertuples(), 2):
        (num1, id1, truth1, vcf1, vcf_ref1) = pair[0]
        (num2, id2, truth2, vcf2, vcf_ref2) = pair[1]
        # get recall
        run_dnadiff(truth1, truth2)
        if os.stat("tmp.dnadiff.snps").st_size == 0:
            print("no snps")
            #results_df[id1][id2] = np.nan
            #results_df[id2][id1] = np.nan
            #results_df[id1][id1] = np.nan
        else:
            comp_name = "%s_%s" %(id1, id2)
            run_compare("tmp.dnadiff.snps", truth1, truth2, vcf1, vcf2, vcf_ref1, comp_name, recall_flank)
            r = minos_stat_recall(comp_name)
            results_df[id1][id2] = r
            #results_df[id2][id1] = np.nan
            #results_df[id1][id1] = np.nan
        # get accuracy
        find_diffs(vcf1, vcf2)
        run_individual(truth1, vcf1.replace(".vcf",".filtered.vcf"), vcf_ref1, id1, precision_flank)
        run_individual(truth2, vcf2.replace(".vcf",".filtered.vcf"), vcf_ref2, id2, precision_flank)
        r = minos_stat_indiv(id1, id2)
        results_df[id2][id1] = r
        #clear up
        files = glob.glob("tmp.*")
        for f in files:
            os.unlink(f)
    results_df.to_csv("results.csv")
    print(results_df)
    
    results_df = results_df*100
    results_df = results_df.round(1)
    mask = results_df.isnull()
    results_df[mask] = 0
    print(results_df)

    fig, ax = plt.subplots()
    fig.set_size_inches(16, 16)
    plt.style.use('seaborn-colorblind')
    sns.heatmap(results_df, mask=mask, cmap='RdYlGn_r', linewidths=0.5, annot=True, square=True, annot_kws={"size":16}, fmt='2g')
    sns.set_context("notebook", font_scale=2.5, rc={"lines.linewidth": 2.5})
    plt.savefig('heatmap.png', transparent=True)

parser = argparse.ArgumentParser(description='Identifies pairs of VCF files called against the same reference in 2 sample directories and runs a minos system call to evaluate how many dnadiff snp differences are captued between each pair of VCFs.')
parser.add_argument('--dir', '-d', type=str,
                    help='Results directory')
parser.add_argument('--recall_flank', '-fr', type=int, default=15,
                    help='Size of flank sequence to use when comparing true alleles to vcf alleles')
parser.add_argument('--precision_flank', '-fp', type=int, default=31,
                    help='Size of flank sequence to use when evaluating vcf alleles')
args = parser.parse_args()

compare_results(args.dir, args.recall_flank, args.precision_flank)
