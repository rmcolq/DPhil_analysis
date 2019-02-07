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

def run_evaluate_recall(truth_vcf, truth_vcf_ref, query_vcf, query_vcf_ref, name, flank, mask):
    '''Runs minos command to evaluate how many variants from truth vcf are called in query vcf'''
    minos_binary = find_binary('minos')
    command = ' '.join([minos_binary, 'check_recall', truth_vcf, truth_vcf_ref, query_vcf, query_vcf_ref, "tmp.recall." + name, "--flank_length", str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--allow_flank_mismatches"])
    if mask:
        command += " --exclude_bed " + mask
    syscall(command)

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

def run_evaluate_precision(truth, vcf, vcf_ref, name, flank, mask, snps):
    '''Runs minos command to check each VCF with truth'''
    minos_binary = find_binary('minos')
    index_ref(truth)
    filter_vcf_by_ref_pos(vcf, vcf_ref, flank)
    filtered_vcf = vcf.replace(".vcf",".good.vcf")
    if snps:
        restrict_to_snps(filtered_vcf)
        filtered_vcf = filtered_vcf.replace(".vcf",".snps.vcf")
    command = ' '.join([minos_binary, 'check_with_ref', filtered_vcf, vcf_ref, truth, "tmp.precision." + name, "--allow_flank_mismatches", "--flank_length",  str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--max_soft_clipped", str(20)])
    if mask:
        command += " --exclude_bed " + mask
    syscall(command)

def minos_to_df(name):
    y_gt = pd.read_table("tmp.recall." + name + ".gt_conf_hist.tsv")
    y_stats = pd.read_table("tmp.recall." + name + ".stats.tsv")

    x_tp = pd.read_table("tmp.precision." + name + '.gt_conf_hist.TP.tsv')
    x_fp = pd.read_table("tmp.precision." + name + '.gt_conf_hist.FP.tsv')
    x_stats = pd.read_table("tmp.precision." + name + '.stats.tsv')

    yscat = []
    xscat = []

    conf_threshold = 0
    c_values = list(y_gt['GT_CONF'].values) + list(x_tp['GT_CONF'].values) + list(x_fp['GT_CONF'].values)
    c_values.sort()

    for confidence in c_values:
        dnadiff_frac = float(sum(y_gt[y_gt['GT_CONF'] >= confidence]['Count']))/float(y_stats['total'].values[0] - y_stats['excluded_vars'].values[0])
        if dnadiff_frac < 0.1:
            continue
        yscat.append(dnadiff_frac)
        sum_fp = sum(x_fp[x_fp['GT_CONF'] >= confidence]['Count'])
        sum_tp = sum(x_tp[x_tp['GT_CONF'] >= confidence]['Count'])
        if sum_fp > 0 and sum_tp > 0:
            xscat.append(sum_fp/float(sum_fp + sum_tp))
        else:
            xscat.append(float(0))

    df = pd.DataFrame()
    df['yscat'] = yscat
    df['xscat'] = xscat
    df['name'] = name
    pk.dump(df, open(name + '.pkl', 'wb'))
    print(df)
    return df

def plot_graphs(items):
    for i in range(5):
        # Define plot
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 12)

        ax.grid(b=True)
        ax.set_axisbelow(b=True)
        plt.style.use('seaborn-colorblind')

        # Label the axes and title the plot
        ax.set_xlabel('Number FPs/Number Genotyped', size=26)
        ax.set_ylabel('Fraction of dnadiff SNPs discoverable from VCFs', size=26)
        # ax.set_title('Precision Recall', size = 30)

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        snippy_i = 0
        nanopolish_i = 0

        # Make a scatter plot
        for x in items:
            print("add df to graph")
            if len(x['name'].values) > 0:
                if x['name'].values[0].startswith("pandora_genotyped_full") or x['name'].values[0].startswith("pandora_recall"):
                    col = colormap_pandora(0)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina") and i > 0:
                    col = colormap_pandora(100)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_30") and i > 1:
                    col = colormap_pandora(200)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("snippy") and i > 2:
                    col = colormap_snippy(snippy_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    snippy_i += 1
                elif x['name'].values[0].startswith("nanopolish") and i > 3:
                    col = colormap_nanopolish(nanopolish_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    nanopolish_i += 1
            else:
                print(x['name'])

        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0 and len(labels) > 0 :
            hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
            handles2, labels2 = zip(*hl)
            ax.legend(handles2, labels2, frameon=False, loc='lower right')
            plt.savefig('roc%d.png' %i, transparent=True)
        else:
            print(handles, labels)


parser = argparse.ArgumentParser(description='Identifies pairs of VCF files called against the same reference in 2 sample directories and runs a minos system call to evaluate how many dnadiff snp differences are captued between each pair of VCFs.')
parser.add_argument('--truth_vcf', type=str,
                    help='VCF of calls to verify')
parser.add_argument('--truth_vcf_ref', type=str,
                    help='Reference FASTA for truth VCF')
parser.add_argument('--mask', type=str, default="",
                    help='BED file of regions of truth reference which are untrustworthy')
parser.add_argument('--sample_vcf', type=str,
                    help='VCF of genotyped')
parser.add_argument('--sample_vcf_ref', type=str,
                    help='Reference FASTA for sample VCF')
parser.add_argument('--recall_flank', '-fr', type=int, default=5,
                    help='Size of flank sequence to use when comparing true alleles to vcf alleles')
parser.add_argument('--precision_flank', '-fp', type=int, default=31,
                    help='Size of flank sequence to use when comparing alleles to truth assembly')
parser.add_argument('--snps', action='store_true',
                    help='Only evaluates SNP VCF records')
args = parser.parse_args()

truth_vcf = args.truth_vcf
truth_vcf_ref = args.truth_vcf_ref
if truth_vcf_ref.endswith(".gz"):
    command = ' '.join(['zcat', truth_vcf_ref, '> truth.fa'])
    syscall(command)
    truth_vcf_ref = "truth.fa"

sample_vcf = args.sample_vcf
sample_vcf_ref = args.sample_vcf_ref
if sample_vcf_ref.endswith(".gz"):
    command = ' '.join(['zcat', sample_vcf_ref, '> sample.fa'])
    syscall(command)
    sample_vcf_ref = "sample.fa"
name = sample_vcf_ref[:-7]

print(name)
if len(glob.glob(name + ".pkl")) == 0:
    run_evaluate_recall(truth_vcf, truth_vcf_ref, sample_vcf, sample_vcf_ref, name, args.recall_flank, args.mask)
    run_evaluate_precision(truth_vcf_ref, sample_vcf, sample_vcf_ref, name, args.precision_flank, args.mask, args.snps)
    df = minos_to_df(name)

