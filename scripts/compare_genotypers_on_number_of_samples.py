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

def run_dnadiff(ref, query):
    '''Runs dnadiff to map ref fasta to query fasta and create a SNPs file of differences.'''
    dnadiff_binary = find_binary('dnadiff')
    command = ' '.join([dnadiff_binary, ref, query, '-p', 'tmp.dnadiff'])
    syscall(command)
    files = glob.iglob('tmp.dnadiff*')
    for f in files:
        if not f.endswith('.snps'):
            os.unlink(f)

def get_vcf_pairs(dir1, dir2):
    '''Compares strings to find pairs of VCF files called against the same VCF reference'''
    paths = glob.glob(dir1 + '/*.ref.fa')
    refs1 = [f.split("/")[-1] for f in paths]

    paths = glob.glob(dir2 + '/*.ref.fa')
    refs2 = [f.split("/")[-1] for f in paths]

    intersect_refs = set(refs1).intersection(refs2)

    results = []
    for r in intersect_refs:
        name = r[:-7]
        vcf_ref = "/".join([dir1, r])
        vcf = name + ".vcf" 
        vcf1 = "/".join([dir1, vcf])
        vcf2 = "/".join([dir2, vcf])
        results.append([name, vcf_ref, vcf1, vcf2])
    return results

def run_compare(snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, name, flank, mask1, mask2):
    '''Runs minos command to compare pair of VCFs to dnadiff snps file'''
    minos_binary = find_binary('minos')
    command = ' '.join([minos_binary, 'check_snps', snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, "tmp.compare." + name, "--flank_length", str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--allow_flank_mismatches"])
    if mask1:
        command += " --exclude_bed1 " + mask1
    if mask2:
        command += " --exclude_bed2 " + mask2
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

def run_individual_minos(truth1, truth2, vcf1, vcf2, vcf_ref, name, flank, mask1, mask2):
    '''Runs minos command to check each VCF with truth'''
    minos_binary = find_binary('minos')
    index_ref(truth1)
    filter_vcf_by_ref_pos(vcf1, vcf_ref, flank)
    filtered_vcf1 = vcf1.replace(".vcf",".good.vcf")
    command = ' '.join([minos_binary, 'check_with_ref', filtered_vcf1, vcf_ref, truth1, "tmp.individual." + name + ".1", "--allow_flank_mismatches", "--flank_length",  str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--max_soft_clipped", str(20)])
    if mask1:
        command += " --exclude_bed " + mask1
    syscall(command)
    index_ref(truth2)
    filter_vcf_by_ref_pos(vcf2, vcf_ref, flank)
    filtered_vcf2 = vcf2.replace(".vcf",".good.vcf")
    command = ' '.join([minos_binary, 'check_with_ref', filtered_vcf2, vcf_ref, truth2, "tmp.individual." + name + ".2", "--allow_flank_mismatches", "--flank_length",  str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--max_soft_clipped", str(20)])
    if mask2:
        command += " --exclude_bed " + mask2
    syscall(command)

def create_or_append(ultimate_file, data_file):
    command = ""
    if os.path.isfile(ultimate_file):
        command = ' '.join(['tail -n+2', data_file, '>>', ultimate_file])
    else:
        command = ' '.join(['cp', data_file, ultimate_file])
    syscall(command)

def append_stats(name):
    create_or_append(name + ".y.gt_conf_hist.tsv", "tmp.compare." + name + ".gt_conf_hist.tsv")
    create_or_append(name + ".y.stats.tsv", "tmp.compare." + name + ".stats.tsv")
    create_or_append(name + ".x.gt_conf_hist.TP.tsv", "tmp.individual." + name + '.1.gt_conf_hist.TP.tsv')
    create_or_append(name + ".x.gt_conf_hist.FP.tsv", "tmp.individual." + name + '.1.gt_conf_hist.FP.tsv')
    create_or_append(name + ".x.stats.tsv", "tmp.individual." + name + '.1.stats.tsv')
    create_or_append(name + ".x.gt_conf_hist.TP.tsv", "tmp.individual." + name + '.2.gt_conf_hist.TP.tsv')
    create_or_append(name + ".x.gt_conf_hist.FP.tsv", "tmp.individual." + name + '.2.gt_conf_hist.FP.tsv')
    create_or_append(name + ".x.stats.tsv", "tmp.individual." + name + '.2.stats.tsv')

def minos_to_df(name):
    y_gt = pd.read_table(name + ".y.gt_conf_hist.tsv")
    y_stats = pd.read_table(name + ".y.stats.tsv")

    x_tp = pd.read_table(name + ".x.gt_conf_hist.TP.tsv")
    x_fp = pd.read_table(name + ".x.gt_conf_hist.FP.tsv")
    x_stats = pd.read_table(name + ".x.stats.tsv")

    yscat = []
    xscat = []

    conf_threshold = 0
    c_values = list(y_gt['GT_CONF'].values) + list(x_tp['GT_CONF'].values) + list(x_fp['GT_CONF'].values)
    c_values.sort()

    ytotal = sum(y_stats['total']) - sum(y_stats['excluded_vars'])
    for confidence in c_values:
        dnadiff_frac = float(sum(y_gt[y_gt['GT_CONF'] >= confidence]['Count']))/float(ytotal)
        if dnadiff_frac < 0.1 or confidence == 0:
            continue
        yscat.append(dnadiff_frac)
        sum_fp = sum(x_fp[x_fp['GT_CONF'] >= confidence]['Count'])
        sum_tp = sum(x_tp[x_tp['GT_CONF'] >= confidence]['Count'])
        if sum_fp > 0 and sum_tp > 0:
            xscat.append(sum_fp/float(sum_fp + sum_tp))
            if len(yscat) > 2 and (yscat[-1] < 0.6 < yscat[-2] or yscat[-2] < 0.6 < yscat[-1]):
                print(name, xscat[-1], yscat[-1])
        else:
            xscat.append(float(0))

    df = pd.DataFrame()
    df['yscat'] = yscat
    df['xscat'] = xscat
    df['name'] = name
    pk.dump(df, open(name + '.pkl', 'wb'))
    print("Max", name, max(yscat))
    #print(df)
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
        plt.xlim(-0.005, 0.02)

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
                if x['name'].values[0].startswith("pandora_genotyped_full") or x['name'].values[0].startswith("pandora_compare"):
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
parser.add_argument('--sample_tsv', '-t1', type=str,
                    help='TSV with a line for each sample, containing ID, TRUTH_FASTA, VCF_DIRECTORY, MASK')
parser.add_argument('--recall_flank', '-fr', type=int, default=5,
                    help='Size of flank sequence to use when comparing true alleles to vcf alleles')
parser.add_argument('--precision_flank', '-fp', type=int, default=31,
                    help='Size of flank sequence to use when comparing alleles to truth assembly')
args = parser.parse_args()

index = pd.read_csv(args.sample_tsv, sep='\t', header=None, names=['id', 'truth', 'vcf_dir', 'mask'])

unzipped_truths = []
for truth in index['truth']:
    if truth.endswith(".gz"):
        new_name = truth.split('/')[-1][:-3]
        command = ' '.join(['zcat', truth, '>', new_name])
        syscall(command)
        unzipped_truths.append(new_name)
    else:
        unzipped_truths.append(truth)
index['unzipped_truth'] = unzipped_truths

names = []
for pair in itertools.combinations(index.itertuples(), 2):
    (num1, id1, zipped_truth1, sample_dir1, mask1, truth1) = pair[0]
    (num2, id2, zipped_truth2, sample_dir2, mask2, truth2) = pair[1] 
    run_dnadiff(truth1, truth2)
    pairs = get_vcf_pairs(sample_dir1, sample_dir2)
    for run in pairs:
        name, vcf_ref, vcf1, vcf2 = run
        print(name, id1, id2)
        names.append(name)
        if len(glob.glob(name +"_" + id1 + "_" + id2 + ".pkl")) == 0:
            run_compare("tmp.dnadiff.snps", truth1, truth2, vcf1, vcf2, vcf_ref, name, args.recall_flank, mask1, mask2)
            run_individual_minos(truth1, truth2, vcf1, vcf2, vcf_ref, name, args.precision_flank, mask1, mask2)
            append_stats(name)
            files = glob.glob("tmp.individual*")
            for f in files:
                if not f.endswith(".tsv"):
                    os.unlink(f)
            files = glob.glob("tmp.compare*")
            for f in files:
                os.unlink(f)
    files = glob.glob(truth1 + ".*")
    for f in files:
        os.unlink(f)
    files = glob.glob(truth2 + ".*")
    for f in files:
        os.unlink(f)
    
dfs = []
for name in names:
    dfs.append(minos_to_df(name))

print("plot graphs")
plot_graphs(dfs)

