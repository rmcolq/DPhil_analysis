#!/usr/bin/env python 

import argparse
import os
import shutil
import glob
import subprocess
import pandas as pd
import numpy as np
import sys
import operator
from Bio import SeqIO
import vcf
import gzip
import itertools
from functools import partial
from multiprocessing import Pool


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

def run_dnadiff(ref, query, name):
    '''Runs dnadiff to map ref fasta to query fasta and create a SNPs file of differences.'''
    dnadiff_binary = find_binary('dnadiff')
    command = ' '.join([dnadiff_binary, ref, query, '-p', 'tmp/dnadiff_%s' %name])
    syscall(command)
    files = glob.glob('tmp/dnadiff*')
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

def get_all_names(dirs):
    '''Compares strings to find pairs of VCF files called against the same VCF reference'''
    refs = []
    for d in dirs:
        paths = glob.glob(d + '/*.ref.fa')
        refs.extend([f.split("/")[-1] for f in paths])

    results = []
    for r in list(set(refs)):
        name = r[:-7]
        results.append(name)
    return results

def get_vcfs(dir1):
    '''Compares strings to find pairs of VCF files called against the same VCF reference'''
    paths = glob.glob(dir1 + '/*.ref.fa')
    refs1 = [f.split("/")[-1] for f in paths]
    results = []
    for r in refs1:
        name = r[:-7]
        vcf_ref = "/".join([dir1, r])
        vcf = name + ".vcf"
        vcf1 = "/".join([dir1, vcf])
        results.append([name, vcf_ref, vcf1])
    return results

def run_compare(snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, name, flank, mask1, mask2):
    '''Runs minos command to compare pair of VCFs to dnadiff snps file'''
    if len(glob.glob("tmp/compare." + name + "*")) == 0:
        minos_binary = find_binary('minos')
        command = ' '.join([minos_binary, 'check_snps', snps_file, truth1, truth2, vcf1, vcf2, vcf_ref, "tmp/compare." + name, "--flank_length", str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--allow_flank_mismatches", "--max_soft_clipped", str(5)])#, "--no_filter_cluster"])
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

def run_individual_minos(truth, vcf, vcf_ref, name, flank, mask):
    '''Runs minos command to check each VCF with truth'''
    if len(glob.glob("tmp/individual." + name + "*")) == 0:
        minos_binary = find_binary('minos')
        index_ref(truth)
        filter_vcf_by_ref_pos(vcf, vcf_ref, flank)
        filtered_vcf = vcf.replace(".vcf",".good.vcf")
        command = ' '.join([minos_binary, 'check_with_ref', filtered_vcf, vcf_ref, truth, "tmp/individual." + name, "--flank_length",  str(flank), "--variant_merge_length",  str(flank), "--include_ref_calls", "--max_soft_clipped", str(20)])
        if mask:
            command += " --exclude_bed " + mask
        syscall(command)

def create_or_append(ultimate_file, data_file):
    command = ""
    if os.path.isfile(ultimate_file):
        command = ' '.join(['tail -n+2', data_file, '>>', ultimate_file])
    else:
        command = ' '.join(['cp', data_file, ultimate_file])
    syscall(command)

def append_stats_x(genotyper_name, run_name):
    create_or_append(genotyper_name + ".x.gt_conf_hist.TP.tsv", "tmp/individual." + run_name + '.gt_conf_hist.TP.tsv')
    create_or_append(genotyper_name + ".x.gt_conf_hist.FP.tsv", "tmp/individual." + run_name + '.gt_conf_hist.FP.tsv')
    create_or_append(genotyper_name + ".x.stats.tsv", "tmp/individual." + run_name + '.stats.tsv')

def append_stats_y(genotyper_name, run_name):
    create_or_append(genotyper_name + ".y.gt_conf_hist.tsv", "tmp/compare." + run_name + ".gt_conf_hist.tsv")
    create_or_append(genotyper_name + ".y.stats.tsv", "tmp/compare." + run_name + ".stats.tsv")

def minos_to_df(name):
    if len(glob.glob(name + ".csv")) != 0:
        print("Already have", name, ".csv")
        return

    print("Construct dataframe for", name)
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
        if dnadiff_frac < 0.001 or confidence == 0:
            continue
        yscat.append(dnadiff_frac)
        sum_fp = sum(x_fp[x_fp['GT_CONF'] >= confidence]['Count'])
        sum_tp = sum(x_tp[x_tp['GT_CONF'] >= confidence]['Count'])
        if sum_fp > 0 and sum_tp > 0:
            xscat.append(sum_fp/float(sum_fp + sum_tp))
            #if len(yscat) > 2 and (yscat[-1] < 0.6 < yscat[-2] or yscat[-2] < 0.6 < yscat[-1]):
            #    print(name, xscat[-1], yscat[-1])
        else:
            xscat.append(float(0))

    df = pd.DataFrame()
    df['yscat'] = yscat
    df['xscat'] = xscat
    df['name'] = name
    print(df)
    print("Save df to file")
    df.to_csv(name + '.csv')
    if len(yscat) > 0:
        print("Max", name, max(yscat))

def unzip_truths(truths):
    unzipped_truths = []
    for truth in truths:
        new_name = ""
        if truth.endswith(".gz"):
            new_name = truth.split('/')[-1][:-3]
            command = ' '.join(['zcat', truth, "| awk '{print $1;}' > tmp.fa"])
            syscall(command)
        else:
            new_name = truth.split('/')[-1]
            command = ' '.join(['cat', truth, "| awk '{print $1;}' > tmp.fa"])
            syscall(command)
        command = ' '.join(['head -n1 tmp.fa >', new_name])
        syscall(command)
        command = ' '.join(["grep -v '>' tmp.fa >>", new_name])
        syscall(command)
        command = ' '.join(['rm tmp.fa'])
        syscall(command)
        unzipped_truths.append(new_name)
    return unzipped_truths

def compare_pair(pair, names, recall_flank):
    (num1, id1, zipped_truth1, sample_dir1, mask1, truth1) = pair[0]
    (num2, id2, zipped_truth2, sample_dir2, mask2, truth2) = pair[1]
    print("Run on sample pair", id1, id2)
    if not os.path.isfile("tmp/dnadiff_%s_%s.snps" %(id1, id2)):
        run_dnadiff(truth1, truth2, id1 + "_" + id2)
    pair_names = names[:]
    pairs = get_vcf_pairs(sample_dir1, sample_dir2)
    for run in pairs:
        name, vcf_ref, vcf1, vcf2 = run
        print("Looking at ref", name, id1, id2)
        pair_names.remove(name)
        comp_name = name + "_" + id1 + "_" + id2
        if len(glob.glob("tmp/compare." + comp_name + "*")) == 0:
            run_compare("tmp/dnadiff_%s_%s.snps" %(id1, id2), truth1, truth2, vcf1, vcf2, vcf_ref, comp_name, recall_flank, mask1, mask2)
            append_stats_y(name, comp_name)
    for name in pair_names:
        print("Add null compare stats for", name, id1, id2)
        comp_name = name + "_" + id1 + "_" + id2
        if len(glob.glob("tmp/compare." + comp_name + "*")) == 0:
            with open("tmp/compare." + comp_name + ".stats.tsv", 'w') as f:
                stats = {x: 0 for x in ['total', 'found_vars', 'missed_vars', 'excluded_vars']}
                with open("tmp/dnadiff_%s_%s.snps" %(id1, id2),'r') as g:
                    for line in g:
                        stats['total'] += 1
                        stats['missed_vars'] += 1
                keys = stats.keys()
                print(*keys, sep='\t', file=f)
                print(*[stats[x] for x in keys], sep='\t', file=f)
            with open("tmp/compare." + comp_name + ".gt_conf_hist.tsv", 'w') as f:
                print('GT_CONF\tCount', file=f)
            append_stats_y(name, comp_name)    

def run_individual(num1, id1, zipped_truth1, sample_dir1, mask1, truth1, precision_flank):
    vcfs = get_vcfs(sample_dir1)
    for run in vcfs:
        name, vcf_ref, vcf1 = run
        print("run individual ", name, vcf1)
        indiv_name = name + "_" + id1
        if len(glob.glob("tmp/individual." + indiv_name + "*")) == 0:
            run_individual_minos(truth1, vcf1, vcf_ref, indiv_name, precision_flank, mask1)
            print("append stats")
            append_stats_x(name, indiv_name)
        #files = glob.glob("tmp/individual*")
        #for f in files:
        #      os.unlink(f)

parser = argparse.ArgumentParser(description='Identifies pairs of VCF files called against the same reference in 2 sample directories and runs a minos system call to evaluate how many dnadiff snp differences are captued between each pair of VCFs.')
parser.add_argument('--sample_tsv', '-t1', type=str,
                    help='TSV with a line for each sample, containing ID, TRUTH_FASTA, VCF_DIRECTORY, MASK')
parser.add_argument('--recall_flank', '-fr', type=int, default=11,
                    help='Size of flank sequence to use when comparing true alleles to vcf alleles')
parser.add_argument('--precision_flank', '-fp', type=int, default=31,
                    help='Size of flank sequence to use when comparing alleles to truth assembly')
parser.add_argument('--num_threads', '-t', type=int, default=1,
                    help='Number of parallel process threads')
args = parser.parse_args()

index = pd.read_csv(args.sample_tsv, sep='\t', header=None, names=['id', 'truth', 'vcf_dir', 'mask'])
index['unzipped_truth'] = unzip_truths(index['truth'])

names = get_all_names(index['vcf_dir'].values)

if not os.path.exists('tmp'):
    os.mkdir("tmp")

# y
pairs = [[[list(p[0]), list(p[1])], names, args.recall_flank] for p in itertools.combinations(index.itertuples(), 2)]
print(pairs)
with Pool(args.num_threads) as pool:
    r = pool.starmap(compare_pair, pairs) 

# x
samples = [list(i) for i in index.itertuples()]
for s in samples:
    s.append(args.precision_flank)
print(samples)
with Pool(args.num_threads) as pool:
    r = pool.starmap(run_individual, samples)

for truth in index['unzipped_truth']:
    files = glob.glob(truth + ".*")
    for f in files:
        os.unlink(f)
    
for name in names:
    print("Try to get df for", name)
    minos_to_df(name)
