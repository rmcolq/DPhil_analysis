#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
import pandas as pd
import numpy as np
import argparse

def get_nums_from_cigar(cigar):
    start = 0
    numbers = []
    for end in range(len(cigar)):
        if cigar[end].isalpha():
            numbers.append(int(cigar[start:end]))
            start = end+1
    return numbers

def plot_sam_dist(sam_file, out_prefix):
    total_gene_bases = 0
    total_mismatch_bases = 0
    num_false_positives = 0
    num_other = 0  
    num_genes = 0
    num_mismatches = []
    lengths = []

    f = open(sam_file, 'r')
    all_names = []
    non_unique_names = []
    for line in f:
        if line != "" and line[0].isalpha():
            name = line.split('\t')[0]
            if name in all_names:
                non_unique_names.append(name)
            else:
                all_names.append(name)
    f.close()

    f = open(sam_file, 'r')
    for line in f:
        if line == "" or not line[0].isalpha() or line.split('\t')[0] in non_unique_names:
            continue
        elif line.split('\t')[1]=="0" or line.split('\t')[1]=="16":
            num = int(line.split('\t')[11].split("NM:i:")[-1])
            num_mismatches.append(num)
            lengths.append(sum(get_nums_from_cigar(line.split('\t')[5])))
            total_mismatch_bases += num
            total_gene_bases += len(line.split('\t')[9])
            num_genes += 1
        elif line.split('\t')[1]!="4":
            num_false_positives += 1
            num_genes += 1
    f.close()

    if total_gene_bases == 0:
        total_gene_bases = 1
    print("False positive gene rate: %d/%d = %f" %(num_false_positives, num_genes, float(num_false_positives)/num_genes*100))
    print("Estimated per base accuracy: 1 - %d/%d = %f" %(total_mismatch_bases, total_gene_bases, (1-(float(total_mismatch_bases)/total_gene_bases))*100))
    print(Counter(num_mismatches))
    plt.rcParams['figure.figsize'] = 10,6
    fig, ax = plt.subplots()
    sns.distplot(num_mismatches, bins=range(0, min(40, max(num_mismatches)+1)), kde=False, norm_hist=True)
    ax.set(xlabel='Number of mismatch bases', ylabel='Frequency')
    plt.grid(b=True, which='major', color='LightGrey', linestyle='-')
    plt.grid(b=True, which='minor', color='GhostWhite', linestyle='-')
    plt.savefig('%s.sam_mismatch_counts.png' %out_prefix, transparent=True)

    fig, ax = plt.subplots()
    plt.scatter(lengths, num_mismatches)
    ax.set(xlabel='Length of sequence', ylabel='Number of mismatch bases')
    plt.grid(b=True, which='major', color='LightGrey', linestyle='-')
    plt.grid(b=True, which='minor', color='GhostWhite', linestyle='-')
    plt.savefig('%s.sam_length_scatter.png' %out_prefix, transparent=True)

    with open('%s_list.txt' %out_prefix, 'w') as f:
        f.write("%s\t" % out_prefix)
        for item in num_mismatches:
            f.write("%s," % item)

parser = argparse.ArgumentParser(description='Plots a histogram of the number of mismatches from a samfile.')
parser.add_argument('--sam', type=str,
                    help='Input SAM')
parser.add_argument('--prefix', type=str,
                    help='Output prefix')
args = parser.parse_args()

plot_sam_dist(args.sam, args.prefix)

