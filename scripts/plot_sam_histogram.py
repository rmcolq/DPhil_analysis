#!/usr/bin/env python3

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter
import pandas as pd
import numpy as np
import argparse

def get_nums_from_cigar(cigar):
    start = 0
    numbers = []
    letters = []
    for end in range(len(cigar)):
        if cigar[end].isalpha():
            numbers.append(int(cigar[start:end]))
            letters.append(cigar[end])
            start = end+1
    return numbers, letters

def plot_sam_dist(sam_file, out_prefix, unique=False):
    total_gene_bases = 0
    total_mismatch_bases = 0
    num_false_positives = 0
    num_other = 0  
    num_genes = 0
    num_mismatches = []
    lengths = []
    indel = []
    largest_indel = []

    f = open(sam_file, 'r')
    all_names = []
    non_unique_names = []
    for line in f:
        if line != "" and line[0]!="@":#line[0].isalpha():
            name = line.split('\t')[0]
            if unique and name in all_names:
                non_unique_names.append(name)
                print(name, [i for i in all_names if i==name])
            else:
                all_names.append(name)
    print("all names", len(all_names))
    print("non unique names", len(non_unique_names))
    f.close()

    f = open(sam_file, 'r')
    for line in f:
        if line == "" or not line[0].isalpha() or line.split('\t')[0] in non_unique_names:
            continue
        elif line.split('\t')[1]=="0" or line.split('\t')[1]=="16":
            num = int(line.split('\t')[11].split("NM:i:")[-1])
            num_mismatches.append(num)
            cigar = line.split('\t')[5]
            nums, letters = get_nums_from_cigar(cigar)
            lengths.append(sum(nums))
            largest = 0
            for i,letter in enumerate(letters):
                if letter == "D" and nums[i] > abs(largest):
                    largest = nums[i]
                elif letter == "I" and nums[i] > abs(largest):
                    largest = nums[i]
            largest_indel.append(largest)
            #if num > 0 and largest > 0:
            #    largest_indel.append(largest/num)
            #else:
            #    largest_indel.append(0)
            if "D" in cigar:
                indel.append("deletion")
            elif "I" in cigar:
                indel.append("insertion")
            else:
                indel.append("match")
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
    df = pd.DataFrame()
    df['lengths']=lengths
    df['num_mismatches']=num_mismatches
    df['indel']=indel
    df['largest_indel'] = largest_indel
    
    uniq=['match','insertion','deletion']
    z = range(1,len(uniq))
    hot = plt.get_cmap('hot')
    cNorm  = colors.Normalize(vmin=0, vmax=len(uniq))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=hot)

    # Plot each scatter
    for i in range(len(uniq)):
        indx = df['indel'] == uniq[i]
        plt.scatter(df['lengths'][indx], df['num_mismatches'][indx], c=scalarMap.to_rgba(i), label=uniq[i])
    ax.set(xlabel='Length of sequence', ylabel='Number of mismatch bases')
    plt.legend(loc='upper left')
    plt.grid(b=True, which='major', color='LightGrey', linestyle='-')
    plt.grid(b=True, which='minor', color='GhostWhite', linestyle='-')
    plt.savefig('%s.sam_length_scatter.png' %out_prefix, transparent=True)

    fig, ax = plt.subplots()
    s = plt.scatter(df['lengths'], df['num_mismatches'], c=df['largest_indel'], cmap='hot')
    ax.set(xlabel='Length of sequence', ylabel='Number of mismatch bases')
    max_diffs = max(df['num_mismatches'].values)
    s.set_clim([0,max_diffs])
    cb = fig.colorbar(s)
    plt.grid(b=True, which='major', color='LightGrey', linestyle='-')
    plt.grid(b=True, which='minor', color='GhostWhite', linestyle='-')
    plt.savefig('%s.sam_length_scatter2.png' %out_prefix, transparent=True)

    with open('%s_list.txt' %out_prefix, 'w') as f:
        f.write("%s\t" % out_prefix)
        for item in num_mismatches:
            f.write("%s," % item)

parser = argparse.ArgumentParser(description='Plots a histogram of the number of mismatches from a samfile.')
parser.add_argument('--sam', type=str,
                    help='Input SAM')
parser.add_argument('--prefix', type=str,
                    help='Output prefix')
parser.add_argument('--unique', action='store_true',
                    help='If used, filters only unique sam matches')
args = parser.parse_args()

plot_sam_dist(args.sam, args.prefix, args.unique)

