#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import seaborn as sns
from collections import Counter
import pandas as pd
import numpy as np

def plot_df_covg(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', 'covg', 'tp', 'fn', 'fp', 'precision', 'recall'])
    print(df.head)
    df['coverage'] = [int(i) for i in df['covg']]
    fixed_names = [s.replace("_ecoli","") for s in df['type']]
    df['type'] = fixed_names
    plt.rcParams['figure.figsize'] = 10,6
    fig, ax = plt.subplots()

    ax.grid(b=True)
    ax.set_axisbelow(b=True)
    plt.style.use('seaborn-colorblind')

    hue_order=list(set(list(df['type'].values)))
    hue_order.sort()
    #hue_order=['illumina_MSv3_150', 'illumina_MSv3_250','illumina_HS25_150', 'nanopore_R9_2D_10000', 'nanopore_R9_2D_50000']
    markers=['v','^','>','o','s','+','x', '<']
    gs = []
    for i in range(len(hue_order)):
        indx = df['type'] == hue_order[i]
        dfs = df[indx]
        g = plt.scatter(dfs['precision'], dfs['recall'], s=dfs['coverage'], marker=markers[i], alpha=0.8)
        gs.append(g)
    ax.set(xlabel="Precision", ylabel="Recall")
    plt.legend(gs, hue_order, loc='lower left', bbox_to_anchor=(1, -0.0))
    
    plt.savefig('gene_finding_covg.png', transparent=True)

parser = argparse.ArgumentParser(description='Plots a scatter for precision and recall of different technologies')
parser.add_argument('--tsv', type=str,
                    help='Input TSV')
args = parser.parse_args()

SMALL_SIZE = 16
MEDIUM_SIZE = 24
BIGGER_SIZE = 30

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plot_df_covg(args.tsv)
