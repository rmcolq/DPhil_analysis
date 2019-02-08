#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import seaborn as sns
from collections import Counter
import pandas as pd
import numpy as np

def plot_df(tsv_file, xvar, xlabel):
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', xvar, 'systime', 'usertime', 'max_mem'])
    df['max_mem_gb'] = df['max_mem']/(1024*1024)
    df['time'] = (df['systime'] + df['usertime'])/(60*60)
    plt.rcParams['figure.figsize'] = 10,6
    fig, ax = plt.subplots()

    ax.grid(b=True)
    ax.set_axisbelow(b=True)
    plt.style.use('seaborn-colorblind')

    g = sns.lmplot( x=xvar, y="time", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
    g = (g.set_axis_labels(xlabel, "Time (h)"))
    plt.legend(loc='lower right')
    plt.savefig('scaling_time.png', transparent=True)

    g = sns.lmplot( x=xvar, y="max_mem_gb", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
    g = (g.set_axis_labels(xlabel, "Max Memory (GB)"))
    plt.legend(loc='lower right')
    plt.savefig('scaling_mem.png', transparent=True)

parser = argparse.ArgumentParser(description='Plots how a variable scales with time and memory.')
parser.add_argument('--tsv_file', type=str,
                    help='TSV file of results')
parser.add_argument('--xvar', type=str,
                    help='Variable name for x')
parser.add_argument('--xlabel', type=str,
                    help='X-axis label')
args = parser.parse_args()

plot_df(args.tsv_file, args.xvar, xlabel)
