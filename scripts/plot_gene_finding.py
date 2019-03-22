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

def plot_df_param1(tsv_file, param1, param2, xlabel, default2):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()

    host_col = 'blue'
    par1_col = 'green'

    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', param1, param2, 'covg','tp', 'fn', 'fp', 'precision', 'recall'])
    df['type'] = [i.replace("_ecoli","") for i in df['type'].values]
    df['type'] = [i.replace("_10000","") for i in df['type'].values]
    df['type'] = [i.replace("_50000","") for i in df['type'].values]
    #order=list(set(list(df['type'].values)))
    order=['illumina_HS25', 'nanopore_R9_1D', 'nanopore_R9_2D']
    order.sort()
    p2s=[default2]
    p2s.sort()
    markers=['o','s','v','^','>','+','x', '<', 'o', 'o', 'o', 'o']
    g_mark = []
    g_host = []
    g_par1 = []
    for i in range(len(order)):
        print(order[i])
        for j in p2s:
            indx = df['type'] == order[i]
            dfs = df[indx]
            indx = dfs[param2] == j
            dfs = dfs[indx]
            dfs = dfs.sort_values(by=[param1])
            g = host.scatter(dfs[param1], dfs['recall'], c='grey', marker=markers[i], alpha=0.8, label=order[i])
            g_mark.append(g)
            g = host.plot(dfs[param1], dfs['recall'], 'o-', c=host_col, marker=markers[i], alpha=0.8, label=order[i])
            g_host.append(g)
            g = par1.plot(dfs[param1], dfs['precision'], 'o-', c=par1_col, marker=markers[i], alpha=0.8, label=order[i])
            g_par1.append(g)
    host.set_xlabel(xlabel)
    host.set_ylabel("Recall")
    par1.set_ylabel("Precision")

    lines = g_mark
    host.legend(lines, [l.get_label() for l in lines], bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)

    for ax in [par1]:
        ax.set_frame_on(True)
        ax.patch.set_visible(False)

        plt.setp(ax.spines.values(), visible=False)
        ax.spines["right"].set_visible(True)

    host.yaxis.label.set_color(host_col)
    par1.yaxis.label.set_color(par1_col)

    host.spines["left"].set_edgecolor(host_col)
    par1.spines["right"].set_edgecolor(par1_col)

    host.tick_params(axis='y', colors=host_col)
    par1.tick_params(axis='y', colors=par1_col)
    plt.savefig('gene_finding_%s.png' %param1, transparent=True)
    
def plot_df_param2(tsv_file, param1, param2, xlabel, default1):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()

    host_col = 'blue'
    par1_col = 'green'

    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', param1, param2, 'covg','tp', 'fn', 'fp', 'precision', 'recall'])
    df['type'] = [i.replace("_ecoli","") for i in df['type'].values]
    df['type'] = [i.replace("_10000","") for i in df['type'].values]
    df['type'] = [i.replace("_50000","") for i in df['type'].values]
    #order=list(set(list(df['type'].values)))
    order=['illumina_HS25', 'nanopore_R9_1D', 'nanopore_R9_2D']
    order.sort()
    p1s=[default1]
    p1s.sort()
    markers=['o','s','v','^','>','+','x', '<', 'o', 'o', 'o', 'o']
    g_mark = []
    g_host = []
    g_par1 = []
    for i in range(len(order)):
        print(order[i])
        for j in p1s:
            indx = df['type'] == order[i]
            dfs = df[indx]
            indx = dfs[param1] == j
            dfs = dfs[indx]
            dfs = dfs.sort_values(by=[param2])
            g = host.scatter(dfs[param2], dfs['recall'], c='grey', marker=markers[i], alpha=0.8, label=order[i])
            g_mark.append(g)
            g = host.plot(dfs[param2], dfs['recall'], 'o-', c=host_col, marker=markers[i], alpha=0.8, label=order[i])
            g_host.append(g)
            g = par1.plot(dfs[param2], dfs['precision'], 'o-', c=par1_col, marker=markers[i], alpha=0.8, label=order[i])
            g_par1.append(g)
    host.set_xlabel(xlabel)
    host.set_ylabel("Recall")
    par1.set_ylabel("Precision")

    lines = g_mark
    host.legend(lines, [l.get_label() for l in lines], bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)

    for ax in [par1]:
        ax.set_frame_on(True)
        ax.patch.set_visible(False)

        plt.setp(ax.spines.values(), visible=False)
        ax.spines["right"].set_visible(True)

    host.yaxis.label.set_color(host_col)
    par1.yaxis.label.set_color(par1_col)

    host.spines["left"].set_edgecolor(host_col)
    par1.spines["right"].set_edgecolor(par1_col)

    host.tick_params(axis='y', colors=host_col)
    par1.tick_params(axis='y', colors=par1_col)
    plt.savefig('gene_finding_%s.png' %param2, transparent=True)


parser = argparse.ArgumentParser(description='Plots a scatter for precision and recall of different technologies')
parser.add_argument('--tsv', type=str,
                    help='Input TSV')
parser.add_argument('--p1', type=str, default="",
                    help='Parameter1')
parser.add_argument('--p2', type=str, default="",
                    help='Parameter2')
parser.add_argument('--l1', type=str, default="",
                    help='Parameter1 label')
parser.add_argument('--l2', type=str, default="",
                    help='Parameter2 label')
parser.add_argument('--d1', type=int, default=0,
                    help='Parameter1 default')
parser.add_argument('--d2', type=int, default=0,
                    help='Parameter2 default')

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

if args.p1=="" and args.p2=="":
    plot_df_covg(args.tsv)
else:
    plot_df_param1(args.tsv, args.p1, args.p2, args.l1, args.d2)
    plot_df_param2(args.tsv, args.p1, args.p2, args.l2, args.d1)
