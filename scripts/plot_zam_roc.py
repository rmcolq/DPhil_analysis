#!/usr/bin/env python 

import argparse
import glob
import operator
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_dfs(directory):
    dfs = []
    files = glob.glob("%s/*.csv" %directory)
    for f in files:
        df = pd.read_csv(f, index_col=0)
        print(df)
        dfs.append(df)
    return dfs

def plot_graphs(items, x_max, y_label):
    for i in range(7):
        # Define plot
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 12)

        ax.grid(b=True)
        ax.set_axisbelow(b=True)
        plt.style.use('seaborn-colorblind')

        # Label the axes and title the plot
        ax.set_xlabel('Number FPs/Number Genotyped')
        ax.set_ylabel(y_label)
        if x_max > 0:
            plt.xlim(-0.0005, x_max)

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        snippy_i = 0
        nanopolish_i = 0

        # Make a scatter plot
        for x in items:
            if len(x['name'].values) > 0:
                if "bin" in x['name'].values[0]:
                    continue
                elif x['name'].values[0].startswith("pandora_genotyped_full"):
                    x['name'] = "pandora_genotyped_nanopore_full"
                elif x['name'].values[0].startswith("pandora_genotyped_30"):
                    x['name'] = "pandora_genotyped_nanopore_30"
                print("add df", x['name'].values[0], "to graph")
                if (x['name'].values[0].startswith("pandora_genotyped_nanopore_full") or x['name'].values[0].startswith("pandora_recall")) and i > 0:
                    col = colormap_pandora(75)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina") and i > 5:
                    col = colormap_pandora(200)
                    if "ref" in x['name'].values[0] or "bin" in x['name'].values[0]:
                        col = colormap_pandora(225)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina_30") and i > 4:
                    col = colormap_pandora(250)
                    if "ref" in x['name'].values[0] or "bin" in x['name'].values[0]:
                        col = colormap_pandora(275)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_nanopore_30"):
                    col = colormap_pandora(0)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_nanopore_discovery") and i > 3:
                    col = colormap_pandora(150)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("snippy") and i > 1:
                    col = colormap_snippy(snippy_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    snippy_i += 1
                elif x['name'].values[0].startswith("nanopolish") and i > 2:
                    col = colormap_nanopolish(nanopolish_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    nanopolish_i += 1
                else:
                    print("failed to add to scatter")
            else:
                print(x['name'].values)

        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0 and len(labels) > 0 :
            hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
            handles2, labels2 = zip(*hl)
            ax.legend(handles2, labels2, frameon=False, loc='lower right')
            plt.savefig('zam_roc%d.png' %i, transparent=True)
        else:
            print(handles, labels)

parser = argparse.ArgumentParser(description='Plots all csv dataframes on axes')
parser.add_argument('--directory', type=str, default=".",
                    help='Directory containing CSV files')
parser.add_argument('--x_max', type=float, default=0,
                    help='Maximum x-axis')
parser.add_argument('--y_label', type=str, default='Fraction of dnadiff SNPs discoverable from VCFs',
                    help='Label for y axis')
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

dfs = load_dfs(args.directory)
plot_graphs(dfs, args.x_max, args.y_label)
