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
        df = pd.read_csv(f)
        print(df)
        dfs.append(df)
    return dfs

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

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        snippy_i = 0
        nanopolish_i = 0

        # Make a scatter plot
        for x in items:
            if len(x['name'].values) > 0:
                print("add df", x['name'].values[0], "to graph")
                if x['name'].values[0].startswith("pandora_genotyped_nanopore_full") or x['name'].values[0].startswith("pandora_recall"):
                    col = colormap_pandora(0)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina") and i > 0:
                    col = colormap_pandora(100)
                    ax.scatter(x['xscat'], x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_nanopore_30") and i > 1:
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
                    print("failed to add to scatter")
            else:
                print(x['name'].values)

        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0 and len(labels) > 0 :
            hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
            handles2, labels2 = zip(*hl)
            ax.legend(handles2, labels2, frameon=False, loc='lower right')
            plt.savefig('roc%d.png' %i, transparent=True)
        else:
            print(handles, labels)

parser = argparse.ArgumentParser(description='Plots all csv dataframes on axes')
parser.add_argument('--directory', type=str, default=".",
                    help='Directory containing CSV files')
args = parser.parse_args()

dfs = load_dfs(args.directory)
plot_graphs(dfs)
