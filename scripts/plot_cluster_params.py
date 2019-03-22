import matplotlib.pyplot as plt
import argparse
import pandas as pd
import numpy as np
import matplotlib

def plot_both_param_cluster(tsv_file, param1, param2, d1, d2, l1, l2):
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', param1, param2, 'covg','tp', 'fn', 'fp', 'precision', 'recall'])
    df['type'] = [i.replace("_ecoli","") for i in df['type'].values]
    df['type'] = [i.replace("_10000","") for i in df['type'].values]
    df['type'] = [i.replace("_50000","") for i in df['type'].values]
    
    plot_param_cluster(df, param1, param2, l1, d2)
    plot_param_cluster(df, param2, param1, l2, d1)
    
def plot_param_cluster(df, param1, param2, xlabel, default=0):
    plt.rcParams['figure.figsize'] = 25,10
    fig, axs = plt.subplots(1, 3)
    plt.style.use('seaborn-colorblind')
    
    # create a horizontal plot
    host_col = 'blue'
    par1_col = 'green'

    paxs = []
    for i in range(3):
        axs[i].set_xlabel(xlabel)
        axs[i].grid(b=True)
        axs[i].set_axisbelow(b=True)
        pax = axs[i].twinx()
        paxs.append(pax)
        axs[i].set_ylabel("Recall")
        paxs[i].set_ylabel("Precision")
        fig.subplots_adjust(right=2.2)

        for ax in [paxs[i]]:
            ax.set_frame_on(True)
            ax.patch.set_visible(False)

            plt.setp(ax.spines.values(), visible=False)
            ax.spines["right"].set_visible(True)

        axs[i].yaxis.label.set_color(host_col)
        paxs[i].yaxis.label.set_color(par1_col)

        axs[i].spines["left"].set_edgecolor(host_col)
        paxs[i].spines["right"].set_edgecolor(par1_col)

        axs[i].tick_params(axis='y', colors=host_col)
        paxs[i].tick_params(axis='y', colors=par1_col)


    order=['illumina_HS25', 'nanopore_R9_1D', 'nanopore_R9_2D']
    order.sort()
    p2s=[default]
    if default == 0:
        p2s=list(set(list(df[param2].values)))
    p2s.sort()
    markers=['o','s','v','^','>','+','x', '<', 'o', 'o', 'o', 'o']

    for i in range(len(order)):
        print(order[i])
        for j in p2s:
            indx = df['type'] == order[i]
            dfs = df[indx]
            indx = dfs[param2] == j
            dfs = dfs[indx]
            dfs = dfs.sort_values(by=[param1])

            axs[i].plot(dfs[param1], dfs['recall'], 'o-', c=host_col, marker=markers[i], alpha=0.8, label=order[i])
            paxs[i].plot(dfs[param1], dfs['precision'], 'o-', c=par1_col, marker=markers[i], alpha=0.8, label=order[i])
            axs[i].set_title(order[i])

    fig.tight_layout()
    plt.savefig('gene_finding_%s.png' %param1, transparent=True, bbox_inches = 'tight')


parser = argparse.ArgumentParser(description='Plots a scatter for precision and recall of different technologies')
parser.add_argument('--tsv', type=str,
                    help='Input TSV')
parser.add_argument('--p1', type=str, default="min_cluster_size",
                    help='Parameter1')
parser.add_argument('--p2', type=str, default="max_cluster_distance",
                    help='Parameter2')
parser.add_argument('--l1', type=str, default="Minimum cluster size",
                    help='Parameter1 label')
parser.add_argument('--l2', type=str, default="Maximum within cluster distance",
                    help='Parameter2 label')
parser.add_argument('--d1', type=int, default=9,
                    help='Parameter1 default')
parser.add_argument('--d2', type=int, default=80,
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

plot_both_param_cluster(args.tsv, args.p1, args.p2, args.d1, args.d2, args.l1, args.l2)
