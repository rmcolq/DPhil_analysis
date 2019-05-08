#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def plot_count_hist(files, prefix="", label_file=""):
    # Load labels if we have them
    label_dict = {}
    if label_file != "":
        with open(label_file, "r") as handle:
            for line in handle:
                i,lab = line.strip().split('\t')
                label_dict[i] = lab
    print(label_dict)

    print("parameters")
    print(files, prefix)
    count_dict = {}
    for filepath in files:
        print(filepath)
        with open(filepath, 'r') as f:
            line = f.readline()
            stype = line.split('\t')[0]
            counts = line.split('\t')[1].split(',')
            icounts = [int(i) for i in counts if len(i) > 0]
            icounts.sort()
            count_dict[stype] = icounts

    plt.rcParams['figure.figsize'] = 10,6
    fig, ax = plt.subplots()
    plt.style.use('seaborn-colorblind')
    ax.grid(b=True)
    ax.set_axisbelow(b=True)
    types = list(count_dict.keys())
    print(count_dict.keys())
    print(types)
    ax.hist([count_dict[t] for t in types],
          density=True,
          label=[label_dict[t] for t in types],
          align="mid",
          bins=range(0, 42, 1))

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel='Number of mismatch bases', ylabel='Frequency')
    plt.savefig(prefix + 'joint.sam_mismatch_counts.png', transparent=True)

    fig, ax = plt.subplots()
    plt.style.use('seaborn-colorblind')
    plt.grid(b=True, which='major', color='LightGrey', linestyle='-')
    plt.grid(b=True, which='minor', color='GhostWhite', linestyle='-')
    ax.set_axisbelow(b=True)
    ms = [max(count_dict[t]) for t in types]
    m = max(ms)
    n, bins, patches = ax.hist([count_dict[t] for t in types], m+1, density=True, histtype='step', lw=1.5, label=[label_dict[t] for t in types], cumulative=True, range=(0,15))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=4, frameon=False)
    ax.set(xlabel='Maximum number of mismatch bases', ylabel='Frequency')
    plt.savefig(prefix + 'joint.step_sam_mismatch_counts.png', transparent=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots a histogram of the number of mismatches from a pair of count files.')
    parser.add_argument('--f1', type=str, default="",
                    help='Count file 1')
    parser.add_argument('--f2', type=str, default="",
                    help='Count file 2')
    parser.add_argument('--f', type=str, default="",
                    help='Count files')
    parser.add_argument('--prefix', type=str, default="",
                    help='String to add to outfile name')
    parser.add_argument('--labels', type=str, default="",
                    help='tsvfile of labels for plot')
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
    
    files = []
    if args.f != "":
        files = args.f.split()
    elif args.f1 != "" and args.f2 != "":
        files.append(args.f1)
        files.append(args.f2)
    if len(files) == 0:
        print("No files provided")
        quit()

    plot_count_hist(files, args.prefix, args.labels)
