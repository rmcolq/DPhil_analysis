#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def plot_count_hist(file1, file2, prefix):
    print("parameters")
    print(file1, file2, prefix)
    count_dict = {}
    for filepath in file1, file2:
        print(filepath)
        with open(filepath, 'r') as f:
            line = f.readline()
            stype = line.split('\t')[0].split('_')[-1]
            counts = line.split('\t')[1].split(',')
            icounts = [int(i) for i in counts if len(i) > 0]
            count_dict[stype] = icounts

    plt.rcParams['figure.figsize'] = 10,6
    fig, ax = plt.subplots()
    plt.style.use('seaborn-deep')
    types = list(count_dict.keys())
    print(count_dict.keys())
    print(types)
    ax.hist([count_dict[types[0]], count_dict[types[1]]],
          density=True,
          label=types,
          align="mid",
          bins=range(0, 42, 1))

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel='Number of mismatch bases', ylabel='Frequency')
    plt.savefig(prefix + 'joint.sam_mismatch_counts.png', transparent=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots a histogram of the number of mismatches from a pair of count files.')
    parser.add_argument('--f1', type=str,
                    help='Count file 1')
    parser.add_argument('--f2', type=str,
                    help='Count file 2')
    parser.add_argument('--prefix', type=str, default="",
                    help='String to add to outfile name')
    args = parser.parse_args()

    plot_count_hist(args.f1, args.f2, args.prefix)
