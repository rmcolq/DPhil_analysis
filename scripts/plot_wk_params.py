import matplotlib.pyplot as plt
import pandas as pd
import argparse

def plot_both_param_wk(tsv_file, param1, param2):
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', param1, param2, 'covg','tp', 'fn', 'fp', 'precision', 'recall'])
    df['type'] = [i.replace("_ecoli","") for i in df['type'].values]
    df['type'] = [i.replace("_10000","") for i in df['type'].values]
    df['type'] = [i.replace("_50000","") for i in df['type'].values]
    
    plot_param_wk(df, param1, param2)

def plot_both_param_wk_times(tsv_file, param1, param2):
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', param1, param2,'utime', 'ctime', 'mem'])
    df['type'] = [i.replace("_ecoli","") for i in df['type'].values]
    df['type'] = [i.replace("_10000","") for i in df['type'].values]
    df['type'] = [i.replace("_50000","") for i in df['type'].values]
    df['time'] = df['utime'] + df['ctime']
    df['time'] = [float(f)/60 for f in df['time'].values]
    df['memory'] = [float(f)/(1024*1024) for f in df['mem'].values]

    plot_param_wk(df, param1, param2, 'time', 'memory', 'Time (m)', 'Max memory (GB)')
    
def plot_param_wk(df, param1, param2, y1='recall', y2='precision', yl1="Recall", yl2="Precision"):
    markers=['o','s','v','^','>','+','x', '<', 'o', 'o', 'o', 'o']
    order=['illumina_HS25', 'nanopore_R9_1D', 'nanopore_R9_2D']
    #order=["Illumina"]
    order.sort()
    p2s=list(set(list(df[param2].values)))
    p2s.sort()
    p1s=list(set(list(df[param1].values)))
    p1s.sort()
    
    for i in range(len(order)):
        plt.rcParams['figure.figsize'] = 15,5
        fig, axs = plt.subplots(1, 2)
        plt.style.use('seaborn-colorblind')
        
        axs[0].set_xlabel('W')
        axs[0].set_ylabel(yl2)
        axs[0].grid(b=True)
        axs[0].set_axisbelow(b=True)
        axs[1].set_xlabel('W')
        axs[1].set_ylabel(yl1)
        axs[1].grid(b=True)
        axs[1].set_axisbelow(b=True)

        indx = df['type'] == order[i]
        dft = df[indx]
        for j in p2s:
            indx = dft[param2] == j
            dfs = dft[indx]
            dfs = dfs.sort_values(by=[param1])        
            axs[0].plot(dfs[param1], dfs[y1], 'o-', marker=markers[i], alpha=0.8, label="%s k=%i" %(order[i], j))
            axs[1].plot(dfs[param1], dfs[y2], 'o-', marker=markers[i], alpha=0.8, label="%s k=%i" %(order[i], j))
            axs[0].set_title(order[i])
            axs[1].set_title(order[i])
            axs[1].legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)
        fig.tight_layout()
        plt.savefig('gene_finding_wk_%s_%s_%s.png' %(y1,y2,order[i]), transparent=True, bbox_inches = 'tight')


parser = argparse.ArgumentParser(description='Plots a scatter for precision and recall of different technologies')
parser.add_argument('--tsv', type=str,
                    help='Input TSV')
parser.add_argument('--p1', type=str, default="w",
                    help='Parameter1')
parser.add_argument('--p2', type=str, default="k",
                    help='Parameter2')
parser.add_argument('--tsv_times', type=str,default="",
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

plot_both_param_wk(args.tsv, args.p1, args.p2)
if args.tsv_times != "":
    plot_both_param_wk_times(args.tsv_times, args.p1, args.p2)
