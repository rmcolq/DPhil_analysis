import pandas as pd
import vcf
import matplotlib.pyplot as plt
import seaborn as sns
import math
import csv
import numpy as np
import pickle
import math
from scipy.optimize import least_squares
import argparse
import os.path

def save_lols(lol, filepath):
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(lol)

def load_lols(filepath):
    lol = []
    with open(filepath, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            lol.append(row)
    return lol

def check_matrix(matrix_file, genes_only=False, min_count=0, max_count=20000):
    matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)
    if genes_only:
        igr = [n for n in list(matrix.index) if n.startswith("Cluster_")]
        matrix = matrix.drop(igr, axis=0)
    num_genes = matrix.sum(axis=0)
    counts = []
    bad_genomes = []
    for name,count in num_genes.iteritems():
        assert(count >= 0)
        if count < min_count or count > max_count:
            print(name, count)
            bad_genomes.append(name)
        counts.append(count)
    return counts, bad_genomes

def get_bin_boundaries(matrix_file, num_bins = 10):
    matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)
    num_samples = len(matrix.columns)
    assert num_bins <= num_samples, "Number of bins %d > number of samples %d" %(num_bins, num_samples)
    bin_size = float(num_samples)/num_bins
    boundaries = [0]
    last_i = 0
    for count in range(num_samples + 2):
        bin_id = math.ceil((count)/bin_size) - 1
        if bin_id > last_i:
            boundaries.append(count)
            last_i = bin_id
    return boundaries
    
def divide_local_graphs_by_frequency(matrix_file, num_bins = 10, genes_only=False, bad_genomes=[]):
    matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)
    if genes_only:
        igr = [n for n in list(matrix.index) if n.startswith("Cluster_")]
        matrix = matrix.drop(igr, axis=0)
    matrix = matrix.drop(bad_genomes, axis=1)
    print("Dividing local graphs by frequency, excluding", str(len(bad_genomes)), "bad genomes")
    num_samples = len(matrix.columns)
    print("Working with", str(num_samples), "samples")
    if num_bins > num_samples:
        num_bins = num_samples
    bin_size = float(num_samples)/num_bins

    partition = []
    for i in range(num_bins):
        partition.append([])

    frequencies = matrix.sum(axis=1)
    for name,count in frequencies.iteritems():
        assert(count >= 0)
        bin_id = math.ceil(count/bin_size) - 1
        assert bin_id < num_bins
        partition[bin_id].append(name)

    return partition

def plot_hist(partition, xlabel="Allele frequency", nbins=20, kde=False, b=[], outfile="allele_frequency_hist.png"):
    plt.rcParams['figure.figsize'] = 10,6
    plt.style.use('seaborn-deep')
    #sns.set_palette(sns.color_palette("Spectral", len(partition)))
    
    fig, ax = plt.subplots()
    for i,p in enumerate(partition):
        l = str(i)
        if len(b) == len(partition) + 1:
            l = "%d <= #G < %d" %(b[i], b[i+1])
        if kde:
            sns.distplot(p,label=l, bins=nbins, kde=False, hist=True)
        else:
            sns.distplot(p,label=l, bins=nbins, kde=True, hist=False)
        #ax.hist(p,
        #  density=True,
        #  label=str(i),
        #  align="mid",
        #  bins=nbins)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel=xlabel, ylabel='Frequency')
    plt.savefig(outfile, transparent=True, bbox_inches='tight')

def plot_hist2(partition, xlabel="Allele frequency", nbins=20, b=[], outfile="allele_frequency_hist2.png"):
    plt.rcParams['figure.figsize'] = 10,6
    plt.style.use('seaborn-deep')
    sns.set_palette(sns.color_palette("Spectral", len(partition)))
    
    for i,p in enumerate(partition):
        fig, ax = plt.subplots()

        l = str(i)
        if len(b) == len(partition) + 1:
            l = "%d <= #G < %d" %(b[i], b[i+1])
        ax.hist(p,
          density=True,
          label=str(i),
          align="mid",
          bins=nbins)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel=xlabel, ylabel='Frequency')
    plt.savefig(outfile, transparent=True, bbox_inches='tight')
    
def plot_bar(partition, outfile="."):
    plt.rcParams['figure.figsize'] = 10,6
    plt.style.use('seaborn-deep')
    
    fig, ax = plt.subplots()
    ax.bar([i+1 for i in range(len(partition))],partition, log=True)
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel="Gene count", ylabel='Frequency')
    plt.savefig(outfile, transparent=True, bbox_inches='tight')
    
def save_obj(obj, filepath):
    with open(filepath, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)

# using https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
def exp_gk(theta,rho,gc,n,k):
    #print("exp_gk:",theta,rho,gc,n,k)
    assert k > 0
    r = float(theta/k)
    for i in range(k):
        r *= (n-i)
        assert (n-i-1+rho) != 0, "exp_gk:(%f,%f,%f,%d,%d): %d-%d-1+%f=%f" %(theta,rho,gc,n,k,n,i,rho,n-i-1+rho)
        r /= (n-i-1+rho)    
    if n == k:
        r += gc
    return r

def fun(x,t,N,y):
    #print("fun:",x,t,n,y)
    assert(len(N) == len(t))
    assert(len(t) == len(y))
    r = []
    for i in range(len(y)):
        try:
            res = y[i] - exp_gk(x[0],x[1],x[2],N[i],t[i])
        except:
            print(i,x[0],x[1],x[2],N[i],t[i],y[i])
        r.append(res)
    assert(len(r) == len(y))
    return r

parser = argparse.ArgumentParser(description='Fits theta, rho and size of core genome to pandora matrix and vcf')
parser.add_argument('--matrix', type=str,
                    help='Input pandora matrix')
parser.add_argument('--vcf', type=str,
                    help='Input pandora multisample VCF')
parser.add_argument('--outdir', type=str,
                    help='Out directory')
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

matrix_file = args.matrix
vcf_file = args.vcf
outdir=args.outdir

if not os.path.exists(outdir):
    os.mkdir(outdir)

if os.path.exists("%s/gene_counts_per_sample.csv" %args.outdir):
    gene_counts_per_sample = load_lols("%s/gene_counts_per_sample.csv" %args.outdir)[0]
    bad_genomes = load_lols("%s/gene_counts_per_sample.csv" %args.outdir)[1]
else:
    gene_counts_per_sample, bad_genomes = check_matrix(matrix_file, genes_only=True, min_count=3500, max_count=6000)
    plot_hist([gene_counts_per_sample], xlabel="Number genes found", nbins=30, kde=True, outfile="%s/gene_counts_per_sample.png" %outdir)
    save_lols([gene_counts_per_sample, bad_genomes], "%s/gene_counts_per_sample.csv" %args.outdir)

if os.path.exists("%s/gene_partition.csv" %args.outdir):
    gene_partition = load_lols("%s/gene_partition.csv" %args.outdir)
else:
    gene_partition = divide_local_graphs_by_frequency(matrix_file, 285, genes_only=True, bad_genomes=bad_genomes)
    plot_bar([len(p) for p in gene_partition], outfile="%s/gene_frequencies.png" %outdir)
    save_lols(gene_partition, "%s/gene_partition.csv" %args.outdir)

n = len(gene_partition)
N = [n for i in range(1,n+1)]
t = [i for i in range(1,n+1)]
y = [len(p) for p in gene_partition]
x0 = [1000.0,1.0,100.0]
bounds = [0,np.inf]

if os.path.exists("%s/res_lsq.pkl" %outdir):
    res_lsq = load_obj("%s/res_lsq.pkl" %outdir)
else:
    res_lsq = least_squares(fun, x0, bounds=bounds, args=(t,N,y))
    print("Linea loss least squares (theta1,rho,gc):",res_lsq.x)
    save_obj(res_lsq, "%s/res_lsq.pkl" %outdir)

if os.path.exists("%s/res_l1.pkl" %outdir):
    res_l1 = load_obj("%s/res_l1.pkl" %outdir)
else:
    res_l1 = least_squares(fun, x0, bounds=bounds, loss='soft_l1', args=(t,N,y))
    print("Soft L1 loss least squares (theta1,rho,gc):",res_l1.x)
    save_obj(res_l1, "%s/res_l1.pkl" %outdir)

if os.path.exists("%s/res_hub.pkl" %outdir):
    res_hub = load_obj("%s/res_hub.pkl" %outdir)
else:
    res_hub = least_squares(fun, x0, bounds=bounds, loss='huber', args=(t,N,y))
    print("Huber loss least squares (theta1,rho,gc):",res_hub.x)
    save_obj(res_hub, "%s/res_hub.pkl" %outdir)

if os.path.exists("%s/res_log.pkl" %outdir):
    res_log = load_obj("%s/res_log.pkl" %outdir)
else:
    res_log = least_squares(fun, x0, bounds=bounds, loss='cauchy', args=(t,N,y))
    print("Cauchy (log) loss least squares (theta1,rho,gc):",res_log.x)
    save_obj(res_log, "%s/res_log.pkl" %outdir)

if os.path.exists("%s/res_at.pkl" %outdir):
    res_at = load_obj("%s/res_at.pkl" %outdir)
else:
    res_at = least_squares(fun, x0, bounds=bounds, loss='arctan', args=(t,N,y))
    print("Arctan loss least squares (theta1,rho,gc):",res_at.x)
    save_obj(res_at, "%s/res_at.pkl" %outdir)

fig, ax = plt.subplots()
gene_frequencies = [len(p) for p in gene_partition]
ax.bar([i+1 for i in range(len(gene_frequencies))],gene_frequencies, log=True)
ax.set(xlabel="Gene count", ylabel='Frequency')
y_lsq = [exp_gk(*res_lsq.x, n, i) for i in t]
y_log = [exp_gk(*res_log.x, n, i) for i in t]
y_l1 = [exp_gk(*res_l1.x, n, i) for i in t]
y_hub = [exp_gk(*res_hub.x, n, i) for i in t]
y_at = [exp_gk(*res_at.x, n, i) for i in t]
plt.plot(t, y_lsq, label='linear loss')
plt.plot(t, y_log, label='cauchy loss')
plt.plot(t, y_l1, label='soft_l1 loss')
plt.plot(t, y_hub, label='huber loss')
plt.plot(t, y_at, label='arctan loss')
plt.xlabel("t")
plt.ylabel("y")
plt.ylim((0.5,10000))
plt.legend()
plt.savefig("%s/gene_frequency_spectrum_with_fits.png" %outdir, transparent=True, bbox_inches='tight')
