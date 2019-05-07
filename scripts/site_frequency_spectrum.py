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

def check_matrix(matrix_file):
    matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)
    num_genes = matrix.sum(axis=0)
    counts = []
    for name,count in num_genes.iteritems():
        assert(count > 0)
        if count < 1000:
            print(name, count)
        counts.append(count)    
    return counts

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
    
def divide_local_graphs_by_frequency(matrix_file, num_bins = 10, genes_only=False):
    matrix = pd.read_csv(matrix_file, sep='\t', header=0, index_col=0)
    num_samples = len(matrix.columns)
    assert num_bins <= num_samples, "Number of bins %d > number of samples %d" %(num_bins, num_samples)
    bin_size = float(num_samples)/num_bins
    
    partition = []
    for i in range(num_bins):
        partition.append([])
        
    frequencies = matrix.sum(axis=1)
    for name,count in frequencies.iteritems():
        if genes_only and name.startswith("Cluster_"):
            continue
        assert(count > 0)
        bin_id = math.ceil(count/bin_size) - 1
        assert bin_id < num_bins
        partition[bin_id].append(name)
    
    return partition

def get_counts_per_gene(vcf_file, max_allele_length=1):
    count_dict = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        if record.CHROM not in count_dict.keys():
            count_dict[record.CHROM] = []
        alleles = [record.REF]
        alleles.extend(record.ALT)
        alleles = [str(a) for a in alleles]
        assert(len(alleles) == 1+len(record.ALT))

        max_len = max([len(a) for a in alleles])
        if max_len > max_allele_length:
            continue
            
        counts = [0 for a in alleles]
        assert(len(counts) == len(alleles))

        # collect counts on each allele from the samples
        for i,sample_record in enumerate(record.samples):
            genotype = sample_record['GT']
            genotypes = genotype.split('/')
            called_alleles = set(genotypes)
            if len(called_alleles) != 1 or '.' in called_alleles:
                continue
            else:
                gt = int(list(called_alleles)[0])
                assert(gt < len(counts))
                counts[gt] += 1
        
        # see if it's a segregating site, and if it is, flip a coin to decide the ancestral allele
        indexes = np.nonzero(counts)
        if len(indexes[0]) > 1:
            j = np.random.randint(0,len(indexes[0]))
            index = indexes[0][j] 
            count_dict[record.CHROM].append(counts[index])
    return count_dict

def plot_hist(partition, xlabel="Allele frequency", nbins=20, kde=False, b=[], outfile="allele_frequency_hist.png"):
    plt.rcParams['figure.figsize'] = 10,6
    plt.style.use('seaborn-deep')
    sns.set_palette(sns.color_palette("Spectral", len(partition)))
    
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
    plt.savefig(outfile, transparent=True)

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
    plt.savefig(outfile, transparent=True)
    
def plot_bar(partition, outfile="."):
    plt.rcParams['figure.figsize'] = 10,6
    plt.style.use('seaborn-deep')
    
    fig, ax = plt.subplots()
    ax.bar([i+1 for i in range(len(partition))],partition, log=True)
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles, labels)
    #plt.legend(labels=types)
    ax.set(xlabel="Gene count", ylabel='Frequency')
    plt.savefig(outfile, transparent=True)
    
def save_obj(obj, filepath):
    with open(filepath, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)

# using https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
def nCr(n,r):
    #only handles integers
    #print("nCr of ", n, " and ", r)
    if r == 0:
        return 1
    if r == 1:
        return n
    f = math.factorial
    ret = f(n) // f(r) // f(n-r)
    return ret

def exp_sfs(rho,theta2,n,k,s):
    assert s <= k
    r = theta2/s
    r *= k
    if s < k:
        r /= n
        r /= nCr(n-1,s)
        t = 0
        for j in range(n-s):
            t += ((j+1)/(j+1+rho)) * nCr(n-j-2,s-1)
        r *= t
    elif s == k:
        r *= s
        t = 0
        for j in range(1,n-s+2):
            u = 1/(j*(j-1+rho))
            u *= nCr(n-k,j-1)/nCr(n,j)
            t += u
        r *= t
    return r

def fun3(x,N,K,S,y):
    assert(len(N) == len(K))
    assert(len(K) == len(S))
    assert(len(S) == len(y))
    r = []
    for i in range(len(y)):
        res = 0
        try:
            res = y[i] - exp_sfs(x[0],x[1],N[i],K[i],S[i])
        except:
            print(i,x[0],x[1],N[i],K[i],S[i],y[i])
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

matrix_file = args.matrix
vcf_file = args.vcf
outdir=args.outdir

gene_counts_per_sample = check_matrix(matrix_file)
plot_hist([gene_counts_per_sample], xlabel="Number genes found", nbins=30, kde=True, outfile="%s/gene_counts_per_sample.png" %outdir)
save_lols([gene_counts_per_sample], "%s/gene_counts_per_sample.csv" %args.outdir)

gene_partition = divide_local_graphs_by_frequency(matrix_file, 285, genes_only=True)
plot_bar([len(p) for p in gene_partition], outfile="%s/gene_frequencies.png" %outdir)
save_lols(gene_partition, "%s/gene_partition.csv" %args.outdir)

count_dict = get_counts_per_gene(vcf_file, max_allele_length=1)
save_obj(count_dict, "%s/count_dict.pkl" %outdir)

N = []
K = []
S = []
y = []

n = len(gene_partition)
for i,p in enumerate(gene_partition):
    k = i+1
    if i == 1 or i == n:
        continue
    
    for gene in p:
        if gene not in count_dict.keys():
            continue
        for s in range(1,k):
            j = sum([1 for c in count_dict[gene] if int(c) == s])
            N.append(n)
            K.append(k)
            S.append(s)
            y.append(j)

x0 = [3.0,0.001]

res_lsq = least_squares(fun3, x0, args=(N,K,S,y))
res_l1 = least_squares(fun3, x0, loss='soft_l1', args=(N,K,S,y))
res_hub = least_squares(fun3, x0, loss='huber', args=(N,K,S,y))
res_log = least_squares(fun3, x0, loss='cauchy', args=(N,K,S,y))
res_at = least_squares(fun3, x0, loss='arctan', args=(N,K,S,y))

for i in [50,100,150,200,250]:
    fig, ax = plt.subplots()
    genes = gene_partition[i-1]
    spectrums = []
    for s in range(1,k):
        spectrums.append([sum([1 for c in count_dict[gene] if int(c) == s]) for gene in genes])
    ax.violinplot(spectrums,showmeans=True)
    ax.set(xlabel="Frequency of SNP", ylabel='Number of SNP sites at this frequency', xticklabels=[s for s in range(1,k)])
    
    y_lsq = [exp_sfs(*res_lsq.x, n, i, t) for t in range(1,i)]
    y_log = [exp_sfs(*res_log.x, n, i, t) for t in range(1,i)]
    y_l1 = [exp_sfs(*res_l1.x, n, i, t) for t in range(1,i)]
    y_hub = [exp_sfs(*res_hub.x, n, i, t) for t in range(1,i)]
    y_at = [exp_sfs(*res_at.x, n, i, t) for t in range(1,i)]
    
    plt.plot([t for t in range(1,i)], y_lsq, label='linear loss')
    plt.plot([t for t in range(1,i)], y_log, label='cauchy loss')
    plt.plot([t for t in range(1,i)], y_l1, label='soft_l1 loss')
    plt.plot([t for t in range(1,i)], y_hub, label='huber loss')
    plt.plot([t for t in range(1,i)], y_at, label='arctan loss')
    
    plt.ylim((0.5,10000))
    plt.legend()
    plt.savefig("%s/site_frequency_spectrum_with_fits_%s.png" %(i,outdir), transparent=True)

print("Linear loss least squares (rho,theta2):", res_lsq.x)
print("Cauchy (log) loss least squares (rho,theta2):",res_log.x)
print("Soft L1 loss least squares (rho,theta2):",res_l1.x)
print("Huber loss least squares (rho,theta2):",res_hub.x)
print("Arctan loss least squares (rho,theta2):",res_at.x)
