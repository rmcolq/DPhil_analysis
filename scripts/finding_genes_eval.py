import argparse
from Bio import SeqIO

def evaluate_gene_finding(prefix, covg, truth, fastq):
    truth_list = []
    with open(truth, 'r') as f:
        ids = f.readline().strip().split(",")
        truth_list = ids

    record_dict = SeqIO.to_dict(SeqIO.parse(fastq, "fastq"))

    tp = [n for n in truth_list if n in record_dict.keys()]
    fn = [n for n in truth_list if n not in record_dict.keys()]
    fp = [n for n in record_dict.keys() if n not in truth_list]

    p = float(len(tp))/(len(tp)+len(fp))
    r = float(len(tp))/(len(tp)+len(fn))

    with open("results.tsv", 'w') as f:
        f.write("%s\t%s\t%d\t%d\t%d\t%f\t%f\n" %(prefix, covg, tp, fn, fn, p, r))

parser = argparse.ArgumentParser(description='Takes pandora FASTQ and a text file list of genes that are true, and evaluates precision and recall')
parser.add_argument('--fastq', type=str,
                    help='Input gene FASTQ')
parser.add_argument('--truth', type=str,
                    help='Input file of true ids')
parser.add_argument('--prefix', type=str,
                    help='Info about run')
parser.add_argument('--covg', type=str,
                    help='Coverage of run')
args = parser.parse_args()

if ".f" in args.prefix:
    prefix = args.prefix.split(".")[0].split("simulated_")[1]
else:
    prefix = args.prefix

evaluate_gene_finding(prefix, args.covg, args.truth, args.fastq)
