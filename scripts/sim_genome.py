from Bio import SeqIO
import argparse
import random
import numpy as np

def rev_comp(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    comp = ''.join(letters)
    return comp[::-1]

def simulate_genome(seq_fasta, num_genes):
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_fasta, "fasta"))

    all_genes = [i for i in record_dict.keys() if i.startswith("GC00")]
    all_intergenic = [i for i in record_dict.keys() if i.startswith("Cluster")]

    genes = random.sample(all_genes, num_genes)
    intergenic = random.sample(all_intergenic, min(num_genes, len(all_intergenic)))
    with open("truth.txt", "w") as f:
        f.write(",".join(genes))
        f.write(",".join(intergenic))

    genome = []
    for i in range(num_genes):
        s = str(record_dict[genes[i]].seq)
        if i < len(intergenic):
            s = s + str(record_dict[intergenic[i]].seq)
        if np.random.binomial(1, 0.01) == 1:
            genome.append(rev_comp(s))
        else:
            genome.append(s)

    with open("sim_genome.fa", "w") as f:
        f.write(">sim_genome\n%s\n" %("".join(genome)))

parser = argparse.ArgumentParser(description='Takes a FASTA of gene sequences and picks randomly with replacement and concatenates to create a simulated reference')
parser.add_argument('--in_fa', type=str,
                    help='Input gene FASTA')
parser.add_argument('--num_genes', type=int,
                    help='Number of genes to concatenate')
args = parser.parse_args()

simulate_genome(args.in_fa, args.num_genes)
