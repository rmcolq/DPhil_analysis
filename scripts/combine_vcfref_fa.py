from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def combine_fa(full_fa, subset_fa, out_fa):
    full_dict = SeqIO.to_dict(SeqIO.parse(full_fa, "fasta"))
    subset_dict = SeqIO.to_dict(SeqIO.parse(subset_fa, "fasta"))
    out_records = []

    for record in full_dict:
        if record in subset_dict:
            out_records.append(subset_dict[record])
        else:
            out_records.append(full_dict[record])

    SeqIO.write(out_records, out_fa, "fasta")


parser = argparse.ArgumentParser(description='Replaces sequences in full with subset where they exist.')
parser.add_argument('--full_fa', type=str,
                    help='Input FASTA')
parser.add_argument('--subset_fa', type=str,
                    help='Input FASTA')
parser.add_argument('--out_fa', type=str,
                    help='Output FASTA')
args = parser.parse_args()

combine_fa(args.full_fa, args.subset_fa, args.out_fa)
