from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import gzip

def remove_flanks(in_fq, out_fa, flank_size):
    in_seqs = []

    print(in_fq, in_fq[-3:])
    if in_fq[-3:] == ".gz":
        with gzip.open(in_fq, 'rt') as input_handle:
            in_seqs = list(SeqIO.parse(input_handle, "fastq"))
    else:
        in_seqs = list(SeqIO.parse(in_fq, "fastq"))
    out_records = []

    for entry in in_seqs:
        new_seq = entry.seq[flank_size:len(entry.seq)-flank_size]
        record = SeqRecord(new_seq, entry.id, '', '')
        out_records.append(record)
    
    SeqIO.write(out_records, out_fa, "fasta")


parser = argparse.ArgumentParser(description='Removes begining and end flank_size of each sequence.')
parser.add_argument('--in_fq', type=str,
                    help='Input FASTQ')
parser.add_argument('--out_fa', type=str,
                    help='Output FASTA')
parser.add_argument('--flank_size', type=int, default=28,
                    help='Size of flank to remove')
args = parser.parse_args()

remove_flanks(args.in_fq, args.out_fa, args.flank_size)

