import argparse
import gzip
import glob

def write_fasta(in_file):
    name = in_file.replace(".gff","").split("/")[-1]
    out_file = in_file.replace(".gff",".fa.gz")
    print("%s\t%s" %(name, out_file))
    fasta = False
    with open(in_file,'r') as in_handle:
        with gzip.open(out_file, 'wt') as out_handle:
            for line in in_handle:
                if fasta:
                    out_handle.write(line)
                if "##FASTA" in line:
                    fasta = True

parser = argparse.ArgumentParser(description='Extracts assembly fasta from gff file')
parser.add_argument('--dir', type=str, default=".",
                    help='Directory of gff files')
args = parser.parse_args()

in_files = glob.glob("%s/*.gff" %args.dir)
for f in in_files:
    write_fasta(f)
