import pysam
import argparse

def get_ref_from_sam(sam_file):
    f = pysam.AlignmentFile(sam_file, "r")
    last_name = None
    for s in f:
        if s.query_name == last_name:
            continue
        print(">" + s.query_name)
        if s.is_unmapped:
            print(s.query_sequence)
        else:
            print(s.get_reference_sequence())
        last_name = s.query_name

parser = argparse.ArgumentParser(description='Uses the equivalent of ref sequence where maps to a query sequence.')
parser.add_argument('--sam', type=str,
                    help='Input SAM')
args = parser.parse_args()

get_ref_from_sam(args.sam)

