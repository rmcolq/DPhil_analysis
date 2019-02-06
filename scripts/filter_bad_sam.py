import pysam
import argparse

def filter_sam(sam_file):
    sam_reader = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    sam_writer = pysam.AlignmentFile("filtered.bam", "wb", template=sam_reader)
    for s in sam_reader.fetch(until_eof=True):
        if s.cigartuples != None:
            sam_writer.write(s)

parser = argparse.ArgumentParser(description='Filters out any lines in SAM which have no cigar, and writes as BAM')
parser.add_argument('--sam', type=str,
                    help='Input SAM')
args = parser.parse_args()

filter_sam(args.sam)
