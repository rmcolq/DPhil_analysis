import vcf
import pysam
from Bio import SeqIO
import numpy as np
import copy
import argparse

def rev_comp(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    comp = ''.join(letters)
    return comp[::-1]

def simulate_ref(vcf_file, vcf_ref, sam_file, truth, out_vcf, var_size_thresh, p):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vcf_writer = vcf.Writer(open(out_vcf, 'w'), vcf_reader)
    record_dict = SeqIO.to_dict(SeqIO.parse(vcf_ref, "fasta"))
    truth_dict = SeqIO.to_dict(SeqIO.parse(truth, "fasta"))
    
    last_record = None
    n, p = 1, .1  # number of trials, probability of each trial
    sum_keep = 0
    sum_throw = 0
    for record in vcf_reader:
        if len(record.REF) < var_size_thresh \
        and (last_record == None \
             or record.CHROM != last_record.CHROM \
             or record.POS > last_record.POS + len(last_record.REF)) \
        and np.random.binomial(n, p) == 1:
            chrom_seq = str(record_dict[record.CHROM].seq)
            ref_seq = chrom_seq[record.POS-1:record.POS+len(record.REF)-1]
            assert(ref_seq == record.REF)
            sam_reader = pysam.AlignmentFile(sam_file, "r", check_sq=False)
            for s in sam_reader.fetch(until_eof=True):
                name = str(s).split("\t")[0]
                if name == record.CHROM:
                    if s.is_unmapped:
                        print("UNMAPPED")
                        break
                    record_copy = copy.deepcopy(record)
                    record.CHROM = s.reference_name
                    orientation = s.is_reverse
                    truth_pos = s.reference_start
                    truth_seq = str(truth_dict[record.CHROM].seq)
                    if orientation == False:
                        truth_local_seq = truth_seq[truth_pos+record.POS-1:truth_pos+record.POS+len(record.REF)-1]
                        record.POS = truth_pos+record.POS-1
                        rand_alt = np.random.randint(len(record.ALT))
                        record.ALT = [record.ALT[int(rand_alt)]]
                    else:
                        match_length = s.infer_query_length()
                        truth_local_seq = rev_comp(truth_seq[truth_pos+match_length-record.POS-len(record.REF)+1:truth_pos+match_length-record.POS+1])
                        record.POS = truth_pos+match_length-record.POS-len(record.REF)+1
                        rand_alt = np.random.randint(len(record.ALT))
                        alt = str(record.ALT[int(rand_alt)])
                        record.ALT[0] = rev_comp(alt)
                        record.ALT = record.ALT[:1]
                        record.REF = rev_comp(record.REF)
                    if(truth_local_seq == ref_seq):
                        vcf_writer.write_record(record)
                        last_record = record_copy
                        break
        
parser = argparse.ArgumentParser(description='Takes a VCF and the VCF reference sequences and randomly selects some. Using a SAM file mapping sequences in the VCF reference to a reference assembly, writes this subset of variants with respect to the reference assembly')
parser.add_argument('--in_vcf', type=str,
                    help='Input VCF')
parser.add_argument('--vcf_ref', type=str,
                    help='Reference sequences for VCF')
parser.add_argument('--sam', type=str,
                    help='SAM file mapping VCF reference sequences to reference assembly')
parser.add_argument('--ref', type=str,
                    help='Reference FASTA with respect to output VCF')
parser.add_argument('--out_vcf', type=str,
                    help='VCF for output')
parser.add_argument('--var_size_thresh', type=int, default=10,
                    help='Maximum length of reference variants to potentially include')
parser.add_argument('--prob', type=float, default=.1,
                    help='Fraction of compatible random variants to include in output VCF')
args = parser.parse_args()

simulate_ref(args.in_vcf, args.vcf_ref, args.sam, args.ref, args.out_vcf, args.var_size_thresh, args.prob)
