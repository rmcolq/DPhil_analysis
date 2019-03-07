#!/usr/bin/env python 

import argparse
import vcf
import numpy as np

def filter_vcf_by_gt_conf(in_vcf, discard_percent):
    percentile = 0
    with open(in_vcf, 'r') as f:
        vcf_reader = vcf.Reader(f)
        gt_confs = []
        covgs = []
        for record in vcf_reader:
            rec = record.samples[0]
            gt = rec['GT'][0]
            if gt not in ['0','1']:
                continue  
            gt = int(gt)
            gt_conf = rec['GT_CONF']
            gt_confs.append(float(gt_conf))
            covg = int(rec['MEAN_FWD_COVG'][gt]) + int(rec['MEAN_REV_COVG'][gt])
            covgs.append(covg)
        percentile_conf = np.percentile(gt_confs, discard_percent)
        percentile_covg = np.percentile(covgs, discard_percent)
        print(percentile_conf, percentile_covg)
        
    vcf_reader = vcf.Reader(open(in_vcf, 'r'))
    out_vcf_good = vcf.Writer(open(in_vcf.replace(".vcf",".%fpercentile.vcf" %discard_percent),'w'), vcf_reader)
    for record in vcf_reader:
        rec = record.samples[0]
        gt = rec['GT'][0]
        if gt not in ['0','1']:
            continue  
        gt = int(gt)
        gt_conf = float(rec['GT_CONF'])
        covg = int(rec['MEAN_FWD_COVG'][gt]) + int(rec['MEAN_REV_COVG'][gt])
        if gt_conf > percentile_conf and covg > percentile_covg:
            out_vcf_good.write_record(record)

parser = argparse.ArgumentParser(description='Filters out the bottom x% of gt_confidences and covgs.')
parser.add_argument('--vcf', type=str,
                    help='VCF of calls to filter')
parser.add_argument('--percentile', type=float,
                    help='Percentile to filter')
args = parser.parse_args()

filter_vcf_by_gt_conf(args.vcf, args.percentile)
