#!/usr/bin/env bash
set -e

if [ $# -ne 3 ]
then
    echo "usage: $0 assembly.fa reads.fq outprefix
Makes a BED file of bad regions, by mapping reads and running bam_to_low_qual_mask.pl.
BWA and samtools faidx indexes the assembly fasta if needed"
    exit
fi

assembly=$1
reads=$2
outprefix=$3
bam=$outprefix.bam
bed=$outprefix.mask.bed

if [ ! -f $assembly.fai ]; then
    samtools faidx $assembly
fi

if [ ! -f $assembly.bwt ]; then
    bwa index $assembly
fi

bwa mem -x intractg $assembly $reads | samtools sort -o $bam
bam_to_low_qual_mask.pl $bam $bed

#belongs to martin
#https://github.com/martinghunt/bioinf-scripts/blob/master/bash/make_low_qual_genome_mask.sh
