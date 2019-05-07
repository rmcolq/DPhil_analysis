mkdir -p /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published 
cd /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published

echo "#!/usr/bin/env bash
set -vex

for accession in $(cat published_accessions)
do
    enaDataGet -f fastq $accession
done" > get_reads.sh

#downloaded accessions from https://www.ncbi.nlm.nih.gov/bioproject/379782 to published_accessions
bsub.py 4 get_reads bash get_reads.sh

echo "for id in $(cat published_accessions2); do echo $id; ls $id; zcat $id/"$id"_1.fastq.gz $id/"$id"_2.fastq.gz > $id/"$id"_both.fastq; gzip $id/"$id"_both.fastq; echo -e "$id\t$(ls $PWD/$id/*both*)" >> INDEX; done" > combine_reads.sh
bsub.py 4 combine_reads bash combine_reads.sh

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/070519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/070519
bsub.py 4 pandora_compare_in_chunks_070519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -l 31 --read_index /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX --chunk_size 2000 --num_samples 266 --max_covg 50 --illumina -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config

