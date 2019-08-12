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

for i in $(cat /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX | cut -f1)
do
echo $i
zcat /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq.gz > /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.45x.fastq -n 1500000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.45x.fastq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.45x.fastq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.45x.fastq" >> /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX_45
rm /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq
done

for i in $(cat /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX | cut -f1)
do
echo $i
zcat /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq.gz > /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.30x.fastq -n 1000000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.30x.fastq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.30x.fastq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/$i.30x.fastq" >> /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX_30
rm /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/$i/"$i"_both.fastq
done

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/070519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/070519
bsub.py 3 pandora_compare_in_chunks_070519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa --w 19 --k 31 --read_index /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX_45 --chunk_size 2000 --num_samples 266 --max_covg 50 --illumina --min_cluster_size 5 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/070519
bsub.py 3 pandora_compare_in_chunks_070519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa --w 19 --k 31 --read_index /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX_30 --chunk_size 2000 --num_samples 266 --max_covg 50 --illumina --min_cluster_size 5 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/280519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/trace/280519
bsub.py 3 pandora_compare_in_chunks_070519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa --w 19 --k 31 --read_index /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina_TRACE_published/INDEX_30 --chunk_size 2000 --num_samples 266 --max_covg 50 --illumina --min_cluster_size 5 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
