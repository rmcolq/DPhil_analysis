# 1. Subsample read files down to 100X
# ILLUMINA
H131800734      /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq
H151080744      /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq
CFT073  /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq
K12_MG1655	/nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz

# NANOPORE
H131800734      /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H131800734_BC03_100x.fastq
H151080744      /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H151080744_BC11_100x.fastq
CFT073  /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299443_100x.fastq
K12_MG1655       /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/loman_k12_pass_porechop.100x.fq

python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq -n 3400000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq -n 3400000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq
zcat /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz > K12_MG1655.fq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py K12_MG1655.fq  /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.100x.fastq -n 1700000
mv /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.100x.fastq.0 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.100x.fastq

python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.30x.fq -n 1000000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.30x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.30x.fq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.30x.fq -n 1000000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.30x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.30x.fq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.30x.fastq -n 1000000
mv /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.30x.fastq.0 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.30x.fastq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py K12_MG1655.fq  /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.30x.fastq -n 570000
mv /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.30x.fastq.0 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.30x.fastq

rm K12_MG1655.fq
# NB nanopore are random so done with head

vim /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv
vim /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv

# 2. Genotype with subsampled files and SNIPPY AND NANOPOLISH
#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734
bsub.py 3.0 genotype_H131800734 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/snippy_nanopolish_genotype.nf --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H131800734_BC03_100x.fastq --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq -resume

#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744
bsub.py 3.0 genotype_H151080744 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/snippy_nanopolish_genotype.nf --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H151080744_BC11_100x.fastq --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq -resume

#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073
bsub.py 3.0 genotype_CFT073 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/snippy_nanopolish_genotype.nf --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq -resume

#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/K12_MG1655
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/K12_MG1655
bsub.py 3.0 genotype_K12_MG1655 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/snippy_nanopolish_genotype.nf --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/K12_MG1655/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.100x.fastq --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --nanopore_reads /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/loman_k12_pass_porechop.100x.fq -resume 

## Indexed PanRG for pandora
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 19 --k 31 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 14 --k 15 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14

## Ran pandora
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/
bsub.py 60 logs/compare_cardio_trace_roc_illumina_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.30x.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/illumina_30 --genotype --max_covg 30 --illumina -w 19 -k 31
bsub.py 100 logs/compare_cardio_trace_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/illumina_100 --genotype --illumina -w 19 -k 31
bsub.py 60 logs/compare_cardio_trace_roc_nanopore_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/nanopore_30 --genotype --max_covg 30
bsub.py 100 logs/compare_cardio_trace_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/nanopore_100 --genotype

# 2 cardio only
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/
bsub.py 30 logs/compare_cardio_only_roc_illumina_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.30x.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/illumina_30 --genotype --max_covg 30 --illumina -w 19 -k 31
bsub.py 50 logs/compare_cardio_only_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/illumina_100 --genotype --illumina -w 19 -k 31
bsub.py 30 logs/compare_cardio_only_roc_nanopore_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/nanopore_30 --genotype --max_covg 30
bsub.py 50 logs/compare_cardio_only_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/nanopore_100 --genotype

bsub.py 50 logs/compare_cardio_only_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/illumina_100_gt1 --genotype --illumina -w 19 -k 31 --genotyping_error_rate 1
bsub.py 50 logs/compare_cardio_only_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/nanopore_100_gt1 --genotype --genotyping_error_rate 1
#--genotyping_error_rate 1 if doesnt work?

## Create single sample vcf and ref.fa files for each sample
for run in illumina_30 illumina_100 nanopore_30 nanopore_100
do
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,12 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,13 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/K12_MG1655/pandora_genotyped_$run.vcf
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/K12_MG1655/pandora_genotyped_$run.ref.fa
done

## and for 2 cardio only
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H131800734
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H151080744
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/
rsync -r cardio_trace_roc/genotyping/H131800734/snippy*.ref.fa cardio_roc/genotyping/H131800734/
rsync -r cardio_trace_roc/genotyping/H131800734/snippy*.vcf cardio_roc/genotyping/H131800734/
rsync -r cardio_trace_roc/genotyping/H151080744/snippy*.vcf cardio_roc/genotyping/H151080744/
rsync -r cardio_trace_roc/genotyping/H151080744/snippy*.ref.fa cardio_roc/genotyping/H151080744/
rsync -r cardio_trace_roc/genotyping/H131800734/nanopolish*.ref.fa cardio_roc/genotyping/H131800734/
rsync -r cardio_trace_roc/genotyping/H131800734/nanopolish*.vcf cardio_roc/genotyping/H131800734/
rsync -r cardio_trace_roc/genotyping/H151080744/nanopolish*.vcf cardio_roc/genotyping/H151080744/
rsync -r cardio_trace_roc/genotyping/H151080744/nanopolish*.ref.fa cardio_roc/genotyping/H151080744/
rm cardio_roc/genotyping/*/pandora*

for run in illumina_30 illumina_100 nanopore_30 nanopore_100
do
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H131800734/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H151080744/pandora_genotyped_$run.vcf
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H131800734/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H151080744/pandora_genotyped_$run.ref.fa
done

#for run in illumina_100_gt1 nanopore_100_gt1
#do
#less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H131800734/pandora_genotyped_$run.vcf
#less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H151080744/pandora_genotyped_$run.vcf
#cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H131800734/pandora_genotyped_$run.ref.fa
#cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyping/H151080744/pandora_genotyped_$run.ref.fa
#done

## make bed of bad regions
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/truths
cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz .
gunzip H151080744_pacbio_assembly.pilon.fa.gz
cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz .
gunzip H131800734_pacbio_assembly.pilon.fa.gz
cp /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta K12_MG1655.assembly.fa

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/masks
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/masks
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh ../truths/H151080744_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq H151080744
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh ../truths/H131800734_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq H131800734
bsub.py 4 build_mask_CF bash make_low_qual_genome_mask.sh ../truths/CFT073_subsampled_pacbio.fasta /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq CFT073
bsub.py 4 build_mask_K12 bash make_low_qual_genome_mask.sh /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz K12_MG1655

## Make DNADIFFs for summary
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/dnadiffs
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/dnadiffs
bsub_big
singularity shell /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img
dnadiff ../truths/H131800734*.fa ../truths/H151080744*.fa -p H131800734_H151080744
dnadiff ../truths/H131800734*.fa ../truths/CFT073*.f*a -p H131800734_CFT073
dnadiff ../truths/H131800734*.fa ../truths/K12_MG1655*.fa -p H131800734_K12_MG1655
dnadiff ../truths/H151080744*.fa ../truths/CFT073*.f*a -p H151080744_CFT073
dnadiff ../truths/H151080744*.fa ../truths/K12_MG1655*.fa -p H151080744_K12_MG1655
dnadiff ../truths/CFT073*.f*a ../truths/K12_MG1655*.fa -p CFT073_K12_MG1655
exit
vim dnadiff.id.tsv
vim dnadiff.alignedbases.tsv
exit
# plotted in python notebook with
#df = pd.read_csv("file.tsv", sep='\t', index_col=0,header=0)
#sns.heatmap(df, cmap='RdYlGn_r', linewidths=0.5, annot=True, fmt='.2f', square=True, annot_kws={"size":12})
#sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5})

## Make ROC
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/genotyped_samples.tsv

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/170519
bsub.py 4.0 2_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/170519
bsub.py 4.0 3_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/170519
bsub.py 6.0 4_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/labels.tsv
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt

## and for 2 cardio
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/genotyped_samples.tsv
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/2_way/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/2_way/170519
bsub.py 4.0 2_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../../genotyped_samples.tsv -t 8

vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/labels.tsv
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_roc/2_way/170519
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt


#######
# and another way around
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/genotyped_samples.tsv

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/170519
bsub.py 4.0 2_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/170519
bsub.py 4.0 3_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/170519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/170519
bsub.py 6.0 4_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/170519
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
#cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/170519
#python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/170519
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
