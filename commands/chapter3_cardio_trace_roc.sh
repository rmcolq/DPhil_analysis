# 1. Subsample read files down to 100X
# ILLUMINA
H131800734      /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq
H151080744      /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq
CFT073  /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq
RHB11-C04       /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351_filtered.fastq

# NANOPORE
H131800734      /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H131800734_BC03_100x.fastq
H151080744      /nfs/leia/research/iqbal/rmcolq/data/cardio/nanopore/H151080744_BC11_100x.fastq
CFT073  /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299443_100x.fastq
RHB11-C04       /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299467_100x.fastq

python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq -n 3400000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.100x.fq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq -n 3400000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.100x.fq
# NB not enough coverage in /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351_filtered.fastq or /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq to need it
# and nanopore are random so done with head

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

#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/RHB11-C04
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/RHB11-C04
bsub.py 3.0 genotype_RHB11-C04 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/snippy_nanopolish_genotype.nf --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/RHB11-C04/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351_filtered.fastq -resume 

## Indexed PanRG for pandora
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 19 --k 31 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 14 --k 15 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14

## Ran pandora
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/
bsub.py 60 logs/compare_cardio_trace_roc_illumina_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/illumina_30 --genotype --max_covg 30 --illumina -w 19 -k 31 --min_cluster_size 5
bsub.py 100 logs/compare_cardio_trace_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/illumina_100 --genotype --illumina -w 19 -k 31 --min_cluster_size 5
bsub.py 60 logs/compare_cardio_trace_roc_nanopore_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/nanopore_30 --genotype --max_covg 30 --min_cluster_size 5
bsub.py 100 logs/compare_cardio_trace_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/nanopore_100 --genotype --min_cluster_size 5

bsub.py 60 logs/compare_cardio_trace_roc_illumina_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/illumina_30 --genotype --max_covg 30 --illumina -w 19 -k 31
bsub.py 100 logs/compare_cardio_trace_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/illumina_100 --genotype --illumina -w 19 -k 31
bsub.py 60 logs/compare_cardio_trace_roc_nanopore_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/nanopore_30 --genotype --max_covg 30
bsub.py 100 logs/compare_cardio_trace_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/160519/nanopore_100 --genotype

# 2 cardio only
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/
bsub.py 30 logs/compare_cardio_only_roc_illumina_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/illumina_30 --genotype --max_covg 30 --illumina -w 19 -k 31
bsub.py 50 logs/compare_cardio_only_roc_illumina_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/illumina_100 --genotype --illumina -w 19 -k 31
bsub.py 30 logs/compare_cardio_only_roc_nanopore_30 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/nanopore_30 --genotype --max_covg 30
bsub.py 50 logs/compare_cardio_only_roc_nanopore_100 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/160519/nanopore_100 --genotype

#--genotyping_error_rate 1 if doesnt work?

## Create single sample vcf and ref.fa files for each sample
for run in illumina_30 illumina_100 nanopore_30 nanopore_100
do
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,12 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample_genotyped.vcf | cut -f1-9,13 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/RHB11-C04/pandora_genotyped_$run.vcf
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/CFT073/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H131800734/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/H151080744/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/090519/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/genotyping/RHB11-C04/pandora_genotyped_$run.ref.fa
done

## make bed of bad regions
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/truths
cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz .
gunzip H151080744_pacbio_assembly.pilon.fa.gz
cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz .
gunzip H131800734_pacbio_assembly.pilon.fa.gz

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/masks
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/masks
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh ../truths/H151080744_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq H151080744
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh ../truths/H131800734_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq H131800734
bsub.py 4 build_mask_RH bash make_low_qual_genome_mask.sh ../truths/RHB11-C04_subsampled_pacbio.fasta /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351_filtered.fastq RHB11-C04
bsub.py 4 build_mask_CF bash make_low_qual_genome_mask.sh ../truths/CFT073_subsampled_pacbio.fasta /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq CFT073

## Make ROC
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/genotyped_samples.tsv

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way/140519
bsub.py 4.0 2_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way/140519
bsub.py 4.0 3_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way/140519
bsub.py 6.0 4_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/labels.tsv
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt

# and another way around
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/genotyped_samples.tsv
vim /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/genotyped_samples.tsv

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/140519
bsub.py 4.0 2_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/140519
bsub.py 4.0 3_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/140519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/140519
bsub.py 6.0 4_way singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv ../genotyped_samples.tsv -t 8

cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/4_way_other/140519
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
#cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/3_way_other/140519
#python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace_roc/2_way_other/140519
python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/plot_roc.py --y_label "Pairwise SNP Recall" --y_min 0.1 --x_max 0.04 --legend_outside --labels ../../labels.tsv &> log_plot_roc.txt
