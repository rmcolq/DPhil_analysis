#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/toy_multisample nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/toy_multisample.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/work2 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/mapping_ani/work
#cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/mapping_ani/work
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/mapping_ani nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/mapping_ani.nf --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --reference_dir /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/subset --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/mapping_ani/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/mapping_ani/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc
#cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/multisample_roc.nf --truth_assembly1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz --vcf_directory1 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping_with_compare/subset/ --truth_assembly2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz --mask1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_mask.bed --vcf_directory2 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping_with_compare/ --mask2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_mask.bed -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config

#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/multisample_roc.nf --truth_assembly1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz --vcf_directory1 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping/chosen_H131800734/ --truth_assembly2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz --mask1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_mask.bed --vcf_directory2 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping/chosen_H151080744/ --mask2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_mask.bed -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc/work5 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config
#cd work7
#bsub.py 20.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc/work7/samples.tsv
#mkdir work8
#cd work8
#bsub.py 20.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/multisample_roc/work7/samples.tsv
