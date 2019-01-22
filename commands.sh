#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/random_path_consensus_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/random_paths_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_consensus_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --reference_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --fasta_for_ref /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.all.fa --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/all
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/all
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/single_sample_roc_all nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/single_sample_roc.nf --truth_assembly /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H152180645/H152180645_pacbio_assembly.pilon.fa --vcf_directory /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/data/ --mask /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H152180645/H152180645_mask.bed --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/all -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/all/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/snps
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/snps
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/single_sample_roc_snps nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/single_sample_roc.nf --truth_assembly /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H152180645/H152180645_pacbio_assembly.pilon.fa --vcf_directory /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/data/ --mask /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H152180645/H152180645_mask.bed --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/snps/ --snps -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/single_sample_roc/snps/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/multisample_roc
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/multisample_roc
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/multisample_roc.nf --truth_assembly1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz --vcf_directory1 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping_with_compare/subset/ --truth_assembly2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz --mask1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_mask.bed --vcf_directory2 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping_with_compare/ --mask2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_mask.bed -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/multisample_roc/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config

