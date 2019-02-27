#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/random_path_consensus_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/random_paths_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/random_path_consensus_evaluation/work4 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config #-resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_consensus_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --reference_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_consensus_evaluation/work4 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config #-resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_evaluation nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work2 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work3 #repeat with a higher number of differences
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work3
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_evaluation3 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work3 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work_pandora_filtered #repeat with a higher number of differences
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work_pandora_filtered
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_evaluation3 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping_evaluation/work_pandora_filtered -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
