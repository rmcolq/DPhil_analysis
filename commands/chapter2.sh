#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_path_consensus_evaluation
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_path_consensus_evaluation
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/perfect_random_path_consensus nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/perfect_random_paths_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_path_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_path_consensus_evaluation/work3 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/noisy_random_path_consensus_evaluation
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/noisy_random_path_consensus_evaluation
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/noisy_random_path_consensus nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/noisy_random_paths_consensus.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/noisy_random_path_consensus_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/noisy_random_path_consensus_evaluation/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_mosaic_evaluation
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_mosaic_evaluation
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_mosaic nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_mosaic_evaluation.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --reference_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_mosaic_evaluation/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_mosaic_evaluation/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_paths_genotyping
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_paths_genotyping
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/perfect_random_paths_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/perfect_random_paths_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/perfect_random_paths_genotyping/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping/work2 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

less /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz | head -n170000 > data/lomank12/loman_k12_pass_porechop.300x.fq
less data/lomank12/loman_k12_pass_porechop.300x.fq | grep -A1 "runid" | grep -v "runid" | wc -c
#1531387765
less /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz | head -n63500 > data/lomank12/loman_k12_pass_porechop.100x.fq
less data/lomank12/loman_k12_pass_porechop.100x.fq | grep -A1 "runid" | grep -v "runid" | wc -c
#504352034
zcat /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz | head -n 6600000 > data/lomank12/Ecoli_S1_L001_R1_001.100x.fq
less data/lomank12/Ecoli_S1_L001_R1_001.100x.fq | grep -A1 "@M0" | grep -v "@M0" | wc -c
#500503821
zcat /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R1_001.fastq.gz | head -n 19800000 > data/lomank12/Ecoli_S1_L001_R1_001.300x.fq
less data/lomank12/Ecoli_S1_L001_R1_001.300x.fq | grep -A1 "@M0" | grep -v "@M0" | wc -c
#1501603318

#mkdir -p /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping/300
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping/300
#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_300
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_300 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --nanopore_reads /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/loman_k12_pass_porechop.300x.fq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/Ecoli_S1_L001_R1_001.300x.fq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_300 -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_300/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

#mkdir -p /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping/100
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/k12_genotyping/100
#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_100
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/k12_genotyping_100 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/k12_genotyping.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/050319_all/ecoli_pangenome_PRG_050319.fa --nanopore_reads /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/loman_k12_pass_porechop.100x.fq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/lomank12/Ecoli_S1_L001_R1_001.100x.fq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/NC_000913.3.fasta --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/logs/sequencing_summary.txt --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_100 -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/k12_genotyping_100/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
