#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_covg/
cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_covg/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/scaling_by_covg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/scaling_by_covg.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_covg/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_covg/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
#mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_num_prg/
#cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_num_prg/
#bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/scaling_by_num_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/scaling_by_num_prg.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Loman_K12/loman_k12_pass_porechop.fastq.gz --illumina_reads /nfs/leia/research/iqbal/rmcolq/data/loman_k12/Ecoli_S1_L001_R2_001.fastq.gz --final_outdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_num_prg/ -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/scaling_by_num_prg/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
