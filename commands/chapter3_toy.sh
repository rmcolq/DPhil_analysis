mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/toy_multisample nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/toy_multisample.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/toy_multisample/work2 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume
