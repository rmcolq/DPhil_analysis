cd /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/multisample_roc
bsub.py 2.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/multisample_roc nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/multisample_roc.nf --truth_assembly1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz --vcf_directory1 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping/chosen_H131800734/ --truth_assembly2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz --mask1 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_mask.bed --vcf_directory2 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping/chosen_H151080744/ --mask2 /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_mask.bed -w /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/multisample_roc/work5 -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config
