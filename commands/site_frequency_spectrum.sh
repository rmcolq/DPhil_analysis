# Have gff assemblies in
/nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files

# Have a klebs PRG in
/nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa

# Want frequency of SNPs by frequency of gene

# To get this, first run pandora with all 265 assemblies
nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa --read_index /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files/INDEX_30 --chunk_size 4000 --num_samples 285 --k 31 --w 19 --max_covg 50 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

bsub.py 3 pandora_compare_in_chunks_100519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa --read_index /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files/INDEX_30 --chunk_size 4000 --num_samples 285 --k 31 --w 19 --max_covg 50 --min_cluster_size 3 --max_diff 100 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/site_frequency_spectrum
bsub.py 10 sfs python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/site_frequency_spectrum.py --matrix /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/250419/pandora_multisample.matrix --vcf /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/250419/pandora_multisample_genotyped.vcf --outdir sfs
