# Have gff assemblies in
/nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files

# Have a klebs PRG in
/nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa

# Want frequency of SNPs by frequency of gene

# To get this, first run pandora with all 265 assemblies
#bsub.py 3 pandora_compare_in_chunks_100519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa --read_index /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files/INDEX_30 --chunk_size 4000 --num_samples 285 --k 31 --w 19 --max_covg 50 --min_cluster_size 3 --max_diff 100 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

bsub.py 2 pandora_compare_in_chunks_140519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa --read_index /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/gff_files/INDEX_30 --chunk_size 4000 --num_samples 285 --k 31 --w 19 --max_covg 50 -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume


# try instead with illumina reads
for i in $(cat /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/INDEX | cut -f1)
do
echo $i
zcat /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/"$i"_both.fastq.gz > /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/"$i"_both.fastq
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/"$i"_both.fastq /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/$i.45x.fastq -n 1500000
mv /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/$i.45x.fastq.0 /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/$i.45x.fastq
echo -e "$i\t/nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/$i.45x.fastq" >> /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/INDEX_45
rm /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/catted_fastq/$i/"$i"_both.fastq
done

bsub.py 3 pandora_compare_in_chunks_100519 nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/pandora_compare_in_chunks.nf --pangenome_prg /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/prg/holt2015/combined_genes_igr10/pangenome_PRG.fa --read_index /nfs/leia/research/iqbal/projects/pandora/klebs/holt_gastro/data/illumina/INDEX_45 --chunk_size 4000 --num_samples 151 --k 31 --w 19 --max_covg 50 --illumina -w pandora_work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume

mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/site_frequency_spectrum/100519
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/site_frequency_spectrum/100519
mkdir gfs
bsub.py 10 gfs python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/gene_frequency_spectrum.py --matrix /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample.matrix --vcf /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample_genotyped.vcf --outdir gfs
mkdir sfs
bsub.py 10 sfs python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/site_frequency_spectrum.py --matrix /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample.matrix --vcf /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample_genotyped.vcf --outdir sfs
mkdir sfs_150519
bsub.py 10 sfs python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/site_frequency_spectrum.py --matrix /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample.matrix --vcf /hps/nobackup2/iqbal/projects/pandora/klebs/holt_gastro/100519/pandora_multisample_genotyped.vcf --outdir sfs_150519
