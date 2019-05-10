#random subsample of fastq
mkdir -p /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples
for i in $(cat /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples.tsv | cut -f1)
do
echo $i
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/$i.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq -n 1500000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq" >> /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random.tsv
done

# run compare on all 16
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/090519
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/090519
bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random.tsv --genotype --max_covg 41 --min_cluster_size 5 --illumina
