#random subsample of fastq
mkdir -p /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples
for i in $(cat /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples.tsv | cut -f1)
do
echo $i
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/$i.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq -n 1500000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.45x.fq" >> /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random.tsv
done

for i in $(cat /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples.tsv | cut -f1)
do
echo $i
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/$i.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.30x.fq -n 1000000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.30x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.30x.fq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.30x.fq" >> /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random30.tsv
done

for i in $(cat /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples.tsv | cut -f1)
do
echo $i
python2 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/subsample_fastq.py /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/$i.fq /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.20x.fq -n 700000
mv /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.20x.fq.0 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.20x.fq
echo -e "$i\t/nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/random_subsamples/$i.20x.fq" >> /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random20.tsv
done

# run compare on all 16
mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/090519
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/090519
bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random.tsv --genotype --max_covg 41 --min_cluster_size 5 --illumina

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/140519_illumina
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/140519_illumina
bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random30.tsv --genotype --max_covg 30 --illumina

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/140519_nanopore
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/140519_nanopore
bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -w 14 -k 15 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/nanopore_samples.tsv --genotype --max_covg 40

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/160519_illumina
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/160519_illumina
bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random30.tsv --genotype --max_covg 30 --illumina

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/230519_illumina
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/230519_illumina
bsub.py 123 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random30.tsv --genotype --max_covg 25 --illumina
#mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/230519_nanopore
#cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/230519_nanopore
#bsub.py 125 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -w 14 -k 15 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/nanopore_samples.tsv --genotype --max_covg 30

## NB, tried to regularize a nanopore vcf with
#java -jar /hps/nobackup/iqbal/rmcolq/apps/bcbio.variation-0.2.6-standalone.jar variant-prep pandora_multisample_genotyped.vcf pandora_multisample.vcf_ref.fa --keep-ref
#but failed with an error about having too many files open...

mkdir -p /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/240519_illumina
cd /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/240519_illumina
bsub.py 123 compare singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-pandora-pandora.img pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -w 19 -k 31 -r /hps/nobackup/iqbal/rmcolq/projects/pandora_compare/analysis/cardio_15/illumina_samples_random30.tsv --genotype --illumina
