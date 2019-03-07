## GENOTYPE SNIPPY AND NANOPOLISH
#nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/chosen/ --final_outdir /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/porechopped/BC03.fastq.gz --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq --truth_assembly /nfs/leia/research/iqbal/mbhall/Projects/Pandora_variation/Data/Cardio/cardio/pacbio/EColi_4346_assembly.fasta.gz -resume
#nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/chosen/ --final_outdir /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/porechopped/BC11.fastq.gz --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq --truth_assembly /nfs/leia/research/iqbal/mbhall/Projects/Pandora_variation/Data/Cardio/cardio/pacbio/EColi_7927_assembly.fasta.gz -resume
#cd /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/logs/CFT073_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/chosen/ --final_outdir /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/genotyping --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/genotyping/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299443.fastq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/trace/CFT073_subsampled_pacbio.fasta -resume
#cd /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/logs/RHB11-C04_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/chosen/ --final_outdir /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/genotyping --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/genotyping/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299467.fastq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351.fastq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/trace/RHB11-C04_subsampled_pacbio.fasta -resume

## Copied over samples from yoda to ebi-cli to run pandora compare on
## FROM YODA - did not work from ebi-cli to yoda
# bsub.py 1 logs/rsync_yoda_to_cli_prg rsync -avcr --no-perms /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all rmcolq@ebi-cli.ebi.ac.uk:/nfs/research1/zi/rmcolq/projects/pandora_compare/prg/ecoli/
#bsub.py 1 logs/rsync_yoda_to_cli_ill rsync -avcr --no-perms /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina rmcolq@ebi-cli.ebi.ac.uk:/nfs/research1/zi/rmcolq/projects/pandora_compare/data/cardio/
#bsub.py 1 logs/rsync_yoda_to_cli_nano rsync -avcr --no-perms /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ rmcolq@ebi-cli.ebi.ac.uk:/nfs/research1/zi/rmcolq/projects/pandora_compare/data/cardio/nanopore/

## ON EBI_CLI
#bsub.py 1 logs/download_simg singularity exec 'shub://rmcolq/pandora:pandora' pandora
#mkdir singularity
#mv rmcolq-pandora-dev-pandora.simg singularity/
#bsub.py --queue 'research-rh74' 100 logs/compare_cardio_trace_nano singularity exec singularity/rmcolq-pandora-dev-pandora.simg pandora compare -p /nfs/research1/zi/rmcolq/projects/pandora_compare/prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa -r /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/read_index_nanopore.tsv -o /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/nanopore_full --genotype --max_covg 100
#bsub.py --queue 'research-rh74' 60 logs/compare_cardio_trace_ill singularity exec singularity/rmcolq-pandora-dev-pandora.simg pandora compare -p /nfs/research1/zi/rmcolq/projects/pandora_compare/prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa -r /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/read_index_illumina.tsv -o /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/illumina --genotype --max_covg 30
#bsub.py --queue 'research-rh74' 30 logs/compare_cardio_trace_nano_30 singularity exec singularity/rmcolq-pandora-dev-pandora.simg pandora compare -p /nfs/research1/zi/rmcolq/projects/pandora_compare/prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa -r /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/read_index_nanopore.tsv -o /nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/nanopore_30 --genotype --max_covg 30

## Copied compare results back to yoda
#bsub.py 1 logs/rsync_cli_to_yoda_ill rsync -avcr --no-perms rmcolq@ebi-cli.ebi.ac.uk:/nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/ /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare 
## Again after updated genotyping on 2nd March 2019
bsub.py 1 logs/rsync_cli_to_yoda rsync -avcr --no-perms rmcolq@ebi-cli.ebi.ac.uk:/nfs/research1/zi/rmcolq/projects/pandora_compare/analysis/compare_trace_and_cardio/ /nfs/lei
a/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/020319/

## Create single sample vcf and ref.fa files for each sample
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/genotyping/snippy* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/genotyping/nanopolish* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H131800734
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping/chosen_H131800734/snippy* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H131800734/
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H131800734/genotyping/chosen_H131800734/nanopolish* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H131800734/
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H151080744
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping/chosen_H151080744/snippy* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H151080744/
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/cardio_H151080744/genotyping/chosen_H151080744/nanopolish* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H151080744/
mkdir /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/genotyping/snippy* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/
rsync /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/genotyping/nanopolish* /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/


for run in illumina nanopore_30 nanopore_full
do
less /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/pandora_genotyped_$run.vcf
less /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H131800734/pandora_genotyped_$run.vcf
less /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample_genotyped.vcf | cut -f1-9,12 > /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H151080744/pandora_genotyped_$run.vcf
less /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample_genotyped.vcf | cut -f1-9,13 > /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/pandora_genotyped_$run.vcf
cp /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample.vcf_ref.fa /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/pandora_genotyped_$run.ref.fa
cp /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample.vcf_ref.fa /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H131800734/pandora_genotyped_$run.ref.fa
cp /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample.vcf_ref.fa /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/H151080744/pandora_genotyped_$run.ref.fa
cp /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/pandora_compare/$run/pandora_multisample.vcf_ref.fa /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/pandora_genotyped_$run.ref.fa
done

## Some snippy vcf files are empty and so these will cause a fail
#for f in /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_*.vcf ; do echo $f; grep -v "#" $f | wc -l; done
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_CP010148.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_CP010163.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_CP010200.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_CP010219.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_CP010237.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_NC_013361.1.vcf
#4
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_NC_017628.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_NZ_AP014857.1.vcf
#3
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_NZ_CM000662.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/CFT073/snippy_NZ_CP012140.1.vcf
#4
#[rmcolq@hh-yoda-08-01 DPhil_analysis]$ for f in /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_*.vcf ; do echo $f; grep -v "#" $f | wc -l; done
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_CP010148.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_CP010163.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_CP010200.1.vcf
#8
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_CP010219.1.vcf
#1
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_CP010237.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_NC_013361.1.vcf
#14
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_NC_017628.1.vcf
#4
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_NZ_AP014857.1.vcf
#0
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_NZ_CM000662.1.vcf
#8
#/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyping/RHB11-C04/snippy_NZ_CP012140.1.vcf
#8
## REMOVED THESE VCFS

## Make tsv of samples
## Fix dots in vcf
for f in analysis/4_way_roc/genotyping/*/*.ref.fa; do cat $f | cut -d "." -f1 > tmp.fa; mv tmp.fa $f; done

## Make ROC
mkdir analysis/4_way_roc/work
cd analysis/4_way_roc/work
bsub.py 20.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/4_way_roc singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/analysis/4_way_roc/genotyped_samples.tsv
