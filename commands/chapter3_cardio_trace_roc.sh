## GENOTYPE SNIPPY AND NANOPOLISH
mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping2/H131800734
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping2/H131800734
bsub.py 3.0 genotype_H131800734 nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/H131800734/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/H131800734/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/porechopped/BC03.fastq.gz --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170616/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq --truth_assembly /nfs/leia/research/iqbal/mbhall/Projects/Pandora_variation/Data/Cardio/cardio/pacbio/EColi_4346_assembly.fasta.gz -resume
#nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/H151080744/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/H151080744/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/porechopped/BC11.fastq.gz --raw_fast5s /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/reads/ --albacore_summary /hps/nobackup/iqbal/mbhall/Pandora_variation/Data/Cardio/ecoli20170412/data/basecalled/sequencing_summary.txt --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq --truth_assembly /nfs/leia/research/iqbal/mbhall/Projects/Pandora_variation/Data/Cardio/cardio/pacbio/EColi_7927_assembly.fasta.gz -resume
#cd /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/CFT073/
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/logs/CFT073_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/CFT073/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/CFT073/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299443.fastq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/trace/CFT073_subsampled_pacbio.fasta -resume
#cd /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/analysis/RHB11-C04/
#bsub.py 3.0 /nfs/leia/research/iqbal/rmcolq/projects/pandora_variation/logs/RHB11-C04_genotyping nextflow run /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/sample_genotype.nf --pangenome_prg /nfs/leia/research/iqbal/rmcolq/projects/pangenome_prg/ecoli/290818_all/ecoli_pangenome_PRG_290818.fa --reference_directory /nfs/leia/research/iqbal/rmcolq/data/references/Escherichia_coli/refseq/4_way_chosen/ --final_outdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/RHB11-C04/ --pipeline_root /nfs/leia/research/iqbal/rmcolq/git/compare_genotypers/nextflow/ -w /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/4_way_roc/genotyping2/RHB11-C04/work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config --nanopore_reads /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5299467.fastq --illumina_reads_1 /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351.fastq --truth_assembly /nfs/leia/research/iqbal/rmcolq/data/trace/RHB11-C04_subsampled_pacbio.fasta -resume

## Indexed PanRG for pandora
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 19 --k 31 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19
bsub.py 2 index_prg nextflow run /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/index_prg.nf --pangenome_prg ecoli_pangenome_PRG_050319.fa --num_prg 37426 --w 14 --k 15 --chunk_size 100 -w work -c /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/nextflow/nextflow.config -resume    ## /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14

## Ran pandora
bsub.py 60 logs/compare_cardio_trace_illumina_30 singularity exec shub://rmcolq/pandora:pandora pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/illumina_30 --genotype --max_covg 30 --illumina --max_diff 32 -w 19 -k 31
bsub.py 100 logs/compare_cardio_trace_illumina_100 singularity exec rmcolq-pandora-dev-pandora.simg pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/illumina_100 --genotype --max_covg 100 --illumina --max_diff 32 -w 19 -k 31
bsub.py 60 logs/compare_cardio_trace_nanopore_30 singularity exec shub://rmcolq/pandora:pandora pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/nanopore_30 --genotype --max_covg 30
bsub.py 100 logs/compare_cardio_trace_nanopore_100 singularity exec rmcolq-pandora-dev-pandora.simg pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k15.w14/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_trace_read_index_nanopore.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/nanopore_100 --genotype --max_covg 100 

bsub.py 60 logs/compare_cardio_illumina_30 singularity exec shub://rmcolq/pandora:pandora pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/100419/illumina_30 --genotype --max_covg 30 --illumina --max_diff 32 -w 19 -k 31  
bsub.py 60 logs/compare_cardio_illumina_100 singularity exec shub://rmcolq/pandora:pandora pandora compare -p /hps/nobackup/iqbal/rmcolq/projects/panrg/ecoli/k31.w19/ecoli_pangenome_PRG_050319.fa -r /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/data/cardio_read_index_illumina.tsv -o /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/2_cardio/100419/illumina_100 --genotype --max_covg 100 --illumina --max_diff 32 -w 19 -k 31 

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


for run in illumina_30 illumina_100 nanopore_30 nanopore_100
do
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample_genotyped.vcf | cut -f1-9,10 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/CFT073/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample_genotyped.vcf | cut -f1-9,11 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/H131800734/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample_genotyped.vcf | cut -f1-9,12 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/H151080744/pandora_genotyped_$run.vcf
less /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample_genotyped.vcf | cut -f1-9,13 > /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/RHB11-C04/pandora_genotyped_$run.vcf
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/CFT073/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/H131800734/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/H151080744/pandora_genotyped_$run.ref.fa
cp /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/pandora_compare/4_way/030419/$run/pandora_multisample.vcf_ref.fa /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/genotyping/RHB11-C04/pandora_genotyped_$run.ref.fa
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

## make bed of bad regions
cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H151080744/H151080744_pacbio_assembly.pilon.fa.gz .
gunzip H151080744_pacbio_assembly.pilon.fa.gz
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh H151080744_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H151080744.fq H151080744

cp /hps/nobackup/iqbal/mbhall/Pandora_variation/Analysis/Cardio/Polishing/data/H131800734/H131800734_pacbio_assembly.pilon.fa.gz .
gunzip H131800734_pacbio_assembly.pilon.fa.gz
bsub.py 4 build_mask_H15 bash make_low_qual_genome_mask.sh H131800734_pacbio_assembly.pilon.fa /nfs/leia/research/iqbal/rmcolq/data/cardio/illumina/H131800734.fq H131800734

bsub.py 4 build_mask_RH bash make_low_qual_genome_mask.sh ../truths/RHB11-C04_subsampled_pacbio.fasta /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287351.fastq RHB11-C04
bsub.py 4 build_mask_CF bash make_low_qual_genome_mask.sh ../truths/CFT073_subsampled_pacbio.fasta /nfs/leia/research/iqbal/rmcolq/data/trace/SRX5287344.fastq CFT073

## Make ROC
#mkdir /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/2_way_roc/080419
cd /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/2_way_roc/080419
bsub.py 1.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/cardio_trace_roc2 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/2_way_roc/genotyped_samples.tsv
mkdir analysis/4_way_roc/work
cd analysis/4_way_roc/work
bsub.py 1.0 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/logs/cardio_trace_roc4 singularity exec /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/singularity/rmcolq-Singularity_recipes-minos.img python3 /nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis/scripts/compare_genotypers_on_number_of_samples.py --sample_tsv /hps/nobackup/iqbal/rmcolq/projects/DPhil_analysis/analysis/cardio_trace/4_way_roc/genotyped_samples.tsv
