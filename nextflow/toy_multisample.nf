params.pangenome_prg = ""
params.number_paths = 2

params.help = false
params.final_outdir = "."
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.max_forks = 100

if (params.help){
    log.info"""
        Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run pandora_compare.nf <arguments>
        Required arguments:
          --pangenome_prg       FILE    PRG file to use as input to pandora

        Optional:
          --final_outdir        DIRECTORY       Where to put final output files
          --max_forks           INT		Number of concurrent jobs, default 100

    """.stripIndent()

    exit 0
}

if (params.pangenome_prg) {
    pangenome_prg = file(params.pangenome_prg).toAbsolutePath()
    if (!pangenome_prg.exists()) {
        exit 1, "Pangenome PRG file not found: ${params.pangenome_prg} -- aborting"
    }
}
else {
    exit 1, "Pangenome PRG file not provided -- aborting"
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}

process pandora_random_paths {
  memory { 7.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input: 
  file prg from pangenome_prg
  val num_paths from params.number_paths

  output:
  file("random_paths_filtered.fa") into random_paths_output

  """
  pandora random_path ${prg} ${num_paths}
  gunzip random_paths.fa.gz
  awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print \$0 "\\n" seq }' random_paths.fa > random_paths_filtered.fa
  """
}

process generate_genomes {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  
  input:
  file(random_paths) from random_paths_output
  
  output:
  file("genomes.fa") into output_genomes
  file("genomes_truth.txt") into truth
  
  """
  python3 ${params.pipeline_root}/scripts/generate_genomes.py --in_fa "${random_paths}"
  """
}
 
path_nano = Channel.create()
path_illumina = Channel.create()
truth_genomes = Channel.create()
output_genomes.splitFasta( file: true ).separate( path_nano, path_illumina, truth_genomes ) { a -> [a, a, a] }

process simulate_nanopore_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8
  time '20m'
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }

  input:
  file(path_fasta) from path_nano

  output:
  file("nanosim*.fa") into sim_reads_nano
  file("index.tsv") into sim_index_nano

  """
  v=\$(head -n1 ${path_fasta})
  name=\${v:1:10}
  echo \$name
  nanosim-h -p ecoli_R9_2D -n 5000 ${path_fasta} --unalign-rate 0 --circular
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  mv simulated.fa nanosim_\$name.fa
  echo -e "\$name\tnanosim_\$name.fa" > index.tsv
  """
}

process simulate_illumina_reads {
  memory { 20.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8
  container {
      'shub://rmcolq/Singularity_recipes:ART'
  }

  input:
  file(path_fasta) from path_illumina

  output:
  file("artsim*.fq") into sim_reads_illumina
  file("index.tsv") into sim_index_illumina

  """
  v=\$(head -n1 ${path_fasta})
  name=\${v:1:10}
  echo \$name
  art_illumina -ss HS25 -i ${path_fasta} -l 150 -f 100 -o simulated
  if [[ -s simulated.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated.fq
  exit 1
  fi
  mv simulated.fq artsim_\$name.fq
  echo -e "\$name\tartsim_\$name.fq" > index.tsv
  """
}

sim_index_nano.collectFile(name: 'nano_index.tsv').set { index_nano }
sim_index_illumina.collectFile(name: 'illumina_index.tsv').set { index_illumina }

pandora_idx = file(pangenome_prg + '.k15.w14.idx')
pandora_kmer_prgs = file(pangenome_prg.parent / 'kmer_prgs')
if (!pandora_idx.exists()) {
    process pandora_index {
        memory { 20.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/pandora:pandora'
        }

        publishDir pangenome_prg.parent, mode: 'copy', overwrite: false

        input:
        file pangenome_prg

        output:
        file "${pangenome_prg}.k15.w14.idx" into pandora_idx
        file "kmer_prgs" into pandora_kmer_prgs

        """
        pandora index -w 14 -k 15 ${pangenome_prg}
        """
    }
}

process pandora_compare_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  file tsv from index_nano
  file reads from sim_reads_nano.collect()
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora_multisample*") into pandora_output_nano
  
  """
  pandora compare -p ${prg} -r ${tsv} -o pandora --max_covg 30 --genome_size 200000 --genotype
  """
} 

process pandora_compare_illumina {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  file prg from pangenome_prg
  file tsv from index_illumina
  file reads from sim_reads_illumina.collect()
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs

  output:
  set file("pandora/pandora_multisample*") into pandora_output_illumina
  set file("pandora_illumina*") into pandora_illumina_vcfs

  """
  pandora compare -p ${prg} -r ${tsv} -o pandora --illumina --max_covg 30 --genome_size 200000 --genotype
  cat pandora/pandora_multisample_genotyped.vcf | cut -f1-9,10 > pandora_illumina_genome1.vcf
  cp pandora/pandora_multisample.vcf_ref.fa pandora_illumina_genome1.vcf_ref.fa
  cat pandora/pandora_multisample_genotyped.vcf | cut -f1-9,11 > pandora_illumina_genome2.vcf
  cp pandora/pandora_multisample.vcf_ref.fa pandora_illumina_genome2.vcf_ref.fa
  cat pandora/pandora_multisample_genotyped.vcf | cut -f1-9,12 > pandora_illumina_genome3.vcf
  cp pandora/pandora_multisample.vcf_ref.fa pandora_illumina_genome3.vcf_ref.fa
  cat pandora/pandora_multisample_genotyped.vcf | cut -f1-9,13 > pandora_illumina_genome4.vcf
  cp pandora/pandora_multisample.vcf_ref.fa pandora_illumina_genome4.vcf_ref.fa
  """
}

process minos_evaluate_recall {
   errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
   maxRetries 1
   memory {1.4.GB * task.attempt}
   container {
       'shub://rmcolq/Singularity_recipes:minos'
   }

   publishDir final_outdir, mode: 'copy', overwrite: true

   input:
   file pandora_files from pandora_illumina_vcfs
   file ref_files from truth_genomes.collect()

   output:

   """
   python3 ${params.pipeline_root}/scripts/compare_toy_genomes.py  --dir "."
   """ 
}
