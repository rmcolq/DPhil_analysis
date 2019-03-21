params.pangenome_prg = ""
params.num_paths = 1
params.num_genes = 5000
params.index_dir = ""
params.w = 14
params.k = 15
params.covg = 50

params.help = false
params.final_outdir = "."
params.max_forks = 100
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"

if (params.help){
    log.info"""
        Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run pandora_compare.nf <arguments>
        Required arguments:
          --pangenome_prg       FILE    PRG file to use as input to pandora

        Optional:
          --index_dir
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
if (params.index_dir) {
    index_dir = file(params.index_dir).toAbsolutePath()
    if (!index_dir.exists()) {
        exit 1, "Index directory not found: ${params.index_dir} -- aborting"
    }   
}   
else {
    index_dir = pangenome_prg.parent
}  

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}

process find_prg_index {
memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    } 
    
    input:
    file index_dir
    val w from params.w
    val k from params.k
    
    output:
    set file("${index_dir}/k${k}.w${w}/*prg.fa"), file("${index_dir}/k${k}.w${w}/*prg.*.idx"), file("${index_dir}/k${k}.w${w}/kmer_prgs"), val("${w}"), val("${k}") into indexes_nano
    set file("${index_dir}/k${k}.w${w}/*prg.fa"), file("${index_dir}/k${k}.w${w}/*prg.*.idx"), file("${index_dir}/k${k}.w${w}/kmer_prgs"), val("${w}"), val("${k}") into indexes_ill
    """
    """
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

  output:
  file("random_paths_filtered.fa") into random_paths

  """
  pandora random_path ${prg} ${params.num_paths}
  gunzip random_paths.fa.gz
  awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print \$0 "\\n" seq }' random_paths.fa > random_paths_filtered.fa
  """
}

process simulate_genome {
  memory { 7.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input: 
  file paths from random_paths

  output:
  file("sim_genome.fa") into nano_genome
  file("sim_genome.fa") into ill_genome
  file("truth.txt") into truth

  """
  python3 ${params.pipeline_root}/scripts/sim_genome.py --in_fa ${paths} --num_genes ${params.num_genes}
  """
}
 
types = Channel.from("ecoli_R9_2D", "ecoli_R9_1D")
process simulate_nanopore_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }

  input:
  file ref_fasta from nano_genome
  val(type) from types

  output:
  set file("simulated*.fa"), val("nanopore_${type}") into nanopore_reads

  """
  nanosim-h -p ${type} -n 50000 ${ref_fasta} --unalign-rate 0 --max-len 10000 --circular
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  mv simulated.fa simulated_nanopore_${type}.fa
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
  file(ref_fasta) from ill_genome

  output:
  set file("simulated*.f*"), val("illumina_HS25") into illumina_reads

  """
  art_illumina -ss HS25 -i ${ref_fasta} -l 150 -f 100 -o simulated_illumina_HS25
  if [[ -s simulated_illumina_HS25.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated*.fq
  if [[ -s simulated_illumina_HS25.aln ]] ; then
  grep -v "#" simulated_illumina_HS25.aln | grep -v "@" > simulated_illumina_HS25.fa
  fi
  if [[ -s simulated_illumina_HS25.fa ]] ; then
  exit 1
  fi
  fi
  """
}

Channel.from(5,7,9,11,13,15).set { min_cluster_size }
Channel.from(10,20,40,80).set { max_cluster_distance }
min_cluster_size.combine(max_cluster_distance).set { thresh }
(thresh_nano, thresh_ill) = thresh.separate(2) { a -> [a, a] }

process pandora_map_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  set val(min_cluster_size), val(max_cluster_distance) from thresh_nano
  set file(reads), val(type) from nanopore_reads
  set file(prg), file(index), file(kmer_prgs), val(w), val(k) from indexes_nano

  output:
  set val("${type}"), file("pandora/pandora.consensus.fq"), val("${min_cluster_size}"), val("${max_cluster_distance}") into pandora_output_nano

  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg ${params.covg} -w ${w} -k ${k} --max_diff ${max_cluster_distance} --min_cluster_size ${min_cluster_size} &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt

  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  gunzip pandora/pandora.consensus.fq.gz
  """
}

process pandora_map_path_illumina {
  memory { 24.GB * task.attempt }
  errorStrategy {task.attempt < 4 ? 'retry' : 'ignore'}
  maxRetries 4
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }   
  
  input:
  set val(min_cluster_size), val(max_cluster_distance) from thresh_ill
  set file(reads), val(type) from illumina_reads
  set file(prg), file(index), file(kmer_prgs), val(w), val(k) from indexes_ill
  
  output:
  set val("${type}"), file("pandora/pandora.consensus.fq"), val("${min_cluster_size}"), val("${max_cluster_distance}") into pandora_output_ill
  
  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --illumina --max_covg ${params.covg} -w ${w} -k ${k} --max_diff ${max_cluster_distance} --min_cluster_size ${min_cluster_size} &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt

  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  gunzip pandora/pandora.consensus.fq.gz
  """
} 

pandora_output_ill.concat(pandora_output_nano).set { consensus }

process evaluate_genes_found {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set val(type), file(pandora_run), val(min_cluster_size), val(max_cluster_distance) from consensus
  file(truth_list) from truth

  output:
  file("results.tsv") into output_tsv

  """
  python3 ${params.pipeline_root}/scripts/finding_genes_eval.py --fastq ${pandora_run} --truth ${truth_list} --prefix "${type}\t${min_cluster_size}\t${max_cluster_distance}" --covg ${params.covg}
  """
}

output_tsv.collectFile(name: "${final_outdir}/gene_finding_by_cluster_params.tsv").set { results }

process make_plot {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true

  input:
  file(tsv) from results

  output:
  file("*.png") into output_plot

  """
  python3 ${params.pipeline_root}/scripts/plot_gene_finding.py --tsv "${tsv}" --p1 'min_cluster_size' --p2 'max_cluster_distance' --l1 "Minimum cluster size" --l2 "Maximum within cluster distance" --d1 11 --d2 80
  """
}
