params.pangenome_prg = ""
params.number_paths = 10
params.number_genes = 1000

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
	  --number_paths	INT		Number of paths to find through each PRG, default 500
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
  awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print \$0 "\\n" seq }' random_paths.fa > tmp.random_paths_filtered.fa
  head -n\$(( 2 * ${params.number_genes} * ${params.number_paths} )) tmp.random_paths_filtered.fa > random_paths_filtered.fa
  
  """
}
 
path_nano = Channel.create()
path_illumina = Channel.create()
random_paths_output.splitFasta( file: true ).separate( path_nano, path_illumina ) { a -> [a, a] }

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

process pandora_map_path_nano {
  memory { 4.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  file read from path_nano
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora_result.fq"), file("${read}") into pandora_output_path_nano
  
  """
  pandora map -p ${prg} -r ${read} --genome_size 1000
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  v=\$(head -n1 ${read})
  zgrep -A3 \$(echo \${v:1:-3}) pandora/pandora.consensus.fq.gz > pandora_result.fq
  if [[ -s pandora_result.fq ]] ; then
  echo "pandora/pandora.consensus.fq.gz has results."
  else
  exit 1
  fi
  
  """
} 

process pandora_map_path_illumina {
  memory { 4.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }   
  
  input:
  file prg from pangenome_prg
  file read from path_illumina
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora_result.fq"), file("${read}") into pandora_output_path_illumina
  
  """
  pandora map -p ${prg} -r ${read} --illumina --genome_size 1000
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  v=\$(head -n1 ${read})
  zgrep -A3 \$(echo \${v:1:-3}) pandora/pandora.consensus.fq.gz > pandora_result.fq
  if [[ -s pandora_result.fq ]] ; then
  echo "pandora/pandora.consensus.fq.gz has results."
  else
  exit 1
  fi
  """
} 

process compare_output_path_to_input_nano {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set file(out_path), file(in_path) from pandora_output_path_nano

  output:
  file("out.sam") into output_sam_nano

  """
  bwa index ${in_path}
  bwa mem ${in_path} ${out_path} > out.sam
  """
}

process compare_output_path_to_input_illumina {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
      
  input:
  set file(out_path), file(in_path) from pandora_output_path_illumina
  
  output:
  file("out.sam") into output_sam_illumina
  
  """
  bwa index ${in_path}
  bwa mem ${in_path} ${out_path} > out.sam
  """ 
} 

output_sam_nano.collectFile(name: final_outdir/'pandora_random_paths_nano.sam').set { full_sam_nano }
output_sam_illumina.collectFile(name: final_outdir/'pandora_random_paths_illumina.sam').set { full_sam_illumina }

process make_plot_nano {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true
  
  input:
  file(samfile) from full_sam_nano
  
  output:
  file("*.png") into output_plot_nano
  file ("*_list.txt") into count_list_nano
  
  """
  python3 ${params.pipeline_root}/scripts/plot_sam_histogram.py --sam "${samfile}" --prefix "Nanopore"
  """
}
process make_plot_illumina {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true

  input:
  file(samfile) from full_sam_illumina

  output:
  file("*.png") into output_plot_ill
  file ("*_list.txt") into count_list_illumina
  
  """
  python3 ${params.pipeline_root}/scripts/plot_sam_histogram.py --sam "${samfile}" --prefix "Illumina"
  """
} 

count_list_nano.concat( count_list_illumina ).set { count_list }

process make_joint_plot {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true

  input:
  file count_files from count_list.collect()

  output:
  file("*.png") into output_joint_plot

  """
  #!/usr/bin/env python3
  import sys
  sys.path.append('${params.pipeline_root}/scripts')
  from plot_joint_histogram import plot_count_hist

  s = "${count_files}".split()
  plot_count_hist(s[0],s[1])
  """
}
