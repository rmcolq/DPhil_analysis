params.pangenome_prg = ""
params.number_paths = 5
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
  set(file("${path_fasta}"),file("filtered_simulated.fa")) into sim_reads_nano

  """
  head -n1 ${path_fasta} > random_path.fa
  for i in {1..250}
  do
  echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" >> random_path.fa
  done
  head -n2 ${path_fasta} | tail -n1 >> random_path.fa
  for i in {1..250}
  do
  echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" >> random_path.fa
  done
  nanosim-h -p ecoli_R9_2D -n 250 random_path.fa --unalign-rate 0 --max-len 10000
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  minimap2 ${path_fasta} simulated.fa > out.mmp
  less out.mmp | cut -f1 > list_ids
  grep -A1 -f list_ids simulated.fa | grep -v \"\\-\\-\" > filtered_simulated.fa
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
  set(file("${path_fasta}"),file("simulated.fq")) into sim_reads_illumina

  """
  head -n1 ${path_fasta} > random_path.fa
  echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" >> random_path.fa
  head -n2 ${path_fasta} | tail -n1 >> random_path.fa
  echo "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" >> random_path.fa
  art_illumina -ss HS25 -i random_path.fa -l 150 -f 100 -o simulated
  if [[ -s simulated.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated.fq
  exit 1
  fi
  """
}

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
  set(file(path), file(reads)) from sim_reads_nano
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora_result.fq"), file("${path}") into pandora_output_path_nano
  
  """
  pandora map -p ${prg} -r ${reads} --genome_size 25000 --max_covg 100
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  v=\$(head -n1 ${path})
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
  set(file(path), file(reads)) from sim_reads_illumina
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora_result.fq"), file("${path}") into pandora_output_path_illumina
  
  """
  pandora map -p ${prg} -r ${reads} --illumina --genome_size 1000
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  v=\$(head -n1 ${path})
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
  file("with_flanks.sam") into output_sam_nano_flanks
  file("without_flanks.sam") into output_sam_nano_none

  """
  bwa index ${in_path}
  bwa mem ${in_path} ${out_path} > with_flanks.sam

  python3 ${params.pipeline_root}/scripts/remove_flank.py --in_fq ${out_path} --out_fa "out.fa" --flank_size 28
  bwa mem ${in_path} out.fa > without_flanks.sam
  python3 ${params.pipeline_root}/scripts/filter_bad_sam.py --sam without_flanks.sam
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
  file("with_flanks.sam") into output_sam_illumina_flanks
  file("without_flanks.sam") into output_sam_illumina_none

  """
  bwa index ${in_path}
  bwa mem ${in_path} ${out_path} > with_flanks.sam

  python3 ${params.pipeline_root}/scripts/remove_flank.py --in_fq ${out_path} --out_fa "out.fa" --flank_size 28
  bwa mem ${in_path} out.fa > without_flanks.sam
  python3 ${params.pipeline_root}/scripts/filter_bad_sam.py --sam without_flanks.sam

  """
}

output_sam_nano_flanks.collectFile(name: final_outdir/'with_flanks_Nanopore.sam').set { full_sam_nano_flanks }
output_sam_nano_none.collectFile(name: final_outdir/'without_flanks_Nanopore.sam').set { full_sam_nano_none }
output_sam_illumina_flanks.collectFile(name: final_outdir/'with_flanks_Illumina.sam').set { full_sam_illumina_flanks }
output_sam_illumina_none.collectFile(name: final_outdir/'without_flanks_Illumina.sam').set { full_sam_illumina_none }

full_sam_nano_flanks.concat( full_sam_nano_none, full_sam_illumina_flanks, full_sam_illumina_none ).set { full_sam }

process make_plot_sam {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true

  input:
  file(samfile) from full_sam

  output:
  file("*.png") into output_plot
  file ("*_list.txt") into count_list

  """
  python3 ${params.pipeline_root}/scripts/plot_sam_histogram.py --sam "${samfile}" --prefix "${samfile.simpleName}"
  """
}

process make_joint_plot {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: true

  input:
  file count_files from count_list.collect()

  output:
  file("*sam_mismatch_counts.png") into output_joint_plot

  """
  python3 ${params.pipeline_root}/scripts/plot_joint_histogram.py --f \"${count_files}\"
  """
}
