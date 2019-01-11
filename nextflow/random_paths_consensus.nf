params.pangenome_prg = ""
params.number_paths = 500

params.help = false
params.final_outdir = "."
params.max_forks = 100

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
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input: 
  file prg from pangenome_prg
  val num_paths from params.number_paths

  output:
  file("random_paths.fa") into random_paths_output

  """
  pandora random_paths ${prg} ${num_paths}
  """
} 

split_paths = Channel
     .fromPath(random_paths_output)
     .splitFasta( file: true )

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

process pandora_map_path {
  memory { 0.8.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks parms.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  file path from split_paths
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora.consensus.fq"), file(${path}) into pandora_output_path
  
  """
  pandora map -p ${prg} -r ${path}
  """
} 

process compare_output_path_to_input {
  memory { 0.8.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks parms.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set file(out_path), file(in_path) from pandora_output_path

  output:
  file("out.sam") into output_sam

  """
  bwa index ${in_path}
  bwa mem ${in_path} ${out_path} > out.sam
  """
}

output_sam.collectFile(name: final_outdir/'pandora_random_paths.sam')
