params.nanopore_reads = ""
params.illumina_reads = ""
params.truth_assembly = ""
params.reference_dir = ""

params.help = false
params.final_outdir = "."
params.max_forks = 100
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"

if (params.help){
    log.info"""
        Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run pandora_compare.nf <arguments>
        Required arguments:
	  --nanopore_reads	FILE
	  --illumina_reads	FILE
	  --truth_assembly	FILE
          --reference_dir	DIRECTORY

        Optional:
          --final_outdir        DIRECTORY       Where to put final output files
          --max_forks           INT		Number of concurrent jobs, default 100

    """.stripIndent()

    exit 0
}

if (params.nanopore_reads) {
    nanopore_reads = file(params.nanopore_reads).toAbsolutePath()
    if (!nanopore_reads.exists()) {
        exit 1, "Nanopore read file not found: ${params.nanopore_reads} -- aborting"
    }
}
else {
    exit 1, "Nanopore read file not provided -- aborting"
}

if (params.illumina_reads) {
    illumina_reads = file(params.illumina_reads).toAbsolutePath()
    if (!illumina_reads.exists()) {
        exit 1, "Illumina read file not found: ${params.illumina_reads} -- aborting"
    }
}   
else {
    exit 1, "Illumina read file not provided -- aborting"
}

if (params.truth_assembly) {
    truth_assembly = file(params.truth_assembly).toAbsolutePath()
    if (!truth_assembly.exists()) {
        exit 1, "Truth assembly file not found: ${params.truth_assembly} -- aborting"
    }
}
else {
    exit 1, "Truth assembly file not provided -- aborting"
}

if (params.reference_dir) {
    reference_dir = file(params.reference_dir).toAbsolutePath()
    if (!reference_dir.exists()) { 
        exit 1, "Reference directory not found: ${params.reference_dir} -- aborting"
    }
}
else {
    exit 1, "Reference directory not provided -- aborting"
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}

process map_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  
  input:
  file reads from nanopore_reads
  file ref from reference_dir
  
  output:
  file("*.sam") into output_nano
  
  """
  bwa mem -x ont2d ${ref} ${reads} > ${ref}.nano.sam
  """
} 

process map_ill {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  file reads from illumina_reads
  file ref from reference_dir

  output:
  file("*.sam") into output_ill

  """
  bwa mem ${ref} ${reads} > ${ref}.ill.sam
  """
}

process fastani {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:fastani'
  }

  input:
  file truth from truth_assembly
  file reference_dir

  output:
  file("out.tsv") into output_fastani

  """
  ls ${reference_dir} > ref_list
  fastANI -q ${truth} -rl ref_list -o out.tsv
  """
}

