params.pangenome_prg = ""
params.nanopore_reads = ""
params.illumina_reads = ""

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
	  --nanopore_reads	FILE
	  --illumina_reads	FILE

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

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}

num_prg = Channel.from( 1..22 ).map { 1000 * it }

process pandora_index {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    }
    time '3s' 
        
    input:
    file pangenome_prg
    val num from num_prg
        
    output:
    set val("${num}"), file("prg.fa"), file("prg.fa.k15.w14.idx"), file("kmer_prgs") into indexed_prgs_nano
    set val("${num}"), file("prg.fa"), file("prg.fa.k15.w14.idx"), file("kmer_prgs") into indexed_prgs_ill
        
    """
    head -n${num} ${pangenome_prg} > prg.fa
    pandora index -w 14 -k 15 prg.fa
    """
}

process pandora_map_nano {
  memory { 8.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  set val(num), file(prg), file(index), file(kmer_prgs) from indexed_prgs_nano
  file reads from nanopore_reads
  
  output:
  set val("Nanopore"), val("${num}"), file("timeinfo.txt") into pandora_output_nano
  
  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg 30 &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
} 

process pandora_map_illumina {
  memory { 12.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  set val(num), file(prg), file(index), file(kmer_prgs) from indexed_prgs_ill
  file reads from illumina_reads

  output:
  set val("Illumina"), val("${num}"), file("timeinfo.txt") into pandora_output_illumina

  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg 30 --illumina &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
}

pandora_output_nano.concat( pandora_output_illumina ).set { pandora_output }

process make_df {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input:
  set val(type), val(num_prgs), file(timeinfo) from pandora_output

  output:
  file("out.tsv") into df_line

  """
  #!/usr/bin/env bash
  sentence=\$(grep "System time (seconds):" ${timeinfo}) 
  stringarray=(\$sentence)
  systime=\${stringarray[3]}
  echo \$systime

  sentence=\$(grep "User time (seconds):" ${timeinfo})
  stringarray=(\$sentence)
  usertime=\${stringarray[3]}
  echo \$usertime

  sentence=\$(grep "Maximum resident set size (kbytes)" ${timeinfo})
  stringarray=(\$sentence)
  maxmem=\${stringarray[5]}
  echo \$maxmem

  echo -e "${type}\t${num_prgs}\t\$systime\t\$usertime\t\$maxmem" > out.tsv
  """
}

df_line.collectFile(name: final_outdir/'scaling.tsv').set { table }

process make_plot {
  memory { 0.1.GB * task.attempt } 
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false
  
  input:
  file(data) from table
  
  output:
  file("*.png") into output_plot

  """
  python3 ${params.pipeline_root}/scripts/plot_scaling.py --tsv_file ${data} --xvar 'num_prgs' --xlabel "Number of Local Graphs in PanRG"
  """
}
