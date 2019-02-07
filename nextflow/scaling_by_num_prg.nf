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

num_prg = Channel.from( 1..30 ).map { 100 * it }

process pandora_index {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }
        
    input:
    file pangenome_prg
    val num from num_prg
        
    output:
    set val("${num}"), file("prg.fa"), file("prg.k15.w14.idx"), file("kmer_prgs") into indexed_prgs
        
    """
    head -n${num} ${pangenome_prg} > prg.fa
    pandora index -w 14 -k 15 prg.fa
    """
}

process pandora_map_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  set val(num), file(prg), file(index), file(kmer_prgs) into indexed_prgs
  file reads from nanopore_reads
  
  output:
  set val("Nanopore"), val("${num}"), file("timeinfo.txt") into pandora_output_nano
  
  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg ${covg} &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
} 

process pandora_map_illumina {
  memory { 55.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  set val(num), file(prg), file(index), file(kmer_prgs) into indexed_prgs
  file reads from illumina_reads

  output:
  set val("Illumina"), val("${num}"), file("timeinfo.txt") into pandora_output_illumina

  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg ${covg} --illumina &> pandora.log" > command.sh
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
  cputime=\${stringarray[3]}
  echo \$cputime

  sentence=\$(grep "Maximum resident set size (kbytes)" ${timeinfo})
  stringarray=(\$sentence)
  maxmem=\${stringarray[5]}
  echo \$maxmem

  echo -e "${type}\t${num_prgs}\t\$cputime\t\$maxmem" > out.tsv
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
  #!/usr/bin/env python3

  import seaborn as sns
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  from collections import Counter
  import pandas as pd
  import numpy as np
 
  def plot_df(tsv_file):
      df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', 'num_prgs', 'time', 'max_mem'])
      df['max_mem_gb'] = df['max_mem']/(1024*1024)
      df['num_prgs'] = df['num_prgs']/2
      plt.rcParams['figure.figsize'] = 10,6
      fig, ax = plt.subplots()

      ax.grid(b=True)
      ax.set_axisbelow(b=True)
      plt.style.use('seaborn-colorblind')

      # Label the axes and title the plot
      ax.set_xlabel('Number of Local Graphs in Index', size = 26)
      ax.set_ylabel('Time(s)', size = 26)
      
      sns.lmplot( x="num_prgs", y="time", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
      plt.legend(loc='lower right')
      plt.savefig('scaling_time.png', transparent=True)

      sns.lmplot( x="num_prgs", y="max_mem", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
      plt.legend(loc='lower right')
      ax.set_ylabel('Max Memory (KB)', size = 26)
      plt.savefig('scaling_mem.png', transparent=True)

  plot_df("${data}")
  """
}

