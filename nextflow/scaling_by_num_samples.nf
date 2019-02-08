params.pangenome_prg = ""
params.nanopore_tsv = ""
params.illumina_tsv = ""

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
	  --nanopore_tsv	FILE
	  --illumina_tsv	FILE

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

if (params.nanopore_tsv) {
    nanopore_tsv = file(params.nanopore_tsv).toAbsolutePath()
    if (!nanopore_tsv.exists()) {
        exit 1, "Nanopore read file not found: ${params.nanopore_tsv} -- aborting"
    }
}
else {
    exit 1, "Nanopore read file not provided -- aborting"
}

if (params.illumina_tsv) {
    illumina_tsv = file(params.illumina_tsv).toAbsolutePath()
    if (!illumina_tsv.exists()) {
        exit 1, "Illumina read file not found: ${params.illumina_tsv} -- aborting"
    }
}   
else {
    exit 1, "Illumina read file not provided -- aborting"
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
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

num_samples_nano = Channel.from( 2..15 )
num_samples_illu = Channel.from( 2..15 )

process pandora_compare_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  val num_samples from num_samples_nano
  file prg from pangenome_prg
  file reads from nanopore_tsv
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set val("Nanopore"), val("${num_samples}"), file("timeinfo.txt") into pandora_output_nano
  
  """
  head -n ${num_samples} ${reads} > subreads.tsv
  echo "pandora compare -p ${prg} -r subreads.tsv --genotype --max_covg 30 &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
} 

process pandora_compare_illumina {
  memory { 55.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  val num_samples from num_samples_illu
  file prg from pangenome_prg
  file reads from illumina_tsv
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs

  output:
  set val("Illumina"), val("${num_samples}"), file("timeinfo.txt") into pandora_output_illumina

  """
  head -n ${num_samples} ${reads} > subreads.tsv
  echo "pandora compare -p ${prg} -r subreads.tsv --genotype --max_covg 30 --illumina &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
}

pandora_output_nano.concat( pandora_output_illumina ).set { pandora_output }

process make_df {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input:
  set val(type), val(num_samples), file(timeinfo) from pandora_output

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

  echo -e "${type}\t${num_samples}\t\$cputime\t\$maxmem" > out.tsv
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
      df = pd.read_csv(tsv_file, sep='\t', header=None, names=['type', 'num_samples', 'time', 'max_mem'])
      df['max_mem_gb'] = df['max_mem']/(1024*1024)
      plt.rcParams['figure.figsize'] = 10,6
      fig, ax = plt.subplots()

      ax.grid(b=True)
      ax.set_axisbelow(b=True)
      plt.style.use('seaborn-colorblind')

      # Label the axes and title the plot
      ax.set_xlabel('Number of Samples', size = 26)
      ax.set_ylabel('Time(s)', size = 26)
      
      g = sns.lmplot( x="num_samples", y="time", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
      g = (g.set_axis_labels("Number of Samples", "Time (s)")
      plt.legend(loc='lower right')
      plt.savefig('scaling_time.png', transparent=True)

      g = sns.lmplot( x="num_samples", y="max_mem_gb", data=df, fit_reg=False, hue='type', legend=False, palette="colorblind")
      g = (g.set_axis_labels("Number of Samples", "Max Memory (GB)")
      plt.legend(loc='lower right')
      ax.set_ylabel('Max Memory (GB)', size = 26)
      plt.savefig('scaling_mem.png', transparent=True)

  plot_df("${data}")
  """
}

