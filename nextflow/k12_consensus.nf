params.pangenome_prg = ""
params.nanopore_reads = ""
params.illumina_reads = ""
params.reference_assembly = ""

params.help = false
params.final_outdir = "."
params.max_forks = 100

if (params.help){
    log.info"""
        Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run pandora_compare.nf <arguments>
        Required arguments:
          --pangenome_prg       FILE    PRG file to use as input to pandora
	  --nanopore_reads	FILE
	  --illumina_reads	FILE
	  --reference_assembly	FILE

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

if (params.reference_assembly) {
    reference_assembly = file(params.reference_assembly).toAbsolutePath()
    if (!reference_assembly.exists()) {
        exit 1, "Reference assembly file not found: ${params.reference_assembly} -- aborting"
    }
}
else {
    exit 1, "Reference assembly file not provided -- aborting"
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

process pandora_map_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  file reads from nanopore_reads
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora.consensus.fq.gz"), val("Nanopore") into pandora_output_nano
  
  """
  pandora map -p ${prg} -r ${reads}
  """
} 

process pandora_map_illumina {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  file prg from pangenome_prg
  file reads from illumina_reads
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs

  output:
  set file("pandora/pandora.consensus.fq.gz"), val("Illumina") into pandora_output_illumina

  """
  pandora map -p ${prg} -r ${reads} --illumina
  """
}

pandora_output_nano.concat( pandora_output_illumina ).set { pandora_output }

process compare_output_to_input {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set file(out_path), val(type) from pandora_output
  file(reference) from reference_assembly

  output:
  set file("out.sam"), val(type) into output_sam

  """
  bwa index ${reference}
  bwa mem ${reference} ${out_path} > out.sam
  """
}

process make_plot {
  memory { 0.1.GB * task.attempt } 
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false
  
  input:
  set file(samfile), val(type) from output_sam
  
  output:
  file("${type}.sam_mismatch_counts.png") into output_plot

  """
  #!/usr/bin/env python3

  import seaborn as sns
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  from collections import Counter
  import pandas as pd
  import numpy as np

  def plot_sam_dist(sam_file, type):
      total_gene_bases = 0
      total_mismatch_bases = 0
      num_false_positives = 0
      num_other = 0
      num_genes = 0
      num_mismatches = []
    
      f = open(sam_file, 'r')
      for line in f:
          if line != "" and line[0].isalpha() and (line.split('\t')[1]=="0" or line.split('\t')[1]=="16"):
              num = int(line.split('\t')[11].split("NM:i:")[-1])
              num_mismatches.append(num)
              total_mismatch_bases += num
              total_gene_bases += len(line.split('\t')[9])
              num_genes += 1
          elif line != "" and line[0].isalpha() and line.split('\t')[1]!="4":
              num_false_positives += 1
              num_genes += 1
              
      f.close()
      if total_gene_bases == 0:
          total_gene_bases = 1
      print("False positive gene rate: %d/%d = %f" %(num_false_positives, num_genes, float(num_false_positives)/num_genes*100))
      print("Estimated per base accuracy: 1 - %d/%d = %f" %(total_mismatch_bases, total_gene_bases, (1-(float(total_mismatch_bases)/total_gene_bases))*100))                
      print(Counter(num_mismatches))
      plt.rcParams['figure.figsize'] = 10,6
      fig, ax = plt.subplots()
      sns.distplot(num_mismatches, bins=range(0, max(num_mismatches)+2, 1), kde=False)
      ax.set(title="Distribution of Number of Mismatch Bases in Consensus", xlabel='Number of mismatch bases', ylabel='Frequency')
      plt.savefig('%s.sam_mismatch_counts.png' %type, transparent=True)
    
  plot_sam_dist("${samfile}", "${type}")
  """
}
