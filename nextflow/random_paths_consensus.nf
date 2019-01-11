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
  memory { 7.GB * task.attempt }
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
  pandora random_path ${prg} ${num_paths}
  gunzip random_paths.fa.gz
  """
} 

random_paths_output.splitFasta( file: true ).set {split_paths}

process simulate_nanopore_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks 8
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }

  input:
  file(path_fasta) from split_paths

  output:
  set(file("${path_fasta}"),file("simulated.fa")) into sim_reads

  """
  nanosim-h -p ecoli_R9_2D -n 100 ${path_fasta}
  #python3 /NanoSim-2.2.0/src/read_analysis.py -i ${training_reads} -r ${training_ref}
  #python3 /NanoSim-2.2.0/src/simulator.py linear -r ${path_fasta} -n 100
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

process pandora_map_path {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  set(file(path), file(reads)) from sim_reads
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora.consensus.fq.gz"), file("${path}") into pandora_output_path
  
  """
  pandora map -p ${prg} -r ${reads}
  """
} 

process compare_output_path_to_input {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks params.max_forks
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

output_sam.collectFile(name: final_outdir/'pandora_random_paths.sam').set { full_sam }

process make_plot {
  memory { 0.1.GB * task.attempt } 
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false
  
  input:
  file(samfile) from full_sam
  
  output:
  file("sam_mismatch_counts.pdf") into output_plot

  """
  #!/usr/bin/env python3

  import seaborn as sns
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  from collections import Counter
  import pandas as pd
  import numpy as np

  def plot_sam_dist(sam_file):
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
      plt.savefig('sam_mismatch_counts.png', transparent=True)
      fig.savefig('sam_mismatch_counts.pdf')
    
  plot_sam_dist("${samfile}")
  """
}
