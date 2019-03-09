params.pangenome_prg = ""
params.nanopore_reads = ""
params.illumina_reads = ""
params.reference_assembly = ""

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
  memory { 55.GB * task.attempt }
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

pandora_output1 = Channel.create()
pandora_output2 = Channel.create()
pandora_output.separate( pandora_output1, pandora_output2 ) { a -> [a, a] }

process compare_output_to_input_with_flanks {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  maxForks params.max_forks

  input:
  set file(out_path), val(type) from pandora_output1
  file(reference) from reference_assembly

  output:
  set file("out.sam"), val(type), val("with_flanks") into output_sam_with_flanks

  """
  bwa index ${reference}
  bwa mem ${reference} ${out_path} > out.sam
  """
}

process compare_output_to_input_without_flanks {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  maxForks params.max_forks

  input:
  set file(out_path), val(type) from pandora_output2
  file(reference) from reference_assembly

  output:
  set file("out.sam"), val(type), val("without_flanks") into output_sam_without_flanks
  set file("filtered.bam"), val(type) into output_sam_for_alignqc

  """
  bwa index ${reference}
  python3 ${params.pipeline_root}/scripts/remove_flank.py --in_fq ${out_path} --out_fa "out.fa" --flank_size 28
  bwa mem ${reference} out.fa > out.sam
  python3 ${params.pipeline_root}/scripts/filter_bad_sam.py --sam out.sam
  """
}

output_sam_with_flanks.concat( output_sam_without_flanks ).set { output_sam }

process make_plot {
  memory { 0.1.GB * task.attempt } 
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false
  
  input:
  set file(samfile), val(type), val(wflanks) from output_sam
  
  output:
  file("${wflanks}_${type}.*.png") into output_plot
  file ("${wflanks}_${type}_list.txt") into count_list

  """
  python3 ${params.pipeline_root}/scripts/plot_sam_histogram.py --sam "${samfile}" --prefix "${wflanks}_${type}" --unique
  """
}

process make_joint_plot {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false

  input:
  file count_files from count_list.collect()

  output:
  file("*joint.*sam_mismatch_counts.png") into output_joint_plot

  """
#!/usr/bin/env python3
import sys
sys.path.append('${params.pipeline_root}/scripts')
from plot_joint_histogram import plot_count_hist

wflank_dict = {}
for file in "${count_files}".split():
    wflank = file.split('_')[0]
    if wflank not in wflank_dict.keys():
        wflank_dict[wflank] = []
    wflank_dict[wflank].append(file)
    
for prefix in wflank_dict.keys():
    if len(wflank_dict[prefix]) == 2:
       print(wflank_dict[prefix][0], wflank_dict[prefix][1], prefix)
       plot_count_hist(wflank_dict[prefix][0], wflank_dict[prefix][1], prefix + "_flanks_")
  """
}

process alignqc_report {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/Singularity_recipes:alignqc'
  }
  publishDir final_outdir, mode: 'copy', overwrite: false

  input:
  set file(bamfile), val(type) from output_sam_for_alignqc
  file(reference_assembly)

  output:
  file("*.xhtml")

  """
  samtools sort -O BAM -o ${bamfile}.sorted.bam ${bamfile}
  alignqc analyze ${bamfile}.sorted.bam -g ${reference_assembly} --no_transcriptome -o ${type}.alignqc_report.xhtml
  """
}
