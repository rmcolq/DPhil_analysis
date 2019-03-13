params.pangenome_prg = ""
params.nanopore_reads = ""
params.illumina_reads = ""

params.help = false
params.final_outdir = "."
params.max_forks = 100
params.covg = 100
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

min_cluster_size_ill = Channel.from(5..15)
max_cluster_distance_ill1 = Channel.from(1..10).map { 5 * it }
max_cluster_distance_ill2 = Channel.from(2..10).map { 50 * it }
max_cluster_distance_ill1.concat ( max_cluster_distance_ill2 ).set { max_cluster_distance_ill }
min_cluster_size_ill.combine(max_cluster_distance_ill).set { thresh_ill }

min_cluster_size_nano = Channel.from(5..15)
max_cluster_distance_nano = Channel.from(1..10).map { 50 * it }
min_cluster_size_nano.combine(max_cluster_distance_nano).set { thresh_nano }

process pandora_map_nano {
  memory { 40.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  set val(min_cluster_size), val(max_cluster_distance) from thresh_nano
  file prg from pangenome_prg
  file reads from nanopore_reads
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set val("Nanopore"), val("${min_cluster_size}"), val("${max_cluster_distance}"), file("timeinfo.txt"), file("pandora.log") into pandora_output_nano
  
  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg ${params.covg} --max_diff ${max_cluster_distance} --min_cluster_size ${min_cluster_size} &> pandora.log" > command.sh
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
  set val(min_cluster_size), val(max_cluster_distance) from thresh_ill
  file prg from pangenome_prg
  file reads from illumina_reads
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs

  output:
  set val("Illumina"), val("${min_cluster_size}"), val("${max_cluster_distance}"), file("timeinfo.txt"), file("pandora.log") into pandora_output_illumina

  """
  echo "pandora map -p ${prg} -r ${reads} --genotype --max_covg ${params.covg} --max_diff ${max_cluster_distance} --min_cluster_size ${min_cluster_size} --illumina &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> timeinfo.txt
  """
}

pandora_output_nano.concat( pandora_output_illumina ).set { pandora_output }

process make_df {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input:
  set val(type), val(min_cluster_size), val(max_cluster_distance), file(timeinfo), file(pandoralog) from pandora_output

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

  sentence=\$(grep "found mean" ${pandoralog} | tail -n1)
  stringarray=(\$sentence)
  mean_covg=\${stringarray[2]}
  echo \$mean_covg

  sentence=\$(grep "found variance" ${pandoralog} | tail -n1)
  stringarray=(\$sentence)
  var_covg=\${stringarray[2]}
  echo \$var_covg

  echo -e "${type}\t${min_cluster_size}\t${max_cluster_distance}\t\$systime\t\$usertime\t\$maxmem\t\$mean_covg\t\$var_covg" > out.tsv
  """
}

df_line.collectFile(name: final_outdir/'parameter_cluster_results.tsv').set { table }

