params.pangenome_prg = ""
params.num_paths = 1
params.num_genes = 5000

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


ks = Channel.from(7,9,11,13,15,17,19,21,23,25,27,29,31)
ws = Channel.from(1..31)
ks.combine(ws).filter { it[1] % 2 == 1 }.filter { it[1] < it[0] }.set { kws }


process index_prg {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 5 ? 'retry' : 'ignore'}
    maxRetries 5
    container {
      'shub://rmcolq/pandora:pandora'
    }

    publishDir pangenome_prg.parent, mode: 'copy', overwrite: false

    input:
    file pangenome_prg
    set val(k), val(w) from kws

    output:
    set file("prg.fa"), file("prg.*.idx"), file("kmer_prgs"), val("${w}"), val("${k}") into indexes_nano
    set file("prg.fa"), file("prg.*.idx"), file("kmer_prgs"),val("${w}"), val("${k}") into indexes_ill
    set val("${w}"), val("${k}"), file("indextimeinfo.txt") into indexes_times

    """
    head -n10000 ${pangenome_prg} > prg.fa
    tail -n10000 ${pangenome_prg} >> prg.fa
    echo "pandora index -w ${w} -k ${k} prg.fa &> pandora.log" > command.sh
    /usr/bin/time -v bash command.sh &> indextimeinfo.txt
    """
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

  output:
  file("random_paths_filtered.fa") into random_paths

  """
  pandora random_path ${prg} ${params.num_paths}
  gunzip random_paths.fa.gz
  awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print \$0 "\\n" seq }' random_paths.fa > random_paths_filtered.fa
  """
}

process simulate_genome {
  memory { 7.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input: 
  file paths from random_paths

  output:
  file("sim_genome.fa") into nano_genome
  file("sim_genome.fa") into ill_genome
  file("truth.txt") into truth

  """
  python3 ${params.pipeline_root}/scripts/sim_genome.py --in_fa ${paths} --num_genes ${params.num_genes}
  """
}
 
process simulate_nanopore_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }

  input:
  file ref_fasta from nano_genome

  output:
  file("simulated*.fa") into sim_reads_nano

  """
  nanosim-h -p ecoli_R9_2D -n 50000 ${ref_fasta} --unalign-rate 0 --max-len 10000
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  mv simulated.fa simulated_nanopore.fa
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
  file(ref_fasta) from ill_genome

  output:
  file("simulated*.fq") into sim_reads_illumina

  """
  art_illumina -ss HS25 -i ${ref_fasta} -l 150 -f 100 -o simulated_illumina
  if [[ -s simulated*.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated*.fq
  exit 1
  fi
  """
}

process pandora_map_path_nano {
  memory { 24.GB * task.attempt }
  errorStrategy {task.attempt < 4 ? 'retry' : 'ignore'}
  maxRetries 4
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file(reads) from sim_reads_nano
  set file(prg), file(index), file(kmer_prgs), val(w), val(k) from indexes_nano
  
  output:
  set val("Nanopore"), file("pandora_result.fq"), val("${w}"), val("${k}") into pandora_output_path_nano
  set val("Nanopore"), val("${w}"), val("${k}"), file("maptimeinfo.txt") into pandora_output_time_nano
  
  """
  echo "pandora map -p ${prg} -r ${reads} --max_covg 100 &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> maptimeinfo.txt

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
  memory { 24.GB * task.attempt }
  errorStrategy {task.attempt < 4 ? 'retry' : 'ignore'}
  maxRetries 4
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }   
  
  input:
  file(reads) from sim_reads_illumina
  set file(prg), file(index), file(kmer_prgs), val(k), val(k) from indexes_ill
  
  output:
  set val("Illumina"), file("pandora_result.fq"), val("${w}"), val("${k}") into pandora_output_path_illumina
  set val("Illumina"), val("${w}"), val("${k}"), file("maptimeinfo.txt") into pandora_output_time_illumina
  
  """
  echo "pandora map -p ${prg} -r ${reads} --illumina --max_covg 100 &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> maptimeinfo.txt

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

pandora_output_path_illumina.concat(pandora_output_path_nano).set { consensus }
pandora_output_time_illumina.concat(pandora_output_time_nano).set { map_times }

process evaluate_genes_found {
  memory { 0.01.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set val(type), file(pandora_run), val(w), val(k) from consensus
  file(truth_list) from truth

  output:
  file("results.tsv") into output_tsv

  """
  python3 ${params.pipeline_root}/scripts/finding_genes_eval.py --fastq ${pandora_run} --truth ${truth_list} --prefix "${type}\t${w}\t${k}" --covg 100
  """
}

output_tsv.collectFile(name: 'gene_finding_by_covg.tsv')

process make_df_index_times {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input:
  set val(w), val(k), file(timeinfo) from indexes_times

  output:
  file("out.tsv") into df_line_indexes

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

  echo -e "${w}\t${k}\t\$systime\t\$usertime\t\$maxmem" > out.tsv
  """
}

df_line_indexes.collectFile(name: final_outdir/'index_parameters.tsv')

process make_df_map_times {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input: 
  set val(type), val(w), val(k), file(timeinfo) from map_times
  
  output:
  file("out.tsv") into df_line_map
  
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
  
  echo -e "${type}\t${w}\t${k}\t\$systime\t\$usertime\t\$maxmem" > out.tsv
  """
} 

df_line_map.collectFile(name: final_outdir/'map_parameters.tsv')
