params.pangenome_prg = ""
params.index_dir = ""
params.num_paths = 1
params.num_genes = 5000
params.covg = 50

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
if (params.index_dir) {
    index_dir = file(params.index_dir).toAbsolutePath()
    if (!index_dir.exists()) {
        exit 1, "Index directory not found: ${params.index_dir} -- aborting"
    }
}
else {
    index_dir = pangenome_prg.parent
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}


ks = Channel.from(7,11,15,19,23,27,31)
ws = Channel.from(1..31)
ks.combine(ws).filter { it[1] % 4 == 2 }.filter { it[1] < it[0] }.set { kws }


/*process index_prg {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 5 ? 'retry' : 'ignore'}
    maxRetries 5
    container {
      'shub://rmcolq/pandora:pandora'
    }

    input:
    file pangenome_prg
    set val(k), val(w) from kws

    output:
    set file("prg.fa"), file("prg.*.idx"), file("kmer_prgs"), val("${w}"), val("${k}") into indexes_nano
    set file("prg.fa"), file("prg.*.idx"), file("kmer_prgs"),val("${w}"), val("${k}") into indexes_ill
    set val("${w}"), val("${k}"), file("indextimeinfo.txt") into indexes_times
    file("${pangenome_prg}") into local_prg

    """
    echo "pandora index -w ${w} -k ${k} ${pangenome_prg} &> pandora.log" > command.sh
    /usr/bin/time -v bash command.sh &> indextimeinfo.txt
    """
}*/

process find_prg_index {
memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    } 
    
    input:
    file index_dir
    set val(k), val(w) from kws
    
    output:
    set file("${index_dir}/k${k}.w${w}/*prg.fa"), file("${index_dir}/k${k}.w${w}/*prg.*.idx"), file("${index_dir}/k${k}.w${w}/kmer_prgs"), val("${w}"), val("${k}") into indexes_nano
    set file("${index_dir}/k${k}.w${w}/*prg.fa"), file("${index_dir}/k${k}.w${w}/*prg.*.idx"), file("${index_dir}/k${k}.w${w}/kmer_prgs"), val("${w}"), val("${k}") into indexes_ill
    """
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
 
types = Channel.from("ecoli_R9_2D", "ecoli_R9_1D")
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
  val(type) from types

  output:
  file("simulated*.fa") into sim_reads_nano

  """
  nanosim-h -p ${type} -n 50000 ${ref_fasta} --unalign-rate 0 --max-len 10000 --circular
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  mv simulated.fa simulated_nanopore_${type}.fa
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
  file("simulated*.f*") into sim_reads_illumina

  """
  art_illumina -ss HS25 -i ${ref_fasta} -l 150 -f 100 -o simulated_illumina_HS25
  if [[ -s simulated_illumina_HS25.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated*.fq
  if [[ -s simulated_illumina_HS25.aln ]] ; then
  grep -v "#" simulated_illumina_HS25.aln | grep -v "@" > simulated_illumina_HS25.fa
  fi
  if [[ -s simulated_illumina_HS25.fa ]] ; then
  exit 1
  fi
  fi
  """
}

process pandora_map_path_nano {
  memory { 50.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file(reads) from sim_reads_nano
  set file(prg), file(index), file(kmer_prgs), val(w), val(k) from indexes_nano
  
  output:
  set val("Nanopore"), file("pandora/pandora.consensus.fq"), val("${w}"), val("${k}") into pandora_output_path_nano
  set val("Nanopore"), val("${w}"), val("${k}"), file("maptimeinfo.txt") into pandora_output_time_nano
  
  """
  echo "pandora map -p ${prg} -r ${reads} --max_covg ${params.covg} -w ${w} -k ${k} &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> maptimeinfo.txt

  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  gunzip pandora/pandora.consensus.fq.gz 
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
  set file(prg), file(index), file(kmer_prgs), val(w), val(k) from indexes_ill
  
  output:
  set val("Illumina"), file("pandora/pandora.consensus.fq"), val("${w}"), val("${k}") into pandora_output_path_illumina
  set val("Illumina"), val("${w}"), val("${k}"), file("maptimeinfo.txt") into pandora_output_time_illumina
  
  """
  echo "pandora map -p ${prg} -r ${reads} --illumina --max_covg ${params.covg} -w ${w} -k ${k} &> pandora.log" > command.sh
  /usr/bin/time -v bash command.sh &> maptimeinfo.txt

  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  gunzip pandora/pandora.consensus.fq.gz
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
  python3 ${params.pipeline_root}/scripts/finding_genes_eval.py --fastq ${pandora_run} --truth ${truth_list} --prefix "${type}\t${w}\t${k}" --covg ${params.covg}
  """
}

output_tsv.collectFile(name: "${final_outdir}/gene_finding_by_wk_params.tsv").set { results }

process make_plot {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }   
  publishDir final_outdir, mode: 'copy', overwrite: true
  
  input:
  file(tsv) from results
  
  output:
  file("*.png") into output_plot
  
  """
  python3 ${params.pipeline_root}/scripts/plot_gene_finding.py --tsv "${tsv}" --p1 'w' --p2 'k' --l1 "w" --l2 "k" --d1 14 --d2 15
  """
} 

/*process make_df_index_times {
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
*/

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

