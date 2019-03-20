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
 
nano_types = Channel.from("ecoli_R9_2D", "ecoli_R9_2D", "ecoli_R9_1D")
nano_lengths = Channel.from(10000, 50000)
nano_types.combine(nano_lengths).set { nano_profile }

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
  set val(type), val(max_len) from nano_profile

  output:
  file("simulated*.fa") into sim_reads_nano

  """
  nanosim-h -p ${type} -n 150000 ${ref_fasta} --unalign-rate 0 --max-len ${max_len} --circular
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi 
  mv simulated.fa simulated_nanopore_${type}_${max_len}.fa
  """
}

ill_types = Channel.from("HS25", "MSv3")
ill_lengths = Channel.from(150, 250)
ill_types.combine(ill_lengths).set { ill_profile }
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
  set val(type), val(length) from ill_profile

  output:
  file("simulated_illumina_${type}_${length}.f*") into sim_reads_illumina

  """
  art_illumina -ss ${type} -i ${ref_fasta} -l ${length} -f 300 -o simulated_illumina_${type}_${length}
  if [[ -s simulated_illumina_${type}_${length}.fq ]] ; then
  echo "simulated.fq has data."
  else
  rm simulated*.fq
  if [[ -s simulated_illumina_${type}_${length}.aln ]] ; then
  grep -v "#" simulated_illumina_${type}_${length}.aln | grep -v "@" > simulated_illumina_${type}_${length}.fa
  fi
  if [[ -s simulated_illumina_${type}_${length}.fa ]] ; then
  exit 1
  fi
  fi
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

max_covg_nano = Channel.from(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300)
sim_reads_nano.combine(max_covg_nano).set { nano_args }
process pandora_map_path_nano {
  memory { 30.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }
  
  input:
  file prg from pangenome_prg
  set file(reads), val(max_covg) from nano_args
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora.consensus.fq"), file("${reads}"), val("${max_covg}") into pandora_output_path_nano
  
  """
  pandora map -p ${prg} -r ${reads} --max_covg ${max_covg}
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi
  
  gunzip pandora/pandora.consensus.fq.gz
  """
} 

max_covg_ill = Channel.from(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300)
sim_reads_illumina.combine(max_covg_ill).set { ill_args }
process pandora_map_path_illumina {
  memory { 30.GB * task.attempt }
  errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
  maxRetries 2
  maxForks params.max_forks
  container {
      'shub://rmcolq/pandora:pandora'
  }   
  
  input:
  file prg from pangenome_prg
  set file(reads), val(max_covg) from ill_args
  file index from pandora_idx
  file kmer_prgs from pandora_kmer_prgs
  
  output:
  set file("pandora/pandora.consensus.fq"), file("${reads}"), val("${max_covg}") into pandora_output_path_illumina
  
  """
  pandora map -p ${prg} -r ${reads} --illumina --max_covg ${max_covg}
  if [[ -f pandora/pandora.consensus.fq.gz ]] ; then
  echo "pandora/pandora.consensus.fq.gz exists"
  else
  exit 1
  fi

  gunzip pandora/pandora.consensus.fq.gz
  """
} 

pandora_output_path_illumina.concat(pandora_output_path_nano).set { consensus }

process evaluate_genes_found {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks params.max_forks
  container {
      'shub://rmcolq/Singularity_recipes:minos'
  }

  input:
  set file(pandora_run), file(sim_reads), val(covg) from consensus
  file(truth_list) from truth

  output:
  file("results.tsv") into output_tsv

  """
  python3 ${params.pipeline_root}/scripts/finding_genes_eval.py --fastq ${pandora_run} --truth ${truth_list} --prefix ${sim_reads} --covg ${covg}
  """
}

output_tsv.collectFile(name: "${final_outdir}/gene_finding_by_covg.tsv").set { results }

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
  python3 ${params.pipeline_root}/scripts/plot_gene_finding.py --tsv "${tsv}"
  """
}
