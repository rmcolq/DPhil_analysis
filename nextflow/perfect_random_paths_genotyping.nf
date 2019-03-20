params.pangenome_prg = ""
params.number_paths = 1
params.number_genes = 5000
params.covg = 30

params.help = false
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.final_outdir = "."
params.max_forks = 100

if (params.help){
    log.info"""
    Pipeline for running pandora and other variant callers and comparing the results.

    Usage: nextflow run compare_genotypers.nf <arguments>
    Required arguments:
      --pangenome_prg    FILE    PRG file to use as input to pandora

    Optional:
      --pipeline_root
      --final_outdir
      --max_forks

    """.stripIndent()

    exit 0
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Output directory not found: ${params.final_outdir} -- aborting"
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

        publishDir pangenome_prg.parent, mode: 'copy', overwrite: true

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

process pandora_random_paths {
  memory { 7.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
  maxRetries 1
  container {
      'shub://rmcolq/pandora:pandora'
  }

  input:
  file prg from pangenome_prg
  val num_paths from params.number_paths

  output:
  file("genome.fa") into truth_assembly

  """
  pandora random_path ${prg} ${num_paths}
  gunzip random_paths.fa.gz
  awk '!/^>/ { next } { getline seq } length(seq) >= 200 { print \$0 "\\n" seq }' random_paths.fa > tmp.random_paths_filtered.fa
  head -n\$(( 2 * ${params.number_genes} * ${params.number_paths} )) tmp.random_paths_filtered.fa > random_paths_filtered.fa
  echo ">genome" > genome.fa
  grep -v ">" random_paths_filtered.fa >> genome.fa

  """
}

process pandora_get_ref_vcf {
    memory { 5.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }

    input:
    file pangenome_prg
    file truth_assembly
    file index from pandora_idx
    file kmer_prgs from pandora_kmer_prgs

    output:
    set(file("pandora/pandora_genotyped.vcf"), file("pandora/pandora_genotyped.ref.fa")) into ref_vcf

    """
    pandora map -p ${pangenome_prg} -r ${truth_assembly} --genotype
    seqtk seq -a pandora/pandora.consensus.fq.gz | awk '{print \$1;}' > pandora/pandora_genotyped.ref.fa
    """
}

process simulate_new_ref {
    memory { 40.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
      'shub://rmcolq/Singularity_recipes:minos'
    }

    input:
    file truth_assembly
    set(file(vcf), file(ref)) from ref_vcf

    output:
    set file("simulated_vars.vcf"), file("truth.fa") into true_variants

    """
    seqtk seq -a ${truth_assembly} | awk '{print \$1;}' | cut -d "." -f1 > truth.fa
    bwa index truth.fa
    bwa mem truth.fa ${ref} > out.sam
    python3 ${params.pipeline_root}/scripts/pick_variants_for_new_ref.py  --in_vcf ${vcf} --vcf_ref ${ref} --ref truth.fa --sam out.sam --out_vcf simulated_vars.vcf --prob .05
    cat tmp.simulated_vars.vcf | grep "#" > simulated_vars.vcf
    cat tmp.simulated_vars.vcf | grep -v "#" | sort -k2 -n >> simulated_vars.vcf
    python3 ${params.pipeline_root}/scripts/filter_overlaps_in_vcf.py --vcf simulated_vars.vcf
    """
}

process simulate_nanopore_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8
  time '20m'
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }

  input:
  file truth_assembly 

  output:
  file("simulated.fa") into nanopore_reads

  """
  nanosim-h --perfect --circular -n 15000 ${truth_assembly} --unalign-rate 0 --max-len 10000
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi
  """
}

process simulate_illumina_reads {
  memory { 10.GB * task.attempt }
  errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
  maxRetries 1
  maxForks 8 
  time '3h' 
  container {
      'shub://rmcolq/Singularity_recipes:nanosimh'
  }   

  input:
  file truth_assembly
  
  output:
  file("simulated.fa") into illumina_reads
  
  """
  nanosim-h --perfect --circular -n 1000000 ${truth_assembly} --unalign-rate 0 --max-len 151 --min-len 150
  if [[ -s simulated.fa ]] ; then
  echo "simulated.fa has data."
  else
  rm simulated.fa
  exit 1
  fi
  """
}

process pandora_genotype_nanopore {
    memory { 40.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }

    publishDir final_outdir, mode: 'copy', overwrite: true

    input:
    file pangenome_prg
    file nanopore_reads
    file index from pandora_idx
    file kmer_prgs from pandora_kmer_prgs

    output:
    set(file("pandora_genotyped_nanopore_${params.covg}.vcf"), file("pandora_genotyped_nanopore_${params.covg}.ref.fa")) into pandora_nanopore_vcf
    set val("Pandora\tNanopore\t${params.covg}"), file("timeinfo.txt") into pandora_nanopore_time

    """
    echo "pandora map -p ${pangenome_prg} -r ${nanopore_reads} --genotype --max_covg ${params.covg} &> pandora.log" > command.sh
    /usr/bin/time -v bash command.sh &> timeinfo.txt
    seqtk seq -a pandora/pandora.consensus.fq.gz | awk '{print \$1;}' > pandora_genotyped_nanopore_${params.covg}.ref.fa
    cp pandora/pandora_genotyped.vcf pandora_genotyped_nanopore_${params.covg}.vcf
    """
}
process pandora_genotype_illumina {
        memory { 65.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/pandora:pandora'
        }

        publishDir final_outdir, mode: 'copy', overwrite: true

        input:
        file pangenome_prg
        file illumina_reads
        file index from pandora_idx
        file kmer_prgs from pandora_kmer_prgs

        output:
        set(file("pandora_genotyped_illumina_${params.covg}.vcf"), file("pandora_genotyped_illumina_${params.covg}.ref.fa")) into pandora_illumina_vcf
        set val("Pandora\tIllumina\t${params.covg}"), file("timeinfo.txt") into pandora_illumina_time

        """
        echo "pandora map -p ${pangenome_prg} -r ${illumina_reads} --genotype --illumina --max_covg ${params.covg} &> pandora.log" > command.sh
        /usr/bin/time -v bash command.sh &> timeinfo.txt
        seqtk seq -a pandora/pandora.consensus.fq.gz | awk '{print \$1;}' > pandora_genotyped_illumina_${params.covg}.ref.fa
        cp pandora/pandora_genotyped.vcf pandora_genotyped_illumina_${params.covg}.vcf
        """
}

pandora_nanopore_vcf.concat(pandora_illumina_vcf).set { all_vcfs }
pandora_nanopore_time.concat(pandora_illumina_time).set { all_times }

process compare_vcfs {
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        memory {1.4.GB * task.attempt}
        container {
              'shub://rmcolq/Singularity_recipes:minos'
            }

        publishDir final_outdir, mode: 'copy', overwrite: true

        input:
        set(file(truth_vcf), file(truth_vcf_ref)) from true_variants
        set(file(vcf), file(vcf_ref)) from all_vcfs

        output:
        file("*.csv") into df

        """
        python3 ${params.pipeline_root}/scripts/compare_genotypers_on_single_sample_vcf.py --truth_vcf ${truth_vcf} --truth_vcf_ref ${truth_vcf_ref} --sample_vcf ${vcf} --sample_vcf_ref ${vcf_ref} --recall_flank 14 --max_var_length 11
        """
}

process make_graph {
    errorStrategy {task.attempt < 2 ? 'retry' : 'fail'}
    maxRetries 2
    memory {0.7.GB * task.attempt}
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    publishDir final_outdir, mode: 'copy', overwrite: true

    input:
    file '*.csv' from df.collect()

    output:
    'roc*.png'

    """
    python3 ${params.pipeline_root}/scripts/plot_roc.py --x_max 0.008 --y_label "Recall"
    """
}

process timeinfo_df {
  memory { 0.1.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3

  input:
  set val(type), file(timeinfo) from all_times

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

  echo -e "${type}\t\$systime\t\$usertime\t\$maxmem" > out.tsv
  """
}

df_line.collectFile(name: final_outdir/'runinfo.tsv')
