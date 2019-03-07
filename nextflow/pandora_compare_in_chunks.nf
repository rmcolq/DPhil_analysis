params.read_index = ""
params.pangenome_prg = ""
params.chunk_size = 750
params.num_samples = 15
params.max_covg = 100

params.help = false
params.testing = false
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
        Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run pandora_compare.nf <arguments>
        Required arguments:
          --read_index              FILE    TAB-delimited file with columns for sample_name and illumina_fq
          --pangenome_prg       FILE    PRG file to use as input to pandora

        Optional:
	  --chunk_size		INT	Size of PRG chunks, default=750
          --num_samples		INT	Number of samples, scales memory requirement, default=15
          --max_covg		INT	Max covg to use, default=100
          --testing             FLAG    Run in test mode on small data which requires less memory
          --final_outdir        DIRECTORY       Where to put final output files
          --max_forks           Not used

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

if (params.read_index) {
    read_tsv = file(params.read_index).toAbsolutePath()
    if (!read_tsv.exists()) {
        exit 1, "Read TSV file not found: ${params.read_tsv} -- aborting"
    }
}   
else {
    exit 1, "Read TSV file not provided -- aborting"
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Final out directory not found: ${params.final_outdir} -- aborting"
}

Channel
    .fromPath(pangenome_prg)
    .splitText(by: 2*params.chunk_size)
    .set{ chunks_ch }

process pandora_index {
  memory { 0.8.GB * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  maxForks params.max_forks

  input: 
  file prg from chunks_ch

  output:
  set(file("${prg}"), file("${prg}.k15.w14.idx"), file("kmer_prgs")) into pandora_idx

  """
  pandora index ${prg}
  """
} 


process pandora_compare_illumina {
  memory { 0.1.GB * params.num_samples * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  maxForks params.max_forks

  input:
  set(file(prg), file(idx), file(kmer_prgs)) from pandora_idx
  file read_tsv
  
  output:
  file("pandora/pandora_multisample_genotyped.vcf") into vcfs
  file("pandora/pandora_multisample.matrix") into matrices
  file("pandora/pandora_multisample.vcf_ref.fa") into vcf_refs
  
  """
  pandora compare -p ${prg} -r ${read_tsv} --genotype --illumina --max_covg ${params.max_covg}
  if [ ! -f pandora/pandora_multisample_genotyped.vcf ]; then
      exit 1
  fi
  if [ ! -f pandora/pandora_multisample_genotyped.matrix ]; then
      exit 1
  fi
  """
} 

/*process pandora_compare_nanopore {
  memory { 0.1.GB * params.num_samples * task.attempt }
  errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
  maxRetries 3
  container {
      'shub://rmcolq/pandora:pandora'
  }
  maxForks params.max_forks

  input:
  set(file(prg), file(idx), file(kmer_prgs)) from pandora_idx
  file read_tsv

  output:
  file("pandora_multisample_genotyped.vcf") into vcfs
  file("vcf_header.txt") into vcf_headers
  file("pandora_multisample.matrix") into matrices
  file("matrix_header.txt") into matrix_headers
  file("pandora/pandora_multisample.vcf_ref.fa") into vcf_refs

  """
  pandora compare -p ${prg} -r ${read_tsv} --genotype --max_covg ${params.max_covg}
  grep -v "#" pandora/pandora_multisample_genotyped.vcf > pandora_multisample_genotyped.vcf
  grep "#" pandora/pandora_multisample_genotyped.vcf > vcf_header.txt
  tail -n+2 pandora/pandora_multisample.matrix > pandora_multisample.matrix
  head -n1 pandora/pandora_multisample.matrix > matrix_header.txt
  """
}*/

vcfs.collectFile(name: final_outdir/'pandora_multisample_genotyped.vcf', keepHeader: true, skip: 12)
vcf_refs.collectFile(name: final_outdir/'pandora_multisample.vcf_ref.fa')
matrices.collectFile(name: final_outdir/'pandora_multisample.matrix', keepHeader:true)

