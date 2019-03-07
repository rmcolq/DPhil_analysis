params.read_index = ""
params.pangenome_prg = ""

params.help = false
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
    Pipeline for running pandora and other variant callers and comparing the results.

    Usage: nextflow run compare_genotypers.nf <arguments>
    Required arguments:
      --read_index    	FILE    Index of read files
      --pangenome_prg   FILE    PRG file to use as input to pandora

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

if (params.read_index) {
    read_index = file(params.read_index).toAbsolutePath()
    if (!read_index.exists()) {
        exit 1, "Read index file not found: ${params.read_index} -- aborting"
    }
}
else {
    exit 1, "Read index file not provided -- aborting"
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

process pandora_genotype_illumina {
    memory { 1.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }

    input:
    file pangenome_prg
    file read_index

    output:
    

    """
    nextflow run ${params.pipeline_root}/pandora_compare_in_chunks.nf --pangenome_prg ${pangenome_prg} --tsv ${read_index} --num_samples 150 --max_covg 50 --max_forks params.max_forks -c ${params.pipeline_root}/nextflow/nextflow.config
    cat pandora_multisample_genotyped.vcf >> vcf_header.txt
    mv vcf_header.txt pandora_multisample_genotyped.vcf
    cat pandora_multisample.matrix >> matrix_header.txt
    mv matrix_header.txt pandora_multisample.matrix
    """
}
