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

illumina_reads = Channel.from(read_index).splitCsv(header: false, sep:'\t')

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
    set(file("pandora_genotyped_illumina.vcf"), file("pandora_genotyped_illumina.ref.fa")) into pandora_illumina_vcf

    """
    pandora map -p ${pangenome_prg} -r ${illumina_reads} --genotype --illumina
    seqtk seq -a pandora/pandora.consensus.fq.gz | awk '{print \$1;}' > pandora_genotyped_illumina.ref.fa
    cp pandora/pandora_genotyped.vcf pandora_genotyped_illumina.vcf
    """
}
