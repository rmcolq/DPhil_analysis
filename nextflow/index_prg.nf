params.pangenome_prg = ""

params.help = false
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.final_outdir = "."
params.max_forks = 50
params.chunk_size = 4000
params.num_prg = 0
params.w = 14
params.k = 15
params.g_offset = 0

if (params.help){
    log.info"""
    Pipeline to index a prg with pandora in parallel.

    Usage: nextflow run index_prg.nf <arguments>
    Required arguments:
      --pangenome_prg    FILE    PRG file to use as input to pandora
      --num_prg		 INT	 Number of PRGs in file

    Optional:
      --pipeline_root
      --final_outdir
      --max_forks	INT	Number of consecutive prgs to index
      --chunk_size	INT	Chunks of PRG file to work on in parallel, needs to be multiple of 4000
      --w		INT
      --k		INT

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

max_num = ((params.num_prg/params.chunk_size).toBigInteger().intValueExact())
nums = Channel.from( 0..max_num ).map { params.chunk_size * it }

process pandora_index {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    }
    maxForks params.max_forks
    time '7d'

    input:
    file pangenome_prg
    val w from params.w
    val k from params.k
    val offset from nums
    val num_prg from params.chunk_size
    val g_offset from params.g_offset

    output:
    file("prg.fa.*.idx") into indexes
    file("kmer_prgs") into kd1
    file("kmer_prgs/*") into kd2
    file("kmer_prgs/*") into kd3

    """
    a=\$(( 2 * ${offset} + 1 ))
    b=\$(( 2 * ${num_prg} ))
    tail -n+\$a ${pangenome_prg} | head -n \$b > prg.fa
    o=\$(( ${offset} + ${g_offset} ))
    pandora index prg.fa --offset \$o -w ${w} -k ${k}
    """
}

process combine {
    memory { 20.GB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    }

    publishDir pangenome_prg.parent, mode: 'copy', overwrite: true

    input:
    file pangenome_prg
    val w from params.w
    val k from params.k
    file("prg.fa.*.idx") from indexes.collect()

    output:
    file("${pangenome_prg}.k*.w*.idx")

    """
    pandora merge_index --outfile ${pangenome_prg}.k${k}.w${w}.idx \$(ls *.idx)
    """
}

process publish_kd1 {
    memory { 1.MB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    } 
    maxForks 1

    publishDir pangenome_prg.parent, mode: 'copy', overwrite: false
    
    input:
    file pangenome_prg
    file(kp) from kd1

    output:
    file("${kp}")
    val("wait") into out_kd1

    """
    """
}

process publish_kd2 {
    memory { 1.MB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    } 
    maxForks 1
    
    publishDir pangenome_prg.parent+"/kmer_prgs", mode: 'copy', overwrite: false
    
    input:
    file pangenome_prg
    file(kp) from kd2
    val(wait) from out_kd1
    
    output:
    file("${kp}")
    val("wait") into out_kd2
    
    """
    """
}

process publish_kgs {
    memory { 1.MB * task.attempt }
    errorStrategy {task.attempt < 1 ? 'retry' : 'fail'}
    maxRetries 1
    container {
      'shub://rmcolq/pandora:pandora'
    }
    maxForks params.max_forks

    publishDir pangenome_prg.parent+"/kmer_prgs/", mode: 'copy', overwrite: true

    input:
    file pangenome_prg
    file(kp) from kd3
    val(wait) from out_kd2

    output:
    file("${kp}/*")

    """
    """
}
