params.truth_assembly = ""
params.illumina_reads_1 = ""
params.illumina_reads_2 = ""
params.nanopore_reads = ""
params.pangenome_prg = ""
params.raw_fast5s = ""
params.albacore_summary = ""

params.help = false
params.testing = false
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
    Pipeline for running pandora and other variant callers and comparing the results.

    Usage: nextflow run compare_genotypers.nf <arguments>
    Required arguments:
      --truth_assembly    FILE    Assembly true sequence for reads
      --pangenome_prg    FILE    PRG file to use as input to pandora
      --fasta_for_ref    FILE    Fasta of annotated sequences to use as reference paths through PRG

    At least one required:
      --nanopore_reads    FILE    Fastaq[gz] file of nanopore reads to make calls from
      --illumina_reads_1    FILE    Fastaq[gz] file of illumina reads to make calls from
      --illumina_reads_2    FILE    Fastaq[gz] file of illumina reads to make calls from if paired

    Optional:
      --raw_fast5s        DIR    Directory where raw fast5 nanopore files are (for use by nanopolish)
      --albacore_summary    FILE    The sequencing_summary.txt file output by albacore during basecalling (speeds up nanopolish)
      --mask        FILE    Mask of regions of truth_assembly to ignore

      --testing
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

if (params.truth_assembly) {
    truth_assembly = file(params.truth_assembly).toAbsolutePath()
    if (!truth_assembly.exists()) {
        exit 1, "Truth assembly file not found: ${params.truth_assembly} -- aborting"
    }
}
else {
    exit 1, "Truth assembly file not provided -- aborting"
}

if (params.fasta_for_ref) {
    fasta_for_ref = file(params.fasta_for_ref).toAbsolutePath()
    if (!fasta_for_ref.exists()) {
        exit 1, "Fasta file of annotated sequences for reference paths not found: ${params.fasta_for_ref} -- aborting"
    }
}       
else {
    exit 1, "Fasta file of annotated sequences for reference paths not provided -- aborting"
}

if (!params.nanopore_reads && !params.illumina_reads_2 && !params.illumina_reads_1) {
    exit 1, "No read file provided -- aborting"
}
if (params.illumina_reads_1) {
    illumina_reads_1 = file(params.illumina_reads_1).toAbsolutePath()
    if (!illumina_reads_1.exists()) {
        exit 1, "Illumina reads file not found: ${params.illumina_reads_1} -- aborting"
    }
}
if (params.illumina_reads_2) {
    illumina_reads_2 = file(params.illumina_reads_2).toAbsolutePath()
    if (!illumina_reads_2.exists()) {
        exit 1, "Illumina reads file not found: ${params.illumina_reads_2} -- aborting"
    }
}
if (params.nanopore_reads) {
    nanopore_reads = file(params.nanopore_reads).toAbsolutePath()
    if (!nanopore_reads.exists()) {
        exit 1, "Nanopore reads file not found: ${params.nanopore_reads} -- aborting"
    }
    if (params.raw_fast5s) {
        raw_fast5s = file(params.raw_fast5s).toAbsolutePath()
        if (!raw_fast5s.exists()) {
            exit 1, "Nanopore raw fast5 reads directory not found: ${params.raw_fast5s} -- aborting"
        }   
    }
}

if (params.albacore_summary) {
        albacore_summary = file(params.albacore_summary).toAbsolutePath()
        if (!albacore_summary.exists()) {
            exit 1, "Albacore sequencing_summary.txt file not found: ${params.albacore_summary} -- aborting"
        }
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

process pandora_get_ref_vcf {
    memory { 40.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
    maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }

    input:
    file pangenome_prg
    file fasta_for_ref
    file truth_assembly
    file index from pandora_idx
    file kmer_prgs from pandora_kmer_prgs

    output:
    set(file("pandora/pandora_genotyped.vcf"), file("${pangenome_prg}.vcf_ref.fa")) into ref_vcf

    """
    pandora get_vcf_ref ${pangenome_prg} ${fasta_for_ref}
    pandora map -p ${pangenome_prg} -r ${truth_assembly} --vcf_refs ${pangenome_prg}.vcf_ref.fa.gz --genotype
    gunzip ${pangenome_prg}.vcf_ref.fa.gz
    """
}

process simulate_new_ref {
    memory { 40.GB * task.attempt }
    errorStrategy {task.attempt < 3 ? 'retry' : 'fail'}
    maxRetries 3
    container {
      'shub://rmcolq/Singularity_recipes:minos'
    }

    input:
    file truth_assembly
    set(file(vcf), file(ref)) from ref_vcf

    output:
    file("simulated_ref.fa") into reference_assemblies_snippy
    file("simulated_ref.fa") into reference_assemblies_nanopolish
    set(file("simulated_vars.vcf"), file("${ref}")) into true_variants

    """
    v=${truth_assembly}
    if [ \${v: -3} == ".gz" ]
    then
      zcat \$v | awk '{print \$1;}' > \${v::-3}
      v=\${v::-3}
    else
      cat \$v | awk '{print \$1;}' > n.\$v
      mv n.\$v \$v
    fi
    bwa index \$v
    bwa mem \$v ${ref} > out.sam
    python3 ${params.pipeline_root}/scripts/pick_variants_for_new_ref.py  --in_vcf ${vcf} --vcf_ref ${ref} --ref \$v --sam out.sam --out_vcf simulated_vars.vcf
    bcftools sort --max_mem 40G simulated_vars.vcf
    bgzip simulated_vars.vcf
    tabix -p vcf simulated_vars.vcf.gz
    cat \$v | vcf-consensus simulated_vars.vcf.gz > simulated_ref.fa
    """
}


/*if (params.nanopore_reads) {
    process pandora_genotype_nanopore {
        memory { 42.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    container {
          'shub://rmcolq/pandora:pandora'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

        input:
        file nanopore_reads
        file pangenome_prg
        file pandora_idx
        file pandora_kmer_prgs

        output:
        set(file("pandora_genotyped_full.vcf"), file("pandora_genotyped_full.ref.fa")) into pandora_full_vcf

        """
    which pandora
        pandora map -p ${pangenome_prg} -r ${nanopore_reads} -w 14 -k 15 --genotype --outdir pandora
    mv pandora/pandora_genotyped.vcf pandora_genotyped_full.vcf
    gunzip -c pandora/pandora.consensus.fq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > pandora_genotyped_full.ref.fa

        """
    }

    process pandora_genotype_30X_nanopore {
        memory { 42.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    container {
      'shub://rmcolq/pandora:pandora'
    }    

        publishDir final_outdir, mode: 'copy', overwrite: false

        input:
        file nanopore_reads
        file pangenome_prg
        file pandora_idx
        file pandora_kmer_prgs

        output:
    set(file("pandora_genotyped_30X.vcf"), file("pandora_genotyped_30X.ref.fa")) into pandora_30X_vcf

        """
    which pandora
        pandora map -p ${pangenome_prg} -r ${nanopore_reads} -w 14 -k 15 --genotype --max_covg 30 --outdir pandora
        mv pandora/pandora_genotyped.vcf pandora_genotyped_30X.vcf
        gunzip -c pandora/pandora.consensus.fq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > pandora_genotyped_30X.ref.fa
        """
    }

    process nanopolish_index {
    container {
          'shub://rmcolq/Singularity_recipes:nanopolish'
        }
        input:
        file nanopore_reads
        file raw_fast5s
        file albacore_summary

        output:
        set(file("${nanopore_reads}.index*")) into nanopolish_index

        """
        nanopolish index -d ${raw_fast5s} -s ${albacore_summary} ${nanopore_reads}

        """
    }

    process nanopolish_genotype_nanopore {
        memory { 64.GB * task.attempt }
        errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
        maxRetries 1
        container {
          'shub://rmcolq/Singularity_recipes:nanopolish'
        }
        cpus 16
    maxForks 5

    publishDir final_outdir, mode: 'copy', overwrite: false
        
        input:
        file nanopore_reads 
    file reference_assembly from reference_assemblies_nanopolish
        file nanopolish_index
        file raw_fast5s
        
        output:
    set(file("nanopolish_*.vcf"), file("nanopolish_*.ref.fa")) into nanopolish_vcf
        
    """
        minimap2 -ax map-ont -t 8 ${reference_assembly} ${nanopore_reads} | samtools sort -o reads.sorted.bam -T reads.tmp
        samtools index reads.sorted.bam
    mkdir -p nanopolish.results/vcf
        python3 /nanopolish/scripts/nanopolish_makerange.py ${reference_assembly} | parallel --results nanopolish.results -P 2 \
    nanopolish variants \
      -t 8 \
      -w {1} \
      --reads ${nanopore_reads} \
      --bam reads.sorted.bam \
          --genome ${reference_assembly} \
      -o nanopolish.results/vcf/nanopolish.{1}.vcf \
      -q dam,dcm \
          --ploidy 1

    v=\$(head -n1 ${reference_assembly})
        ref_id=\${v:1:\${#v}}

    cp \$(ls nanopolish.results/vcf/nanopolish.*.vcf | head -n1) nanopolish_\$ref_id.vcf
    for f in \$(ls nanopolish.results/vcf/nanopolish.*.vcf | tail -n+2)
    do
    cat \$f | grep -v "#" >> nanopolish_\$ref_id.vcf
    done
    cp ${reference_assembly} nanopolish_\$ref_id.ref.fa
        """
    }
}

if (params.illumina_reads_1 && params.illumina_reads_2) {
    process combine_illumina_reads {
    input:
        file illumina_reads_1
        file illumina_reads_2
        
    output:
    file ${illumina_reads_1} into illumina_reads
    
    """
    cat ${illumina_reads_2} >> ${illumina_reads_1}
    """
    }

    process snippy_genotype_pe_illumina {
        memory { 24.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/Singularity_recipes:snippy'
        }
    cpus 8
        maxForks 5

    publishDir final_outdir, mode: 'copy', overwrite: false

        input:
        file reference_assembly from reference_assemblies_snippy
        file illumina_reads_1
        file illumina_reads_2

        output:
        set(file("snippy_*.vcf"), file("snippy_*.ref.fa")) into snippy_vcf

        """
        snippy --cpus 8 --outdir snippy_outdir --reference ${reference_assembly} --pe1 ${illumina_reads_1} --pe2 ${illumina_reads_2}

    v=\$(head -n1 ${reference_assembly})
    ref_id=\${v:1:\${#v}}
        cp snippy_outdir/snps.filt.vcf snippy_\$ref_id.vcf
    cp ${reference_assembly} snippy_\$ref_id.ref.fa
        """
    }
}
else if (params.illumina_reads_1) {
    illumina_reads = illumina_reads_1

    process snippy_genotype_se_illumina {
        memory { 24.GB * task.attempt } 
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/Singularity_recipes:snippy'
        }
    cpus 8

    publishDir final_outdir, mode: 'copy', overwrite: false
        
        input:
        file reference_assembly from reference_assemblies_snippy
        file illumina_reads_1
        
        output:
        set(file("snippy_*.vcf"), file("snippy_*.ref.fa")) into snippy_vcf
        
    """
        snippy --cpus 8 --outdir snippy_outdir --reference ${reference_assembly} --se ${illumina_reads_1}

    v=\$(head -n1 ${reference_assembly})
        ref_id=\${v:1:\${#v}}
        cp snippy_outdir/snps.filt.vcf snippy_\$ref_id.vcf
        cp ${reference_assembly} snippy_\$ref_id.ref.fa
        """
    }
}

if (params.illumina_reads_1) {
    process pandora_genotype_illumina {
        memory { 64.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
    container {
          'shub://rmcolq/pandora:pandora'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

        input:
        file illumina_reads
        file pangenome_prg
        file pandora_idx
    file pandora_kmer_prgs

        output:
    set(file("pandora_genotyped_illumina.vcf"), file("pandora_genotyped_illumina.ref.fa")) into pandora_illumina_vcf

        """
    
        pandora map -p ${pangenome_prg} -r ${illumina_reads} -w 14 -k 15 --genotype --illumina --outdir pandora
    mv pandora/pandora_genotyped.vcf pandora_genotyped_illumina.vcf
        gunzip -c pandora/pandora.consensus.fq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > pandora_genotyped_illumina.ref.fa
        """
    }
}

if (params.illumina_reads_1) {
    pandora_illumina_vcf.concat(snippy_vcf).set { illumina_channels }
} else {
    illumina_channels = Channel.from()
}
if (params.nanopore_reads) {
    pandora_full_vcf.concat(pandora_30X_vcf, nanopolish_vcf).set { nanopore_channels }
} else {
    nanopore_channels = Channel.from()
}
nanopore_channels.concat(illumina_channels).set { all_vcfs }

if (params.mask) {
    process compare_vcfs {
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        memory {1.4.GB * task.attempt}
        container {
              'shub://rmcolq/Singularity_recipes:minos'
            }

        publishDir final_outdir, mode: 'copy', overwrite: false

        input:
        set(file(truth_vcf), file(truth_vcf_ref)) from true_variants
        set(file(vcf), file(vcf_ref)) from all_vcfs
        file mask

        output:
        file('*.png') into output_graphs

        """
        python3 ${params.pipeline_root}/scripts/compare_genotypers_on_sample.py --truth_vcf ${truth_vcf} --truth_vcf_ref ${truth_vcf_ref} --directory ${final_outdir} --mask ${mask} --flank 7
        """
    }
}
else {
    process compare_vcfs {
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        memory {1.4.GB * task.attempt}
        container {
              'shub://rmcolq/Singularity_recipes:minos'
            }

        publishDir final_outdir, mode: 'copy', overwrite: false
    
        input:
        set(file(truth_vcf), file(truth_vcf_ref)) from true_variants
        set(file(vcf), file(vcf_ref)) from all_vcfs
              
        output:
        file('*.png') into output_graphs
        
        """
        python3 ${params.pipeline_root}/scripts/compare_genotypers_on_sample.py --truth_vcf ${truth_vcf} --truth_vcf_ref ${truth_vcf_ref} --directory ${final_outdir} --flank 7
        """
    }   
}       
*/
