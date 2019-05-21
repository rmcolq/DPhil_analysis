params.reference_directory = ""
params.illumina_reads_1 = ""
params.illumina_reads_2 = ""
params.nanopore_reads = ""
params.raw_fast5s = ""
params.albacore_summary = ""

params.help = false
params.testing = false
params.pipeline_root = ""
params.final_outdir = "."
params.max_forks = 10

if (params.help){
    log.info"""
	Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run compare_genotypers.nf <arguments>
        Required arguments:
          --reference_directory	DIRECTORY	Directory containing different references to use when calling with e.g. snippy

	At least one required:
	  --nanopore_reads	FILE	Fastaq[gz] file of nanopore reads to make calls from
	  --illumina_reads_1	FILE	Fastaq[gz] file of illumina reads to make calls from
	  --illumina_reads_2	FILE	Fastaq[gz] file of illumina reads to make calls from if paired

	Optional:
	  --raw_fast5s	DIR	Directory where raw fast5 nanopore files are (for use by nanopolish)
	  --albacore_summary	FILE	The sequencing_summary.txt file output by albacore during basecalling (speeds up nanopolish)

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

if (params.reference_directory) {
    reference_assemblies_in = Channel.fromPath("${params.reference_directory}/*.{fa,fa.gz,fasta,fasta.gz,fna,fna.gz}")
}
else {
    exit 1, "Reference assembly directory not provided -- aborting"
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

process unzip_reference_assembly {
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    file reference_assembly from reference_assemblies_in

    output:
    file "*.{fa,fna,fasta}" into reference_assemblies_snippy
    file "*.{fa,fna,fasta}" into reference_assemblies_nanopolish

    """
    v=${reference_assembly}
    if [ \${v: -3} == ".gz" ]
    then
        zcat ${reference_assembly} | awk '{print \$1;}' > \${v::-3}
    else
        cat ${reference_assembly} | awk '{print \$1;}' > n.\$v
    fi
    """
}

if (params.nanopore_reads) {

    process nanopolish_index {
	container {
          'shub://rmcolq/Singularity_recipes:nanopolish'
        }
        input:
        file nanopore_reads
        file raw_fast5s
        file albacore_summary

        output:
        set file("${nanopore_reads}.index*") into nanopolish_index

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

	publishDir final_outdir, mode: 'copy', overwrite: true
        
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

        echo "Nanopolish run successful"

        wait(10)

        v=\$(head -n1 ${reference_assembly})
        ref_id=\${v:1:\${#v}}
        echo "Found id \$ref_id"

        cp \$(ls nanopolish.results/vcf/nanopolish.*.vcf | head -n1) nanopolish_\$ref_id.vcf
        for f in \$(ls nanopolish.results/vcf/nanopolish.*.vcf | tail -n+2)
        do
        cat \$f | grep -v "#" >> nanopolish_\$ref_id.vcf
        done
        echo "Combined output files"

        cp ${reference_assembly} nanopolish_\$ref_id.ref.fa
        echo "Copied ref assembly"

        echo "Finish"
        """
    }
}

if (params.illumina_reads_1 && params.illumina_reads_2) {

    process snippy_genotype_pe_illumina {
        memory { 24.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/Singularity_recipes:snippy'
        }
	cpus 8
        maxForks 5

	publishDir final_outdir, mode: 'copy', overwrite: true

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

    process snippy_genotype_se_illumina {
        memory { 24.GB * task.attempt } 
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
        container {
          'shub://rmcolq/Singularity_recipes:snippy'
        }
	cpus 8

	publishDir final_outdir, mode: 'copy', overwrite: true
        
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

