params.reference_directory = ""
params.illumina_reads_1 = ""
params.illumina_reads_2 = ""
params.nanopore_reads = ""
params.pangenome_prg = ""
params.vcf_ref = ""
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
          --pangenome_prg	FILE	PRG file to use as input to pandora

	At least one required:
	  --nanopore_reads	FILE	Fastaq[gz] file of nanopore reads to make calls from
	  --illumina_reads_1	FILE	Fastaq[gz] file of illumina reads to make calls from
	  --illumina_reads_2	FILE	Fastaq[gz] file of illumina reads to make calls from if paired

	Optional:
	  --vcf_ref	FILE	Fastaq[gz] file specifying reference sequence for each PRG (sequences must exist in PRG)
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

if (!params.vcf_ref) {
    process pandora_get_vcf_ref {
	memory 0.7.GB
	container {
          'shub://rmcolq/pandora:pandora'
        }

	input:
	file pangenome_prg

	output:
        file "${pangenome_prg}.vcf_ref.fa.gz" into vcf_ref

	"""
	pandora get_vcf_ref ${pangenome_prg}
	"""
    }
}
else {
    vcf_ref = file(params.vcf_ref).toAbsolutePath()
    if (!vcf_ref.exists()) {
        exit 1, "VCF ref file not found: ${params.vcf_ref} -- aborting"
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
    /*process pandora_genotype_nanopore {
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
	file vcf_ref

    	output:
    	set(file("pandora_genotyped_full.vcf"), file("pandora_genotyped_full.ref.fa"), file("pandora_consensus_full.fq.gz")) into pandora_full_vcf

    	"""
	which pandora
    	pandora map -p ${pangenome_prg} -r ${nanopore_reads} -w 14 -k 15 --vcf_refs ${vcf_ref} --genotype --outdir pandora
	mv pandora/pandora.consensus.fq.gz pandora_consensus_full.fq.gz
	mv pandora/pandora_genotyped.vcf pandora_genotyped_full.vcf

	v=${vcf_ref}
	if [ \${v: -3} == ".gz" ]
        then
        zcat ${vcf_ref} | awk '{print \$1;}' > pandora_genotyped_full.ref.fa
	else
	mv ${vcf_ref} pandora_genotyped_full.ref.fa
	fi

    	"""
    }

    process pandora_genotype_30X_nanopore {
        memory { 42.GB * task.attempt }
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
	container {
	  'shub://rmcolq/pandora:pandora'
	}	

        publishDir final_outdir, mode: 'copy', overwrite: true

        input:
        file nanopore_reads
        file pangenome_prg
        file pandora_idx
        file pandora_kmer_prgs
        file vcf_ref

        output:
	set(file("pandora_genotyped_30X.vcf"), file("pandora_genotyped_30X.ref.fa"), file("pandora_consensus_30X.fq.gz")) into pandora_30X_vcf

        """
	which pandora
        pandora map -p ${pangenome_prg} -r ${nanopore_reads} -w 14 -k 15 --vcf_refs ${vcf_ref} --genotype --max_covg 30 --outdir pandora
	mv pandora/pandora.consensus.fq.gz pandora_consensus_30X.fq.gz
	mv pandora/pandora_genotyped.vcf pandora_genotyped_30X.vcf

	v=${vcf_ref}
        if [ \${v: -3} == ".gz" ]
        then
        zcat ${vcf_ref} | awk '{print \$1;}' > pandora_genotyped_30X.ref.fa
        else
        mv ${vcf_ref} pandora_genotyped_30X.ref.fa
        fi
        """
    }*/

/*    process canu_assembly_nanopore {
	memory { 68.GB }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
	container {
	  'shub://rmcolq/Singularity_recipes:canu'
	}
	cpus 16

	input:
        file nanopore_reads

	output:
	file "canu_assembly.fa" into canu_assembly

	"""
	which canu
        canu -p assembly -d canu_dir genomeSize=5m -nanopore-raw ${nanopore_reads} useGrid=false maxThreads=16 maxMemory=64g
	mv canu_dir/assembly.contigs.fasta canu_assembly.fa
        """
    }*/

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

/*    process racon_polish_assembly {
        memory { 10.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        container {
          'shub://rmcolq/Singularity_recipes:racon'
        }
        cpus 8

        input:
	file nanopore_reads
        file draft_assembly from canu_assembly

        output:
        file "racon_assembly.fa" into racon_assembly

        """
	minimap2 -t 8 ${draft_assembly} ${nanopore_reads} > draft.minimap1.paf
	racon -t 8 ${nanopore_reads} draft.minimap1.paf ${draft_assembly} > draft.racon1.fasta

	minimap2 -t 8 draft.racon1.fasta ${nanopore_reads} > draft.minimap2.paf
        racon -t 8 ${nanopore_reads} draft.minimap2.paf draft.racon1.fasta > racon_assembly.fa
        """
    }

    process nanopolish_polish_assembly {
        memory { 24.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        container {
          'shub://rmcolq/Singularity_recipes:nanopolish'
        }
        cpus 8

        input:
        file nanopore_reads
        file draft_assembly from racon_assembly
        file nanopolish_index
        file raw_fast5s

        output:
        file "nanopolish_assembly.fa" into nanopolish_assembly

        """
	minimap2 -ax map-ont -t 8 ${draft_assembly} ${nanopore_reads} | samtools sort -o reads.sorted.bam -T reads.tmp
        samtools index reads.sorted.bam
	mkdir -p nanopolish.results/consensus
        python3 /nanopolish/scripts/nanopolish_makerange.py ${draft_assembly} | parallel --results nanopolish.results -P 8 \
	  nanopolish variants --consensus \
	  -o nanopolish.results/consensus/nanopolish.{1}.vcf \
	  -w {1} \
	  --reads ${nanopore_reads} \
	  --bam reads.sorted.bam \
	  --genome ${draft_assembly} \
	  --max-haplotypes=2000 \
	  -q dam,dcm
        nanopolish vcf2fasta -g ${draft_assembly} nanopolish.results/consensus/nanopolish.*.vcf > nanopolish_assembly.fa

        """
    }*/

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

	v=\$(head -n1 ${reference_assembly})
        ref_id=\${v:1:\${#v}}
        echo "Found id $ref_id"

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
/*
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
    illumina_reads = illumina_reads_1

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

if (params.illumina_reads_1) {
    process pandora_genotype_illumina {
        memory { 64.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
	container {
          'shub://rmcolq/pandora:pandora'
        }

	publishDir final_outdir, mode: 'copy', overwrite: true

        input:
        file illumina_reads
        file pangenome_prg
        file pandora_idx
	file pandora_kmer_prgs
        file vcf_ref

        output:
	set(file("pandora_genotyped_illumina.vcf"), file("pandora_genotyped_illumina.ref.fa"), file("pandora_consensus_illumina.fq.gz")) into pandora_illumina_vcf

        """
	
        pandora map -p ${pangenome_prg} -r ${illumina_reads} -w 14 -k 15 --vcf_refs ${vcf_ref} --genotype --illumina --outdir pandora
        mv pandora/pandora.consensus.fq.gz pandora_consensus_illumina.fq.gz
	mv pandora/pandora_genotyped.vcf pandora_genotyped_illumina.vcf

	v=${vcf_ref}
        if [ \${v: -3} == ".gz" ]
        then
        zcat ${vcf_ref} | awk '{print \$1;}' > pandora_genotyped_illumina.ref.fa
        else
        mv ${vcf_ref} pandora_genotyped_illumina.ref.fa
        fi
        """
    }
}*/

if (params.illumina_reads_1) {
    pandora_illumina_vcf.set { illumina_channels }
} else {
    illumina_channels = Channel.from()
}
if (params.nanopore_reads) {
    pandora_full_vcf.concat(pandora_30X_vcf).set { nanopore_channels }
} else {
    nanopore_channels = Channel.from()
}
nanopore_channels.concat(illumina_channels).set { pandora_vcfs }

/*process filter_pandora_vcfs {
    memory { 2.GB * task.attempt }
    errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
    maxRetries 2
    container {
      'shub://rmcolq/Singularity_recipes:minos'
    }

    publishDir final_outdir, mode: 'copy', overwrite: false

    input:
    set(file(pandora_vcf), file(pandora_ref), file(pandora_consensus)) from pandora_vcfs

    output:
    set(file("pandora_genotyped_*.vcf"), file(pandora_ref), file(pandora_consensus)) into filtered_pandora_vcfs

    """
    python /nfs/leia/research/iqbal/rmcolq/scripts/filter_vcf_by_pos_in_ref.py --vcf ${pandora_vcf} --ref ${pandora_ref} -n 31
    v=${pandora_ref}
    file_id=\${v:0:\${#v}-7}
    mv *.good.vcf \$file_id.vcf
    """
}*/

