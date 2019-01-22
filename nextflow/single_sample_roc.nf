params.truth_assembly = ""
params.vcf_directory = ""

params.help = false
params.testing = false
params.pipeline_root = "/nfs/leia/research/iqbal/rmcolq/git/DPhil_analysis"
params.final_outdir = "."
params.max_forks = 10
params.snps = false

if (params.help){
    log.info"""
	Pipeline for running pandora and other variant callers and comparing the results.

        Usage: nextflow run compare_genotypers.nf <arguments>
        Required arguments:
	  --truth_assembly	FILE	Assembly against which to compare the calls
          --vcf_directory	DIRECTORY	Directory containing vcfs to compare

	Optional arguments:
	  --mask		FILE	BED file of regions to mask out of analyses
	  --snps		FLAG    Restrict to SNPs?
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

if (params.mask) {
    mask = file(params.mask).toAbsolutePath()
    if (!mask.exists()) {
        exit 1, "Mask BED file not found: ${params.mask} -- aborting"
    }
}

vcfs_to_check = Channel.fromFilePairs("${params.vcf_directory}/*.{vcf,ref.fa}", flat:true)
pandora_consensus = Channel.fromPath("${params.vcf_directory}/pandora*consensus*")

if (params.snps) {
    process filter_vcfs {
        memory { 1.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        container {
            'shub://rmcolq/Singularity_recipes:minos'
        }
    
        input:
        set(val(name), file(vcf_reference), file(vcf)) from vcfs_to_check

        output:
        set (val("${name}"), file("${vcf_reference}"), file("*.filtered.vcf")) into filtered_vcfs

        """
        python3 ${params.pipeline_root}/scripts/filter_vcf.py --vcf ${vcf} --vcf_ref ${vcf_reference} --flank 41 --snps
        """
    }
} 
else {
    process filter_vcfs {
        memory { 1.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        container {
            'shub://rmcolq/Singularity_recipes:minos'
        }
            
        input:
        set(val(name), file(vcf_reference), file(vcf)) from vcfs_to_check
        
        output:
        set (val("${name}"), file("${vcf_reference}"), file("*.filtered.vcf")) into filtered_vcfs
        
        """
        python3 ${params.pipeline_root}/scripts/filter_vcf.py --vcf ${vcf} --vcf_ref ${vcf_reference} --flank 41
        """
    }
}

if (params.mask) {
    process check_vcf {
        memory { 1.GB * task.attempt }
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        container {
            'shub://rmcolq/Singularity_recipes:minos'
        }

        input:
        set(val(name), file(vcf_reference), file(vcf)) from filtered_vcfs
        file truth_assembly
        file mask

        output:
        set file("minos.*"), val("${name}"), file("${vcf_reference}"), file("${truth_assembly}") into minos_output

        """
        echo "${vcf}"
        echo "${vcf_reference}"
        v=${vcf_reference}
        if [ \${v: -3} == ".gz" ]
        then
          zcat \$v | awk '{print \$1;}' > \${v::-3}
          v=\${v::-3}
        else
          cat \$v | awk '{print \$1;}' > n.\$v
          mv n.\$v \$v
        fi
        bwa index ${truth_assembly}
        minos check_with_ref ${vcf} \$v ${truth_assembly} minos --allow_flank_mismatches --flank_length 41 --variant_merge_length 41 --include_ref_calls --exclude_bed ${mask}
        """
    }
} 
else {
process check_vcf {
    memory { 1.GB * task.attempt }
    errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        } 
    input:
    set(val(name), file(vcf_reference), file(vcf)) from filtered_vcfs
    file truth_assembly
    file mask
    
    output:
    set file("minos.*"), val("${name}"), file("${vcf_reference}"), file("${truth_assembly}") into minos_output
    
    """
    echo "${vcf}"
    echo "${vcf_reference}"
    v=${vcf_reference}
    if [ \${v: -3} == ".gz" ]
    then
      zcat \$v | awk '{print \$1;}' > \${v::-3}
      v=\${v::-3}
    else
      cat \$v | awk '{print \$1;}' > n.\$v
      mv n.\$v \$v
    fi
    bwa index ${truth_assembly}
    minos check_with_ref ${vcf} \$v ${truth_assembly} minos --allow_flank_mismatches --flank_length 41 --variant_merge_length 41 --include_ref_calls
    """
}
}
   

process minos_to_df {
    errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    input:
    set file(minos_files), val(name), file(vcf_reference), file(truth_assembly) from minos_output

    output:
    file 'sample.pkl' into pickled_dfs
    set file('minos.vcf'), stdout, val(name), file("${vcf_reference}"), file("${truth_assembly}") into filter_input
    
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import _pickle as pk

    tp = pd.read_table('minos.gt_conf_hist.TP.tsv')
    fp = pd.read_table('minos.gt_conf_hist.FP.tsv')
    stats = pd.read_table('minos.stats.tsv')

    precision = []
    recall = []
    yscat = []
    xscat = []
    sensitivity = []
    specificity_p = []
    conf_threshold = 0
    c = list(tp['GT_CONF'].values) + list(fp['GT_CONF'].values)
    c.sort()
    old_total = stats['total'].values[0] - stats['gt_excluded'].values[0]
    total_tp = sum(tp['Count'])
    total_fp = sum(fp['Count'])
    total = total_tp + total_fp
    print(old_total, total, total_tp, total_fp)
    for confidence in c:
        num_tp = sum(tp[tp['GT_CONF'] >= confidence]['Count'])
        num_fp = sum(fp[fp['GT_CONF'] >= confidence]['Count'])
        genotyped = num_tp + num_fp
        if genotyped > 0.2*total:
            precision.append(num_tp/float(genotyped))
            recall.append(num_tp/float(total))
            yscat.append(genotyped/float(total))
            xscat.append(num_fp/float(genotyped))
            sensitivity.append(num_tp/float(stats['gt_correct'].values[0]))
            specificity_p.append(num_fp/float(stats['gt_wrong'].values[0]))
            if recall[-1] < 0.8 and conf_threshold == 0:
                conf_threshold = confidence

    df = pd.DataFrame()
    df['precision'] = precision
    df['recall'] = recall
    df['yscat'] = yscat
    df['xscat'] = xscat
    df['sensitivity'] = sensitivity
    df['specificity_p'] = specificity_p
    df['name'] = "${name}"
    df.to_csv("df.tsv", sep='\t')
    pk.dump(df, open('sample.pkl', 'wb'))
    print(conf_threshold)
    """
}

dfs = pickled_dfs.collectFile(name: 'all.pkl')
(pr_dfs, roc_dfs, fpc_dfs) = dfs.separate(3) { a -> [a, a, a] }


process make_pr {
    errorStrategy {task.attempt < 2 ? 'retry' : 'fail'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

    input:
    file 'all.pkl' from pr_dfs
   
    output:
    'pr.png'

    """
    #!/usr/bin/env python3
    import pandas as pd
    import _pickle as pickle
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import operator

    def loadall(filename):
        with open(filename, "rb") as f:
            while True:
                try:
                    yield pickle.load(f)
                except EOFError:
                    break

    items = loadall('all.pkl')

    for i in range(5):
        # Define plot
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 12)

        ax.grid(b=True)
        ax.set_axisbelow(b=True)
        plt.style.use('seaborn-colorblind')

        # Label the axes and title the plot
        ax.set_xlabel('Number true positives/Number genotyped', size = 26)
        ax.set_ylabel('Number true positives/Total number of sites', size = 26)
        #ax.set_title('Precision Recall', size = 30)

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        colormap_other = plt.cm.summer
        snippy_i = 0
        nanopolish_i = 0
        other_i = 0
    
        # Make a scatter plot
        items = loadall('all.pkl')
        for x in items:
            if len(x['name'].values) > 0:
                if x['name'].values[0].startswith("pandora_genotyped_full"):
                    col = colormap_pandora(0)
                    ax.scatter(x['precision'],x['recall'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina") and i > 0:
                    col = colormap_pandora(100)
                    ax.scatter(x['precision'],x['recall'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_30") and i > 1:
                    col = colormap_pandora(200)
                    ax.scatter(x['precision'],x['recall'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("snippy") and i > 2:
                    col = colormap_snippy(snippy_i*20)
                    ax.scatter(x['precision'],x['recall'], label=x['name'].values[0], c=col, alpha=0.5)
                    snippy_i += 1
                elif x['name'].values[0].startswith("nanopolish") and i > 3:
                    col = colormap_nanopolish(nanopolish_i*20)
                    ax.scatter(x['precision'],x['recall'], label=x['name'].values[0], c=col, alpha=0.5)
                    nanopolish_i += 1
                elif x['name'].values[0].startswith("pandora_genotyped_"):
                    col = colormap_other(other_i*100)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    other_i += 1
            else:
                print(x['name'])

        handles, labels = ax.get_legend_handles_labels()
        hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
        handles2, labels2 = zip(*hl)
    
        ax.legend(handles2, labels2, frameon=False, loc='lower right')
        plt.savefig('pr%d.png' %i, transparent=True)
    """
}

process make_fpc {
    errorStrategy {task.attempt < 2 ? 'retry' : 'fail'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

    input:
    file 'all.pkl' from fpc_dfs

    output:
    'fpc4.png'

    """
    #!/usr/bin/env python3
    import pandas as pd
    import _pickle as pickle
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import operator

    def loadall(filename):
        with open(filename, "rb") as f:
            while True:
                try:
                    yield pickle.load(f)
                except EOFError:
                    break

    items = loadall('all.pkl')

    for i in range(5):
        # Define plot
        fig, ax = plt.subplots()
        fig.set_size_inches(16, 12)

        ax.grid(b=True)
        ax.set_axisbelow(b=True)
        plt.style.use('seaborn-colorblind')

        # Label the axes and title the plot
        ax.set_xlabel('Number false positives/Number genotyped', size = 26)
        ax.set_ylabel('Number genotyped/Total number of sites', size = 26)
        #ax.set_title('Precision Recall', size = 30)

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        colormap_other = plt.cm.summer
        snippy_i = 0
        nanopolish_i = 0
        other_i = 0

        # Make a scatter plot
        items = loadall('all.pkl')
        for x in items:
            if len(x['name'].values) > 0:
                if x['name'].values[0].startswith("pandora_genotyped_full"):
                    col = colormap_pandora(0)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_illumina") and i > 0:
                    col = colormap_pandora(100)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("pandora_genotyped_30") and i > 1:
                    col = colormap_pandora(200)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                elif x['name'].values[0].startswith("snippy") and i > 2:
                    col = colormap_snippy(snippy_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    snippy_i += 1
                elif x['name'].values[0].startswith("nanopolish") and i > 3:
                    col = colormap_nanopolish(nanopolish_i*20)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    nanopolish_i += 1
                elif x['name'].values[0].startswith("pandora_genotyped_"):
                    col = colormap_other(other_i*100)
                    ax.scatter(x['xscat'],x['yscat'], label=x['name'].values[0], c=col, alpha=0.5)
                    other_i += 1
            else:
                print(x['name'])

        handles, labels = ax.get_legend_handles_labels()
        hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
        handles2, labels2 = zip(*hl)

        ax.legend(handles2, labels2, frameon=False, loc='lower right')
        plt.savefig('fpc%d.png' %i, transparent=True)
    """
}

process make_roc {
    errorStrategy {task.attempt < 2 ? 'retry' : 'fail'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

    input:
    file 'all.pkl' from roc_dfs

    output:
    'roc.png'

    """
    #!/usr/bin/env python3
    import pandas as pd
    import _pickle as pickle
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import operator
    
    def loadall(filename):
        with open(filename, "rb") as f:
            while True:
                try:
                    yield pickle.load(f)
                except EOFError:
                    break
                    
    items = loadall('all.pkl')
    
    # Define plot
    fig, ax = plt.subplots()
    fig.set_size_inches(16, 12)
    
    ax.grid(b=True)
    ax.set_axisbelow(b=True)
    plt.style.use('seaborn-colorblind')
    
    # Label the axes and title the plot
    ax.set_ylabel('Sensitivity', size = 26)
    ax.set_xlabel('1 - Specificity', size = 26)
    ax.set_title('ROC', size = 30)
    
    # Set colormaps
    colormap_pandora = plt.cm.autumn
    colormap_snippy = plt.cm.winter
    colormap_nanopolish = plt.cm.cool
    pandora_i = 0
    snippy_i = 0
    nanopolish_i = 0

    # Make a scatter plot
    items = loadall('all.pkl')
    for x in items:
        if len(x['name'].values) > 0:
            if x['name'].values[0].startswith("pandora"):
                col = colormap_pandora(pandora_i*100)
                ax.scatter(x['specificity_p'],x['sensitivity'], label=x['name'].values[0], c=col, alpha=0.5)
                pandora_i += 1
            elif x['name'].values[0].startswith("snippy"):
                col = colormap_snippy(snippy_i*20)
                ax.scatter(x['specificity_p'],x['sensitivity'], label=x['name'].values[0], c=col, alpha=0.5)
                snippy_i += 1
            elif x['name'].values[0].startswith("nanopolish"):
                col = colormap_nanopolish(nanopolish_i*20)
                ax.scatter(x['specificity_p'],x['sensitivity'], label=x['name'].values[0], c=col, alpha=0.5)
                nanopolish_i += 1
            else:
                ax.scatter(x['specificity_p'],x['sensitivity'], label=x['name'].values[0])
        else:
            print(x['name'])
            
    handles, labels = ax.get_legend_handles_labels()
    hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)

    ax.legend(handles2, labels2, frameon=False, loc='lower right')
    plt.savefig('roc.png', transparent=True)
    """
}

/*process filter_vcf {
    errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
    maxRetries 2
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }
    
    input:
    set file(minos_vcf), val(conf_threshold), val(name), file(vcf_reference), file(truth_assembly) from filter_input
    
    output:
    set file("out.vcf"), val(conf_threshold), val(name), file(vcf_reference), file(truth_assembly) into filter_output
    
    """
    #!/usr/bin/env python3
    import os
    import vcf

    vcf_reader = vcf.Reader(open("${minos_vcf}", 'r'))
    vcf_writer = vcf.Writer(open('out.vcf', 'w'), vcf_reader)
    thresh = ${conf_threshold}
    for record in vcf_reader:
        if 'GT' in record.FORMAT and [sample['GT'] for sample in record.samples][0] == '1/1' and \
           'GT_CONF' in record.FORMAT and float([sample['GT_CONF'] for sample in record.samples][0]) > float(thresh):
            vcf_writer.write_record(record) 
    """
}

process evaluate_genome_covered {
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        } 
        
    input:
    set file(minos_vcf), val(conf_threshold), val(name), file(vcf_reference), file(truth_assembly) from filter_output
    
    output:
    file('summary.txt') into vcf_summary
    
    """   
    bgzip ${minos_vcf}
    tabix -p vcf ${minos_vcf}.gz
    v=${vcf_reference}
    if [ \${v: -3} == ".gz" ]
    then
        zcat ${vcf_reference} | vcf-consensus ${minos_vcf}.gz > corrected.fa
    else
        cat ${vcf_reference} | vcf-consensus ${minos_vcf}.gz > corrected.fa
    fi
    zcat ${truth_assembly} > truth_assembly.fa
    dnadiff truth_assembly.fa corrected.fa
    echo "${name}\t\$(grep "AlignedBases" out.report | awk '{print \$2}')\t\$(grep "TotalSNPs" out.report | awk '{print \$2}')" > summary.txt
    """
}  

pandora_consensus_in = Channel.from(pandora_consensus).map { file -> tuple("${file}.baseName", file) }

process evaluate_genome_covered_consensus {
    errorStrategy {task.attempt < 1 ? 'retry' : 'ignore'}
    maxRetries 1
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        } 
        
    input:
    set name, file(consensus_seqs) from pandora_consensus_in
    file truth_assembly

    output:
    file('summary.txt') into consensus_summary

    """   
    zcat ${consensus_seqs} > corrected.fa
    zcat ${truth_assembly} > truth_assembly.fa
    dnadiff truth_assembly.fa corrected.fa
    echo "${name}\t\$(grep "AlignedBases" out.report | awk '{print \$2}')\t\$(grep "TotalSNPs" out.report | awk '{print \$2}')" > summary.txt
    """
}

vcf_summary.concat(consensus_summary).set { summary }

summary.collectFile(name: final_outdir/'summary.txt')
*/
