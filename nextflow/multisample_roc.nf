params.truth_assembly1 = ""
params.vcf_directory1 = ""
params.truth_assembly2 = ""
params.vcf_directory2 = ""

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
          --truth_assembly1	FILE	Assembly against which to compare the calls for sample 1
          --vcf_directory1	DIRECTORY	Directory containing vcfs to compare for sample 1
          --mask1		FILE	BED file of regions to exclude from analysis for sample 1
          --truth_assembly2	FILE	Assembly against which to compare the calls for sample 2
          --vcf_directory2	DIRECTORY	Directory containing vcfs to compare for sample 2
          --mask2               FILE    BED file of regions to exclude from analysis for sample 2

        Optional arguments:
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

if (params.truth_assembly1) {
    truth_assembly1 = file(params.truth_assembly1).toAbsolutePath()
    if (!truth_assembly1.exists()) {
        exit 1, "Truth assembly file not found: ${params.truth_assembly1} -- aborting"
    }
}
else {
    exit 1, "Truth assembly file 1 not provided -- aborting"
}

if (params.truth_assembly2) {
    truth_assembly2 = file(params.truth_assembly2).toAbsolutePath()
    if (!truth_assembly2.exists()) {
        exit 1, "Truth assembly file not found: ${params.truth_assembly2} -- aborting"
    } 
}
else {
    exit 1, "Truth assembly file 2 not provided -- aborting"
}

if (params.mask1) {
    mask1 = file(params.mask1).toAbsolutePath()
    if (!mask1.exists()) {
        exit 1, "Mask file not found: ${params.mask1} -- aborting"
    }   
}   
if (params.mask2) {
    mask2 = file(params.mask2).toAbsolutePath()
    if (!mask2.exists()) {
        exit 1, "Mask file not found: ${params.mask2} -- aborting"
    }   
}

if (params.vcf_directory1) {
    vcf_directory1 = file(params.vcf_directory1).toAbsolutePath()
    if (!vcf_directory1.exists()) {
        exit 1, "Truth assembly file not found: ${params.vcf_directory1} -- aborting"
    }
}
else {
    exit 1, "VCF directory 1 not provided -- aborting"
}

if (params.vcf_directory2) {
    vcf_directory2 = file(params.vcf_directory2).toAbsolutePath()
    if (!vcf_directory2.exists()) {
        exit 1, "Truth assembly file not found: ${params.vcf_directory2} -- aborting"
    }
}   
else {    
    exit 1, "VCF directory 2 not provided -- aborting"
}

if (params.mask1 && params.mask2) {
    process compare_vcfs {
        errorStrategy {task.attempt < 2 ? 'retry' : 'ignore'}
        maxRetries 2
        memory {1.4.GB * task.attempt}
        container {
              'shub://rmcolq/Singularity_recipes:minos'
            }

        input:
        file truth_assembly1
        file truth_assembly2
        file vcf_directory1
        file vcf_directory2
        file mask1
        file mask2
        file params.pipeline_root

        output:
        file('*.pkl') into pickled_dfs

        """
        python3 ${params.pipeline_root}/scripts/compare_genotypers_on_pair_of_samples.py --truth1 ${truth_assembly1} --truth2 ${truth_assembly2} --sample_dir1 ${vcf_directory1} --sample_dir2 ${vcf_directory2} --mask1 ${mask1} --mask2 ${mask2}
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

        input:
        file truth_assembly1
        file truth_assembly2
        file vcf_directory1
        file vcf_directory2
        file params.pipeline_root

        output:
        file('*.pkl') into pickled_dfs
    
        """
        python3 ${params.pipeline_root}/scripts/compare_genotypers_on_pair_of_samples.py --truth1 ${truth_assembly1} --truth2 ${truth_assembly2} --sample_dir1 ${vcf_directory1} --sample_dir2 ${vcf_directory2}
        """
    }
}

/*dfs = pickled_dfs.collectFile(name: 'all.pkl')

process make_graph {
    errorStrategy {task.attempt < 2 ? 'retry' : 'fail'}
    maxRetries 2
    memory {0.7.GB * task.attempt}
    container {
          'shub://rmcolq/Singularity_recipes:minos'
        }

    publishDir final_outdir, mode: 'copy', overwrite: false

    input:
    file 'all.pkl' from dfs

    output:
    'roc*.png'

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
        ax.set_xlabel('Number FPs/Number Genotyped', size = 26)
        ax.set_ylabel('Fraction of dnadiff SNPs discoverable from VCFs', size = 26)
        #ax.set_title('Precision Recall', size = 30)
        ax.set_xlim(right=0.1)

        # Set colormaps
        colormap_pandora = plt.cm.autumn
        colormap_snippy = plt.cm.winter
        colormap_nanopolish = plt.cm.cool
        snippy_i = 0
        nanopolish_i = 0

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
            else:
                print(x['name'])

        handles, labels = ax.get_legend_handles_labels()
        hl = sorted(zip(handles, labels),key=operator.itemgetter(1))
        handles2, labels2 = zip(*hl)

        ax.legend(handles2, labels2, frameon=False, loc='lower right')
        plt.savefig('roc%d.png' %i, transparent=True)
    """
}
*/
