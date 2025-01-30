rule download_gene_annot:
    output:
        human_annot_file = "data/genome/hg38_gene_annot.bed",
        mouse_annot_file = "data/genome/mm10_gene_annot.bed"
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/analysis/download_gene_annot.py"


rule download_hg19_to_hg38:
    output:
        chain_file = "data/genome/hg19ToHg38.over.chain.gz"
    shell:
        """
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O {output.chain_file}
        # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O {output.chain_file}
        
        """


rule map_PCHic_regions_from_hg19_to_hg38:
    input:
        chain_file = "data/genome/hg19ToHg38.over.chain.gz",
        PCHiC = "data/PCHiC/PBMC_hg19_PCHiC.tsv"
    output:
        hg19_enhancers = "data/PCHiC/hg19_enhancers.bed",
        hg19_promoters = "data/PCHiC/hg19_promoters.bed",
        hg38_enhancers = "data/PCHiC/hg38_enhancers.bed",
        hg38_promoters = "data/PCHiC/hg38_promoters.bed",
        PCHiC_hg38 = "data/PCHiC/PBMC_hg38_PCHiC.tsv"
    #singularity:
    #    "envs/circe.sif"
    shell:
        """
        python scripts/analysis/PCHiC_hg19_regions_to_bed.py \
            --PC_HiC {input.PCHiC} \
            --bed_enhancers {output.hg19_enhancers} \
            --bed_promoters {output.hg19_promoters}

        head {output.hg19_enhancers}

        module load UCSC-tools/v390
        liftOver {output.hg19_enhancers} {input.chain_file} {output.hg38_enhancers} {output.hg38_enhancers}.unmapped
        liftOver {output.hg19_promoters} {input.chain_file} {output.hg38_promoters} {output.hg38_promoters}.unmapped

        python scripts/analysis/PCHiC_to_hg38.py \
            --PC_HiC {input.PCHiC} \
            --hg19_enhancers {output.hg19_enhancers} \
            --hg19_promoters {output.hg19_promoters} \
            --hg38_enhancers {output.hg38_enhancers} \
            --hg38_promoters {output.hg38_promoters} \
            --hg38_PCHiC {output.PCHiC_hg38}
        """

rule download_promoter_capture_HiC:
    params:
        url = "https://www.cell.com/cms/10.1016/j.cell.2016.09.037/attachment/5bc79f6f-1b69-4192-8cb8-4247cc2e0f39/mmc4.zip"
    output:
        human_pchic = "data/PCHiC/PBMC_hg19_PCHiC.tsv",
    shell:
        """
        # Dataset has to be downloaded and extract manually (PCHiC_meak_matrix_cutoff5.tsv)

        cp data/PCHiC/PCHiC_peak_matrix_cutoff5.tsv {output.human_pchic}
        """


rule auroc_PCHiC_overlap:
    input:
        PCHiC = "data/PCHiC/PBMC_hg38_PCHiC.tsv",
        atac = "data/datasets/{sample}/{sample}.h5ad",
        predictions = "results/{sample}/{method}/{method}_{clustering}.tsv"
    output:
        results = "results/{sample}/{method}/{method}_{clustering}_PCHiC_overlap.tsv"
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/analysis/PCHiC_overlap.py"


rule auroc_table:
    input:
        auroc_files = expand("results/{sample}/{method}_PCHiC_overlap.tsv",
            sample=["pbmc10x"],
            method=[
                "circe/circe_singlecell",
                "circe/circe_bin_singlecell",
                "circe/circe_pseudocell",
                "circe/circe_bin_pseudocell",
                "cicero/cicero_singlecell",
                "cicero/cicero_pseudocell",
                ])
    output:
        auroc_table = "results/{sample}/auroc_table.tsv"
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/analysis/compute_auroc_table.py"


rule correlation_table:
    input:
        circe_singlecell = "results/{sample}/circe/circe_singlecell_from_cicero.tsv",
        cicero_singlecell = "results/{sample}/cicero/cicero_singlecell.tsv",
        cicero_pseudocell = "results/{sample}/cicero/cicero_pseudocell.tsv",
        circe_pseudocell = "results/{sample}/circe/circe_pseudocell_from_cicero.tsv",
    output:
        correlation_table_pearson = "results/{sample}/correlation_table_pearson.tsv",
        correlation_table_spearman = "results/{sample}/correlation_table_spearman.tsv"
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/analysis/compute_correlation_table.py"
