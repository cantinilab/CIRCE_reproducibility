rule download_genomesizes:
    output:
        hg38='data/genome/hg38.txt',
        mm10='data/genome/mm10.txt',
    shell:
        """
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes' -O {output.hg38}
        wget 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes' -O {output.mm10}
        """

rule run_singlecell_cicero:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac_path = "data/datasets/{sample}/{sample}.h5ad",
        human_genome = "data/genome/hg38.txt",
        mouse_genome = "data/genome/mm10.txt"
    output:
        cicero_results = "results/{sample}/cicero/cicero_singlecell.tsv",
    params:
        organism = lambda wildcards: config[wildcards.sample]["organism"],
        distance_threshold = 500_000,
        sample_num = 100
    benchmark:
        "benchmark/{sample}/singlecell_cicero_time.txt"
    singularity:
        "envs/cicero.sif"
    script:
        "../scripts/cicero/run_cicero_singlecell.R"

rule run_pseudocell_cicero:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac_path = "data/datasets/{sample}/{sample}.h5ad",
        human_genome = "data/genome/hg38.txt",
        mouse_genome = "data/genome/mm10.txt"
    output:
        cicero_results = "results/{sample}/cicero/cicero_pseudocell.tsv",
    params:
        number_cells_per_clusters = 50,
        organism = lambda wildcards: config[wildcards.sample]["organism"],
        distance_threshold = 500_000,
        sample_num = 100
    benchmark:
        "benchmark/{sample}/pseudocell_cicero_time.txt"
    singularity:
        "envs/cicero.sif"
    script:
        "../scripts/cicero/run_cicero_pseudocell.R"


# Extract the same cells than the ones used in the cicero run
# To run Circe on the same ones.
# Could be extracted from run_{clustering}_circe, but would have to rerun all
rule extract_singlecell_cicero:
    input:
        atac_path = "data/datasets/{sample}/{sample}.h5ad",
    output:
        cicero_cds_tsv = "results/{sample}/cicero/cicero_singlecell_cells.tsv",
    params:
        organism = lambda wildcards: config[wildcards.sample]["organism"],
        distance_threshold = 500_000,
        sample_num = 100
    singularity:
        "envs/cicero.sif"
    script:
        "../scripts/cicero/extract_cicero_singlecell.R"


rule extract_pseudocell_cicero:
    input:
        atac_path = "data/datasets/{sample}/{sample}.h5ad",
    output:
        cicero_cds_tsv = "results/{sample}/cicero/cicero_pseudocell_cells.tsv",
    params:
        number_cells_per_clusters = 50,
        organism = lambda wildcards: config[wildcards.sample]["organism"],
        distance_threshold = 500_000,
        sample_num = 100
    singularity:
        "envs/cicero.sif"
    script:
        "../scripts/cicero/extract_cicero_pseudocell.R"