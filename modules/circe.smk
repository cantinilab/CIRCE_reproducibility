rule run_singlecell_circe_bin:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac = "data/datasets/{sample}/{sample}.h5ad"
    params:
        binarize = True,
        distance_threshold = 500_000,
        n_samples = 100,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 20
    benchmark:
        "benchmark/{sample}/singlecell_circe_time_bin.txt"
    singularity:
        "envs/circe.sif"
    output:
        circe_results = "results/{sample}/circe/circe_bin_singlecell.tsv"
    script:
        "../scripts/circe/run_circe_singlecell.py"

rule run_singlecell_circe:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac = "data/datasets/{sample}/{sample}.h5ad"
    params:
        binarize = False,
        distance_threshold = 500_000,
        n_samples = 100,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 20
    benchmark:
        "benchmark/{sample}/singlecell_circe_time.txt"
    singularity:
        "envs/circe.sif"
    output:
        circe_results = "results/{sample}/circe/circe_singlecell.tsv"
    script:
        "../scripts/circe/run_circe_singlecell.py"


rule run_pseudocell_circe_bin:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac = "data/datasets/{sample}/{sample}.h5ad"
    params:
        binarize = True,
        distance_threshold = 500_000,
        number_cells_per_clusters = 50,
        n_samples = 100,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 20
    benchmark:
        "benchmark/{sample}/pseudocell_circe_time_bin.txt"
    singularity:
        "envs/circe.sif"
    output:
        circe_results = "results/{sample}/circe/circe_bin_pseudocell.tsv"
    script:
        "../scripts/circe/run_circe_pseudocell.py"

rule run_pseudocell_circe:
    resources:
        cores=config["resources"]["cores"],
        mem_mb=config["resources"]["mem_mb"]
    input:
        atac = "data/datasets/{sample}/{sample}.h5ad"
    params:
        binarize = False,
        distance_threshold = 500_000,
        number_cells_per_clusters = 50,
        n_samples = 100,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 20
    benchmark:
        "benchmark/{sample}/pseudocell_circe_time.txt"
    singularity:
        "envs/circe.sif"
    output:
        circe_results = "results/{sample}/circe/circe_pseudocell.tsv"
    script:
        "../scripts/circe/run_circe_pseudocell.py"


rule circe_from_cicero:
    input:
        cicero_cds = "results/{sample}/cicero/cicero_singlecell_cells.tsv",
        human_genome = "data/genome/hg38.txt",
        mouse_genome = "data/genome/mm10.txt"
    params:
        organism = lambda wildcards: config[wildcards.sample]["organism"],
        distance_threshold = 500_000,
        n_samples = 500,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 20
    singularity:
        "envs/circe.sif"
    output:
        circe_results = "results/{sample}/circe/circe_{clustering}_from_cicero.tsv"
    script:
        "../scripts/circe/circe_from_cicero.py"
