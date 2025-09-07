
# preprocess human fetal atlas data
rule download_human_fetal_atlas:
    output:
        anndata = "data/datasets/human_fetal_atlas/pre_human_fetal_atlas.h5ad",
    params:
        url = "http://download.gao-lab.org/GLUE/dataset/Domcke-2020.h5ad"
    shell:
        """
        wget {params.url} -O {output.anndata}
        """
        # gunzip {output.anndata}

rule preprocess_human_fetal_atlas:
    input:
        anndata = "data/datasets/human_fetal_atlas/pre_human_fetal_atlas.h5ad"
    output:
        anndata = "data/datasets/human_fetal_atlas/human_fetal_atlas.h5ad",
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/preprocess_human_fetal_atlas.py"

rule run_circe_atlas:
    input:
        atac = "data/datasets/human_fetal_atlas/human_fetal_atlas.h5ad",
    singularity:
        "envs/circe.sif"
    params:
        binarize = False,
        distance_threshold = 500_000,
        n_samples = 100,
        n_maxtry = 500,
        max_n_iter_alpha = 100,
        n_jobs = 70
    benchmark:
        "benchmark/human_fetal_atlas/singlecell_circe_time.txt"
    output:
        circe_results = "results/human_fetal_atlas/circe_singlecell.tsv"
    script:
        "../scripts/circe/run_circe_singlecell.py"