# preprocess smallpbmc data
rule download_smallpbmc_data:
    output:
        mudata = "data/datasets/smallpbmc/smallpbmc.h5mu",
        #"data/{sample}_ground_truth.tsv"
    params:
        url = "https://figshare.com/ndownloader/files/44639476?private_link=b0840d90e42e37fa165f"
    shell:
        """
        wget {params.url} -O {output.mudata}
        """

rule preprocess_full_small_pbmc_data:
    input:
        mudata = "data/datasets/smallpbmc/smallpbmc.h5mu"
    output:
        anndata = "data/datasets/smallpbmc/smallpbmc.h5ad",
#        csv = "data/smallpbmc.csv.gz"
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/preprocess_mudata.py"


# preprocess celllines data
rule download_celllines_data:
    output:
        csv = "data/datasets/celllines/celllines.csv",
    params:
        url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-08205-7/MediaObjects/41467_2018_8205_MOESM6_ESM.xls"
    shell:
        """
        wget {params.url} -O {output.csv}
        """

rule save_h5mu_celllines_data:
    input:
        csv = "data/datasets/celllines/celllines.csv"
    output:
        mudata = "data/datasets/celllines/celllines.h5mu",
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/csv_to_mdata.py"

rule preprocess_full_celllines_data:
    input:
        mudata = "data/datasets/celllines/celllines.h5mu"
    output:
        anndata = "data/datasets/celllines/celllines.h5ad",
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/preprocess_mudata.py"


# preprocess pbmc10x data
rule download_bmmcneurips:
    output:
        anndata = "data/datasets/bmmcneurips/multiome_bmmcneurips.h5ad",
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122%5Fopenproblems%5Fneurips2021%5Fmultiome%5FBMMC%5Fprocessed%2Eh5ad%2Egz"
    shell:
        """
        wget {params.url} -O {output.anndata}.gz
        gunzip {output.anndata}.gz
        """

rule preprocess_neuripsbmmc:
    input:
        anndata = "data/datasets/bmmcneurips/multiome_bmmcneurips.h5ad"
    output:
        anndata = "data/datasets/bmmcneurips/bmmcneurips.h5ad",
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/preprocess_neuripsbmmc.py"