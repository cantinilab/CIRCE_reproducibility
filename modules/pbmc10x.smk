# Already processed
rule zenodo_download_public:
    output: "data/datasets/pbmc10x/pbmc10x.h5mu"
    params:
        record="17661273",
        filename="pbmc10x.h5mu"
    shell:
        r"""
        wget -O {output} \
          "https://zenodo.org/api/records/{params.record}/files/{params.filename}/content"
        """

###################################################################################
# Download and preprocessing of PBMC10X dataset from GRETA benchmark repository:  #
#      - https://github.com/saezlab/greta_benchmark/                              #
###################################################################################


#rule download_pbmc10x:
#    output:
#        mudata = "data/pbmc_neurips/pbmc_neurips.h5mu",
#    params:
#        url_matrix: 'https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5'
#        url_peak_annot: 'https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv'
#        url_atac_frags: 'https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
#        url_atac_index: 'https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi'
#
#    shell:
#        """
#        wget {params.url} -O {output.mudata}
#        """
#
#
#rule prcannot_pbmc10x:
#    output:
#        tmp=temp(directory(local('datasets/pbmc10x/tmp'))),
#        annot=temp(local('datasets/pbmc10x/annot.csv')),
#    shell:
#        """
#        python workflow/scripts/datasets/pbmc10x/prc_annot.py \
#        -t {output.tmp} \
#        -a {output.annot}
#        """
#
#rule callpeaks_pbmc10x:
#    input:
#        frags='datasets/pbmc10x/smpl.frags.tsv.gz',
#        annot='datasets/pbmc10x/annot.csv',
#    singularity:
#        'workflow/envs/gretabench.sif'
#    output:
#        tmp=temp(directory(local('datasets/pbmc10x/tmp_peaks'))),
#        peaks=temp(local('datasets/pbmc10x/peaks.h5ad'))
#    resources:
#        mem_mb=64000,
#    threads: 16
#    shell:
#        """
#        python workflow/scripts/datasets/callpeaks.py \
#        -f {input.frags} \
#        -a {input.annot} \
#        -t {output.tmp} \
#        -o {output.peaks}
#        """
#
#rule annotate_pbmc10x:
#    input:
#        annot='datasets/pbmc10x/annot.csv',
#        g='gdata/geneids',
#        peaks='datasets/pbmc10x/peaks.h5ad',
#    singularity:
#        'workflow/envs/gretabench.sif'
#    output:
#        tmp=temp(directory(local('datasets/pbmc10x/tmp_annot'))),
#        out='datasets/pbmc10x/pbmc10x.h5mu'
#    params:
#        organism=config['datasets']['pbmc10x']['organism'],
#    resources:
#        mem_mb=32000,
#    shell:
#        """
#        python workflow/scripts/datasets/pbmc10x/pbmc10x.py \
#        -a {output.tmp} \
#        -b {input.annot} \
#        -c {input.g} \
#        -d {params.organism} \
#        -e {input.peaks} \
#        -f {output.out}
#        """


rule preprocess_pbmc10x:
    input:
        mudata = "data/datasets/pbmc10x/pbmc10x.h5mu"
    output:
        anndata = "data/datasets/pbmc10x/pbmc10x.h5ad",
    singularity:
        "envs/circe.sif"
    script:
        "../scripts/preprocessing/preprocess_mudata.py"