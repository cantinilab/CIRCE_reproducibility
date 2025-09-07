# Benchmarking  Cicero / atac_networks (python)
from snakemake.utils import min_version
min_version("7.0.0")

configfile: "config/data_config.yaml"
samples = ["bmmcneurips"]  #["smallpbmc", "pbmc10x"]
method = ["circe", "cicero"]
clustering = ["singlecell", "pseudocell"]


# Target rule to define the desired final output
rule all:
    input:
        # Cicero runs : for both single-cells and metacells
#        expand("results/{sample}/cicero/cicero{bin}_{clustering}.tsv",
#            sample=samples, method=method, clustering=clustering, bin=[""]),
        # CIRCE runs : for both single-cells and metacells
        expand("results/{sample}/circe/circe{bin}_{clustering}.tsv",
            sample=samples, clustering=clustering, bin=["", "_bin"]),
        expand("results/{sample}/circe/circe_{clustering}_from_cicero.tsv",
            sample=samples, clustering=clustering),

        # For correlation between Cicero and Circe

        # For evaluation on PCHiC - pbmc10x
        expand("results/pbmc10x/cicero/cicero{bin}_{clustering}_PCHiC_overlap.tsv",
            clustering=clustering, bin=[""]),        # Add the binarized version for Circe
        expand("results/pbmc10x/circe/circe{bin}_{clustering}_PCHiC_overlap.tsv",
            clustering=clustering, bin=["", "_bin"]),     
        # Table auroc
        expand("results/{sample}/auroc_table.tsv", sample=["pbmc10x"]),

# 5. Run CIRCE on the human fetal atlas
module human_fetal_atlas:
    snakefile: "modules/human_fetal_atlas.smk"
    config: config
use rule * from human_fetal_atlas

# 4. Compare Cicero and atac_networks
module enhancers_analysis:
    snakefile: "modules/enhancers_analysis.smk"
    config: config
use rule * from enhancers_analysis

## 3. Module to run Cicero
#module cicero:
#    snakefile: "modules/cicero.smk"
#    config: config
#use rule * from cicero

# 2. Module to run circe
module circe:
    snakefile: "modules/circe.smk"
    config: config
use rule * from circe

# 1. Generate/Preprocess the fake scATAC-seq data
module data_generation:
    snakefile: "modules/data_processing.smk"
    config: config
use rule * from data_generation

# 1. Generate/Preprocess the pbmc10x data
module preprocess_pbmc10x:
    snakefile: "modules/pbmc10x.smk"
    config: config
use rule * from preprocess_pbmc10x
