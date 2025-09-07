This repository contains the code to replicate the experiment of CIRCE's manuscript, accessible [here]()

<p align="left">
    <img alt="Circe overview figure" src=circe_figure.png width="600">
</p>

## Organisation
The code is a snakemake workflow, from data downloading to network inference.
<br>The figures are then generated in the notebook [Figures_notebook](Figures_notebook.ipynb).

You will find here:
1. The short benchmark of the packages [CIRCE](https://github.com/cantinilab/Circe), in python, and [Cicero](https://github.com/cole-trapnell-lab/cicero-release), in R.
    - This comparison is done on 2 datasets (PBMC_10X and BMMC_neurips2021) but can easily be extended to other datasets.
2. The scripts to process the [human fetal scATAC atlas](https://www.science.org/doi/10.1126/science.aba7612), used to evaluate CIRCE on very large datasets.

## Execution
### 1. Build singularities
You will first need to install singularity/apptainer and to re-create the containers of CIRCE and Cicero:

```singularity build envs/circe.sif envs/circe.def```

### 2. Run jobs
You can now produce the output of interest for you with:

```
snakemake -c {number_of_cores} {output_requested} --use-singularity
# Replacing {number_of_cores} and {output_requested} by your values.
```

Prior to this command, you can check the jobs that will be launched to produce this output with:
```
snakemake -n {output_requested}
```
! The scripts calling Cicero might run in more than a day, hence I would advice you to run only one/few job(s) at a time.

### 3. Draw figure panels
You can then reproduce the figure panels which have all necessary outputs available. 





