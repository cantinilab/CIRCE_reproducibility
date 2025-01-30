import numpy as np
import scanpy as sc
import circe as ci

print('Running circe')
# Parameters
cicero_cds = snakemake.input['cicero_cds']
distance_threshold = snakemake.params['distance_threshold']
circe_results_path = snakemake.output['circe_results']
n_samples = snakemake.params['n_samples']
n_maxtry = snakemake.params['n_maxtry']
max_n_iter_alpha = snakemake.params['max_n_iter_alpha']
n_jobs = snakemake.params['n_jobs']

# Load data
atac = sc.read_csv(cicero_cds, delimiter=',').transpose()
sc.pp.filter_genes(atac, min_cells=1)
# sc.pp.filter_cells(atac, min_genes=1)

print(atac.var_names)
print(atac.obs_names)
print(np.isnan(atac.X).sum())
atac = ci.add_region_infos(atac, sep=('-', '-'))

print("clustering", snakemake.wildcards["clustering"])
if snakemake.wildcards["clustering"]=="singlecell":
    pass
elif snakemake.wildcards["clustering"]=="pseudocell":
    atac = ci.metacells.compute_metacells(atac)
else:
    raise ValueError("snakemake.wildcards['clustering'] is not 'singlecell' not 'pseudocell': {}".format(snakemake.wildcards["clustering"])) 
# Compute network
ci.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=n_samples,
    n_samples_maxtry=n_maxtry,
    max_alpha_iteration=max_n_iter_alpha,
    njobs=n_jobs,
    seed=2024
)

# Extract and save results
peak_layer = ci.extract_atac_links(
    atac,
    columns=['Peak1', 'Peak2', 'coaccess'])
peak_layer.to_csv(circe_results_path, sep='\t', index=False)
