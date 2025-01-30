import scipy as sp
import scanpy as sc
import circe as ci

print('Running circe')
# Parameters
atac_path = snakemake.input['atac']
distance_threshold = snakemake.params['distance_threshold']
num_cells = snakemake.params['number_cells_per_clusters']
circe_results_path = snakemake.output['circe_results']
n_samples = snakemake.params['n_samples']
n_maxtry = snakemake.params['n_maxtry']
max_n_iter_alpha = snakemake.params['max_n_iter_alpha']
n_jobs = snakemake.params['n_jobs']
binarize = snakemake.params["binarize"]

# Load data
atac = sc.read_h5ad(atac_path)
atac = ci.add_region_infos(atac, sep=('_', '_'))

# Make data sparse
atac.X = sp.sparse.csr_matrix(atac.X)

# Pseudocells
if binarize:
    atac.X = atac.X > 0


atac = ci.metacells.compute_metacells(
    atac,
    max_metacells = 5000,
    k=50,
    dim_reduction='lsi'
)

import numpy as np
print(np.isnan(atac.X).sum())

print("""'atac data': {}""".format(atac))

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
)

# Extract and save results
peak_layer = ci.extract_atac_links(
    atac,
    columns=['Peak1', 'Peak2', 'coaccess'])

peak_layer.to_csv(circe_results_path, sep='\t', index=False)
