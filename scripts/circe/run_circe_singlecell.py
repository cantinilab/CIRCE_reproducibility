import scanpy as sc
import circe as ci

print('Running circe')
# Parameters
atac_path = snakemake.input['atac']
distance_threshold = snakemake.params['distance_threshold']
circe_results_path = snakemake.output['circe_results']
n_samples = snakemake.params['n_samples']
n_maxtry = snakemake.params['n_maxtry']
max_n_iter_alpha = snakemake.params['max_n_iter_alpha']
n_jobs = snakemake.params['n_jobs']
binarize = snakemake.params['binarize']

# Load data
atac = sc.read_h5ad(atac_path)
print(atac.var_names)
atac = ci.add_region_infos(atac, sep=('_', '_'))

if binarize:
    atac.X = atac.X>0
# Compute network
ci.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=n_samples,
    n_samples_maxtry=n_maxtry,
    max_alpha_iteration=max_n_iter_alpha,
    njobs=n_jobs
)

# Extract and save results
peak_layer = ci.extract_atac_links(
    atac,
    columns=['Peak1', 'Peak2', 'coaccess'])
peak_layer.to_csv(circe_results_path, sep='\t', index=False)
