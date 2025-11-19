import scipy as sp
import scanpy as sc
import circe as ci
import numpy as np

print('Running circe')

if __name__ == '__main__':
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
    print("""'atac raw data': {}""".format(atac))
    atac = ci.add_region_infos(atac, sep=('_', '_'))

    # Pseudocells
    if binarize:
        atac.X = atac.X > 0
    elif 'bmmcneurips' in atac_path:
        print('Using counts layer for bmmcneurips')
        atac.X  = atac.layers['counts']  # use counts layer for bmmcneurips
                    # Since the main data matrix has been previously binarized

    atac = ci.metacells.compute_metacells(
        atac,
        max_metacells = 5000,
        k=50,
        dim_reduction='lsi'
    )

    print("""'atac metacell data': {}""".format(atac))

    # Make data sparse
    if not sp.sparse.issparse(atac.X):
        atac.X = sp.sparse.csr_matrix(atac.X)

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
