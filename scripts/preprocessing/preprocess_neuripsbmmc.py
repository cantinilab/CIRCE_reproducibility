import scanpy as sc

atac = sc.read_h5ad(snakemake.input['anndata'])

atac = atac[:, atac.var['feature_types'] == 'ATAC']
atac = atac[atac.obs['batch'] == "s3d7"]

#sc.pp.subsample(atac, n_obs=500, random_state=0)
#atac = atac[:, :10000]
atac = atac[atac.X.sum(1) > 0, :]
atac = atac[:, atac.X.sum(0) > 0]

atac.var_names = atac.var_names.str.replace(':', '_', regex=False)
atac.var_names = atac.var_names.str.replace('-', '_', regex=False)

print(atac.var_names)
print(atac.obs_names)
print(atac)

atac.X = atac.X.toarray()
atac.write_h5ad(snakemake.output['anndata'])