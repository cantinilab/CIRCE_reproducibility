import anndata as ad
import scanpy as sc
import muon as mu
import scipy as sp

# Open the CSV file
adata = sc.read_csv(snakemake.input["csv"], delimiter="\t").T
print(adata)

# Filter out low quality cells
adata.X = adata.X[adata.X.sum(axis=1) > 0,]
adata.X = adata.X[:, adata.X.sum(axis=0) > 0]

adata.X = sp.sparse.csr_matrix(adata.X)
# Save the AnnData object as h5ad
mdata = mu.MuData(
    {"atac": adata},
)
mdata.write_h5mu(snakemake.output["mudata"])
