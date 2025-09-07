import scanpy as sc

atac_input_file = snakemake.input.anndata
atac_output_file = snakemake.output.anndata

atac = sc.read_h5ad(atac_input_file)
atac.var_names = atac.var_names.str.replace('-', '_')

# Change the format to sparse matrix csc_matrix
atac.X = atac.X.tocsc()

# Save the preprocessed AnnData object
atac.write_h5ad(atac_output_file)