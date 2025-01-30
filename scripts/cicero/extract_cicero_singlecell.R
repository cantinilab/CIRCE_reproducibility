library(rhdf5)
library(cicero)

# Parameters
reduction_method <- "UMAP"
seed <- 2
number_cells_per_clusters <- snakemake@params$number_cells_per_clusters
atac_path <- snakemake@input$atac_path
cicero_results_path <- snakemake@output$cicero_results


# Read peaks and atac
print("Open Data")
indata <- H5Fopen(atac_path, flags='H5F_ACC_RDONLY')
atac <- as.matrix(indata$X)
colnames(atac) <- indata$obs$`_index`
rownames(atac) <- indata$var$`_index`
rownames(atac) <- gsub(":", "-", rownames(atac))
rownames(atac) <- gsub("_", "-", rownames(atac))
#atac <- atac[1:1000, 1:1000]

# Matrix to edgelist
print("Reshape to edge list")
acc <- reshape2::melt(atac)
colnames(acc) <- c("V1", "V2", "V3")

# Prepare directly cicero input
print("prepare Cicero input")
set.seed(seed)
cicero_cds <- cicero::make_atac_cds(acc, binarize = FALSE)

row.names(fData(cicero_cds))

print("Save Cicero cells")
#data.table::fwrite(
 #   as.matrix(cicero_cds@assays@data$counts),
  #  snakemake@output$cicero_cds_tsv
   # )
# mat = apply(as.matrix(cicero_cds@assays@data$counts), c(1,2), function(x) ifelse (x==0, "", x))
mat = as.matrix(cicero_cds@assays@data$counts)
colnames(mat) = colnames(atac)
rownames(mat) = row.names(fData(cicero_cds))
#rownames(mat)
write.csv(
    mat,
    quote=FALSE,
    row.names=TRUE,
    snakemake@output$cicero_cds_tsv
    )



