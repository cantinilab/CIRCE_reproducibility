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

# Matrix to edgelist
print("Reshape to edge list")
acc <- reshape2::melt(atac)
colnames(acc) <- c("V1", "V2", "V3")

# Prepare cicero input
print("prepare Cicero input")
set.seed(seed)

input_cds <- cicero::make_atac_cds(acc, binarize = TRUE)

# Monocle3 preprocessing
if (length(which(colSums(as.matrix(monocle3::exprs(input_cds))) == 0)) == 0
) {
  print("Run Monocle3 dimensionality reduction")
  input_cds <- monocle3::estimate_size_factors(
    input_cds)
  input_cds <- monocle3::preprocess_cds(
    input_cds,
    method = "LSI")  # Preprocessing using LS
  input_cds <- monocle3::reduce_dimension(
    input_cds,
    reduction_method = reduction_method,
    preprocess_method = "LSI")  # Dimensionality reduction using UMAP/other key
} else {
  print("Error: there is at least one cell with no signal.")
}

# Compute pseudocells
print("Compute pseudo-cells")
umap_coords <- SingleCellExperiment::reducedDims(input_cds)$UMAP
cicero_cds <- cicero::make_cicero_cds(
  input_cds,  # Create a Cicero CDS object
  reduced_coordinates = umap_coords,
  k = number_cells_per_clusters,  #number neighbors/ Default = 50
  summary_stats = NULL,
  size_factor_normalize = TRUE,
  silent = FALSE)


# Run Cicero
print("Run Cicero")
if (snakemake@params$organism == "mm10"){ # obtain chromosome sizes
  chromosome_sizes <- read.table(snakemake@input$mouse_genome)
} else if (snakemake@params$organism == "hg38") {
  chromosome_sizes <- read.table(snakemake@input$human_genome)
} else {
  stop("Species not supported")
}

print('Run cicero')
cicero <- run_cicero(
  cds = cicero_cds, # Infer peak-links
  genomic_coords = chromosome_sizes,
  window = snakemake@params$distance_threshold,             # Default = 5e+05
  silent = FALSE,             # Default = FALSE
  sample_num = snakemake@params$sample_num) # Default = 100

# Format cicero output
print("Format cicero output")
cicero <- cicero[which(!is.na(cicero$coaccess)), ]  # Remove NAs
my_cols <- which(as.character(cicero$Peak1) <= as.character(cicero$Peak2))
cicero <- cicero[my_cols, ] # Remove double edges
cicero$Peak1 <- gsub("_", "-", cicero$Peak1)
cicero$Peak2 <- gsub("_", "-", cicero$Peak2)

# Write cicero output
print("Write cicero output")
write.table(cicero, cicero_results_path, sep = "\t", quote = FALSE)

