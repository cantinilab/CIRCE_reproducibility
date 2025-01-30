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

# Prepare directly cicero input
print("prepare Cicero input")
set.seed(seed)
cicero_cds <- cicero::make_atac_cds(acc, binarize = FALSE)

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
  window = snakemake@params$distance_threshold,  # Default = 5e+05
  silent = FALSE,                              # Default = FALSE
  sample_num = snakemake@params$sample_num)  # Default = 100

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

