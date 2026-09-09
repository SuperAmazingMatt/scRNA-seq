# Single-sample Seurat analysis
# Open in RStudio, set the working directory to the repository root,
# and review each section before running it.
# Each run writes to a new results folder.

# 1. Check the RStudio session and installed packages -------------------------
if (!interactive() || Sys.getenv("RSTUDIO") != "1") {
  stop("Open this script in RStudio before running it.")
}
required_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       ". Make these packages available in RStudio before continuing.")
}
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
if (packageVersion("Seurat") < "5.0.0") {
  stop("This script requires Seurat 5 or later.")
}
options(Seurat.object.assay.version = "v5")
set.seed(42)

# 2. Select the input and create a new output folder -----------------
project_root <- getwd()
sample_id <- "single_sample"
matrix_dir <- file.path(project_root, "data", sample_id,
                        "filtered_feature_bc_matrix")
if (!dir.exists(matrix_dir)) {
  stop("Cannot find the input matrix: ", matrix_dir)
}
matrix_files <- file.path(matrix_dir, c("matrix.mtx", "features.tsv", "barcodes.tsv"))
if (any(!file.exists(matrix_files) & !file.exists(paste0(matrix_files, ".gz")))) {
  stop("The matrix folder is missing a matrix, feature, or barcode file.")
}

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(project_root, "results", "single_sample",
                     paste0("run_", run_stamp))
if (dir.exists(run_dir)) {
  stop("This output folder already exists. Run this section again for a new timestamp.")
}
if (!dir.create(run_dir, recursive = TRUE)) stop("Cannot create output folder: ", run_dir)
message("Input matrix: ", matrix_dir)
message("This run's outputs: ", run_dir)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), "")
)
utils::write.csv(package_versions, file.path(run_dir, "package_versions.csv"), row.names = FALSE)
writeLines(c(
  paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("RStudio session:", Sys.getenv("RSTUDIO")),
  paste("Input:", matrix_dir),
  "Initial settings: min.cells=3; min.features=200; nFeature_RNA > 200.",
  "LogNormalize scale.factor=10000; vst nfeatures=2000; scale all retained genes.",
  "PCA uses variable features; UMAP uses PCs 1:10; random seed=42.",
  "No mitochondrial filter, clustering, or cell annotation is applied."
), file.path(run_dir, "run_details.txt"))

# Print plots explicitly so RStudio Source displays them, and save a PNG copy.
save_plot <- function(plot, filename, width = 8, height = 5) {
  print(plot)
  ggplot2::ggsave(file.path(run_dir, filename), plot = plot,
                  width = width, height = height, units = "in", dpi = 300)
}

# 3. Read the Cell Ranger matrix and create the Seurat object -----------------
counts <- Seurat::Read10X(data.dir = matrix_dir)
if (is.list(counts)) {
  if (!"Gene Expression" %in% names(counts)) stop("Gene Expression matrix not found.")
  counts <- counts[["Gene Expression"]]
}
# Check the supplied matrix without assuming dataset-specific dimensions.
if (nrow(counts) == 0L || ncol(counts) == 0L) {
  stop("The input matrix has no features or barcodes.")
}
cell_counts <- data.frame(
  stage = "Cell Ranger filtered matrix",
  features = nrow(counts), cells = ncol(counts)
)
seurat_object <- Seurat::CreateSeuratObject(
  counts = counts, project = sample_id, min.cells = 3, min.features = 200
)
seurat_object$sample_id <- sample_id
cell_counts <- rbind(cell_counts, data.frame(
  stage = "CreateSeuratObject min.cells 3 and min.features 200",
  features = nrow(seurat_object), cells = ncol(seurat_object)
))
print(seurat_object)

# 4. Inspect cell quality and apply the initial feature filter --------
p_qc_violin <- Seurat::VlnPlot(
  seurat_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2
)
save_plot(p_qc_violin, "qc_violin_before_filter.png")
p_qc_scatter <- Seurat::FeatureScatter(
  seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
)
save_plot(p_qc_scatter, "qc_scatter_before_filter.png")

# This example retains cells with more than 200 detected genes.
# Select additional thresholds only after inspecting the data.
# No mitochondrial or upper-count filter is applied.
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200)
cell_counts <- rbind(cell_counts, data.frame(
  stage = "After nFeature_RNA greater than 200",
  features = nrow(seurat_object), cells = ncol(seurat_object)
))
print(cell_counts)
utils::write.csv(cell_counts, file.path(run_dir, "cell_counts_by_stage.csv"), row.names = FALSE)
utils::write.csv(seurat_object[[]], file.path(run_dir, "cell_metadata_after_filter.csv"))

# 5. Normalise counts and identify 2,000 variable features --------------------
seurat_object <- Seurat::NormalizeData(
  seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4
)
seurat_object <- Seurat::FindVariableFeatures(
  seurat_object, selection.method = "vst", nfeatures = 2000
)
top_vf <- head(Seurat::VariableFeatures(seurat_object), 10)
print(top_vf)
utils::write.csv(data.frame(gene = top_vf),
                 file.path(run_dir, "top10_variable_features.csv"), row.names = FALSE)
p_variable <- Seurat::VariableFeaturePlot(seurat_object)
p_variable_labels <- Seurat::LabelPoints(
  plot = p_variable, points = top_vf, repel = TRUE, xnudge = 0, ynudge = 0
)
save_plot(p_variable, "variable_features.png")
save_plot(p_variable_labels, "variable_features_top10.png")

# 6. Scale all retained genes and perform PCA on variable features ------------
# Scaling every retained gene uses more memory than
# scaling only the variable genes. PCA still uses the selected variable genes.
all.genes <- rownames(seurat_object)
seurat_object <- Seurat::ScaleData(seurat_object, features = all.genes)
seurat_object <- Seurat::RunPCA(
  seurat_object, features = Seurat::VariableFeatures(object = seurat_object),
  seed.use = 42
)
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
writeLines(capture.output(print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)),
           file.path(run_dir, "pca_top_features.txt"))

p_pca_loadings <- Seurat::VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
save_plot(p_pca_loadings, "pca_loadings_dims1_2.png", width = 10, height = 6)
p_pca <- Seurat::DimPlot(seurat_object, reduction = "pca")
save_plot(p_pca, "pca_cells.png")
# fast=FALSE returns a plot that can be printed and saved with ggsave.
# Display PC 1 with 500 cells and balanced feature selection.
p_pca_heatmap <- Seurat::DimHeatmap(
  seurat_object, dims = 1, cells = 500, balanced = TRUE, fast = FALSE
)
save_plot(p_pca_heatmap, "pca_heatmap_dim1.png", width = 10, height = 7)

# 7. Calculate and display UMAP using PCs 1-10 -------------------
# The seed improves repeatability within the recorded software environment.
# This is an exploratory embedding; no clusters or cell types are assigned.
seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:10, seed.use = 42)
p_umap <- Seurat::DimPlot(seurat_object, reduction = "umap")
save_plot(p_umap, "umap_cells.png")

# 8. Record the session and save the object
# Run this final section only after the earlier sections complete successfully.
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
rds_path <- file.path(run_dir, "single_sample_seurat_object.rds")
message("Saving the Seurat object to: ", rds_path)
saveRDS(seurat_object, file = rds_path)
# STOP HERE
