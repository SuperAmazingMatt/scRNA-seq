# Three-sample Seurat analysis
# Open in RStudio, set the working directory to the repository root,
# and review each section before running it.
# Inputs are filtered feature-barcode matrices containing raw UMI counts.
# merge() combines the samples without batch correction or integration.

# Setup -------------------------------------------------------------------
if (!interactive() || Sys.getenv("RSTUDIO") != "1") {
  stop("Open this script in RStudio and run it there, not with Rscript.")
}
required_packages <- c("Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "),
       ". Make these packages available in RStudio before continuing.")
}
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
if (packageVersion("Seurat") < "5.0.0") {
  stop("This merged workflow requires Seurat 5 or later.")
}
options(Seurat.object.assay.version = "v5")
set.seed(42)

project_dir <- getwd()
sample_info <- data.frame(
  sample_id = c("sample_1", "sample_2", "sample_3"),
  sample_label = c("sample_1", "sample_2", "sample_3"),
  stringsAsFactors = FALSE
)
sample_info$matrix_dir <- file.path(
  project_dir, "data", sample_info$sample_id, "filtered_feature_bc_matrix"
)
required_files <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
for (matrix_dir in sample_info$matrix_dir) {
  if (!all(file.exists(file.path(matrix_dir, required_files)))) {
    stop("A required matrix file is missing from: ", matrix_dir)
  }
}

# Each run has its own output folder so previous results are preserved.
run_id <- paste0("merged_", format(Sys.time(), "%Y%m%d_%H%M%S"))
output_dir <- file.path(project_dir, "results", "merged", run_id)
if (dir.exists(output_dir)) stop("Output folder already exists: ", output_dir)
if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create: ", output_dir)
write.csv(sample_info, file.path(output_dir, "sample_inputs.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_start.txt"))
message("Results will be saved in: ", output_dir)

# Read matrices and create one object per sample ---------------------------
# Initial filters: genes detected in >=3 cells within each sample, and cells
# with >=200 detected genes at object creation. We record those removals.
sample_objects <- setNames(vector("list", nrow(sample_info)), sample_info$sample_id)
cell_counts <- data.frame(sample_id = sample_info$sample_id,
                          cellranger_barcodes = NA_integer_,
                          after_creation = NA_integer_)
for (i in seq_len(nrow(sample_info))) {
  counts <- Read10X(data.dir = sample_info$matrix_dir[i])
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) stop("Gene Expression counts missing.")
    counts <- counts[["Gene Expression"]]
  }
  cell_counts$cellranger_barcodes[i] <- ncol(counts)
  sample_object <- CreateSeuratObject(
    counts = counts, project = sample_info$sample_id[i],
    min.cells = 3, min.features = 200
  )
  if (ncol(sample_object) == 0L) stop("No cells retained for ", sample_info$sample_id[i])
  sample_object$sample_id <- sample_info$sample_id[i]
  sample_object$sample_label <- sample_info$sample_label[i]
  cell_counts$after_creation[i] <- ncol(sample_object)
  sample_objects[[i]] <- sample_object
}

# Merge all three samples before the common analysis ----------------------
# Prefix barcodes to keep them unique and preserve each cell's sample labels.
seurat_merged <- merge(
  x = sample_objects[["sample_1"]],
  y = list(sample_objects[["sample_2"]], sample_objects[["sample_3"]]),
  add.cell.ids = sample_info$sample_id,
  project = "merged_samples", merge.data = FALSE
)
stopifnot(ncol(seurat_merged) == sum(cell_counts$after_creation),
          !anyDuplicated(colnames(seurat_merged)))

# Seurat 5 preserves separate count layers when merging. Join counts into one
# layer so variable-feature selection uses the pooled dataset rather than a
# consensus of per-sample variable features.
# API: https://satijalab.github.io/seurat-object/reference/SplitLayers.html
layers_before <- SeuratObject::Layers(seurat_merged[["RNA"]])
seurat_merged <- SeuratObject::JoinLayers(seurat_merged, assay = "RNA")
layers_after <- SeuratObject::Layers(seurat_merged[["RNA"]])
stopifnot(identical(layers_after, "counts"))
capture.output(list(before_join = layers_before, after_join = layers_after),
               file = file.path(output_dir, "merged_layers.txt"))
seurat_merged$sample_label <- factor(
  seurat_merged$sample_label, levels = c("sample_1", "sample_2", "sample_3")
)
Idents(seurat_merged) <- "sample_label"
rm(counts, sample_object, sample_objects)

# QC plots and the initial feature filter -----------------------------------
qc_violin <- VlnPlot(
  seurat_merged, features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "sample_label", ncol = 2
)
print(qc_violin)
ggsave(file.path(output_dir, "qc_violin_before_filter.png"), qc_violin,
       width = 10, height = 5, dpi = 300)
qc_scatter <- FeatureScatter(
  seurat_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
  group.by = "sample_label"
)
print(qc_scatter)
ggsave(file.path(output_dir, "qc_scatter_before_filter.png"), qc_scatter,
       width = 7, height = 5, dpi = 300)
write.csv(seurat_merged[[]], file.path(output_dir, "cell_metadata_before_filter.csv"))

# This example retains cells with more than 200 detected genes.
# No mitochondrial or upper-count threshold is applied.
# Choose suitable thresholds after inspecting the data.
seurat_merged <- subset(seurat_merged, subset = nFeature_RNA > 200)
cell_counts$after_feature_filter <- as.integer(table(factor(
  seurat_merged$sample_id, levels = sample_info$sample_id
)))
print(cell_counts)
write.csv(cell_counts, file.path(output_dir, "cell_counts_by_stage.csv"), row.names = FALSE)
write.csv(seurat_merged[[]], file.path(output_dir, "cell_metadata_after_filter.csv"))
if (any(cell_counts$after_feature_filter == 0L)) {
  stop("At least one sample has no cells after filtering. Review QC before continuing.")
}

# Normalise and identify variable genes -----------------------------------
seurat_merged <- NormalizeData(
  seurat_merged, normalization.method = "LogNormalize", scale.factor = 1e4
)
seurat_merged <- FindVariableFeatures(
  seurat_merged, selection.method = "vst", nfeatures = 2000
)
top_vf <- head(VariableFeatures(seurat_merged), 10)
print(top_vf)
write.csv(data.frame(gene = VariableFeatures(seurat_merged)),
          file.path(output_dir, "variable_genes.csv"), row.names = FALSE)
plot1 <- VariableFeaturePlot(seurat_merged)
plot2 <- LabelPoints(plot1, points = top_vf, repel = TRUE, xnudge = 0, ynudge = 0)
print(plot1)
print(plot2)
ggsave(file.path(output_dir, "variable_genes_top10.png"), plot2,
       width = 8, height = 6, dpi = 300)

# Scale all retained genes and run PCA ------------------------------------
# Scaling all retained genes uses more memory than scaling variable genes only.
all.genes <- rownames(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, features = all.genes)
seurat_merged <- RunPCA(
  seurat_merged, features = VariableFeatures(seurat_merged), seed.use = 42
)
# [["pca"]] accesses the stored PCA reduction.
print(seurat_merged[["pca"]], dims = 1:5, nfeatures = 5)
capture.output(print(seurat_merged[["pca"]], dims = 1:5, nfeatures = 5),
               file = file.path(output_dir, "pca_top_loadings.txt"))
pca_loadings <- VizDimLoadings(seurat_merged, dims = 1:2, reduction = "pca")
print(pca_loadings)
ggsave(file.path(output_dir, "pca_loadings.png"), pca_loadings,
       width = 10, height = 5, dpi = 300)
pca_plot <- DimPlot(seurat_merged, reduction = "pca", group.by = "sample_label")
print(pca_plot)
ggsave(file.path(output_dir, "pca_by_sample.png"), pca_plot,
       width = 7, height = 5, dpi = 300)
# fast=FALSE changes the drawing method so the same heatmap can be saved.
pca_heatmap <- DimHeatmap(
  seurat_merged, dims = 1, cells = 500, balanced = TRUE, fast = FALSE
)
print(pca_heatmap)
ggsave(file.path(output_dir, "pca_heatmap.png"), pca_heatmap,
       width = 10, height = 7, dpi = 300)

# UMAP and the final saved object -----------------------------------------
# PCs 1-10 are an initial setting to review for the supplied data.
seurat_merged <- RunUMAP(seurat_merged, dims = 1:10, seed.use = 42)
umap_plot <- DimPlot(seurat_merged, reduction = "umap", group.by = "sample_label")
print(umap_plot)
ggsave(file.path(output_dir, "umap_by_sample.png"), umap_plot,
       width = 7, height = 5, dpi = 300)
write.csv(seurat_merged[[]], file.path(output_dir, "cell_metadata_final.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo_final.txt"))
message("Saving the merged object. No analysis after STOP HERE is included.")
saveRDS(seurat_merged, file = file.path(output_dir, "seurat_merged_saved.rds"))
# STOP HERE
