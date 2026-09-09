# scRNA-seq

Two R scripts for exploratory single-cell RNA sequencing analysis with Seurat 5: one processes a single sample, and one merges three samples while preserving their labels. Both stop after saving a Seurat object.

The scripts are prepared examples reviewed as text only. They have not been executed or validated on a supplied dataset. This repository contains code and a general explanation; it includes no input data or analysis outputs.

## Files

| File | Purpose |
|---|---|
| [single_sample_02_seurat.R](scripts/single_sample_02_seurat.R) | Process one independent sample. |
| [merged_02_seurat.R](scripts/merged_02_seurat.R) | Combine three samples for pooled processing and plots coloured by sample. |

## Prepare the inputs

Use an interactive RStudio session with R packages `Seurat` (version 5 or later), `SeuratObject`, `dplyr`, `ggplot2` and `patchwork` available. The scripts check required packages and do not install them.

Set RStudio's working directory to the repository root before running a script. The scripts use `getwd()` to locate the following input folders:

```text
data/
  single_sample/filtered_feature_bc_matrix/
  sample_1/filtered_feature_bc_matrix/
  sample_2/filtered_feature_bc_matrix/
  sample_3/filtered_feature_bc_matrix/
```

The single-sample script uses only `single_sample`. The merged script uses only `sample_1`, `sample_2` and `sample_3`. These names are placeholders that you can change in the configuration section. The single-sample input is not automatically added to the merge.

Each input is a Cell Ranger-style filtered matrix of raw UMI counts, with `matrix.mtx.gz`, `features.tsv.gz` and `barcodes.tsv.gz`. The single-sample script also accepts uncompressed versions. Provide the matching feature and barcode files from the same matrix. No particular organism, cell type, matrix dimensions or biological groups are assumed.

## How the code works

1. **Check the environment and inputs.** The code requires an interactive RStudio session, checks package availability and matrix files, and creates a new results folder for the run.
2. **Create Seurat objects.** `Read10X()` loads the matrix. `CreateSeuratObject()` initially retains features detected in at least three cells and cells with at least 200 detected features. Cell counts are recorded at each filtering stage.
3. **Merge samples when requested.** The merged script creates one object per sample, prefixes barcodes with sample IDs, and stores sample labels in cell metadata. `JoinLayers()` combines the count layers for pooled processing. Merging does not perform batch correction or establish biological comparability.
4. **Inspect quality and filter cells.** Violin and scatter plots show detected features (`nFeature_RNA`) and total counts (`nCount_RNA`). The subsequent filter retains cells with strictly more than 200 detected features. No mitochondrial or upper-count threshold is applied.
5. **Normalise and select variable features.** `NormalizeData()` uses `LogNormalize` with scale factor 10,000. `FindVariableFeatures()` selects up to 2,000 variable features with the `vst` method.
6. **Scale and reduce dimensions.** All retained genes are scaled, then PCA uses the selected variable features. The scripts save loading plots, a PCA plot and a PC 1 heatmap. UMAP uses PCs 1-10 with random seed 42.
7. **Save results and stop.** Plots, cell metadata, count summaries and session information are written under `results/`. Each script ends immediately after `saveRDS()` at `# STOP HERE`.

Review the QC plots before continuing through the later sections. The numeric settings are starting examples, not validated thresholds for every dataset. Choose filters and the number of principal components using the supplied data and experimental design. The input needs enough cells and variable features for the requested PCA and UMAP calculations. Scaling all genes can require substantial memory.

Clustering, marker detection, cell-type annotation and integration are outside these scripts. An embedding alone does not identify cell types or establish differences between groups.

## Outputs and version control

Outputs are created in a new timestamped folder under `results/single_sample/` or `results/merged/`. Generated files stay local and are ignored by Git, together with input data, R session files and environment configuration. Review any generated metadata or session information separately before sharing it.

## Method documentation

- [Seurat introductory workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Merging Seurat objects](https://satijalab.github.io/seurat-object/reference/merge.Seurat.html)
- [Splitting and joining assay layers](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
