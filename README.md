# scRNA-seq

A terminal script runs Cell Ranger on one sample in a Slurm job. Two R scripts then support exploratory analysis with Seurat 5: one processes a single sample, and one merges three samples while preserving their labels. Both R scripts stop after saving a Seurat object.

These generalised examples have not been run on a supplied dataset. The shell script has been checked with `bash -n`; the R scripts have been reviewed as text only. The code examples include no input data. A separate [poster portfolio](posters/README.md) shows original bulk RNA-seq graphs with study-specific findings obscured.

## Files

| File | Purpose |
|---|---|
| [seurat_workflows.Rmd](seurat_workflows.Rmd) | R Markdown with both R workflows, step explanations and analysis disabled when knitting. |
| [sample_01_cellranger.sh](scripts/sample_01_cellranger.sh) | Run Cell Ranger for one sample through Slurm, using generic input paths. |
| [single_sample_02_seurat.R](scripts/single_sample_02_seurat.R) | Process one independent sample. |
| [merged_02_seurat.R](scripts/merged_02_seurat.R) | Combine three samples for pooled processing and plots coloured by sample. |

## Read the R Markdown guide

Open [seurat_workflows.Rmd](seurat_workflows.Rmd) in RStudio. It contains the full
code from both R scripts in separate workflows, with explanations of inputs,
parameters, plots and saved outputs. Choose the single-sample or three-sample
workflow to match the inputs.

With `rmarkdown`, `knitr` and Pandoc available, **Knit** creates an HTML guide.
All analysis chunks have `eval=FALSE`, so knitting displays the code without
running Seurat or reading data. Only the document-options chunk runs. See the
[knitr option reference](https://yihui.org/knitr/options/).

Interactive RStudio Run controls do execute code. Use them only for an intended
analysis, in a dedicated session with the repository root as the working
directory, and follow the chosen workflow in order. Both workflows stop after
saving the object. If changing code, keep the matching R script and R Markdown
chunks consistent. Generated HTML stays local and is ignored by Git.

## Run Cell Ranger from the terminal

The shell script is a reusable version of the recorded single-sample job workflow. It needs Bash, Slurm, Cell Ranger 8.0.1, paired FASTQs from a compatible gene-expression library, and a matching prepared Cell Ranger reference. Use a Linux cluster terminal. If matrices are already available, continue with the R input section below.

Arrange the inputs relative to the repository root:

```text
data/
  reference/                 # Prepared reference root containing reference.json
  sample_1/fastq/             # Matching R1/R2 FASTQs from every lane for this sample
  sample_2/fastq/
  sample_3/fastq/
```

The `reference` directory may be a symbolic link to an existing reference. Preserve the original FASTQ names. For example, if a file is named `library_A_S1_L001_R1_001.fastq.gz`, its FASTQ prefix is `library_A`; the script's sample ID can still be `sample_1`. Each job processes one sample, with its lanes kept together.

Review the script's `module load cellranger/8.0.1` line for your cluster. For a standalone installation, replace that line with the appropriate PATH setup. The script requests one task, eight CPUs, 64 GB of memory and two hours; Cell Ranger uses the allocated CPU count and a 60 GB memory limit. These are example resources to review for the dataset. BAM generation is enabled.

Create the log directory **before** submission. From the repository root, replace `YOUR_ACCOUNT` and `YOUR_PARTITION` with your cluster settings, and supply the actual FASTQ prefix as the final argument:

```bash
mkdir -p logs/01_cellranger_logs
sbatch --parsable --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/sample_01_cellranger.sh sample_1 library_A
```

The second argument is optional when the FASTQ prefix matches the sample ID. For three independent inputs whose prefixes are `sample_1`, `sample_2` and `sample_3`, the submission commands are:

```bash
for sample in sample_1 sample_2 sample_3; do
  sbatch --parsable --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
    scripts/sample_01_cellranger.sh "$sample"
done
```

Use either the single-job example or the loop as appropriate. A returned job ID confirms submission. Check scheduling and the final state with `squeue -j JOB_ID` and `sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS`, replacing `JOB_ID` with the returned number. Inspect the Cell Ranger log, `web_summary.html` and `metrics_summary.csv` before interpreting the results. The script refuses an existing sample results path to avoid an accidental rerun; an intentional restart needs separate review.

Each run writes to `results/01_cellranger_results/<sample_id>/outs/`. After successful completion and output checks, link its filtered matrix to the corresponding R input path. For example, from the repository root:

```bash
ln -s ../../results/01_cellranger_results/sample_1/outs/filtered_feature_bc_matrix \
  data/sample_1/filtered_feature_bc_matrix
```

Repeat the link for `sample_2` and `sample_3` when using the merged R script. For the independent R script, use the sample ID `single_sample` and its matching directory throughout. Existing links or directories must be reviewed before replacement. No job is submitted merely by downloading these files.

## Prepare the R inputs

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

## How the R code works

1. **Check the environment and inputs.** The code requires an interactive RStudio session, checks package availability and matrix files, and creates a new results folder for the run.
2. **Create Seurat objects.** `Read10X()` loads the matrix. `CreateSeuratObject()` initially retains features detected in at least three cells and cells with at least 200 detected features. Cell counts are recorded at each filtering stage.
3. **Merge samples when requested.** The merged script creates one object per sample, prefixes barcodes with sample IDs, and stores sample labels in cell metadata. `JoinLayers()` combines the count layers for pooled processing. Merging does not perform batch correction or establish biological comparability.
4. **Inspect quality and filter cells.** Violin and scatter plots show detected features (`nFeature_RNA`) and total counts (`nCount_RNA`). The subsequent filter retains cells with strictly more than 200 detected features. No mitochondrial or upper-count threshold is applied.
5. **Normalise and select variable features.** `NormalizeData()` uses `LogNormalize` with scale factor 10,000. `FindVariableFeatures()` selects up to 2,000 variable features with the `vst` method.
6. **Scale and reduce dimensions.** All retained genes are scaled, then PCA uses the selected variable features. The scripts save loading plots, a PCA plot and a PC 1 heatmap. UMAP uses PCs 1-10 with random seed 42.
7. **Save results and stop.** Plots, cell metadata, count summaries and session information are written under `results/`. Each script ends immediately after `saveRDS()` at `# STOP HERE`.

Review the QC plots before continuing through the later sections. The numeric settings are starting examples, not validated thresholds for every dataset. Choose filters and the number of principal components using the supplied data and experimental design. The input needs enough cells and variable features for the requested PCA and UMAP calculations. Scaling all genes can require substantial memory.

Clustering, marker detection, cell-type annotation and integration are outside the R scripts. An embedding alone does not identify cell types or establish differences between groups.

## Outputs and version control

Cell Ranger results are stored under `results/01_cellranger_results/`, with Slurm logs under `logs/01_cellranger_logs/`. R outputs are created in a new timestamped folder under `results/single_sample/` or `results/merged/`. Generated files stay local and are ignored by Git, together with input data, R session files and environment configuration. Review any generated metadata or session information separately before sharing it.

## Method documentation

- [Running Cell Ranger count, version 8.0](https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/running-pipelines/cr-gex-count)
- [Slurm sbatch reference](https://slurm.schedmd.com/sbatch.html)
- [Seurat introductory workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Merging Seurat objects](https://satijalab.github.io/seurat-object/reference/merge.Seurat.html)
- [Splitting and joining assay layers](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
