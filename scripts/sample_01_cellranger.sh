#!/usr/bin/env bash
#SBATCH --job-name=cellranger
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/01_cellranger_logs/%x_%j.out
#SBATCH --error=logs/01_cellranger_logs/%x_%j.err

# Generalised single-sample Cell Ranger job for a Slurm cluster.
# Submit from the repository root after creating the log directory.
# Supply the cluster account and partition to sbatch, not in this shared file.
# Usage: sbatch [cluster options] scripts/sample_01_cellranger.sh SAMPLE_ID [FASTQ_PREFIX]
# This version is prepared for reuse; it has not been run on a supplied dataset.

set -euo pipefail

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }

[[ $# -ge 1 && $# -le 2 ]] || fail 'Provide SAMPLE_ID and optionally FASTQ_PREFIX.'
sample_id="$1"
fastq_prefix="${2:-$sample_id}"
[[ "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || fail 'Use letters, numbers, underscores or hyphens for SAMPLE_ID.'
[[ "$fastq_prefix" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail 'FASTQ_PREFIX must be a single filename prefix without spaces, slashes or wildcards.'
[[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]] || fail 'Submit this script through sbatch from the repository root.'
[[ "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]] || fail 'The job needs a positive Slurm CPU allocation.'

project_root="$SLURM_SUBMIT_DIR"
fastq_dir="$project_root/data/$sample_id/fastq"
reference_dir="$project_root/data/reference"
results_root="$project_root/results/01_cellranger_results"

[[ -d "$fastq_dir" ]] || fail "Missing FASTQ directory: $fastq_dir"
[[ -s "$reference_dir/reference.json" ]] || fail 'Place a prepared Cell Ranger reference at data/reference, or link that directory to one.'
[[ ! -e "$results_root/$sample_id" && ! -L "$results_root/$sample_id" ]] || fail 'This sample already has a results path. Review the existing run before any intentional restart.'

# Accept paired reads across one or more lanes, compressed or uncompressed.
# Cell Ranger performs the full input and reference validation.
shopt -s nullglob
read1=("$fastq_dir/${fastq_prefix}_"*_R1_*.fastq "$fastq_dir/${fastq_prefix}_"*_R1_*.fastq.gz)
read2=("$fastq_dir/${fastq_prefix}_"*_R2_*.fastq "$fastq_dir/${fastq_prefix}_"*_R2_*.fastq.gz)
[[ ${#read1[@]} -gt 0 && ${#read2[@]} -gt 0 ]] || fail 'No matching R1/R2 FASTQs were found for FASTQ_PREFIX.'
for input_file in "${read1[@]}" "${read2[@]}"; do
    [[ -s "$input_file" ]] || fail "Empty FASTQ: $input_file"
done

# Edit this module line to match the software installation on your cluster.
# With a standalone installation, make cellranger available on PATH instead.
module load cellranger/8.0.1
command -v cellranger >/dev/null 2>&1 || fail 'Cell Ranger is not available.'
cellranger --version

mkdir -p "$results_root"
cd "$results_root"
printf 'Run ID: %s\nFASTQ prefix: %s\n' "$sample_id" "$fastq_prefix"

cellranger count \
    --id="$sample_id" \
    --transcriptome="$reference_dir" \
    --fastqs="$fastq_dir" \
    --sample="$fastq_prefix" \
    --create-bam=true \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem=60

# Accounting queried inside the batch job may still show RUNNING.
# Check the final job state and exit code again after the job has ended.
if command -v sacct >/dev/null 2>&1; then
    sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS || true
fi
