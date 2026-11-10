#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=dofasta_and_mask
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00               # hh:mm:ss, per sample
#SBATCH --output=logs/%x-%A_%a.out    # %A = array job ID, %a = task number
#SBATCH --error=logs/%x-%A_%a.err

# Make a masked, genome-wide consensus sequence for each sample with
# dofasta_and_mask.sh (README section "Make masked sequences with ANGSD").
#
# This is an ARRAY job with one task per sample: each task uses its number
# ($SLURM_ARRAY_TASK_ID) to pick "the Nth BAM file" in the BAM directory.
#
# Submit from the 3_phylogeny directory, with one task per BAM file:
#   cd 3_phylogeny
#   mkdir -p logs
#   n=$(ls ../data/map_and_call_out_standard/02_bamfiles/*.dedup.bam | wc -l)
#   sbatch --array=1-"$n" ../slurmscripts/3_dofasta_and_mask_array.sh

# stop the script as soon as any command fails
set -euo pipefail

# ---- Software ---------------------------------------------------------------
# The programs can come from the tutorial container (recommended) or from conda.
#
# Container: set POPGEN_SIF to the path of popgen_tutorial.sif, either by
# running this in your terminal before you submit (the job inherits it):
#   export POPGEN_SIF=/path/to/popgen_tutorial.sif
# or by writing the path directly here instead, e.g. sif=/path/to/popgen_tutorial.sif
sif=${POPGEN_SIF:-}
#
# Conda: leave sif empty. If you submit the job with popgen_tutorial_env
# activated, the programs are usually found. Otherwise uncomment:
# ml PDC miniconda3                    # Dardel only: load conda as a module
# source activate popgen_tutorial_env  # or: conda activate popgen_tutorial_env

# $run is written in front of every program below. With the container it is
# "apptainer exec <sif>", which runs the program inside the container. With
# conda it is empty, so the program runs as usual.
run=""
if [ -n "$sif" ]; then
    run="apptainer exec $sif"
fi

# check that all programs are found before we start
for tool in angsd samtools bedtools; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
bamdir=../data/map_and_call_out_standard/02_bamfiles
maskdir=../data/map_and_call_out_standard/03_genotypes/02_maskfiles

# pick the BAM file that belongs to this array task (line number N of the list)
bam=$(ls "$bamdir"/*.dedup.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(basename "$bam" .dedup.bam)
mask="$maskdir/${sample}_snp_mask.bed.gz"

# ---- Run --------------------------------------------------------------------
echo "Task $SLURM_ARRAY_TASK_ID: sample $sample"
$run bash dofasta_and_mask.sh "$bam" "$mask"
