#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=cram_to_bam
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00               # hh:mm:ss, per CRAM file
#SBATCH --output=logs/%x-%A_%a.out    # %A = array job ID, %a = task number
#SBATCH --error=logs/%x-%A_%a.err

# Convert CRAM files to BAM files (README section "Using angsd").
#
# This is an ARRAY job: slurm starts one task per CRAM file, and every task
# gets its own number in $SLURM_ARRAY_TASK_ID (1, 2, 3, ...). We use that
# number to pick "the Nth CRAM file" from the directory, so all files are
# converted at the same time instead of one after another.
#
# Submit from the 1.2_pca_with_angsd directory, with one task per CRAM file:
#   cd 1.2_pca_with_angsd
#   mkdir -p logs
#   n=$(ls ../data/map_and_call_out_standard/02_bamfiles/*.dedup.cram | wc -l)
#   sbatch --array=1-"$n" ../slurmscripts/1.2_cram_to_bam.sh

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
for tool in samtools; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
bamdir=../data/map_and_call_out_standard/02_bamfiles
reference=../data/map_and_call_out_standard/00_input_data/00_reference_genome/reference_genome.fa.gz

# pick the CRAM file that belongs to this array task (line number N of the list)
cram=$(ls "$bamdir"/*.dedup.cram | sed -n "${SLURM_ARRAY_TASK_ID}p")
bam="$bamdir/$(basename "$cram" .dedup.cram).dedup.bam"

# ---- Run --------------------------------------------------------------------
echo "Task $SLURM_ARRAY_TASK_ID: converting $cram -> $bam"
$run samtools view -b -o "$bam" -T "$reference" "$cram"
$run samtools index "$bam"
