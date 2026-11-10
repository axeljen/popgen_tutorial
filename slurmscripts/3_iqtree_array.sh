#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=iqtree
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00               # hh:mm:ss, per window
#SBATCH --output=logs/%x-%A_%a.out    # %A = array job ID, %a = task number
#SBATCH --error=logs/%x-%A_%a.err

# Estimate one tree per window with IQ-TREE 2 (README section "Estimate a
# tree for each window"). Requires the filtered alignments from
# 3_extract_filter_windows.sh.
#
# This is an ARRAY job with one task per window: each task uses its number
# ($SLURM_ARRAY_TASK_ID) to pick "the Nth alignment" in windows_filtered/.
#
# Submit from the 3_phylogeny directory, with one task per alignment:
#   cd 3_phylogeny
#   mkdir -p logs
#   n=$(ls windows_filtered/*.fa | wc -l)
#   sbatch --array=1-"$n" ../slurmscripts/3_iqtree_array.sh
#
# Tip: with thousands of windows, add a limit on how many tasks run at the same
# time, e.g. --array=1-"$n"%50 runs at most 50 at once.

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
for tool in iqtree2; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
# pick the alignment that belongs to this array task (line number N of the list)
alignment=$(ls windows_filtered/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
window=$(basename "$alignment" .fa)

# ---- Run --------------------------------------------------------------------
mkdir -p trees

echo "Task $SLURM_ARRAY_TASK_ID: window $window"
$run iqtree2 \
    -s "$alignment" \
    -m MFP \
    -B 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix "trees/$window"
