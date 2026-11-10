#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=admixture
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00               # hh:mm:ss, per value of K
#SBATCH --array=2-6                   # one task per value of K
#SBATCH --output=logs/%x-%A_%a.out    # %A = array job ID, %a = task number (= K)
#SBATCH --error=logs/%x-%A_%a.err

# Run ADMIXTURE for K = 2 to 6 (README section "Run ADMIXTURE for several
# values of K").
#
# This is an ARRAY job: the line "#SBATCH --array=2-6" makes slurm start five
# tasks, numbered 2, 3, 4, 5 and 6. Each task gets its number in
# $SLURM_ARRAY_TASK_ID, which we simply use as K. All five values of K then
# run at the same time instead of one after another.
#
# Submit from the 2.1_admixture directory, after preparing the input with
# prepare_admix_input.sh:
#   cd 2.1_admixture
#   mkdir -p logs
#   sbatch ../slurmscripts/2.1_admixture_array.sh

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
for tool in admixture; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
# change this to admix_bp-space_5k/admix_bp-space_5k.bed for the thinned data set
input=admix_stdfilt/admix_stdfilt.bed
K=$SLURM_ARRAY_TASK_ID

# ---- Run --------------------------------------------------------------------
echo "Running ADMIXTURE on $input with K=$K"
$run bash run_admixture.sh "$input" "$K"
