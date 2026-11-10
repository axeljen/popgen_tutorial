#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=ngsadmix
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00               # hh:mm:ss, per value of K
#SBATCH --array=2-6                   # one task per value of K
#SBATCH --output=logs/%x-%A_%a.out    # %A = array job ID, %a = task number (= K)
#SBATCH --error=logs/%x-%A_%a.err

# Run NGSadmix for K = 2 to 6 (README section "With NGSadmix"). Requires the
# thinned likelihood file from 2.2_angsd_ngsadmix_input.sh.
#
# This is an ARRAY job: "#SBATCH --array=2-6" starts five tasks numbered
# 2 to 6, and we use each task's number ($SLURM_ARRAY_TASK_ID) as K.
#
# Submit from the 2.2_ngsadmix directory:
#   cd 2.2_ngsadmix
#   mkdir -p logs
#   sbatch ../slurmscripts/2.2_ngsadmix_array.sh

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

# check that all programs are found before we start (NGSadmix comes with angsd)
for tool in NGSadmix; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
likes=angsd_ngsadmix.thinned5kb.beagle.gz
K=$SLURM_ARRAY_TASK_ID

# ---- Run --------------------------------------------------------------------
echo "Running NGSadmix on $likes with K=$K"
$run NGSadmix \
    -likes "$likes" \
    -K "$K" \
    -minMaf 0.05 \
    -seed 42 \
    -P "$SLURM_CPUS_PER_TASK" \
    -o "ngsadmix_K${K}"
