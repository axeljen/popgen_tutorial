#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=astral
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00               # hh:mm:ss
#SBATCH --output=logs/%x-%j.out       # %x = job name, %j = job ID
#SBATCH --error=logs/%x-%j.err

# Collapse poorly supported branches in the window trees and summarize them
# into a population tree with ASTRAL-IV (README sections "Estimate a tree for
# each window" and "Combine the window trees with ASTRAL-IV"). Requires the
# window trees from 3_iqtree_array.sh.
#
# Submit from the 3_phylogeny directory:
#   cd 3_phylogeny
#   mkdir -p logs
#   sbatch ../slurmscripts/3_collapse_astral.sh

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
for tool in nw_ed astral4; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Run --------------------------------------------------------------------
# collapse internal branches with bootstrap support below 20.
# With many windows this loop calls nw_ed thousands of times. Starting the
# container takes a fraction of a second each time, so we start it only once
# and run the whole loop inside it: everything between "<<'EOF'" and "EOF" is
# passed to bash, which runs inside the container if you use one.
$run bash <<'EOF'
set -euo pipefail
for tree in trees/*.treefile; do
    nw_ed "$tree" 'i & b < 20' o > "${tree%.treefile}.collapsed.nwk"
done
EOF

cat trees/*.collapsed.nwk > window_trees.nwk
echo "Number of window trees: $(wc -l < window_trees.nwk)"

# population tree with ASTRAL-IV
$run astral4 \
    -t "$SLURM_CPUS_PER_TASK" \
    -o population_tree.nwk \
    window_trees.nwk

cat population_tree.nwk
