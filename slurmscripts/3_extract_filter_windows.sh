#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=extract_windows
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00               # hh:mm:ss, increase if you have many windows
#SBATCH --output=logs/%x-%j.out       # %x = job name, %j = job ID
#SBATCH --error=logs/%x-%j.err

# Extract one alignment per window from the masked consensus sequences, and
# filter them for missing data with fastaparser (README section "Make masked
# sequences with ANGSD"). Requires windows.bed and the consensus sequences
# from 3_dofasta_and_mask_array.sh.
#
# Submit from the 3_phylogeny directory:
#   cd 3_phylogeny
#   mkdir -p logs
#   sbatch ../slurmscripts/3_extract_filter_windows.sh

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
for tool in samtools fastaparser; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Run --------------------------------------------------------------------
mkdir -p windows_raw windows_filtered logs

# The loops below call samtools and fastaparser thousands of times. Starting the
# container takes a fraction of a second each time, so instead of writing $run
# in front of every call, we start it only once and run the whole block inside
# it: everything between "<<'EOF'" and "EOF" is passed to bash, which runs
# inside the container if you use one.
$run bash <<'EOF'
set -euo pipefail

# 1) extract each window from every sample's consensus sequence
while IFS=$'\t' read -r chrom start end; do
    region="${chrom}:$((start + 1))-${end}"
    > "windows_raw/${chrom}_${start}_${end}.fa"
    for consensus in consensus/*.masked.fa; do
        sample=$(basename "$consensus" .masked.fa)
        samtools faidx "$consensus" "$region" | sed "s/^>.*$/>$sample/" \
            >> "windows_raw/${chrom}_${start}_${end}.fa"
    done
done < windows.bed

# 2) filter each window alignment for missing data
for alignment in windows_raw/*.fa; do
    window=$(basename "$alignment" .fa)
    if fastaparser filter \
        -i "$alignment" \
        --max-missing-seq 0.30 \
        --max-missing-site 0.20 \
        -o "windows_filtered/${window}.fa" 2> "logs/${window}.filter.log"; then
        :
    elif grep -Eq 'No sequences remain after filtering|No sites remain after filtering' "logs/${window}.filter.log"; then
        rm -f "windows_filtered/${window}.fa"
        echo "Discarded empty alignment: $window"
    else
        cat "logs/${window}.filter.log" >&2
        exit 1
    fi
done
EOF

echo "Done. $(ls windows_filtered/*.fa | wc -l) windows passed filtering."
