#!/bin/bash
#SBATCH --account=<SLURM_ACCOUNT>     # your compute project / allocation
#SBATCH --partition=<PARTITION>       # e.g. "shared" on Dardel
#SBATCH --job-name=angsd_ngsadmix
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00               # hh:mm:ss, increase for big data sets
#SBATCH --output=logs/%x-%j.out       # %x = job name, %j = job ID
#SBATCH --error=logs/%x-%j.err

# Calculate genotype likelihoods with ANGSD and thin them to one site per 5 kb,
# as input for NGSadmix (README section "With NGSadmix").
#
# Submit from the 2.2_ngsadmix directory:
#   cd 2.2_ngsadmix
#   mkdir -p logs
#   sbatch ../slurmscripts/2.2_angsd_ngsadmix_input.sh

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
for tool in angsd python; do
    if ! $run bash -c "command -v $tool" > /dev/null; then
        echo "ERROR: '$tool' was not found. Check the software settings at the top of this script." >&2
        exit 1
    fi
done

# ---- Input and output -------------------------------------------------------
bamdir=../data/map_and_call_out_standard/02_bamfiles
popfile=../data/sample_popfile.tsv
min_distance=5000   # minimum distance (bp) between retained sites

# ---- Run --------------------------------------------------------------------
# make the BAM list in the same order as the population file
> bamlist.txt
while read -r sample pop; do
    echo "$bamdir/$sample.dedup.bam" >> bamlist.txt
done < "$popfile"

# genotype likelihoods in beagle format
# $SLURM_CPUS_PER_TASK is the number of CPUs requested above
$run angsd \
    -b bamlist.txt \
    -out angsd_ngsadmix \
    -GL 1 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -minMapQ 30 \
    -minQ 20 \
    -SNP_pval 1e-6 \
    -nThreads "$SLURM_CPUS_PER_TASK"

# thin the likelihoods to reduce linkage between sites
$run python ../bin/thin_beagle.py angsd_ngsadmix.beagle.gz "$min_distance" | gzip > angsd_ngsadmix.thinned5kb.beagle.gz

echo "Sites before thinning: $(zcat angsd_ngsadmix.beagle.gz | wc -l)"
echo "Sites after thinning:  $(zcat angsd_ngsadmix.thinned5kb.beagle.gz | wc -l)"
