#!/bin/bash
# Check that all software for the popgen tutorial is available and starts.
# Exits with an error if anything is missing.

failed=0

check() {
    # $1 = name to print, the rest = command to run
    local name=$1
    shift
    if "$@" > /dev/null 2>&1; then
        echo "OK       $name"
    else
        echo "MISSING  $name"
        failed=1
    fi
}

for tool in plink2 admixture angsd NGSadmix pcangsd samtools bgzip tabix bedtools \
            iqtree2 nw_ed astral4 Rscript python fastaparser tree; do
    check "$tool" command -v "$tool"
done

# a few deeper checks that the programs actually run
check "plink2 runs"        plink2 --version
check "samtools runs"      samtools --version
check "bedtools runs"      bedtools --version
check "iqtree2 runs"       iqtree2 --version
check "pcangsd runs"       pcangsd --help
check "fastaparser runs"   fastaparser --help
check "R tidyverse loads"  Rscript -e 'suppressPackageStartupMessages(library(tidyverse))'

if [ "$failed" -ne 0 ]; then
    echo "Some software is missing or broken." >&2
    exit 1
fi
echo "All software found."
