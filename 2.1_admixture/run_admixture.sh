#!/bin/bash
input=$1

# K
K=$2

# go to the directory of the input
cd "$(dirname "$input")"

# change input to the basename only
input=$(basename "$input")

# run ADMIXTURE and capture the standard output a log file
admixture --cv=10 --seed=42 "$input" "$K" > "${input/.bed/}"."$K".log 2>&1

# return to the original directory
cd -