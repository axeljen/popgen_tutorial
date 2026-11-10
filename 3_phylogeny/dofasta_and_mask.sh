#!/bin/bash

# directories for storing inverse (negative) masks and consensus sequences
mkdir -p masks
mkdir -p consensus

# path to reference index
ref=../data/map_and_call_out_standard/00_input_data/00_reference_genome/reference_genome.fa.gz


# bam file for the sample
bam=$1
# mask file for the sample
mask=$2

# extract the sample name from the BAM file path (e.g. ind10.dedup.bam -> ind10)
sample=$(basename "$bam" .dedup.bam)


# first, we'll run bedtools complement to create the exclusion intervals for masking
bedtools complement \
    -i "$mask" \
    -g "$ref.fai" > "masks/${sample}.exclude.bed"

# if the bam isn't indexed, create an index
if [ ! -f "$bam.bai" ]; then
    samtools index "$bam"
fi

# next run angsd to generate a genome-wide consensus sequence from the mapped BAM file
angsd \
    -i "$bam" \
    -ref "$ref" \
    -doCounts 1 \
    -doFasta 1 \
    -minMapQ 30 \
    -minQ 20 \
    -out "consensus/$sample"

# bedtools maskfasta can't read gzipped fasta files, so decompress the angsd output first
gunzip -c "consensus/${sample}.fa.gz" > "consensus/${sample}.fa"

# apply the negative mask to the consensus sequence
bedtools maskfasta \
    -fi "consensus/${sample}.fa" \
    -bed "masks/${sample}.exclude.bed" \
    -fo "consensus/${sample}.masked.fa"

# index the finished masked consensus sequence
samtools faidx "consensus/${sample}.masked.fa"

# and remove the intermediate files to keep the workspace clean
rm "consensus/${sample}.fa.gz" "consensus/${sample}.fa" "masks/${sample}.exclude.bed"