

# take the vcf file as first argument
vcf=$1

# and a prefix name for the plink files as second argument
outname=$2

## filtering parameters
missing_ind=0.10 # minimum individual missingness threshold, i.e., individuals with more than 10% missing genotypes will be excluded
missing_snp=0.10 # minimum SNP missingness threshold, i.e., SNPs with more than 10% missing genotypes will be excluded
minor_allele_freq=0.05 # minimum minor allele frequency threshold, i.e., SNPs with a minor allele frequency below 5% will be excluded
linkage_pruning_threshold=0.2 # threshold for linkage disequilibrium pruning, i.e., SNPs with an r^2 value above this threshold will be pruned
plink2 \
    --vcf $vcf \
    --allow-extra-chr \
    --set-all-var-ids '@:#:$r:$a' \
    --snps-only just-acgt \
    --max-alleles 2 \
    --mind $missing_ind \
    --geno $missing_snp \
    --maf $minor_allele_freq \
    --make-pgen \
    --indep-pairwise 50 5 $linkage_pruning_threshold \
    --bad-ld \
    --rm-dup exclude-all \
    --freq \
    --out $outname