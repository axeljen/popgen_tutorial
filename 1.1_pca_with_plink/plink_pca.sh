
prefix=${1}

plink2 \
    --pfile ${prefix} \
    --extract ${prefix}.prune.in \
    --allow-extra-chr \
    --read-freq ${prefix}.afreq \
    --pca 20 \
    --out ${prefix}_pca