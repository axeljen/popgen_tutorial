
prefix=${1}

plink2 \
    --pfile ${prefix} \
    --allow-extra-chr \
    --read-freq ${prefix}.afreq \
    --pca 20 \
    --out ${prefix}_pca