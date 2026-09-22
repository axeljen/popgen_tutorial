import sys
import gzip

# Thin a Beagle-format genotype likelihood file to reduce linkage disequilibrium.
# Usage: python thin_beagle.py input.beagle.gz min_distance > output.beagle.gz
# where `min_distance` is the minimum distance between retained sites – i.e., keep only sites that are at least `min_distance` bases apart.

input_file = sys.argv[1]
min_distance = int(sys.argv[2])

outfile = sys.stdout

with gzip.open(input_file, 'rt') as f:
    header = next(f)
    outfile.write(header)
    chrom = None
    last_chrom = None
    last_pos = -min_distance
    for i, line in enumerate(f):
        # get current chrom and pos
        chrom = "_".join(line.split()[0].split("_")[0:-1])
        pos = int(line.split()[0].split("_")[-1])
        if chrom != last_chrom or pos - last_pos >= min_distance:
            outfile.write(line)
            last_chrom = chrom
            last_pos = pos

