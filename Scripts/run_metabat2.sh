#!/bin/sh
# run_metabat2.sh

THREADS="64"
SAMPLES="metagenome_samples.txt"

mkdir -p metabat2_bins && cd metabat2_bins

while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]
do
    # Skip empty lines
    if [ -z "$SAMPLE" ]; then
        continue
    fi

    echo "Binning sample: $SAMPLE"

    # Command to run for each sample
    runMetaBat.sh -m 2500 --maxP 95 --minS 60 --minCV 1 --minClsSize 200000 --numThreads $THREADS ../fasta/${SAMPLE}-contigs-prefix-formatted-only.fasta ../04_MAPPING/${SAMPLE}/*.bam
    echo "Done ${SAMPLE}"
done < ../"$SAMPLES"

