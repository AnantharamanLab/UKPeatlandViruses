#!/bin/sh
# run_checkm.sh

THREADS="64"
SAMPLES="metagenome_samples.txt"

mkdir -p checkm

while IFS= read -r SAMPLE || [ -n "$SAMPLE" ]
do
    # Skip empty lines
    if [ -z "$SAMPLE" ]; then
        continue
    fi

    # Find the path to the metaBAT2 folder for the current sample
    # Using 'find' with wildcards '*' because folder names can vary
    MAG_DIR=$(find metabat2_bins -name ${SAMPLE}*.metabat-bins*)

    echo "Running CheckM on sample $SAMPLE"

    # Command to run for each sample
    checkm lineage_wf \
        --extension fa \
        --threads "$THREADS" \
        --file "checkm/${SAMPLE}/${SAMPLE}_quality.txt" \
	    --tab_table \
        "$MAG_DIR" \
        "checkm/${SAMPLE}"
    echo "Done ${SAMPLE}"
done < "$SAMPLES"
