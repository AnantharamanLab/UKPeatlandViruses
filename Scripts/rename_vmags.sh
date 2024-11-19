#!/bin/bash

# Define the destination directory
destination_dir="vMAGs"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Loop through each subdirectory in the ViWrap_results directory
for sample_dir in ./ViWrap_results/*/08_ViWrap_summary_outdir/Virus_genomes_files/; do
  # Get the sample name by extracting the directory name directly under ViWrap_results
  sample_name=$(basename $(dirname $(dirname "$sample_dir")))

  # Loop through each file in the current sample directory
  for file in "$sample_dir"/*; do
    # Get the filename without the directory path
    filename=$(basename "$file")

    # Construct the new filename with the sample name prefix
    new_filename="${sample_name}_${filename}"

    # Copy the file to the destination directory with the new name
    cp "$file" "$destination_dir/$new_filename"
  done
done

echo "Files copied and renamed successfully."