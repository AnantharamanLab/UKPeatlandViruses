#!/bin/bash

# Directory containing the fastq.gz files
input_dir="./01_QC"

# Output CSV file
output_file="seq_depth.csv"

# Header for the CSV file
echo "Sample,Pair,File,Number.of.Reads,Hundred.Millions.Reads" > $output_file

# Loop through each fastq.gz file in the directory
for file in "$input_dir"/*_R1.fastq.gz; do
  # Extract filename without directory
  filename=$(basename "$file")
  
  # Extract Sample and Pair components
  sample=$(echo "$filename" | cut -d'-' -f1)
  pair=$(echo "$filename" | grep -oP '_R\d+' | cut -d'_' -f2)
  
  # Count the number of reads (assuming one read per line in the file)
  num_reads=$(zcat "$file" | echo $((`wc -l`/4)))
  
  # Calculate Hundred.Millions.Reads
  hundred_millions_reads=$(echo "scale=8; $num_reads / 100000000" | bc)
  
  # Append the results to the output file
  echo "$sample,$pair,$filename,$num_reads,$hundred_millions_reads" >> $output_file
done
