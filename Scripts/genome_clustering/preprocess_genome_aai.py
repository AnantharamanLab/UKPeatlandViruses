#!/usr/bin/env python3

import os
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

import parse_faa

def gene_to_genome(input_faa_files):
    """
    Obtains gene-to-genome mapping from prodigal-formatted input
    faa(s)

    input_files: path(s) to amino-acid fasta files to be parsed.
    
    Returns a Polars DataFrame with genes mapped to the genomes
    that encode them.
    """

    # Accumulate all faa data
    faa_data = []

    # Assumes the first file is not vMAG and subsequent files are
    for i, file in enumerate(input_faa_files):
        is_vMAG = i != 0  # First file is not vMAG, others are
        data = parse_faa.parse_faa_file(file, is_vMAG=is_vMAG, vMAG_name=os.path.basename(file).rsplit('.', 1)[0])
        faa_data.extend(data)  # Aggregate data

    # Convert the list of dictionaries (faa_data) to a Polars DataFrame
    faa_dataframe = pl.DataFrame(faa_data)
    
    # Subset the Polars DataFrame to create the gene-to-genome table
    gene_to_genome = faa_dataframe.select(pl.col(["member","genome"]))
    gene_to_genome = gene_to_genome.rename({"member": "ptn"})

    return gene_to_genome

def main():
    # Snakemake integration
    input_files = snakemake.params.input_files.split()
    wdir = snakemake.params.wdir
    out_g2g = snakemake.output.gene_map_file
    num_cpus = snakemake.threads
    
    gene_map = gene_to_genome(input_files)
    os.makedirs(wdir, exist_ok=True)
    gene_map.write_csv(out_g2g, separator="\t")
    
if __name__ == "__main__":
    main()