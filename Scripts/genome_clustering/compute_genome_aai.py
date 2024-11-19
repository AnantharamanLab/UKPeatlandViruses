#!/usr/bin/env python3

import os
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

from pathlib import Path
from typing import Literal

### The following functions are courtesy of Cody Martin from Anantharaman Lab at the University of Wisconsin-Madison (2024)

FilePath = str | Path
ClusterTaxaLevel = Literal["family", "genus", "species"]

def read_gene_map(file: FilePath) -> pl.DataFrame:
    return pl.read_csv(file, separator="\t")

def read_alignments(file, gene_map: pl.DataFrame) -> pl.DataFrame:
    cols = {
        "query": pl.String,
        "target": pl.String,
        "fident": pl.Float32,
        "alnlen": pl.Int16,
        "mismatch": pl.Int16,
        "gapopen": pl.Int16,
        "qstart": pl.Int16,
        "qend": pl.Int16,
        "tstart": pl.Int16,
        "tend": pl.Int16,
        "evalue": pl.Float32,
        "bits": pl.Float32,
        "qlen": pl.Int16,
        "tlen": pl.Int16,
        "qcov": pl.Float32,
        "tcov": pl.Float32,
    }

    gene_map = gene_map.lazy()

    aln_df = pl.read_csv(
        file,
        separator='\t',
        has_header=False,
        new_columns=list(cols.keys()),  # Define new column names
        # dtypes=list(cols.values())  # Define data types for the columns
        schema_overrides=list(cols.values())  # Define data types for the columns
    )

    # Now, let's prepare gene_map for the joins
    # Ensure that 'ptn' and 'genome' columns are named properly for joining
    gene_map = gene_map.rename({"ptn": "query_ptn", "genome": "query_genome"})  # Prepare for first join
    gene_map_target = gene_map.rename({"query_ptn": "target_ptn", "query_genome": "target_genome"})  # Prepare for second join

    # Using lazy evaluation for efficient processing
    aln_df = aln_df.lazy() \
        .drop(["mismatch", "gapopen"])  # Drop these columns as they are not needed
    
    # Joining gene_map to alignment dataframe for 'query' and 'target'
    aln_df = aln_df \
        .join(gene_map, left_on="query", right_on="query_ptn") \
        .rename({"query_genome": "query_genome"})  # Rename to avoid confusion, though not strictly necessary
    
    aln_df = aln_df \
        .join(gene_map_target, left_on="target", right_on="target_ptn") \
        .rename({"target_genome": "target_genome"})  # Rename to use in filter

    # Apply filter and cast categorical transformations
    aln_df = aln_df \
        .filter(pl.col("query_genome") != pl.col("target_genome")) \
        .with_columns(
            [pl.col(column_name).cast(pl.Categorical) for column_name in ["query", "target", "query_genome", "target_genome"]]
        ) \
        .collect()  # Collect results to realize the computation
        
    return aln_df

def preprocess(aln_df: pl.DataFrame) -> pl.DataFrame:
    # for each pair of genomes, only want to keep the most significant matches per proteins
    # ie suppose we have genomes A and B with 3 ptns each
    # we could have ptn A1 with matches to B1 AND B2
    # but suppose that A1-B2 is the more significant match -- that's the one we want to keep
    # basically AAI is only calculated for the most significant RECIPROCAL matches
    # ie ptn A1 is most related to B2 out of all 3 ptns in genome B AND
    # ptn B2 is most related to A1 out of all 3 ptns in genome A
    processed_df = (
        aln_df.sort("bits", descending=True)
        .lazy()
        .unique(subset=["query", "query_genome", "target_genome"], keep="first")
        .unique(subset=["target", "query_genome", "target_genome"], keep="first")
        .collect()
    )
    return processed_df

def aai(aln_df: pl.DataFrame, gene_map: pl.DataFrame) -> pl.DataFrame:
    ptn_counts = gene_map.group_by("genome").agg(n_ptns=pl.len()).lazy().with_columns([pl.col("genome").cast(pl.Categorical)])

    aai_df = (
        aln_df.lazy()
        .group_by("query_genome", "target_genome")
        .agg(
            aai=pl.mean("fident"),
            shared_genes=pl.len(),
        )
        .join(ptn_counts, left_on="query_genome", right_on="genome")
        .rename({"n_ptns": "query_n_ptns"})
        .join(ptn_counts, left_on="target_genome", right_on="genome")
        .rename({"n_ptns": "target_n_ptns"})
        .with_columns(
            query_shared=(pl.col("shared_genes") / pl.col("query_n_ptns")).cast(
                pl.Float32
            ),
            target_shared=(pl.col("shared_genes") / pl.col("target_n_ptns")).cast(
                pl.Float32
            ),
        )
        .collect()
    )

    return aai_df

###

def main():
    # Snakemake integration
    alignments_file = snakemake.input.search_output
    gene_map_file = snakemake.input.gene_map_file
    aai_file = snakemake.output.aai
    
    gene_map = read_gene_map(file=gene_map_file)
    aln_df = read_alignments(file=alignments_file, gene_map=gene_map)
    processed_df = preprocess(aln_df=aln_df)
    aai_df = aai(aln_df=processed_df, gene_map=gene_map)
    
    aai_df.write_csv(aai_file, separator='\t')

if __name__ == "__main__":
    main()