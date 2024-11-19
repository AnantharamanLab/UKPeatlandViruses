#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
from typing import Literal
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

FilePath = str | Path
ClusterTaxaLevel = Literal["class", "family", "genus", "species"]

def read_mcl_clusters(file: Path) -> pl.DataFrame:
    clusters = []
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            for genome in line.strip().split('\t'):
                clusters.append((genome, i + 1))
    
    return pl.DataFrame(clusters, schema=["genome", "cluster"])

def mcl(
    mcl_abc_file: FilePath,
    inflation: float,
    outdir: Path,
    mcl_out_file: FilePath,
    threads: int,
):
    mcl_bin_file = os.path.join(outdir, "mcl.mci")
    mcl_tab_file = os.path.join(outdir, "mcl.mcxtab")
    mcl_log_file = os.path.join(outdir, "mcl.log")

    mcxload_cmd = (
        f"mcxload -abc {mcl_abc_file} -o {mcl_bin_file} -write-tab {mcl_tab_file}"
    )

    mcl_cmd = f"mcl {mcl_bin_file} -I {inflation} -te {threads} -o {mcl_out_file} -use-tab {mcl_tab_file}"
    with Path(mcl_log_file).open("w") as log:
        subprocess.run(mcxload_cmd.split(), stdout=log, stderr=log)
        subprocess.run(mcl_cmd.split(), stdout=log, stderr=log)

def cluster(
    aai_df: pl.DataFrame,
    level: ClusterTaxaLevel,
    outdir: Path,
    threads: int,
    mcl_out_file: FilePath,
    inflation: float = 2.0,
):
    if level == "class":
        min_genes = 1
        min_shared = 0.05
        min_aai = 0.1
    elif level == "family":
        min_genes = 8
        min_shared = 0.35
        min_aai = 0.35
    elif level == "genus":
        min_genes = 16
        min_shared = 0.20
        min_aai = 0.40
    elif level == "species":
        min_genes = 20
        min_shared = 0.50
        min_aai = 0.50
    else:
        raise ValueError(f"Invalid taxa level: {level}")

    mcl_abc_file = os.path.join(outdir, f"mcl_{level}.abc")
    filtered_df = (
        aai_df.filter(
            (
                (pl.col("shared_genes") >= min_genes)
                | (
                    (pl.col("query_shared") >= min_shared)
                    & (pl.col("target_shared") >= min_shared)
                )
            )
            & (pl.col("aai") >= min_aai)
        )
        .with_columns(
            score=pl.col("aai") * pl.min_horizontal("query_shared", "target_shared")
        )
        .select("query_genome", "target_genome", "score")
    )

    filtered_df.write_csv(mcl_abc_file, include_header=False, separator="\t")
    mcl(
        mcl_abc_file=mcl_abc_file,
        inflation=inflation,
        outdir=outdir,
        threads=threads,
        mcl_out_file=mcl_out_file
    )

def hierarchical_clustering(aai_df: pl.DataFrame, all_genomes_list_df: pl.DataFrame, outdir: Path, outfile: Path, threads: int, cluster_ranks: list):
    levels = {
        # class-level inflation factor is arbitrary
        "class": 0.8,
        # family and genus-level inflation factors are derived from
        # Nayfach et al. (2021) Nat. Microbiol. https://doi.org/10.1038/s41564-021-00928-6
        "family": 1.2,
        "genus": 2.0,
        
        # species-level inflation factor is derived from
        # Roux et al. (2016) eLife https://doi.org/10.7554/eLife.08490 and
        # Lima-Mendez et al. (2008) Molecular Biology and Evolution  https://doi.org/10.1093/molbev/msn023
        "species": 4.0 # Roux et al. value
        # "species": 6.0 # Lima-Mendez et al. value
    }

    def get_cluster_assignments(aai_df, level, inflation, outdir, threads):
        mcl_out_file = os.path.join(outdir, f"genome_clusters_{level}.tsv")
        cluster(
            aai_df=aai_df,
            level=level,
            outdir=outdir,
            threads=threads,
            mcl_out_file=mcl_out_file,
            inflation=inflation
        )
        return read_mcl_clusters(mcl_out_file).rename({"cluster": f"{level}_cluster"})

    rank_dfs = {}
    max_cluster_ids = {level: 0 for level in cluster_ranks}
    cluster_names = {level: set() for level in cluster_ranks}

    for level in cluster_ranks:
        rank_df = get_cluster_assignments(aai_df, level, levels[level], outdir, threads)
        if not rank_df.is_empty():
            max_cluster_num = rank_df[f"{level}_cluster"].max()
            if max_cluster_num is not None:
                rank_df = rank_df.with_columns(
                    (pl.col(f"{level}_cluster") + max_cluster_ids[level])
                    .cast(pl.Utf8)
                    .map_elements(lambda x: f"{level}_{x}")
                    .alias(f"{level}_cluster")
                )
                max_cluster_ids[level] += max_cluster_num
                cluster_names[level].update(rank_df[f"{level}_cluster"].unique().to_list())

        rank_dfs[level] = rank_df

    # Perform an outer join on all rank dataframes
    final_df = rank_dfs[cluster_ranks[0]]
    for level in cluster_ranks[1:]:
        df = rank_dfs[level]
        # Ensure the 'genome' column is of type 'str' even if it contains null values
        df = df.with_columns(
            pl.col("genome").cast(pl.Utf8).alias("genome")
        )
        final_df = final_df.join(df, on="genome", how="outer")
        final_df = final_df.select(["genome"] + [col for col in final_df.columns if col.endswith("_cluster")])
    
    # Add any remaining unclustered input genomes that didn't cluster at any level
    existing_genomes = final_df.select("genome")
    missing_genomes = all_genomes_list_df.join(existing_genomes, on="genome", how="anti")

    # Create a DataFrame with missing genomes and null values for other columns
    missing_genomes_dict = {col: [None] * len(missing_genomes) for col in final_df.columns if col != "genome"}
    missing_genomes_dict["genome"] = missing_genomes["genome"]
    missing_genomes_df = pl.DataFrame(missing_genomes_dict)

    # Concatenate the missing genomes DataFrame to the existing DataFrame
    final_df = pl.concat([final_df, missing_genomes_df], how = "diagonal")
    
    # Fill null values in the cluster columns (i.e. singleton genomes) with pseudo-cluster assignments
    for level in cluster_ranks:
        max_cluster_num = max_cluster_ids[level]
        def fill_cluster(col_name):
            nonlocal max_cluster_num
            new_col_values = []
            for val in final_df[col_name]:
                if val is None:
                    max_cluster_num += 1
                    new_value = f"{level}_{max_cluster_num}"
                    new_col_values.append(new_value)
                else:
                    new_col_values.append(val)
            filled_col = pl.Series(new_col_values)
            return filled_col.alias(col_name)
        final_df = final_df.with_columns(fill_cluster(f"{level}_cluster"))
        
        # Update cluster names and counts after filling clusters
        cluster_names[level].update(final_df[f"{level}_cluster"].unique().to_list())
        max_cluster_ids[level] = max_cluster_num
    
    print(final_df)
        
    # Ensure unique hierarchical names
    def ensure_unique_hierarchical_names(df, levels):
        for i, level in enumerate(levels[1:], 1):
            parent_level = levels[i - 1]
            df = df.with_columns(
                (pl.col(f"{parent_level}_cluster") + "-" + pl.col(f"{level}_cluster")).alias(f"{level}_cluster")
            )
        return df

    final_df = ensure_unique_hierarchical_names(final_df, cluster_ranks)

    # Custom sorting by extracting numeric parts
    def numeric_sort(col):
        return col.str.extract(r'(\d+)$').cast(pl.Int32)
    
    final_df = final_df.sort(
        [numeric_sort(pl.col(f"{rank}_cluster")) for rank in cluster_ranks] + ["genome"], 
        nulls_last=True
    )

    # Add numeric columns for sorting
    for level in cluster_ranks:
        final_df = final_df.with_columns(
            numeric_sort(pl.col(f"{level}_cluster")).alias(f"{level}_numeric")
        )

    # Sort by numeric columns
    final_df = final_df.sort(
        [pl.col(f"{rank}_numeric") for rank in cluster_ranks] + ["genome"]
    )

    # Rename clusters within each rank after sorting
    for level in cluster_ranks[1:]:
        parent_level = cluster_ranks[cluster_ranks.index(level) - 1]
        final_df = final_df.sort([pl.col(f"{parent_level}_numeric"), pl.col(f"{level}_numeric")])
        
        cluster_column = final_df[f"{level}_cluster"]
        unique_clusters = cluster_column.unique().to_list()
        cluster_mapping = {old_name: f"{level}_{i+1}" for i, old_name in enumerate(unique_clusters)}
        final_df = final_df.with_columns(
            pl.col(f"{level}_cluster").replace(cluster_mapping).alias(f"{level}_cluster")
        )

    # Drop numeric columns
    final_df = final_df.drop([f"{level}_numeric" for level in cluster_ranks])

    # Add the final_cluster column using only the ranks that exist
    final_df = final_df.with_columns(
        pl.concat_str([pl.col(f"{rank}_cluster").cast(pl.Utf8) for rank in cluster_ranks], separator="-").alias("final_cluster")
    )
    
    # Rearrange columns to have final_cluster after genome
    final_df = final_df.select(["genome", "final_cluster"] + [col for col in final_df.columns if col not in ["genome", "final_cluster"]])

    # Remove rows with all nulls
    final_df = final_df.filter(~pl.col("genome").is_null())

    # Write final cluster table
    final_df.write_csv(outfile, separator="\t")
    
    # Print the number of clusters for each rank by counting unique values
    for level in cluster_ranks:
        unique_clusters = final_df.select(pl.col(f"{level}_cluster")).unique()
        print(f"Number of clusters for {level}: {unique_clusters.height}")

def main():
    # Snakemake integration
    aai = snakemake.input.aai
    gene_map = snakemake.params.gene_map
    output_dir = snakemake.params.wdir
    outfile = snakemake.output.genome_clusters
    num_cpus = snakemake.threads
    cluster_ranks = snakemake.params.cluster_taxa_levels

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    aai_df = pl.read_csv(aai, separator='\t')
    all_genomes_list_df = pl.read_csv(gene_map, separator='\t').select("genome").unique()

    hierarchical_clustering(
        aai_df=aai_df,
        all_genomes_list_df=all_genomes_list_df,
        outdir=output_dir,
        outfile=outfile,
        threads=num_cpus,
        cluster_ranks=cluster_ranks
    )

if __name__ == "__main__":
    main()