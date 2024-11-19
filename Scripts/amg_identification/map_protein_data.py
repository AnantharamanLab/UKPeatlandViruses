#!/usr/bin/env python3

import os
import gc
from pathlib import Path
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from parse_faa import parse_faa_file

def process_data(tsv, single_contig_prots, vmag_prots, cluster_files, ranks):
    """
    Obtains gene-level data from prodigal-formatted input faa(s)

    tsv: input tsv of mmseqs protein clusters.
    input_files: path(s) to amino-acid fasta files to be parsed.
    cluster_files: path to the folder containing the mmseqs
    protein cluster files.
    
    Returns a Polars DataFrame with gene- and genome-level
    data mapped to each protein in the mmseqs tsv.
    """

    # Accumulate all faa data in chunks to avoid large memory usage
    faa_data = []

    if vmag_prots:
        for file in vmag_prots:
            data = parse_faa_file(file, is_vMAG=True, vMAG_name=os.path.basename(file).rsplit('.', 1)[0])
            faa_data.extend(data)
            
            # Explicitly trigger garbage collection
            gc.collect()
    
    if single_contig_prots:
        data = parse_faa_file(single_contig_prots, is_vMAG=False, vMAG_name=None)
        faa_data.extend(data)
        
        # Explicitly trigger garbage collection
        gc.collect()

    # Convert the list of dictionaries (faa_data) to a Polars DataFrame
    faa_dataframe = pl.DataFrame(faa_data)

    if cluster_files and tsv and ranks:
        # Read the tsv and other data using Polars
        mmseqs = pl.read_csv(tsv, separator='\t', infer_schema_length=0)
        mmseqs.head()

        # Process cluster files from each rank directory
        all_cluster_data = []
        for rank in ranks:
            # Optimize data types
            mmseqs = mmseqs.with_columns([
                mmseqs[f"{rank}_genome_cluster"].cast(pl.Categorical),
                mmseqs[f"{rank}_protein_cluster"].cast(pl.Categorical)
            ])
            
            rank_dir = os.path.join(cluster_files, rank)
            rank_files = []
            for root, dirs, files in os.walk(rank_dir):
                for file in files:
                    if file.endswith('.tsv'):
                        rank_files.append(os.path.join(root, file))
            
            for file in rank_files:
                cluster_id = os.path.basename(file).split('.')[0].split('_')[-1]
                genome_cluster = str(cluster_id)
                cluster_df = pl.read_csv(file, separator='\t', infer_schema_length=0, new_columns=["protein_cluster_rep", "protein_internal"])
                lookup_path = os.path.join(Path(file).parent, f"{os.path.basename(Path(file).parent)}.db.lookup")
                lookup_df = pl.read_csv(lookup_path, separator='\t', infer_schema_length=0, new_columns=["protein_internal", "protein", "_"]).drop(["_"])
                cluster_df = cluster_df.join(lookup_df, on="protein_internal", how="inner").drop(["protein_internal"])
                cluster_df = cluster_df.with_columns([
                    pl.lit(genome_cluster).alias("genome_cluster"),
                    pl.lit(rank).alias("rank")
                ])
                all_cluster_data.append(cluster_df)

        # Concatenate all cluster data
        all_clusters_df = pl.concat(all_cluster_data)

        # Merge with mmseqs and faa_dataframe
        merged_df = mmseqs.join(all_clusters_df, left_on = "member", right_on="protein", how="left")
        merged_df = merged_df.join(faa_dataframe, left_on="member", right_on="member", how="left")
        
        return merged_df
    else:
        return faa_dataframe

def main():
    # Snakemake integration
    tsv = snakemake.params.tsv
    db_files = snakemake.params.db
    all_genes = snakemake.params.database
    ranks = snakemake.params.cluster_taxa_levels
    out_parent = snakemake.params.out_parent
    num_cpus = snakemake.threads
    
    protein_dir = snakemake.params.protein_dir
    vmag_proteins_subdir = snakemake.params.vmag_proteins_subdir
    
    if not os.path.exists(out_parent):
        os.makedirs(out_parent, exist_ok=True)
    
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        single_contig_prots = os.path.join(protein_dir, 'single_contig_proteins.faa')
    else:
        single_contig_prots = None
    
    if os.path.exists(vmag_proteins_subdir) and os.path.isdir(vmag_proteins_subdir):
        vmag_prots = [os.path.join(vmag_proteins_subdir, f) for f in os.listdir(vmag_proteins_subdir) if f.endswith('.faa')]
    else:
        vmag_prots = None
    
    processed_data = process_data(tsv, single_contig_prots, vmag_prots, db_files, ranks)
    processed_data.write_csv(all_genes, separator="\t")

if __name__ == "__main__":
    main()