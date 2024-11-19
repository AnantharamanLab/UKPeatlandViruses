#!/usr/bin/env python3

import os
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def prefix_columns(dataframe, prefix):
    # Select 'sequence' only once, and prefix other columns as needed
    cols_to_select = [pl.col("sequence")]
    cols_to_select.extend(
        pl.col(col).alias(f"{prefix}_{col}") for col in dataframe.columns if col not in ["sequence", "db"]
    )
    return dataframe.select(cols_to_select)

def widen_hmm_results(hmm_results_df):
    # Create a new column for each combination of `db` and the associated scores/ids/V-scores
    hmm_results_df = hmm_results_df.with_columns([
        pl.when(pl.col("db") == db).then(pl.col(col)).otherwise(None).alias(f"{db}_{col}")
        for db in hmm_results_df.select("db").unique().to_series().to_list()
        for col in ["score", "hmm_id"]
    ])

    # Remove the original score, hmm_id, V-score, and db columns as they are now redundant
    hmm_results_df = hmm_results_df.drop(["score", "hmm_id", "db"])

    # Pivot the DataFrame to widen it, aggregating by sequence and filling missing values as needed
    hmm_results_df_wide = hmm_results_df.group_by("sequence").agg([
        pl.max(col).alias(col) for col in hmm_results_df.columns if col != "sequence"
    ])

    return hmm_results_df_wide

def assign_db(db_path):
    if "KEGG" in db_path or "kegg" in db_path or "kofam" in db_path:
        return "KEGG"
    elif "Pfam" in db_path or "pfam" in db_path:
        return "Pfam"
    elif "dbcan" in db_path or "dbCAN" in db_path or "dbCan" in db_path:
        return "dbCAN"
    elif "METABOLIC_custom" in db_path or "metabolic_custom" in db_path:
        return "METABOLIC"
    elif "VOG" in db_path or "vog" in db_path:
        return "VOG"
    elif "eggNOG" in db_path or "eggnog" in db_path:
        return "eggNOG"
    elif "PHROG" in db_path or "phrog" in db_path:
        return "PHROG"
    else:
        return None 
    
def main():
    # Integration with Snakemake
    tsv = snakemake.params.gene_index
    vscores = snakemake.input.vscores
    all_hmm_results = snakemake.input.all_hmm_results
    db_dir = snakemake.params.db_dir
    output = snakemake.params.gene_index_annotated
    num_cpus = snakemake.threads
    
    logging.info("Processing of V-scores/L-scores and HMM results starting...")
    
    # Load input dataframes
    logging.info(f"Loading input dataframes {tsv}, {vscores}, and {all_hmm_results}")
    tsv_df = pl.read_csv(tsv, separator="\t")
    vscores_df = pl.read_csv(vscores, separator="\t")
    hmm_df = pl.read_csv(all_hmm_results, separator="\t")
    hmm_df_wide = widen_hmm_results(hmm_df)

    merged_df = tsv_df.join(hmm_df_wide, left_on='member', right_on='sequence', how='left')
    
    cols_to_remove = [
                    'eggNOG_L-score',
                    'KEGG_L-score',
                    'VOG_L-score',
                    'Pfam_L-score',
                    'dbCAN_L-score',
                    'METABOLIC_L-score',
                    'PHROG_L-score',
                    'eggNOG_V-score',
                    'KEGG_V-score',
                    'VOG_V-score',
                    'Pfam_V-score',
                    'dbCAN_V-score',
                    'METABOLIC_V-score',
                    'PHROG_V-score',
                    'ID',
                    'rank', # these three are critical, since there are multiple entries for genes (one per rank)
                    'protein_cluster_rep',
                    'genome_cluster'
                    ]
    cols_to_remove += [col for col in merged_df.columns if col.endswith('_right')]
    for col in cols_to_remove:
        if col in merged_df.columns:
            merged_df = merged_df.drop(col)
    merged_df = merged_df.unique()
    
    # Split the DataFrame by 'db' value
    df_pfam = vscores_df.filter(pl.col("db") == "Pfam")
    df_kegg = vscores_df.filter(pl.col("db") == "KEGG")

    # Add prefixes
    pfam_prefixed = prefix_columns(df_pfam, "Pfam")
    kegg_prefixed = prefix_columns(df_kegg, "KEGG")

    # Join on 'sequence'
    wide_df = pfam_prefixed.join(kegg_prefixed, on="sequence", how="full")
    cols_to_remove = [col for col in wide_df.columns if col.endswith('_right')]
    for col in cols_to_remove:
        if col in wide_df.columns:
            wide_df = wide_df.drop(col)
    
    # Merge V-scores and L-scores with input TSV, and save results
    logging.info(f"Merging V-scores/L-scores and annotations with input {tsv} and writing results to {output}")
    merged_df = merged_df.join(wide_df, left_on='member', right_on='sequence', how='left')
    
    # Remove other columns ending with '_right', resulting from the merges
    cols_to_remove = [col for col in merged_df.columns if (col.endswith('_right') and col.endswith('_score_right') == False)]
    for col in cols_to_remove:
        if col in merged_df.columns:
            merged_df = merged_df.drop(col)
    
    merged_df = merged_df.drop(col for col in merged_df.columns if col.endswith("_right"))        
    
    # Ensure all DBs are represented in columns, even if there were no hits
    for db_path in os.listdir(db_dir):
        db_name = assign_db(db_path)
        col_found = False
        for col in merged_df.columns:
            if db_name in col:
                col_found = True
                break
        if not col_found:
            logging.debug(f"No {db_name} found in the merged DataFrame ({merged_df.columns}); adding columns")
            merged_df = merged_df.with_columns(
                pl.lit(None).alias(f"{db_name}_score"),
                pl.lit(None).alias(f"{db_name}_hmm_id")
                )
    
    merged_df.write_csv(output, separator="\t")
    logging.info("Processing of V-scores/L-scores and HMM results completed.")

if __name__ == "__main__":
    main()