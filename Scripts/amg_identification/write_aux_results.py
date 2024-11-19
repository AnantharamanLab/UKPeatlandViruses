#!/usr/bin/env python3

import os
import load_prot_paths
from Bio import SeqIO
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def summarize_aux_table(input_table, hmm_descriptions):
    """
    Summarizes the auxiliary table by selecting relevant columns,
    filtering by genes with metabolic annotations, merging with HMM descriptions,
    identifying viral flanking genes, and adding other information.
    
    Args:
        input_table (str): Path to the input TSV file containing the auxiliary table.
        hmm_descriptions (DataFrame): Polars DataFrame containing HMM descriptions.
        metabolic_kos (DataFrame): Polars DataFrame containing metabolic KO identifiers.
    
    Returns:
        DataFrame: A summarized Polars DataFrame.
    """
    table = pl.read_csv(input_table, separator="\t")
    table = table.select([
        "member",
        "contig",
        "circular_contig",
        "genome",
        "KEGG_hmm_id",
        "Pfam_hmm_id",
        "dbCAN_hmm_id",
        "METABOLIC_hmm_id",
        "PHROG_hmm_id",
        "KEGG_score",
        "Pfam_score",
        "dbCAN_score",
        "METABOLIC_score",
        "PHROG_score",
        "window_avg_KEGG_L-score_viral",
        "window_avg_Pfam_L-score_viral",
        "KEGG_verified_flank",
        "Pfam_verified_flank"
    ])
    table = table.rename({"member": "Protein"})
    
    # Join with KEGG descriptions
    table = table.join(hmm_descriptions, left_on="KEGG_hmm_id", right_on="id", how="left").rename({"name": "KEGG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # Join with Pfam descriptions
    table = table.join(hmm_descriptions, left_on="Pfam_hmm_id", right_on="id", how="left").rename({"name": "Pfam_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # Custom join with dbCAN descriptions to handle IDs with underscores
    table = table.with_columns(pl.col("dbCAN_hmm_id").str.replace(r'_(.*)', '', literal=False).alias("dbCAN_hmm_id_no_underscore"))
    table = table.join(hmm_descriptions, left_on="dbCAN_hmm_id_no_underscore", right_on="id", how="left").rename({"name": "dbCAN_Description"})
    table = table.drop("dbCAN_hmm_id_no_underscore")
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # Join with METABOLIC descriptions
    table = table.join(hmm_descriptions, left_on="METABOLIC_hmm_id", right_on="id", how="left").rename({"name": "METABOLIC_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
        
    # Join with PHROG descriptions
    table = table.join(hmm_descriptions, left_on="PHROG_hmm_id", right_on="id", how="left").rename({"name": "PHROG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
            
    # Mark genes with viral flanking genes
    table = table.with_columns(
        pl.when(pl.col("KEGG_verified_flank") | pl.col("Pfam_verified_flank"))
        .then(True)
        .otherwise(False)
        .alias("Viral_Flanking_Genes")
    )
    
    # Ensure that all score columns are treated as floats
    table = table.with_columns([
        pl.col("KEGG_score").cast(pl.Float64),
        pl.col("Pfam_score").cast(pl.Float64),
        pl.col("dbCAN_score").cast(pl.Float64),
        pl.col("METABOLIC_score").cast(pl.Float64),
        pl.col("PHROG_score").cast(pl.Float64),
    ])
    
    # Initialize new columns as empty
    table = table.with_columns([
        pl.lit(None).alias("top_hit_hmm_id"),
        pl.lit(None).alias("top_hit_description"),
        pl.lit(None).alias("top_hit_db")
    ])
    
    # Compare scores and update top_hit columns
    table = table.with_columns([
        pl.when(pl.col("KEGG_score").is_not_null() & 
                (pl.col("KEGG_score") >= pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("KEGG_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("KEGG_hmm_id")).otherwise(pl.col("top_hit_hmm_id")).alias("top_hit_hmm_id")
    ])
    table = table.with_columns([
        pl.when(pl.col("KEGG_score").is_not_null() & 
                (pl.col("KEGG_score") >= pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("KEGG_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("KEGG_Description")).otherwise(pl.col("top_hit_description")).alias("top_hit_description")
    ])
    table = table.with_columns([
        pl.when(pl.col("KEGG_score").is_not_null() & 
                (pl.col("KEGG_score") >= pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("KEGG_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("KEGG_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.lit("KEGG")).otherwise(pl.col("top_hit_db")).alias("top_hit_db")
    ])

    table = table.with_columns([
        pl.when(pl.col("Pfam_score").is_not_null() & 
                (pl.col("Pfam_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("Pfam_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("Pfam_hmm_id")).otherwise(pl.col("top_hit_hmm_id")).alias("top_hit_hmm_id")
    ])
    table = table.with_columns([
        pl.when(pl.col("Pfam_score").is_not_null() & 
                (pl.col("Pfam_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("Pfam_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("Pfam_Description")).otherwise(pl.col("top_hit_description")).alias("top_hit_description")
    ])
    table = table.with_columns([
        pl.when(pl.col("Pfam_score").is_not_null() & 
                (pl.col("Pfam_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("dbCAN_score").fill_null(float('-inf'))) & 
                (pl.col("Pfam_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("Pfam_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.lit("Pfam")).otherwise(pl.col("top_hit_db")).alias("top_hit_db")
    ])

    table = table.with_columns([
        pl.when(pl.col("dbCAN_score").is_not_null() & 
                (pl.col("dbCAN_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("dbCAN_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("dbCAN_hmm_id")).otherwise(pl.col("top_hit_hmm_id")).alias("top_hit_hmm_id")
    ])
    table = table.with_columns([
        pl.when(pl.col("dbCAN_score").is_not_null() & 
                (pl.col("dbCAN_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("dbCAN_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("dbCAN_Description")).otherwise(pl.col("top_hit_description")).alias("top_hit_description")
    ])
    table = table.with_columns([
        pl.when(pl.col("dbCAN_score").is_not_null() & 
                (pl.col("dbCAN_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("dbCAN_score") >= pl.col("METABOLIC_score").fill_null(float('-inf'))) &
                (pl.col("dbCAN_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.lit("dbCAN")).otherwise(pl.col("top_hit_db")).alias("top_hit_db")
    ])

    table = table.with_columns([
        pl.when(pl.col("METABOLIC_score").is_not_null() & 
                (pl.col("METABOLIC_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("METABOLIC_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("METABOLIC_hmm_id")).otherwise(pl.col("top_hit_hmm_id")).alias("top_hit_hmm_id")
    ])
    table = table.with_columns([
        pl.when(pl.col("METABOLIC_score").is_not_null() & 
                (pl.col("METABOLIC_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("METABOLIC_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.col("METABOLIC_Description")).otherwise(pl.col("top_hit_description")).alias("top_hit_description")
    ])
    table = table.with_columns([
        pl.when(pl.col("METABOLIC_score").is_not_null() & 
                (pl.col("METABOLIC_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("METABOLIC_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("METABOLIC_score") >= pl.col("PHROG_score").fill_null(float('-inf'))))
        .then(pl.lit("METABOLIC")).otherwise(pl.col("top_hit_db")).alias("top_hit_db")
    ])

    table = table.with_columns([
        pl.when(pl.col("PHROG_score").is_not_null() & 
                (pl.col("PHROG_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("PHROG_score") > pl.col("METABOLIC_score").fill_null(float('-inf'))))
        .then(pl.col("PHROG_hmm_id")).otherwise(pl.col("top_hit_hmm_id")).alias("top_hit_hmm_id")
    ])
    table = table.with_columns([
        pl.when(pl.col("PHROG_score").is_not_null() & 
                (pl.col("PHROG_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("PHROG_score") > pl.col("METABOLIC_score").fill_null(float('-inf'))))
        .then(pl.col("PHROG_Description")).otherwise(pl.col("top_hit_description")).alias("top_hit_description")
    ])
    table = table.with_columns([
        pl.when(pl.col("PHROG_score").is_not_null() & 
                (pl.col("PHROG_score") > pl.col("KEGG_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("Pfam_score").fill_null(float('-inf'))) & 
                (pl.col("PHROG_score") > pl.col("dbCAN_score").fill_null(float('-inf'))) &
                (pl.col("PHROG_score") > pl.col("METABOLIC_score").fill_null(float('-inf'))))
        .then(pl.lit("PHROG")).otherwise(pl.col("top_hit_db")).alias("top_hit_db")
    ])

    # Select only relevant columns for output
    table = table.select([
        "Protein",
        "contig",
        "circular_contig",
        "genome",
        "KEGG_hmm_id",
        "KEGG_Description",
        "Pfam_hmm_id",
        "Pfam_Description",
        "dbCAN_hmm_id",
        "dbCAN_Description",
        "METABOLIC_hmm_id",
        "METABOLIC_Description",
        "PHROG_hmm_id",
        "PHROG_Description",
        "top_hit_hmm_id",
        "top_hit_description",
        "top_hit_db",
        "window_avg_KEGG_L-score_viral",
        "window_avg_Pfam_L-score_viral",
        "Viral_Flanking_Genes"
    ])
    table = table.rename({
        "contig": "Contig",
        "genome": "Genome"
    })
    
    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])

def filter_annots(table, metabolic_kos):
    """
    Filter results to exclude false AMGs.
    """
    # Filter results to only include genes in virus-like windows
    table = table.filter(
        (pl.col("window_avg_KEGG_L-score_viral") == True) |
        (pl.col("window_avg_Pfam_L-score_viral") == True)
        )
    
    false_substrings_desc = [
        "EC:3.1.21", "EC:3.1.26", "EC:3.1.27", "nuclease",
        "EC 3.1.21", "EC 3.1.26", "EC 3.1.27",
        "EC:3.2", "EC 3.2", "glycosylhydrolase", 
        "EC:3.2.1.17", "EC 3.2.1.17", "lysozyme", "muramidase",
        "EC:2.4", "EC 2.4", "glycosyltransferase",
        "EC:2.1", "EC 2.1", "methyltransferase", "methylase", "methylation",
        "adenylyltransferase", "peptidase",
        "transcriptional", "tRNA synthetase", "tRNA",
        "phage tail", "phage head",
        "tail fiber", "tail peptide",
        "baseplate", "portal",
        "Ribonucleotide Reductase", "Ribonucleotide reductase","ribonucleotide reductase", "ribonucleoside-diphosphate reductase", "Ribonucleoside-diphosphate reductase",
        "Methyltransferase", "methyltransferase",
        "Methylase", "methylase",
        "Glycosyltransferase", "glycosyltransferase", "Glycosyl transferase", "glycosyl transferase",
        "Glycosyl hydrolase", "glycosyl hydrolase", "Glycoside hydrolase", "glycoside hydrolase",
        "UDP-glucose/GDP-mannose dehydrogenase",
        "dTDP", "dUDP", "dADP", "dCDP", "dGDP", "dCMP",
        "TDP", "UDP", "ADP", "CDP", "GDP", "CMP",
        "SPFH domain / Band 7 family",
        "PLD-like domain",
        "deoxyribose", "Deoxyribose",
        "Nucleotidyl transferase", "nucleotidyl transferase", "nucleotidyltransferase", "Nucleotidyltransferase",
        "Sigma-70", "sigma-70", "sigma70", "Sigma70",
        "Concanavalin A-like lectin/glucanases superfamily", "Concanavalin A-like lectin", # glycosyl hydrolase in disguise
        "Histidine phosphatase", "histidine phosphatase", "histidinol phosphatase", "Histidinol phosphatase", "histidinol-phosphatase", "Histidinol-phosphatase", # Nucleases in disguise
        "Cytidylate kinase",  "cytidylate kinase", "CMP kinase", "CMPK", # Nucleotide metabolism in disguise
        "4Fe-4S ferredoxin", # Nuclease/nucleotide metabolism in disguise
        "Radical SAM superfamily", "anaerobic magnesium-protoporphyrin IX monomethyl ester cyclase", # DNA modification in disguise (also contains 4Fe-4S clusters)
        # From Cody Martin et al. AMG perspective:
        "QueC", "Queuosine", "queuosine", # tRNA maturation in disguise
        "DsrC like protein", # tRNA modification in disguise
        "cysteine", "Cysteine", "methionine", "Methionine", # DNA modification in disguise
        "chitinase", "Chitinase" # lysozyme in disguise
    ]
    
    for substring in false_substrings_desc:
        table = table.filter(
            ~(
                pl.col("KEGG_Description").str.contains(substring).fill_null(False) |
                pl.col("Pfam_Description").str.contains(substring).fill_null(False) |
                pl.col("dbCAN_Description").str.contains(substring).fill_null(False) |
                pl.col("METABOLIC_Description").str.contains(substring).fill_null(False) |
                pl.col("PHROG_Description").str.contains(substring).fill_null(False)
            )
        )
    
    # dcCAN/CAZyme-spescific false AMGs 
    false_substrings_dbcan_id = [
        "CBM", # Carbohydrate-binding modules usually match to tail peptides used for receptor binding
        "GT", # Glycosyltransferases
        "GH", # Glycoside hydrolases
    ]   
    for substring in false_substrings_dbcan_id:
        table = table.filter(~(pl.col("dbCAN_hmm_id").str.contains(substring).fill_null(False)))

    
    # Remove .X or .XX suffixes from top_hit_hmm_id for proper matching of Pfam hits
    table = table.with_columns(
        pl.col("top_hit_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("top_hit_hmm_id_clean")
    )
    
    # Filter results to only include top hits that are metabolic
    # Genes are considered metabolic if:
    # 1. they are in the METABOLIC custom database or dbCAN, or
    # 2. they are present in the VIBRANT AMG list of KOs, or
    # 3. they are present in the DRAM-V list of Pfam IDs and KOs
    # Even if the DRAM-V list contains likely non-AMGs,
    # they should not be present in the dataframe because of fitlering steps above.
    # The VIBRANT AMG list should already have been pre-filtered upon creation
    # to exclude non-AMGs based on KEGG module, descriptions, and ECs.

    # Filter results as outlined above
    table = table.filter(
        (
            ((pl.col("top_hit_db") == "METABOLIC") | (pl.col("top_hit_db") == "dbCAN")) &
            (pl.col("top_hit_hmm_id") != "None")
        ) |
        (
            (pl.col("top_hit_db") == "KEGG") &
            (pl.col("top_hit_hmm_id").is_in(metabolic_kos["KO"]))
        ) |
        (
            (pl.col("top_hit_db") == "PHROG") &
            (
                (
                    pl.col("KEGG_hmm_id").is_in(metabolic_kos["KO"])
                ) |
                (
                    (pl.col("METABOLIC_hmm_id").is_not_null()) &
                    (pl.col("METABOLIC_hmm_id") != "None")
                ) |
                (
                    (pl.col("dbCAN_hmm_id").is_not_null()) &
                    (pl.col("dbCAN_hmm_id") != "None")
                )
            )
        )
    )
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")
    
    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])

def write_fasta_files(aux_table, prot_paths, fasta_outdir):
    """
    Writes FASTA files based on the auxiliary table and protein paths. Proteins are written to different
    files based on the presence of viral flanking genes and metabolic status.
    
    Args:
        aux_table (DataFrame): Polars DataFrame containing the summarized auxiliary table.
        prot_paths (list): List of strings containing paths to protein FASTA files.
        fasta_outdir (str): Output directory for the FASTA files.
    """
    
    aux_no_flank_path = os.path.join(fasta_outdir, "AMGs_no_viral_flank.faa")
    aux_with_flank_path = os.path.join(fasta_outdir, "AMGs_with_viral_flank.faa")
    aux_all_path = os.path.join(fasta_outdir, "AMGs_all.faa")
    
    no_flank_records = []
    with_flank_records = []
    all_records = []
    
    # Load all protein sequences
    prot_records = {}
    for prot_path in prot_paths:
        for record in SeqIO.parse(prot_path, "fasta"):
            prot_records[record.id] = record

    # Match proteins in aux_table with sequences and categorize them
    for protein in aux_table.iter_rows(named=True):
        record = prot_records.get(protein['Protein'])
        if record:
            # Modify the sequence header
            if protein['top_hit_hmm_id'] and protein['top_hit_description']:
                record.description = f'"{protein["top_hit_hmm_id"]} {protein["top_hit_description"]}"'
            record.id = record.id.split(" ")[0]  # Keep only the sequence name
            all_records.append(record)
            
            if protein['circular_contig']:
                with_flank_records.append(record)
            else:
                if protein['Viral_Flanking_Genes']:
                    with_flank_records.append(record)
                else:
                    no_flank_records.append(record)

    # Write the sequences to the respective FASTA files
    SeqIO.write(no_flank_records, aux_no_flank_path, "fasta")
    SeqIO.write(with_flank_records, aux_with_flank_path, "fasta")
    SeqIO.write(all_records, aux_all_path, "fasta")

def main():
    input_table  = snakemake.input.context_table
    hmm_ref = snakemake.params.hmm_ref
    amg_ref = snakemake.params.vibrant_amgs
    input_prots_dir = snakemake.params.protein_dir
    out_table = snakemake.params.aux_table
    all_annot_out_table = snakemake.params.all_annot_out_table
    fasta_outdir = snakemake.params.aux_fasta_dir

    logging.info("Starting the curation of AMG results...")
    
    hmm_descriptions = pl.read_csv(hmm_ref,
                                   schema={"id": pl.Utf8,
                                           "db": pl.Utf8,
                                           "name": pl.Utf8}).select(["id", "db", "name"]).filter((pl.col("db") == "KEGG") | (pl.col("db") == "Pfam") | (pl.col("db") == "dbCAN") | (pl.col("db") == "METABOLIC") | (pl.col("db") == "PHROG"))
    metabolic_kos = pl.read_csv(amg_ref, separator="\t", schema={"KO": pl.Utf8})
        
    annot_table = summarize_aux_table(input_table, hmm_descriptions)    
    aux_table = filter_annots(annot_table, metabolic_kos)
    
    annot_table.drop("window_avg_KEGG_L-score_viral", "window_avg_Pfam_L-score_viral").write_csv(all_annot_out_table, separator="\t")
    aux_table.write_csv(out_table, separator="\t")
    
    prot_paths = load_prot_paths.load_prots(input_prots_dir)    
    write_fasta_files(aux_table, prot_paths, fasta_outdir)
    
    logging.info("AMG results curation completed.")
    logging.info(f"Total number of genes analyzed: {annot_table.shape[0]}")
    logging.info(f"Number of identified AMGs: {aux_table.shape[0]}")
    logging.info(f"AMGs with viral flanking genes or encoded on circular contigs: {aux_table.filter((aux_table['Viral_Flanking_Genes'] == True) | (aux_table['circular_contig'] == True)).shape[0]}")
    logging.info(f"AMGs without viral flanking genes encoded on non-circular contigs: {aux_table.filter((aux_table['Viral_Flanking_Genes'] == False) & (aux_table['circular_contig'] == False)).shape[0]}")

if __name__ == "__main__":
    main()