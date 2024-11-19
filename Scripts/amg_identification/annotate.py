#!/usr/bin/env python3

import os
import resource
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import collections
import pyhmmer
from pyhmmer import easel, plan7
from Bio import SeqIO
import load_prot_paths

# Function to set memory limit in GB
def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024  # Convert GB to bytes
    resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_hmms(hmmdbs):
    """Load HMM profiles from given HMM database files."""
    hmms = {}
    for db_path in hmmdbs:
        with plan7.HMMFile(db_path) as hmm_file:
            hmms[db_path] = list(hmm_file)
    return hmms

def filter_results(results):
    """
    Filter the search results to keep only the best scoring hit for each sequence within each database.

    Args:
        results (list of namedtuples): The list of search results, each entry is a namedtuple
                                        with fields including 'sequence', 'score', and 'db_path'.

    Returns:
        list: A list of filtered results with only the best scoring hit per sequence per database.
    """
    best_results = {}
    for result in results:
        # Use a tuple of sequence and db_path as the key to ensure uniqueness per sequence per database
        key = (result.sequence, result.db_path)
        if key not in best_results or result.score > best_results[key].score:
            best_results[key] = result
    return list(best_results.values())

def get_hmm_coverage(domain):
    """
    Calculate the alignment coverage for a given domain.
    """
    n_aligned_positions = domain.alignment.hmm_to - domain.alignment.hmm_from + 1
    return n_aligned_positions / domain.alignment.hmm_length

def run_hmmsearch_on_sequences(sequence_files, hmms, e_value_threshold=1e-5, num_cpus=1, cov_fraction=0.5):
    """Perform an HMM search against a list of sequences using HMM profiles from each database."""
    Result = collections.namedtuple("Result", ["hmm_id", "sequence", "score", "db_path", "hmm_coverage"])

    all_results = []
    for sequence_file in sequence_files:
        # Validate the sequence file format using biopython
        try:
            with open(sequence_file, 'r') as f:
                sequences = list(SeqIO.parse(f, "fasta"))
                if not sequences:
                    raise ValueError(f"No sequences found in file {sequence_file}.")
        except Exception as e:
            raise ValueError(f"Error reading sequence file {sequence_file}: {e}")
        
        # Create the amino acid alphabet
        aa_alphabet = easel.Alphabet.amino()
        
        # Read the sequences using pyhmmer.easel after validation
        with easel.SequenceFile(sequence_file, format="fasta", digital=True, alphabet=aa_alphabet) as seqs_file:
            proteins = list(seqs_file)

        for db_path, hmm_list in hmms.items():
            results = []
            
            for hits in pyhmmer.hmmsearch(queries=hmm_list, sequences=proteins, cpus=num_cpus, E=e_value_threshold):
                for hit in hits:
                    for domain in hit.domains.included:
                        hmm_coverage = get_hmm_coverage(domain)
                        if hmm_coverage >= cov_fraction:  # Enforce alignment coverage threshold
                            if "Pfam" in db_path:
                                results.append(Result(hits.query_accession.decode(), hit.name.decode(), hit.score, db_path, hmm_coverage))
                            elif "eggNOG" in db_path:
                                results.append(Result(hits.query_name.decode().split(".")[0], hit.name.decode(), hit.score, db_path, hmm_coverage))
                            else:
                                if ".wlink.txt.mafft" in Result(hits.query_name.decode(), hit.name.decode(), hit.score, db_path, hmm_coverage)[0]:
                                    results.append(Result(hits.query_name.decode().split(".")[1], hit.name.decode(), hit.score, db_path, hmm_coverage))
                                else:
                                    results.append(Result(hits.query_name.decode().replace("_alignment", "").replace(".mafft", "").replace(".txt", "").replace(".hmm", "").replace("_protein.alignment", ""), hit.name.decode(), hit.score, db_path, hmm_coverage))
            # Filter results for this database and add to all results
            all_results.extend(filter_results(results))
    
    return all_results

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
    genome_or_build = snakemake.params.genome_or_build
    hmm_vscores = snakemake.params.hmm_vscores
    cov_fraction = snakemake.params.cov_fraction
    db_dir = snakemake.params.db_dir
    output = snakemake.output.vscores
    all_hmm_results = snakemake.output.all_hmm_results
    num_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    
    logging.info("Protein HMM alignments starting...")
    logging.info(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    set_memory_limit(mem_limit)
    
    prots = load_prot_paths.load_prots(protein_dir)
    
    hmmdbs = [os.path.join(db_dir, db) for db in os.listdir(db_dir) if db.endswith('.H3M') or db.endswith('.h3m')]
    logging.info(f"Loading HMM profiles at {hmmdbs}.")
    hmms = load_hmms(hmmdbs)
    logging.info("Running HMM search on sequences.")
    all_results = run_hmmsearch_on_sequences(prots, hmms, 1E-5, num_cpus, cov_fraction)
    
    # Prepare data for DataFrame creation
    logging.info("Formatting and writing results.")
    data = [{"hmm_id": r.hmm_id, "sequence": r.sequence, "score": r.score, "db_path": r.db_path} for r in all_results]
    df = pl.DataFrame(data)
    df = df.rename({"hmm_id": "id"})
    df = df.with_columns(pl.col("db_path").map_elements(assign_db, return_dtype=pl.Utf8).alias("db"))
    df = df.drop('db_path')
    df = df.rename({"id": "hmm_id"})
    df = df.drop("hmm_id").insert_column(1, df.get_column("hmm_id"))
    df = df.sort(['sequence', 'score', 'db', 'hmm_id'])

    # Write all hmm results to file
    df.write_csv(all_hmm_results, separator="\t")

    # Load V-scores CSV
    vscores_df = pl.read_csv(hmm_vscores, schema_overrides={"id": pl.Utf8, "V-score": pl.Float64, "L-score": pl.Float64, "db": pl.Categorical, "name": pl.Utf8})

    merged_df = df.join(vscores_df, left_on='hmm_id', right_on='id', how='left').filter(pl.col("V-score").is_not_null())
    df = df.drop("hmm_id").insert_column(1, df.get_column("hmm_id"))
    merged_df = merged_df.drop('name', 'db_right')
    merged_df = merged_df.sort(['sequence', 'score', 'V-score', 'db'])
    merged_df.write_csv(output, separator="\t")
    logging.info("Protein HMM alignments completed.")

if __name__ == "__main__":
    main()