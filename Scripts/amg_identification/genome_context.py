import os
import resource
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import numpy as np
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor
import logging

def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024  # Convert GB to bytes
    resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_gene_lengths(data):
    """
    Calculate gene lengths and protein lengths in amino acids.

    Parameters:
    data (pl.DataFrame): DataFrame containing gene data with start and end positions.

    Returns:
    pl.DataFrame: DataFrame with added gene length columns.
    """
    data = data.with_columns([
        (pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1).alias('gene_length_bases'),
        ((pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1) / 3).cast(pl.Int32, strict=False).alias('prot_length_AAs')
    ])
    
    # Cast gene_length_bases to Int64, handling any non-integer results
    data = data.with_columns([
        pl.col('gene_length_bases').cast(pl.Int64, strict=False).alias('gene_length_bases')
    ])
    
    # Identify rows where casting failed
    invalid_gene_lengths = data.filter(pl.col('gene_length_bases').is_null())
    if invalid_gene_lengths.shape[0] > 0:
        logging.warning(f"Found {invalid_gene_lengths.shape[0]} rows with invalid gene_length_bases. Skipping these rows.")
        data = data.filter(pl.col('gene_length_bases').is_not_null())
    
    return data

def calculate_contig_statistics(data, circular_contigs):
    """
    Calculate the contig statistics including mean for V-scores and L-scores, 
    and check if contig is in the complete contigs set.

    Parameters:
    data (pl.DataFrame): DataFrame containing contig data.
    circular_contigs (set): Set of contig names previously checked to be circular.

    Returns:
    pl.DataFrame: DataFrame with added contig mean columns and circular_contig status.
    """
    # Calculate statistics for V-scores and L-scores
    stats = data.group_by("contig").agg([
        pl.col("KEGG_V-score").mean().alias("contig_avg_KEGG_V-score"),
        pl.col("Pfam_V-score").mean().alias("contig_avg_Pfam_V-score"),
        pl.col("KEGG_L-score").mean().alias("contig_avg_KEGG_L-score"),
        pl.col("Pfam_L-score").mean().alias("contig_avg_Pfam_L-score")
    ])
    
    # Join statistics back to original data
    result = data.join(stats, on="contig")
    
    # Add a column 'circular_contig' that checks if contig is in circular_contigs
    result = result.with_columns(pl.col("contig").is_in(circular_contigs).alias("circular_contig"))
    
    return result

def window_avg(scores, lengths, index, window_size, minimum_percentage):
    """
    Calculate the average score within a variable-length window surrounding each gene.

    Parameters:
    scores (list): List of scores (either L-scores or V-scores) for each gene.
    lengths (list): List of gene lengths corresponding to each score.
    index (int): Index of the current gene for which the window average is being calculated.
    window_size (int): Length of the window in bases.
    minimum_percentage (float): Minimum percentage of valid gene length required for the window.

    Returns:
    float or NaN: The average score within the window if the valid length percentage criteria is met; otherwise, NaN.
    """   
    left_length, right_length = 0, 0
    left_scores, right_scores = [], []

    # Check left side of the window
    for i in range(index - 1, -1, -1):
        try:
            length = int(lengths[i])
        except (ValueError, TypeError):
            logging.debug(f"Invalid length at index {i}: {lengths[i]}")
            length = np.nan  # Skip if length is not valid
        
        left_length += length
        if left_length >= window_size:
            break
        left_scores.append(scores[i])

    # Check right side of the window
    for i in range(index + 1, len(scores)):
        try:
            length = int(lengths[i])
        except (ValueError, TypeError):
            logging.debug(f"Invalid length at index {i}: {lengths[i]}")
            length = np.nan  # Skip if length is not valid

        right_length += length
        if right_length >= window_size:
            break
        right_scores.append(scores[i])

    total_length = left_length + right_length + lengths[index]
    valid_scores = [score for score in left_scores + [scores[index]] + right_scores if score is not None]
    
    try:
        valid_length = sum([
            length if score is not None else 0
            for score, length in zip(
                left_scores + [scores[index]] + right_scores,
                [lengths[i] for i in range(index - 1, -1, -1)] + [lengths[index]] + [lengths[i] for i in range(index + 1, len(scores))]
            )
        ])
    except (ValueError, TypeError):
        logging.debug(f"Invalid length at index {index}: {lengths[index]}")
        valid_length = 0

    if total_length > 0 and (valid_length / total_length) * 100 >= minimum_percentage:
        return sum(valid_scores) / len(valid_scores) if valid_scores else float('nan')
    else:
        return float('nan')

def calculate_window_statistics(data, window_size, minimum_percentage):
    """
    Calculate the average L-scores and V-scores within a variable-length window surrounding each gene.

    Parameters:
    data (pl.DataFrame): DataFrame containing contig data.
    window_size (int): Length of the window in bases.
    minimum_percentage (float): Minimum percentage of valid gene length required for the window.

    Returns:
    pl.DataFrame: DataFrame with added columns for average L-scores and V-scores within the window.
    """
    # Log the initial structure of the dataframe
    logging.debug(f"Initial data structure: {data.dtypes}")
    
    # Ensure that the relevant columns are cast to the correct types
    data = data.with_columns([
        pl.col("KEGG_L-score").cast(pl.Float64, strict=False).alias("KEGG_L-score"),
        pl.col("Pfam_L-score").cast(pl.Float64, strict=False).alias("Pfam_L-score"),
        pl.col("KEGG_V-score").cast(pl.Float64, strict=False).alias("KEGG_V-score"),
        pl.col("Pfam_V-score").cast(pl.Float64, strict=False).alias("Pfam_V-score"),
        pl.col("gene_length_bases").cast(pl.Int64, strict=False).alias("gene_length_bases")
    ])

    # Log the casting result to check if columns were successfully cast
    logging.debug(f"Data after casting to numeric types: {data.dtypes}")

    # Sort the data
    data = data.sort(["contig", "contig_pos_start"])
    
    def calculate_window_avg(df, score_col, length_col):
        scores = df[score_col].to_list()
        lengths = df[length_col].to_list()

        # Log the scores and lengths for debugging
        logging.debug(f"Scores ({score_col}): {scores[:10]}")  # Show only the first 10 for brevity
        logging.debug(f"Lengths ({length_col}): {lengths[:10]}")

        indices = list(range(len(scores)))
        avg_scores = [window_avg(scores, lengths, i, window_size, minimum_percentage) for i in indices]
        return pl.DataFrame({f"window_avg_{score_col}": avg_scores})

    # Calculate window averages for all score columns grouped by 'contig'
    data = data.group_by("contig").map_groups(lambda df: df.with_columns(
        pl.Series(calculate_window_avg(df, "KEGG_L-score", "gene_length_bases")),
        pl.Series(calculate_window_avg(df, "Pfam_L-score", "gene_length_bases")),
        pl.Series(calculate_window_avg(df, "KEGG_V-score", "gene_length_bases")),
        pl.Series(calculate_window_avg(df, "Pfam_V-score", "gene_length_bases"))
    ))

    logging.debug(f"Final data structure after window statistics: {data.dtypes}")
    return data

def check_window_avg_lscore(data, min_window_avg_lscore):
    """
    Check whether the window average L-scores for KEGG and Pfam exceed a given threshold.

    Parameters:
    data (pl.DataFrame): DataFrame containing contig data with window averages.
    min_window_avg_lscore (float): Minimum threshold for window average L-scores.

    Returns:
    pl.DataFrame: DataFrame with flags indicating whether the window averages exceed the threshold.
    """
    data = data.with_columns([
        pl.when(((pl.col('window_avg_KEGG_L-score') > min_window_avg_lscore) | (pl.col('window_avg_Pfam_L-score') > min_window_avg_lscore))).then(True).otherwise(False).alias('window_avg_KEGG_L-score_viral'),
        pl.when(((pl.col('window_avg_KEGG_L-score') > min_window_avg_lscore) | (pl.col('window_avg_Pfam_L-score') > min_window_avg_lscore))).then(True).otherwise(False).alias('window_avg_Pfam_L-score_viral')
    ])
    
    return data

def verify_flanking_vscores(contig_data, minimum_vscore, max_flank_length):
    """
    Verify whether flanking genes have V-scores above a given threshold within a max flank length.

    Parameters:
    contig_data (pl.DataFrame): DataFrame containing contig data.
    minimum_vscore (float): Minimum V-score threshold for flanking genes.
    max_flank_length (int): Maximum allowable distance (in bases) for flanking genes.

    Returns:
    pl.DataFrame: DataFrame with flags for verified flanking genes based on V-scores.
    """
    def flank_verification(vscores, lengths, minimum_vscore, max_flank_length):
        flanks = []
        
        for i in range(len(vscores)):
            left_length = right_length = 0
            left_flank = right_flank = False
            
            # Check left flank
            for j in range(i - 1, -1, -1):
                try:
                    length = int(lengths[j])
                except (ValueError, TypeError):
                    logging.debug(f"Invalid length at index {j}: {lengths[j]}")
                    length = np.nan  # Skip if length is not valid
                
                left_length += length
                if vscores[j] is not None and vscores[j] >= minimum_vscore:
                    left_flank = True
                if left_length >= max_flank_length:
                    break
            
            # Check right flank
            for j in range(i + 1, len(vscores)):
                try:
                    length = int(lengths[j])
                except (ValueError, TypeError):
                    logging.debug(f"Invalid length at index {j}: {lengths[j]}")
                    length = np.nan  # Skip if length is not valid

                right_length += length
                if vscores[j] is not None and vscores[j] >= minimum_vscore:
                    right_flank = True
                if right_length >= max_flank_length:
                    break
            
            flanks.append(left_flank and right_flank)
        
        return flanks


    gene_lengths = contig_data['gene_length_bases']
    kegg_vscore_flanks = flank_verification(contig_data['KEGG_V-score'], gene_lengths, minimum_vscore, max_flank_length)
    pfam_vscore_flanks = flank_verification(contig_data['Pfam_V-score'], gene_lengths, minimum_vscore, max_flank_length)

    contig_data = contig_data.with_columns([
        pl.Series('KEGG_verified_flank', kegg_vscore_flanks).alias('KEGG_verified_flank'),
        pl.Series('Pfam_verified_flank', pfam_vscore_flanks).alias('Pfam_verified_flank')
    ])

    return contig_data

def process_genomes(data, circular_contigs, chunk_no, minimum_percentage, min_window_avg_lscore, window_size, max_flank_length, minimum_vscore):
    """
    Process the genome data to calculate statistics and verify gene contexts using V-scores.

    Parameters:
    data (pl.DataFrame): DataFrame containing the genome data.
    circular_contigs: set of circular contigs.
    minimum_percentage (float): Minimum percentage of valid gene length for window calculations.
    min_window_avg_lscore (float): Minimum L-score threshold for window averages.
    window_size (int): Size of the window in bases for averaging.
    max_flank_length (int): Maximum flank length for verification.
    minimum_vscore (float): Minimum V-score for flanking gene verification.

    Returns:
    pl.DataFrame: Processed DataFrame (eager) with gene context verification results.
    """
    logging.debug(f"Calculating gene lengths for chunk {chunk_no}.")
    data = calculate_gene_lengths(data)
    
    logging.debug(f"Calculating contig statistics for chunk {chunk_no}.")
    data = calculate_contig_statistics(data, circular_contigs)

    logging.debug(f"Calculating window statistics for chunk {chunk_no}.")
    data = calculate_window_statistics(data, window_size, minimum_percentage)
    data = check_window_avg_lscore(data, min_window_avg_lscore)
    data = data.unique()

    data = verify_flanking_vscores(data, minimum_vscore, max_flank_length)

    data = data.unique()
    
    # If 'circular_contig' exists, rearrange columns to place it after 'member' and 'contig'
    if 'circular_contig' in data.columns:
        columns = data.columns
        # Place 'member', 'contig', and 'circular_contig' in order, then all other columns
        reordered_columns = (
            [col for col in columns if col == 'member'] +
            [col for col in columns if col == 'contig'] +
            ['circular_contig'] +
            [col for col in columns if col not in {'member', 'contig', 'circular_contig'}]
        )
        data = data.select(reordered_columns)
    
    return data

def split_dataframe_by_num_chunks(data, num_chunks):
    """
    Efficiently split the dataframe into a specified number of chunks by genome, ensuring no genome is split across chunks.

    Parameters:
    data (pl.DataFrame): The input dataframe to split.
    num_chunks (int): The number of chunks to split the data into.

    Returns:
    list: A list of smaller DataFrames split by genome.
    """
    # Ensure the data is sorted by the genome column
    data = data.sort("genome")

    # Get unique genome IDs and their corresponding row counts
    genome_boundaries = data.group_by("genome").agg([pl.len().alias("genome_count")])

    # Calculate number of genomes per chunk
    total_genomes = genome_boundaries.shape[0]
    genomes_per_chunk = total_genomes // num_chunks

    chunks = []
    current_chunk = []
    current_genome_count = 0

    # Iterate through each genome and assign them to chunks
    for genome, genome_count in zip(genome_boundaries["genome"].to_list(), genome_boundaries["genome_count"].to_list()):
        current_chunk.append(genome)
        current_genome_count += 1

        if current_genome_count >= genomes_per_chunk and len(chunks) < num_chunks - 1:
            chunk_data = data.filter(pl.col("genome").is_in(current_chunk))
            chunks.append(chunk_data)
            current_chunk = []
            current_genome_count = 0

    # Add remaining genomes to the last chunk
    if current_chunk:
        chunk_data = data.filter(pl.col("genome").is_in(current_chunk))
        chunks.append(chunk_data)

    return chunks

def process_genomes_from_temp(temp_filepath, circular_contigs, chunk_no, minimum_percentage, min_window_avg_lscore, window_size, max_flank_length, minimum_vscore):
    """
    Load a temporary file, process the genome data, and save the results to a temporary file.

    Parameters:
    temp_filepath (str): Path to the temporary file to process.
    """
    logging.debug(f"Processing temporary file: {temp_filepath}")
    data = pl.read_csv(temp_filepath, separator='\t')

    # Clean and process the data
    data = process_genomes(data, circular_contigs, chunk_no, minimum_percentage, min_window_avg_lscore, window_size, max_flank_length, minimum_vscore)
    temp_output_filepath = temp_filepath.replace('.tsv', '_processed.tsv')
    
    logging.debug(f"Writing processed data to: {temp_output_filepath}")
    data.write_csv(temp_output_filepath, separator='\t')

def combine_processed_files(temp_folder, output_file):
    """
    Combine all processed files in the temporary folder into a single output file.

    Parameters:
    temp_folder (str): Path to the temporary folder containing processed files.
    output_file (str): Path to the final output file.
    """
    temp_files = [os.path.join(temp_folder, f) for f in os.listdir(temp_folder) if f.endswith('_processed.tsv')]
    
    combined_data = None

    for temp_file in temp_files:
        logging.debug(f"Combining file: {temp_file}")
        temp_data = pl.read_csv(temp_file, separator='\t')

        # Ensure consistent schema between files
        if combined_data is None:
            combined_data = temp_data
        else:
            # Align schema: ensure columns that should be float are cast to float
            for col in combined_data.columns:
                if combined_data[col].dtype == pl.Float64 and temp_data[col].dtype != pl.Float64:
                    temp_data = temp_data.with_columns(pl.col(col).cast(pl.Float64, strict=False))

            combined_data = pl.concat([combined_data, temp_data], how="vertical").unique()

    logging.debug(f"Writing final combined data to: {output_file}")
    combined_data.write_csv(output_file, separator='\t')

def clean_data(data):
    """
    Clean the data by removing rows with null contig positions and logging a warning.
    
    Parameters:
    data (pl.DataFrame): The input data containing gene annotations.
    
    Returns:
    pl.DataFrame: Cleaned DataFrame with valid contig positions.
    """
    logging.info(f"Loaded data with {data.shape[0]} rows and {data.shape[1]} columns.")
    
    # Remove rows with null contig positions
    missing_positions = data.filter(pl.col('contig_pos_start').is_null() | pl.col('contig_pos_end').is_null())
    if missing_positions.shape[0] > 0:
        logging.warning(f"Found {missing_positions.shape[0]} rows with missing contig positions. Skipping these rows.")
    cleaned_data = data.filter(pl.col('contig_pos_start').is_not_null() & pl.col('contig_pos_end').is_not_null())
    
    # Attempt to cast contig_pos_start and contig_pos_end to integers
    cleaned_data = cleaned_data.with_columns([
        pl.col('contig_pos_start').cast(pl.Int64, strict=False).alias('contig_pos_start'),
        pl.col('contig_pos_end').cast(pl.Int64, strict=False).alias('contig_pos_end')
    ])
    
    # Identify rows where casting failed (resulted in nulls)
    casting_issues = cleaned_data.filter(pl.col('contig_pos_start').is_null() | pl.col('contig_pos_end').is_null())
    if casting_issues.shape[0] > 0:
        logging.warning(f"Found {casting_issues.shape[0]} rows with non-integer contig positions. Skipping these rows.")
    
    # Remove rows where casting failed
    cleaned_data = cleaned_data.filter(pl.col('contig_pos_start').is_not_null() & pl.col('contig_pos_end').is_not_null())
    
    logging.debug(f"Data after cleaning contains {cleaned_data.shape[0]} rows and {cleaned_data.shape[1]} columns.")
    
    return cleaned_data

def main():
    input_file = snakemake.params.gene_index_annotated
    output_file = snakemake.output.context_table
    circular_contigs_file = snakemake.params.circular_contigs
    minimum_percentage = snakemake.params.annotation_percent_threshold
    min_window_avg_lscore = snakemake.params.min_window_avg_lscore
    window_size = snakemake.params.window_size
    minimum_vscore = snakemake.params.minimum_flank_vscore
    max_flank_length = snakemake.params.max_flank_length
    outparent = snakemake.params.outparent
    n_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    tmp_dir = snakemake.params.tmp_dir  # Variable for the temp folder location

    logging.info("Starting genome context analysis...")
    logging.info(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    set_memory_limit(mem_limit)
    
    os.makedirs(outparent, exist_ok=True)

    logging.info(f"Reading input files: {input_file} and {circular_contigs_file}")
    data = pl.read_csv(input_file, separator='\t')  # Eager loading
    data = clean_data(data)
    circular_contigs = set(pl.read_csv(circular_contigs_file, separator='\t')['contig'].to_list())

    total_rows = data.shape[0]
    
    baseline_multiplier = 1e-5
    if mem_limit > 1000:
        multiplier = baseline_multiplier
    elif mem_limit > 100:
        multiplier = baseline_multiplier * 10
    elif mem_limit > 10:
        multiplier = baseline_multiplier * 100
    else:
        multiplier = baseline_multiplier * 1000
        
    n_chunks = max(1, total_rows * multiplier)
    logging.debug(f"Total rows: {total_rows}, Number of chunks: {n_chunks}")
    if n_chunks > 1e4:
        logging.warning(f"Number of chunks needed to split the input table {input_file} is very high: {n_chunks}!")
        logging.warning(f"Reducing to 10000, but the amount of memory required to process each chunk may exceed the maximum you specified ({mem_limit} GB) and cause this script to terminate with an error!")
        n_chunks = 1e4

    # Split data into smaller chunks by genome
    chunks = split_dataframe_by_num_chunks(data, n_chunks)

    logging.debug(f"The input table {input_file} has been split into {len(chunks)} temporary chunks.")

    # Create a temporary folder beneath the directory specified by tmp_dir
    temp_folder = tempfile.mkdtemp(dir=tmp_dir)
    logging.debug(f"Using temporary folder: {temp_folder}")
    
    try:
        # If the dataframe is smaller than n_chunks, process it directly
        if n_chunks == 1:
            temp_filepath = os.path.join(temp_folder, "temp_chunk_0.tsv")
            data.write_csv(temp_filepath, separator='\t')
            logging.debug(f"Data written to temporary file: {temp_filepath}")
            process_genomes_from_temp(temp_filepath, circular_contigs, n_chunks, minimum_percentage, min_window_avg_lscore, window_size, max_flank_length, minimum_vscore)
        else:
            del data # Remove original data from memory
            # Process each chunk one at a time
            for idx, chunk in enumerate(chunks):
                temp_filepath = os.path.join(temp_folder, f"temp_chunk_{idx}.tsv")
                chunk.write_csv(temp_filepath, separator='\t')
                logging.debug(f"Processing chunk {idx} with {chunk.shape[0]} rows. Temporary file: {temp_filepath}")
                process_genomes_from_temp(temp_filepath, circular_contigs, idx, minimum_percentage, min_window_avg_lscore, window_size, max_flank_length, minimum_vscore)
        
        # Combine all processed chunks into the final output file
        logging.info("Combining processed files into final output.")
        combine_processed_files(temp_folder, output_file)
        logging.info(f"Final output written to: {output_file}")

    finally:
        # Clean up the temporary folder after processing
        logging.debug(f"Removing temporary folder: {temp_folder}")
        shutil.rmtree(temp_folder)

    logging.info("Genome context analysis completed.")

if __name__ == "__main__":
    main()