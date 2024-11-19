#!/usr/bin/env python3

import os
import logging
import resource
import subprocess
import re
from Bio import SeqIO

# Function to set memory limit in GB
def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024  # Convert GB to bytes
    resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_pyrodigal(input_fasta, output_faa, threads, log_file):
    cmd = [
        "pyrodigal-gv",
        "-i", input_fasta,
        "-a", output_faa,
        "-j", str(threads)
    ]
    with open(log_file, "a") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.PIPE, text=True)
        
    # Check if the output faa file is empty
    if os.path.getsize(output_faa) == 0:
        logging.error(f"The output .faa file {output_faa} is empty. Did Pyrodigal-GV terminate because of memory usage?")
        raise RuntimeError(f"Error: The output .faa file {output_faa} is empty.")

def run_pyrodigal_on_vmag(vmag_fasta_files, vmag_proteins_subdir, threads, log_file):
    # Run pyrodigal-gv on each vMAG file separately to produce their own output files
    for vmag_fasta in vmag_fasta_files:
        if vmag_fasta.endswith(('.fna', '.fasta')):
            vmag_name = os.path.splitext(os.path.basename(vmag_fasta))[0]
            output_faa = os.path.join(vmag_proteins_subdir, f"{vmag_name}.faa")
            logging.debug(f"Running Pyrodigal-GV on vMAG {vmag_name}...")
            run_pyrodigal(vmag_fasta, output_faa, threads, log_file)
        else:
            raise ValueError(f"File {vmag_fasta} is not a valid FASTA file.")

def write_gene_to_genome(all_results, gene_to_genome_output_path):
    with open(gene_to_genome_output_path, 'w') as gene_to_genome_output:
        gene_to_genome_output.write("gene\tcontig\tgenome\tgene_number\n")  # Header
        for seq_id, predicted_genes, is_vmag, vmag_name in all_results:
            for gene_i, gene in enumerate(predicted_genes, 1):
                gene_name = f"{seq_id}_{gene_i}"
                genome_name = vmag_name if is_vmag else seq_id
                gene_to_genome_output.write(f"{gene_name}\t{seq_id}\t{genome_name}\t{gene_i}\n")

def parse_faa_file(faa_file, is_vmag, vmag_name=None):
    parsed_results = []
    for record in SeqIO.parse(faa_file, "fasta"):
        # Parse gene information from the header
        header = record.description
        gene_info = re.search(
            r">(.+) # (\d+) # (\d+) # ([+-]) # ID=(.+?);partial=(\d)(\d);start_type=(.+?);",
            header
        )

        if gene_info:
            seq_id = gene_info.group(1)
            begin = int(gene_info.group(2))
            end = int(gene_info.group(3))
            strand = gene_info.group(4)
            gene_number = gene_info.group(5)
            partial_begin = bool(int(gene_info.group(6)))
            partial_end = bool(int(gene_info.group(7)))
            start_type = gene_info.group(8)

            # Create a mock `Gene` object to hold this information, as pyrodigal-gv might have its own `Gene` class
            gene = {
                "seq_id": seq_id,
                "begin": begin,
                "end": end,
                "strand": strand,
                "gene_number": gene_number,
                "partial_begin": partial_begin,
                "partial_end": partial_end,
                "start_type": start_type,
                "translation": str(record.seq)
            }

            parsed_results.append(gene)

    return parsed_results

def parse_all_faa_files(output_single_contig_prots, vmag_proteins_subdir):
    all_results = []

    # Parse single contig results if present
    if os.path.exists(output_single_contig_prots):
        single_contig_results = parse_faa_file(output_single_contig_prots, is_vmag=False)
        for gene in single_contig_results:
            all_results.append((gene["seq_id"], [gene], False, None))
    
    # Parse vMAG results
    if os.path.exists(vmag_proteins_subdir):
        for faa_file in os.listdir(vmag_proteins_subdir):
            if faa_file.endswith(".faa"):
                vmag_name = os.path.splitext(faa_file)[0]
                vmag_file_path = os.path.join(vmag_proteins_subdir, faa_file)
                vmag_results = parse_faa_file(vmag_file_path, is_vmag=True, vmag_name=vmag_name)
                for gene in vmag_results:
                    all_results.append((gene["seq_id"], [gene], True, vmag_name))

    return all_results

def main():
    input_single_contig_genomes = snakemake.params.input_single_contig_genomes
    input_vmag_fastas = snakemake.params.input_vmag_fastas
    wdir = snakemake.params.wdir
    output_dir = snakemake.params.output_dir
    output_single_contig_prots = snakemake.params.single_contig_prots
    vmag_proteins_subdir = snakemake.params.vmag_proteins_subdir
    gene_to_genome = snakemake.params.gene_to_genome
    n_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    log_file = os.path.join(wdir, "pyrodigal-gv.log")

    logging.info("Pyrodigal-GV run starting...")
    logging.info(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    if not os.path.exists(wdir):
        os.makedirs(wdir, exist_ok=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    single_contig_fasta = input_single_contig_genomes
    vmag_fasta_files = input_vmag_fastas.split()

    # Run pyrodigal-gv on single contig genomes
    if single_contig_fasta:
        logging.info(f"Running Pyrodigal-GV on single contig genomes with {n_cpus} CPUs")
        run_pyrodigal(single_contig_fasta, output_single_contig_prots, n_cpus, log_file)
    
    # Run pyrodigal-gv on vMAGs if provided
    if vmag_fasta_files:
        os.makedirs(vmag_proteins_subdir, exist_ok=True)
        logging.info(f"Predicting and translating genes from vMAGs in {vmag_proteins_subdir}")
        run_pyrodigal_on_vmag(vmag_fasta_files, vmag_proteins_subdir, n_cpus, log_file)
    
    # Collect gene to genome mapping
    logging.info("Writing gene to genome mappings.")
    all_results = parse_all_faa_files(output_single_contig_prots, vmag_proteins_subdir)
    write_gene_to_genome(all_results, gene_to_genome)

    logging.info("Pyrodigal-GV run completed.")

if __name__ == "__main__":
    main()