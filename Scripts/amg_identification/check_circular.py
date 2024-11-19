#!/usr/bin/env python3

# This script is modified from portions of the tool CheckV, originally developed by
# S. Nayfach et al. (2021) in their publication "CheckV: assessing the quality of metagenome-assembled viral genomes" (Nature Biotechnology, https://doi.org/10.1038/s41587-020-00774-7). 
# The original source code can be found at: 
# https://bitbucket.org/berkeleylab/checkv/src/master
# 
# Specifically, code modifications were made based on:
# 1. utility.py - https://bitbucket.org/berkeleylab/checkv/src/master/checkv/utility.py
# 2. complete_genomes.py - https://bitbucket.org/berkeleylab/checkv/src/master/checkv/modules/complete_genomes.py
#
# These modifications allow for assessing sequence completeness based on direct terminal repeats (DTRs)
# without having to run CheckV in its entirety. Thus, this script just checks for circularity
# and inverted terminal repeats (ITRs) without reliance on CheckV-generated intermediate files.

import os
import sys
import csv
import time
import logging
import re
import Bio.SeqIO
import collections
import multiprocessing as mp
import resource

# Function to set memory limit in GB
def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024  # Convert GB to bytes
    resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Genome class
class Genome:
    def __init__(self):
        self.seq = None
        self.id = None
        self.length = None
        self.dtr = None

# TR class
class TR:
    def __init__(self):
        pass

def fetch_dtr(fullseq, min_length=20):
    startseq = fullseq[0:min_length]
    matches = [
        m.start() for m in re.finditer("(?={0})".format(re.escape(startseq)), fullseq)
    ]
    matches = [_ for _ in matches if _ >= len(fullseq) / 2]
    for matchpos in matches:
        endseq = fullseq[matchpos:]
        if fullseq[0 : len(endseq)] == endseq:
            return endseq
    return ""

def reverse_complement(seq):
    trans = str.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)

def fetch_itr(seq, min_len, max_len=1000):
    rev = reverse_complement(seq)
    if seq[:min_len] == rev[:min_len]:
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    return ""

def calculate_kmer_frequency(genome_seq, k = 15):
    """Calculate kmer frequency for a genome sequence."""
    kmer_counts = collections.Counter()
    for i in range(len(genome_seq) - k + 1):
        kmer = genome_seq[i:i+k]
        kmer_counts[kmer] += 1
    avg_kmer_freq = sum(kmer_counts.values()) / len(kmer_counts)
    return avg_kmer_freq

def main():
    input_fasta = snakemake.params.input_single_contig_genomes
    input_vmag_fastas = snakemake.params.input_vmag_fastas
    output_path = snakemake.params.circularity_tbl
    tr_min_len = snakemake.params.tr_min_len
    tr_max_len = snakemake.params.tr_max_len
    tr_max_count = snakemake.params.tr_max_count
    tr_max_ambig = snakemake.params.tr_max_ambig
    tr_max_basefreq = snakemake.params.tr_max_basefreq
    kmer_max_freq = snakemake.params.kmer_max_freq
    k = snakemake.params.k
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logging.info("Sequence circularity check starting...")
    single_contig_fasta = input_fasta
    vmag_fasta_files = input_vmag_fastas.split()
    
    genomes = {}
    # Parse input single-contig sequences
    if single_contig_fasta:
        logging.info(f"Reading single-contig sequences from {single_contig_fasta}")
        for record in Bio.SeqIO.parse(input_fasta, "fasta"):
            if len(record.seq) == 0:
                continue
            genome = Genome()
            genome.id = record.id
            genome.seq = str(record.seq).upper()
            genome.length = len(genome.seq)
            genome.kmer_freq = calculate_kmer_frequency(genome.seq, k)
            genomes[genome.id] = genome
    
    # Parse input vMAG sequences
    if vmag_fasta_files:
        logging.info(f"Reading vMAG sequences from {os.path.dirname(vmag_fasta_files[0])}")
        for fasta in vmag_fasta_files:
            # Below, each contig/scaffold in a vMAG is named as a 'genome' for consistency
            for record in Bio.SeqIO.parse(fasta, "fasta"):
                if len(record.seq) == 0:
                    continue
                genome = Genome()
                genome.id = record.id
                genome.seq = str(record.seq).upper()
                genome.length = len(genome.seq)
                genome.kmer_freq = calculate_kmer_frequency(genome.seq, k)
                genomes[genome.id] = genome

    # Find DTRs and ITRs
    for genome in genomes.values():
        genome.tr = TR()
        dtr = fetch_dtr(genome.seq, tr_min_len)
        itr = fetch_itr(genome.seq, tr_min_len, tr_max_len)
        if len(dtr) < tr_min_len and len(itr) < tr_min_len:
            genome.tr.type = None
        elif len(dtr) >= len(itr):
            genome.tr.type = "DTR"
            genome.tr.seq = dtr
            genome.tr.length = len(dtr)
        else:
            genome.tr.type = "ITR"
            genome.tr.seq = itr
            genome.tr.length = len(itr)

    # Filter terminal repeats and calculate flags
    for genome in genomes.values():
        if genome.tr.type is not None:
            mode_base, mode_count = collections.Counter(genome.tr.seq).most_common(1)[0]
            genome.tr.mode_freq = 1.0 * mode_count / len(genome.tr.seq)
            genome.tr.n_freq = 1.0 * genome.tr.seq.count("N") / len(genome.tr.seq)
            genome.tr.count = genome.seq.count(genome.tr.seq)
            
            # Flag based on user-defined thresholds
            if genome.tr.n_freq > tr_max_ambig:
                genome.flagged = True
                genome.reason = "Too many ambiguous bases in TR"
            elif genome.tr.count > tr_max_count:
                genome.flagged = True
                genome.reason = "Repetitive TR sequence"
            elif genome.tr.mode_freq > tr_max_basefreq:
                genome.flagged = True
                genome.reason = "Low complexity TR"
            else:
                genome.flagged = False
        
        # Check for kmer frequency
        if genome.kmer_freq > kmer_max_freq:
            genome.flagged = True
            genome.reason = "Multiple genome copies detected"

    # Write results to output
    with open(output_path, "w") as out:
        header = [
            "contig", "contig_length", "kmer_freq", "prediction_type",
            "repeat_length", "repeat_count", "repeat_n_freq", 
            "repeat_mode_base_freq", "repeat_seq"
        ]
        out.write("\t".join(header) + "\n")
        for genome in genomes.values():
            if genome.tr.type is not None:
                row = [
                    genome.id, genome.length, genome.kmer_freq, genome.tr.type,
                    genome.tr.length, genome.tr.count, genome.tr.n_freq, 
                    genome.tr.mode_freq, genome.tr.seq
                ]
                out.write("\t".join(map(str, row)) + "\n")
    
    logging.info(f"Number of sequences checked: {len(genomes)}")
    logging.info(f"Number of circular sequences detected: {len([g for g in genomes.values() if g.tr.type is not None])}")
    logging.info("Circularity check completed.")

if __name__ == "__main__":
    main()