import os

# Define rules to run by default
rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "write_results.done")

# Check circulatiry of user genomes
rule check_circular:
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "check_circular.done"))
    params:
        input_single_contig_genomes = config["input_single_contig_genomes"],
        input_vmag_fastas = config["input_vmag_fastas"],
        tr_min_len = 20,
        tr_max_len = 1000,
        tr_max_count = 8,
        tr_max_ambig = 0.2,
        tr_max_basefreq = 0.70,
        kmer_max_freq = 1.5,
        k = 15,
        circularity_tbl = os.path.join(config["paths"]["output_dir"], "wdir", "circular_contigs.tsv")
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Checking the circularity of input sequences and writing results to {params.circularity_tbl}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "check_circular.py")

# Annotate user genomes
rule run_pyrodigal_gv:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "check_circular.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "run_pyrodigal_gv.done"))
    params:
        input_single_contig_genomes = config["input_single_contig_genomes"],
        input_vmag_fastas = config["input_vmag_fastas"],
        wdir = os.path.join(config["paths"]["output_dir"], "wdir"),
        output_dir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv"),
        single_contig_prots = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv", "single_contig_proteins.faa"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv", "vMAG_proteins")),
        gene_to_genome = os.path.join(config["paths"]["output_dir"], "wdir", "gene_to_genome.txt")
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Predicting genes in input genomes with pyrodigal-gv & translating"
    script:
        os.path.join(config["paths"]["scripts_dir"], "run_pyrodigal_parallel.py")

# Obtain gene information from input (prodigal-formatted) .faa and genome information from...
rule index_genes:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "run_pyrodigal_gv.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "index_genes.done"))
    params:
        cluster_taxa_levels = None,
        database = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        db = "",
        tsv = "",
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv", "vMAG_proteins")),
        out_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
    threads:
        config["threads"]
    message:
        "Adding gene- and genome-level data from pyrodigal-gv and writing to {params.database}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "map_protein_data.py")

# Assign V-scores and L-scores to the proteins in the database
rule assign_vscores:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "index_genes.done")
    output:
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv")
    params:
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv"),
        hmm_vscores = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        db_dir = config["paths"]["db_dir"],
        cov_fraction = config["cov_fraction"],
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Assigning V-scores and L-scores to proteins in {params.protein_dir} using an HMMsearch of the annotations in {params.hmm_vscores}"
    script:
        os.path.join(config["paths"]["scripts_dir"], "annotate.py")

# Merge the calcualted v-scores and L-scoe with the protein database
rule add_vscores:
    input:
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "add_annots.done"))
    params:
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        db_dir = config["paths"]["db_dir"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Adding calculated v-scores in {input.vscores} to {params.gene_index} and writing to {params.gene_index_annotated}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "add_annots.py")

# Analyze the genomic context of V-scores/L-scores
rule genome_context:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "add_annots.done")
    output:
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv")
    params:
        outparent = os.path.join(config["paths"]["output_dir"], "results"),
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        circular_contigs = os.path.join(config["paths"]["output_dir"], "wdir", "circular_contigs.tsv"),
        annotation_percent_threshold = config["annotation_percent_threshold"],
        min_window_avg_lscore = config["min_window_avg_lscore"],
        window_size = config["window_size"],
        minimum_flank_vscore = config["minimum_flank_vscore"],
        max_flank_length = config["max_flank_length"],
        vscore_ref = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        tmp_dir = os.path.join(config["paths"]["output_dir"], "wdir")
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Analyzing the genomic context of V-scores and L-scores in {params.gene_index_annotated} and writing results to {output.context_table}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "genome_context.py")

# Write final results files
rule write_results:
    input:
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "write_results.done"))
    params:
        aux_table = os.path.join(config["paths"]["output_dir"], "results", "aux_results.tsv"),
        all_annot_out_table = os.path.join(config["paths"]["output_dir"], "results", "all_annotations.tsv"),
        aux_fasta_dir = os.path.join(config["paths"]["output_dir"], "results"),
        hmm_ref = os.path.join(config["paths"]["files_dir"], "hmm_id_to_name.csv"),
        vibrant_amgs = os.path.join(config["paths"]["files_dir"], "VIBRANT_AMGs_filtered.tsv"),
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv")
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Writing the final auxiliary protein results to {params.aux_table} and writing auxiliary sequences to {params.aux_fasta_dir}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "write_aux_results.py")
