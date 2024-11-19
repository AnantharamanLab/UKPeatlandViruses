import os

# Define rules to run by default
rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "build_hmm_profiles.done")

# Create database from input proteins
rule create_db:
    # TODO: Direct mmseqs output to null or specific snakemake .log file unless in debug mode for logger
    output:
        protein_db = os.path.join(config["paths"]["output_dir"], "protein_db", "protein_db")
    params:
        input_files = config["input_protein_files"],
        out_db_dir = os.path.join(config["paths"]["output_dir"], "protein_db")
    message:
        "Executing command: mmseqs createdb {params.input_files} {output.protein_db}"
    shell:
        """
        mkdir -p {params.out_db_dir}
        mmseqs createdb {params.input_files} {output.protein_db}
        """

# Run an mmseqs search on the database against itself
rule search_db:
    # TODO: Direct mmseqs output to null or specific snakemake .log file unless in debug mode for logger
    input:
        protein_db = os.path.join(config["paths"]["output_dir"], "protein_db", "protein_db")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "protein_search_db.done"))
    threads:
        config["threads"]
    params:
        wdir_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "search"),
        search_parent = os.path.join(config["paths"]["output_dir"], "protein_search_db"),
        search_db = os.path.join(config["paths"]["output_dir"], "protein_search_db", "protein_search_db"),
        search_sensitivity = config["mmseqs_params"]["sensitivity"],
        min_seq_id = 0, # Keep minimum sequence identity at 0.0 since we want to sample the entire range of sequence identity
        cov_fraction = config["mmseqs_params"]["cov_fraction"]
    message:
        "Executing command: mmseqs search {input.protein_db} {input.protein_db} {params.search_db} {params.wdir} -s {params.search_sensitivity} --min-seq-id {params.min_seq_id} -c {params.cov_fraction} --threads {threads}"
    shell:
        """
        mkdir -p {params.wdir_parent}
        mkdir -p {params.search_parent}
        mmseqs search {input.protein_db} {input.protein_db} {params.search_db} {params.wdir} -s {params.search_sensitivity} --min-seq-id {params.min_seq_id} -c {params.cov_fraction} --threads {threads}
        """

# Extract the output table from mmseqs search
rule extract_search_output:
    # TODO: Direct mmseqs output to null or specific snakemake .log file unless in debug mode for logger
    input:
        search_db = os.path.join(config["paths"]["output_dir"], "wdir", "protein_search_db.done"),
        protein_db = os.path.join(config["paths"]["output_dir"], "protein_db", "protein_db")
    output:
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.tsv")
    threads:
        config["threads"]
    params:
        search_db = os.path.join(config["paths"]["output_dir"], "protein_search_db", "protein_search_db")
    message:
        "Executing command: mmseqs convertalis {input.protein_db} {input.protein_db} {params.search_db} {output.search_output} --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov --threads {threads}"
    shell:
        """
        mmseqs convertalis {input.protein_db} {input.protein_db} {params.search_db} {output.search_output} --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov --threads {threads}
        """

# Preprocess search output for AAI calculations
rule preprocess_aai:
    input:
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.tsv")
    output:
        gene_map_file = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv")
    threads:
        config["threads"]
    params:
        input_files = config["input_protein_files"],
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai")
    message:
        "Preparing gene-to-genome map for AAI calculations."
    script:
        os.path.join(config["paths"]["scripts_dir"], "preprocess_genome_aai.py")

# Compute the amino-acid identity of all proteins in the database
rule compute_aai:
    # TODO: Direct mmseqs output to null or specific snakemake .log file unless in debug mode for logger
    input:
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.tsv"),
        gene_map_file = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv")
    output:
        aai = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "genome_aai.tsv")
    threads:
        config["threads"]
    params:
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai")
    message:
        "Computing pairwise amino-acid identities of all proteins in the database."
    script:
        os.path.join(config["paths"]["scripts_dir"], "compute_genome_aai.py")

# Cluster genomes using amino-acid identities
rule cluster_aai:
    # TODO: Direct mmseqs output to null or specific snakemake .log file unless in debug mode for logger
    input:
        aai = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "genome_aai.tsv")
    output:
        genome_clusters = os.path.join(config["paths"]["output_dir"], "wdir", "genome_clusters.tsv")
    threads:
        config["threads"]
    params:
        gene_map = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv"),
        cluster_taxa_levels = config["cluster_ranks"],
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai")
    message:
        "Clustering genomes based on AAI results in {input.aai} at the taxonomic rank(s): {params.cluster_taxa_levels}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "cluster_genome_aai.py")
