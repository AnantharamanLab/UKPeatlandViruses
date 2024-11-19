#!/usr/bin/env python3

def parse_faa_file(file_path, is_vMAG=False, vMAG_name=None):
    """
    Obtains gene-level data from prodigal-formatted input faa(s)

    file_path: input amino-acid fasta filepath
    is_vMAG: treats input file as a vMAG with multiple contigs
            to a genome. Defaults to False.
    vMAG_name: Genome name, if vMAG. Defaults to None.
    
    Returns a generator that yields gene-level data encoded in
    the sequence headers for the input file.
    """

    sequence_data = ''
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if sequence_data:  # Yield previous entry
                    entry['sequence'] = sequence_data.replace('\n', '')
                    yield entry
                    sequence_data = ''

                parts = line.strip().split(' # ')
                contig_and_gene = parts[0].split()[0].lstrip('>')  # Remove '>'
                
                contig = "_".join(contig_and_gene.rsplit('_', 1)[:-1])
                gene_number = contig_and_gene.rsplit('_', 1)[-1]
                genome = vMAG_name if is_vMAG else contig
                try:
                    contig_pos_start, contig_pos_end, frame = parts[1:4]
                    attributes = dict(item.split('=') for item in parts[4].split(';') if item)

                    entry = {
                        'member': contig_and_gene,
                        'contig': contig,
                        'gene_number': gene_number,
                        'genome': genome,
                        'contig_pos_start': contig_pos_start,
                        'contig_pos_end': contig_pos_end,
                        'frame': frame,
                        'is_vMAG': is_vMAG,
                    }
                    entry.update(attributes)
                except ValueError:
                    print(f"WARNING! Did not find gene start, stop, and frame for record: {contig_and_gene}. These values will be left null, which may cause errors down the line.")
                    entry = {
                        'member': contig_and_gene,
                        'contig': contig,
                        'gene_number': gene_number,
                        'genome': genome,
                        'contig_pos_start': None,
                        'contig_pos_end': None,
                        'frame': None,
                        'is_vMAG': is_vMAG,
                    }
            else:
                sequence_data += line.strip()
    
    if sequence_data:  # Don't forget the last entry
        entry['sequence'] = sequence_data.replace('\n', '')
        yield entry