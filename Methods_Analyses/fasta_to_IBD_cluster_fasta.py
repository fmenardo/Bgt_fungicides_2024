from Bio import SeqIO

def filter_fasta_by_individuals(input_fasta, output_fasta, individuals_list):
    """
    Filters sequences from a fasta file based on a list of individuals.

    :param input_fasta: Path to input fasta file.
    :param output_fasta: Path to output fasta file.
    :param individuals_list: List of individual names to be extracted.
    """
    # Read the fasta file and filter sequences
    filtered_sequences = []
    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # If the sequence id is in the individuals list, add it to the filtered sequences
            if record.id in individuals_list:
                filtered_sequences.append(record)
    
    # Write the filtered sequences to a new fasta file
    with open(output_fasta, "w") as output_file:
        SeqIO.write(filtered_sequences, output_file, "fasta")
    
    print(f"Filtered sequences saved to {output_fasta}")

def read_individuals_from_file(individuals_file):
    """
    Reads individual names from a file.

    :param individuals_file: Path to the file containing the individual names (one name per line).
    :return: List of individual names.
    """
    with open(individuals_file, "r") as file:
        individuals_list = [line.strip() for line in file]
    return individuals_list

# Example usage:
input_fasta = "/shares/menardo.bgt.uzh/nikos/fungicide/resistance/0_data/alignments/het_as_X/Bgt_Eur+_cyp51_p.fa"  # Replace with your input fasta file path
output_fasta = "/shares/menardo.bgt.uzh/nikos/fungicide/resistance/0_data/alignments/het_as_X/Bgt_Eur+_cyp51_p_IBDcluster8.fa"  # Replace with the desired output file path
individuals_file = "/shares/menardo.bgt.uzh/nikos/fungicide/resistance/2_output/cyp51_IBD_cluster8_indlist.txt"  # Replace with your file containing the individual names

# Read the individuals from the file
individuals_list = read_individuals_from_file(individuals_file)

# Filter and save the fasta sequences
filter_fasta_by_individuals(input_fasta, output_fasta, individuals_list)
