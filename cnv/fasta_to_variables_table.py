import argparse
import csv
from Bio import SeqIO

# Initialize the argument parser
parser = argparse.ArgumentParser(description="Parse protein and nucleotide FASTA files and extract values at specified positions into CSV files")

# Add the input protein FASTA file argument
parser.add_argument('-amino', '--input_amino_fasta', metavar='', default='', help='input protein fasta file', nargs=1, type=str)
# Add the input nucleotide FASTA file argument
parser.add_argument('-nucl', '--input_nucl_fasta', metavar='', default='', help='input nucleotide fasta file', nargs=1, type=str)
# Add the output protein CSV file argument
parser.add_argument('-o', '--output_csv', metavar='', default='output.csv', help='output protein CSV file', nargs=1, type=str)
# Add the output nucleotide CSV file argument
parser.add_argument('-co', '--codon_output_csv', metavar='', default='codon_output.csv', help='output codon CSV file', nargs=1, type=str)
# Add the positions argument (1-based indexing)
parser.add_argument('-p', '--positions', metavar='', help='positions to extract (1-based indexing)', nargs='+', type=int)

# Parse the arguments
arguments = parser.parse_args()

# Extract the input files and output CSV files from arguments
input_amino_fasta = arguments.input_amino_fasta[0]
input_nucl_fasta = arguments.input_nucl_fasta[0]
output_csv = arguments.output_csv[0]
codon_output_csv = arguments.codon_output_csv[0]

# Convert positions to 0-based indexing
positions = [pos - 1 for pos in arguments.positions]

# Prepare the CSV headers
header = ['ID'] + [f'Position_{pos + 1}' for pos in positions]
codon_header = ['ID'] + [f'Codon_Position_{pos + 1}' for pos in positions]

# Dictionary to store the nucleotide sequences
nucl_sequences = {record.id: record.seq for record in SeqIO.parse(input_nucl_fasta, "fasta")}

# Open the output CSV files
with open(output_csv, 'w', newline='') as csvfile, open(codon_output_csv, 'w', newline='') as codon_csvfile:
    csvwriter = csv.writer(csvfile)
    codon_csvwriter = csv.writer(codon_csvfile)
    
    csvwriter.writerow(header)  # Write the header row for amino acids
    codon_csvwriter.writerow(codon_header)  # Write the header row for codons
    
    # Parse the protein FASTA file and extract values at specified positions
    for record in SeqIO.parse(input_amino_fasta, "fasta"):
        sequence = record.seq
        name = record.id
        
        # Write the amino acid information
        aa_row = [name]
        for pos in positions:
            if len(sequence) > pos:
                aa_row.append(sequence[pos])
            else:
                aa_row.append("N/A")  # Handle cases where the sequence is shorter than the position
        csvwriter.writerow(aa_row)  # Write the amino acid data row
        
        # Write the codon information
        codon_row = [name]
        nucl_seq = nucl_sequences.get(name, "")
        for pos in positions:
            nucl_start = pos * 3
            nucl_end =  nucl_start +3
            if len(nucl_seq) >= nucl_end:
                codon_row.append(nucl_seq[nucl_start:nucl_end])
            else:
                codon_row.append("N/A")  # Handle cases where the nucleotide sequence is shorter than needed
        codon_csvwriter.writerow(codon_row)  # Write the codon data row
