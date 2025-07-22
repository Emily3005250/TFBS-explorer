
import argparse
import re
from Bio.Seq import Seq
import pandas as pd

def scan_exact_motif(input_file, motif, output_file):

    # Prepare the regex patterns for forward and reverse motifs
    motif = motif.upper()
    forward_pattern = re.compile(motif)
    reverse_motif = str(Seq(motif).reverse_complement())
    reverse_pattern = re.compile(reverse_motif)

    # Prepare the output list to store results
    # Each entry will be a list: [gene_name, gene_id, strand, start, end, sequence_type]
    results = []

    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Iterate through each row in the DataFrame
    for idx, row in df.iterrows():
        gene_name = row['gene_name']
        gene_id = row['ensembl_id']

        # Extract sequences from the row
        upstream = row['Upstream']
        gene = row['Gene_seq']
        downstream = row['Downstream']
        utr5 = row['5_UTR']
        cds = row['CDS']
        utr3 = row['3_UTR']

        # Combine sequences to form genomic and transcript sequences
        genomic_sequence = upstream + gene + downstream # Contain Introns
        transcript_sequence = utr5 + cds + utr3 # Without Introns

        ## Find matches in Genomic
        if not genomic_sequence:  # Check if genomic_sequence is empty
            print(f"Warning: Empty genomic sequence for {gene_name} ({gene_id})")
            continue
        if forward_pattern.search(genomic_sequence) or reverse_pattern.search(genomic_sequence):
            # Find all matches in forward
            for match in forward_pattern.finditer(genomic_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Genomic'])

            # Find all matches in reverse
            for match in reverse_pattern.finditer(genomic_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Genomic'])

        ## Find matches in Transcript(mRNA)
        # Forward pattern matches
        if not transcript_sequence:  # Check if transcript_sequence is empty
            print(f"Warning: Empty transcript sequence for {gene_name} ({gene_id})")
            continue
        if forward_pattern.search(transcript_sequence) or reverse_pattern.search(transcript_sequence):
            
            # Find all matches for forward
            for match in forward_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Transcript'])
            
            # Find all matches for reverse
            for match in reverse_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Transcript'])
        
    # Create a DataFrame from the results
    columns = ['Gene_Name', 'Gene_ID', 'Strand', 'Start', 'End', 'Sequence_Type']
    results_df = pd.DataFrame(results, columns=columns)

    # Save to CSV
    results_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Scan exact match in DNA sequence')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file - DNA sequence data for scanning')
    parser.add_argument('--motif', '-m', required=True, help='Input motif sequence')
    parser.add_argument('--output', '-o', required=True, help='Assign the output file name')

    args = parser.parse_args()

    scan_exact_motif(args.input, args.motif, args.output)

if __name__ == '__main__':
    main()