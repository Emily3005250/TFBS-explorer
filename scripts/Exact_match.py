
import argparse
import re
from Bio.Seq import Seq
import pandas as pd

def scan_exact_motif(input_file, motif_file, output_file):

    # Load motifs for all TFs from the CSV file
    motif_df = pd.read_csv(motif_file)
    if motif_df.empty:
        print("Error: Motif file is empty.")
        return

    # Prepare the output list to store results
    # Each entry will be a list: [gene_name, gene_id, strand, start, end, sequence_type]
    results = []

    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Iterate through each row in the DataFrame
    for idx, row in df.iterrows():
        gene_name = row['gene_name']
        gene_id = row['ensembl_id']
        species = row['species']
        transcript_id = row['Transcript_ID'] if 'Transcript_ID' in row else None

        # Extract sequences from the row
        upstream = row['Upstream'] if pd.notna(row['Upstream']) else ''
        gene = row['Gene_seq'] if pd.notna(row['Gene_seq']) else ''
        downstream = row['Downstream'] if pd.notna(row['Downstream']) else ''
        # Extract UTR and CDS sequences
        utr5 = row['5_UTR'] if pd.notna(row['5_UTR']) else ''
        cds = row['CDS'] if pd.notna(row['CDS']) else ''
        utr3 = row['3_UTR'] if pd.notna(row['3_UTR']) else ''

        # Combine sequences to form genomic and transcript sequences
        genomic_sequence = upstream + gene + downstream # Contain Introns
        transcript_sequence = utr5 + cds + utr3 # Without Introns

        ## Find matches in Genomic
        if not genomic_sequence:  # Check if genomic_sequence is empty
            print(f"Warning: Empty genomic sequence for {gene_name} ({gene_id})")
            continue
        
        for _, motif_row in motif_df.iterrows():
            tf = motif_row['Transcription_Factor']
            forward_motif = motif_row['exact_forward']
            reverse_motif = motif_row['exact_reverse']

            # Compile regex patterns for forward and reverse motifs
            forward_pattern = re.compile(forward_motif)
            reverse_pattern = re.compile(reverse_motif)

            if forward_pattern.search(genomic_sequence) or reverse_pattern.search(genomic_sequence):
                # Scan the matches forward pattern
                for match in forward_pattern.finditer(genomic_sequence):
                    start, end = match.start(), match.end()
                    results.append([species, gene_name, gene_id, tf, '+', start, end, 'Genomic'])

                # Scan the matches reverse pattern
                for match in reverse_pattern.finditer(genomic_sequence):
                    start, end = match.start(), match.end()
                    results.append([species, gene_name, gene_id, tf, '-', start, end, 'Genomic'])

        if forward_pattern.search(genomic_sequence) or reverse_pattern.search(genomic_sequence):
            # Find all matches in forward
            for match in forward_pattern.finditer(genomic_sequence):
                start = match.start()
                end = match.end()
                results.append([species, gene_name, gene_id, tf, '+', start, end, 'Genomic'])

            # Find all matches in reverse
            for match in reverse_pattern.finditer(genomic_sequence):
                start = match.start()
                end = match.end()
                results.append([species, gene_name, gene_id, tf, '-', start, end, 'Genomic'])

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
                results.append([species, gene_name, transcript_id, tf, '+', start, end, 'Transcript'])
            
            # Find all matches for reverse
            for match in reverse_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([species, gene_name, transcript_id, tf, '-', start, end, 'Transcript'])
        
    # Create a DataFrame from the results
    columns = ['Species','Gene_Name', 'ID', 'Transcription_Factor', 'Strand', 'Start', 'End', 'Sequence_Type']
    results_df = pd.DataFrame(results, columns=columns)
    results_df = results_df.drop_duplicates()  # Remove duplicates if any

    # Save to CSV
    results_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Scan exact match in DNA sequence')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file - DNA sequence data for scanning')
    parser.add_argument('--motif', '-m', required=True, help='Input motif sequence CSV file to scan')
    parser.add_argument('--output', '-o', required=True, help='Assign the output file name')

    args = parser.parse_args()

    scan_exact_motif(args.input, args.motif, args.output)

if __name__ == '__main__':
    main()