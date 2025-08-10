
import argparse
import re
import pandas as pd
from Bio.Seq import Seq


def scan_with_regex(input_file, motif_file, output_file):
    #Load the motif from the CSV file
    motif_df = pd.read_csv(motif_file)
    if motif_df.empty:
        print("Error: Motif file is empty.")
        return

    # Prepare the output list to store results
    # Each entry will be a list: [species, gene_name, ID, Transcription_Factor, strand, start, end, sequence_type]
    results = []

    # Read the CSV file
    df = pd.read_csv(input_file)

    # Iterate through each row in the DataFrame
    for idx, row in df.iterrows():
        gene_name = row['gene_name']
        gene_id = row['ensembl_id']
        species = row['species']

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
        transcript_sequence = utr5 + cds + utr3 # With out Introns - mRNA

        ## Find matches in Genomic
        if not genomic_sequence:  # Check if genomic_sequence is empty
            print(f"Warning: Empty genomic sequence for {gene_name} ({gene_id})")
            continue

        for _, motif_row in motif_df.iterrows():
            tf = motif_row['Transcription_Factor']
            forward_motif = motif_row['regex_forward']
            reverse_motif = motif_row['regex_reverse']

            # Compile regex patterns for forward and reverse motifs
            forward_pattern = re.compile(forward_motif)
            reverse_pattern = re.compile(reverse_motif)

            if forward_pattern.search(genomic_sequence) or reverse_pattern.search(genomic_sequence): # Check if either motif is present it might save time
                # Scan the matches forward pattern
                for match in forward_pattern.finditer(genomic_sequence):
                    start, end = match.start(), match.end()
                    results.append([species, gene_name, gene_id, tf, '+', start, end, 'Genomic'])

                # Scan the matches reverse pattern
                for match in reverse_pattern.finditer(genomic_sequence):
                    start, end = match.start(), match.end()
                    results.append([species, gene_name, gene_id, tf, '-', start, end, 'Genomic'])
            
            ## Find matches in Transcript(mRNA)
            if forward_pattern.search(transcript_sequence) or reverse_pattern.search(transcript_sequence):
                # Scan the matches forward pattern
                for match in forward_pattern.finditer(transcript_sequence):
                    start = match.start()
                    end = match.end()
                    results.append([species, gene_name, gene_id, tf, '+', start, end, 'Transcript'])

                # Scan the matches reverse pattern
                for match in reverse_pattern.finditer(transcript_sequence):
                    start, end = match.start(), match.end()
                    results.append([species, gene_name, gene_id, tf, '-', start, end, 'Transcript'])

    # Create a DataFrame from the results
    columns = ['Species','Gene_Name', 'ID', 'Transcription_Factor','Strand', 'Start', 'End', 'Sequence_Type']
    results_df = pd.DataFrame(results, columns=columns)
    results_df = results_df.drop_duplicates()  # Remove duplicates if any

    # Save to CSV
    results_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Scan TF binding site using regex (forward + reverse pattern).')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file from Ensembl')
    parser.add_argument('--motif', '-m', required=True, help='CSV file with motif sequences')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file name after scanning')

    args = parser.parse_args()

    scan_with_regex(args.input, args.motif, args.output)

if __name__ == '__main__':
    main()