
import argparse
import re
import pandas as pd


def scan_with_regex(input_file, output_file):
    # Define the regex motif patterns
    forward_pattern = re.compile(r'GAAA[ATGC]{2}GAAA')
    reverse_pattern = re.compile(r'TTTC[ATGC]{2}TTTC')

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
        transcript_sequence = utr5 + cds + utr3 # With out Introns - mRNA

        ## Find matches in Genomic
        if not genomic_sequence:  # Check if genomic_sequence is empty
            print(f"Warning: Empty genomic sequence for {gene_name} ({gene_id})")
            continue
        if forward_pattern.search(genomic_sequence) or reverse_pattern.search(genomic_sequence):
            # Scan the matches forward pattern
            for match in forward_pattern.finditer(genomic_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Genomic'])

            # Scan the matches reverse pattern
            for match in reverse_pattern.finditer(genomic_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Genomic'])
        
        ## Find matches in Transcript(mRNA)
        if not transcript_sequence:  # Check if transcript_sequence is empty
            print(f"Warning: Empty transcript sequence for {gene_name} ({gene_id})")
            continue
        if forward_pattern.search(transcript_sequence) or reverse_pattern.search(transcript_sequence):
            # Scan the matches forward pattern
            for match in forward_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Transcript'])

            # Scan the matches reverse pattern
            for match in reverse_pattern.finditer(transcript_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Transcript'])

    # Create a DataFrame from the results
    columns = ['Gene_Name', 'Gene_ID', 'Strand', 'Start', 'End', 'Sequence_Type']
    results_df = pd.DataFrame(results, columns=columns)

    # Save to CSV
    results_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Scan TF binding site using regex (forward + reverse pattern).')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file from Ensembl')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file name after scanning')

    args = parser.parse_args()

    scan_with_regex(args.input, args.output)

if __name__ == '__main__':
    main()