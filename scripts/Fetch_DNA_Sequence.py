
#### Fetch DNA sequences from Ensembl database ####

# Import the libraries
import requests
import argparse
import pandas as pd


# Function to retrieve the Transcript_ID
def fetch_transcript_id(gene_id):
    server = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{gene_id}?expand=1"
    headers = {"Accept": "application/json"}

    r = requests.get(server + endpoint, headers=headers)
    if not r.ok:
        print(f"Error : Cannot fetch transcript ID for {gene_id}")
        return None
    
    data = r.json()
    transcripts = data.get("Transcript", [])

    for tx in transcripts:
        if tx.get("is_canonical"):
            return tx["id"]
    print(f"Warning : No canonical transcript found for {gene_id}") # logger for further investigation
    return None

# Main function to make a CSV file (contains gene meta data)
def fetch_sequences(input_file, output_file ,upstream, downstream):
    
    df = pd.read_csv(input_file) # Read the input file with pandas
    modified_df = df[['Species', 'Gene', 'ENSEMBL ID']] # Select the columns we need
    modified_df = modified_df.rename(columns={'Species': 'species', 'Gene': 'gene_name', 'ENSEMBL ID': 'ensembl_id'}) # Rename the columns
    
    # Add new columns to the DataFrame for storing sequences and transcript IDs
    modified_df['Upstream'] = ''
    modified_df['Gene_seq'] = ''
    modified_df['Downstream'] = ''
    modified_df['Transcript_ID'] = ''
    modified_df['5_UTR'] = ''
    modified_df['CDS'] = ''
    modified_df["3_UTR"] = ''

    # Set the server database to retrieve
    server = "https://rest.ensembl.org"
    headers = {"Accept": "application/json"}

    for idx, row in modified_df.iterrows():
        gene_id = row['ensembl_id']  # Get the Ensembl ID from the row
        gene_name = row['gene_name']  # Get the gene name from the row

        url = f"{server}/sequence/id/{gene_id}?expand_5prime={upstream}&expand_3prime={downstream}"
        r = requests.get(url, headers=headers)
        if not r.ok:
            print(f"Error while genomic sequence retrieving for {gene_id}: {r.status_code}")
            continue
        data = r.json()  # Store extract data in json format
        sequence = data.get('seq', '')  # extract the sequence from data in json format
        # Trimming the sequence by upstream (promoter), gene, downstream & Add to csv file
        modified_df.at[idx, 'Upstream'] = upstream_seq = sequence[:upstream]
        modified_df.at[idx, 'Downstream'] = downstream_seq = sequence[-downstream:]
        modified_df.at[idx, 'Gene_seq'] = gene_seq = sequence[upstream:-downstream]

        # Confirm the sequence length before writing CSV file
        assert len(upstream_seq) == upstream, f"Promoter length mismatch for {gene_id}"
        assert len(downstream_seq) == downstream, f"Downstream length mismatch for {gene_id}"
        assert len(gene_seq) == (len(sequence) - upstream - downstream), f"Gene sequence length mismatch for {gene_id}"
        
        # Get the Transcript_ID using Gene_ID
        transcript_id = fetch_transcript_id(gene_id)
        if not transcript_id:
            print(f"Warning : No transcript found for {gene_id} : {gene_name}")
            continue    
        modified_df.at[idx, 'Transcript_ID'] = transcript_id  # Store the transcript ID in the DataFrame
        
        # Set the URL for cDNA and CDS using Transcript ID
        cdna_url = f'{server}/sequence/id/{transcript_id}?type=cdna'
        cds_url = f"{server}/sequence/id/{transcript_id}?type=cds"
        
        # Request to Ensembl to get the data for cDNA and CDS
        cdna_data = requests.get(cdna_url, headers=headers)
        cds_data = requests.get(cds_url, headers=headers)
        if not cdna_data.ok or not cds_data.ok:
            print(f"Warning : Failed to fetch sequences for {transcript_id}")
            continue
        
        # Extract the sequence data
        cdna = cdna_data.json().get("seq", "")
        cds = cds_data.json().get("seq", "")
        
        # Split the cDNA using CDS
        if cds and cds in cdna:
            utr5 = cdna.split(cds)[0]
            utr3 = cdna.split(cds)[-1]
        else:
            utr5, utr3 = "", ""
            print(f"Warning : CDS not found inside cDNA for {gene_id} - {transcript_id}")
        
        # Store the UTR and CDS sequences in the DataFrame
        modified_df.at[idx, '5_UTR'] = utr5
        modified_df.at[idx, 'CDS'] = cds
        modified_df.at[idx, '3_UTR'] = utr3
    
    # Save the modified DataFrame to the output file
    modified_df.to_csv(output_file, index=False)  # Save the modified DataFrame to the output file
    print(f"Sequences fetched and saved to {output_file}")  # Print the message when the sequences are fetched and saved

## Initiate function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch the DNA sequence from Ensembl')

    ## Add argumentn to parser
    parser.add_argument('-i', '--input_file', required=True, help='Input CSV file with gene information') # Provide the input file name
    parser.add_argument('-o', '--output_file', required=True, help='Set the output file name') # Provide the output file name
    parser.add_argument('-u', '--upstream', type=int, default=500, help='Set the length how long') # Provide the length we want to fetch
    parser.add_argument('-d', '--downstream', type=int, default=500, help='Set the length of downstream of gene') # Provide the length we want to fetch

    # Set the args to use the argument with order (index)
    args = parser.parse_args()

    fetch_sequences(args.input_file, args.output_file, args.upstream, args.downstream)