
#### Scanning Exact Matches in Promoter Region ####

# Import the libraries
import requests
import csv
import argparse
import os

def fetch_sequences(input_file, output_file ,upstream, downstream, species):
    # Set the dict to store gene_id and gene_name
    gene_ids = {}

    # Read input file with 'r' - not change the original file
    with open(input_file, "r") as file:
        next(file) # Skip the header row
        for line in file:
            data = line.split(",") # Split the row by comma
            gene_id = data[1] # Extract the gene_id from file
            gene_name = data[2] # Extract the gene_name from the file
            gene_ids[gene_name] = gene_id # Store the gene name and gene id to dict

    # Print this message to confirm fectching
    print(f"Loaded {len(gene_ids)} genes from {input_file}.")

    # Set the server database to retrieve
    server = "https://rest.ensembl.org"
    headers = {"Accept": "application/json"}

    # write a csv file for promoter sequences
    with open("output_file", "w", newline = "") as csvfile:
        writer = csv.writer(csvfile) 
        writer.writerow(['Species', 'Gene_Name', 'Ensembl_ID', 'Sequence']) # Set the header row for CSV file
        
        # Iterate through dict to get the sequece data from Ensembl
        for gene_name, gene_id in gene_ids.items():
            # Set the endpoint to get 5' and 3' end expanded sequences
            endpoint = f"/sequence/id/{gene_id}?expand_5prime={upstream}&expand_3prime={downstream}"
            url = server + endpoint

            r = requests.get(url, headers=headers)

            if not r.ok:
                print(f"Error while retrieve {gene_id}: {r.status_code}") 
                continue
            
            data = r.json() # Store extract data in json format
            sequence = data.get('seq', '') # extract the sequence from data in json format
            
            # Trimming the sequence by upstream (promoter), gene, downstream
            promoter_seq = sequence[ :upstream]
            downstream_seq = sequence[-downstream:]
            gene_seq = sequence[upstream:-downstream]
            
            # Confirm the sequence length before writing CSV file
            assert len(promoter_seq) == upstream, f"Promoter length mismatch for {gene_id}"
            assert len(downstream_seq) == downstream, f"Downstream length mismatch for {gene_id}"
            assert len(gene_seq) == (len(sequence) - upstream - downstream), f"Gene sequence length mismatch for {gene_id}"

            # Write the csv file with data
            writer.writerow([species, gene_name, gene_id, promoter_seq, gene_seq ,downstream_seq])
    
    # Print to confirm
    print(f'Save sequence to {output_file}')

## Initiate function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch the DNA sequence from Ensembl')

    ## Add argumentn to parser
    parser.add_argument('-i', '--input_file', required=True, help='Input CSV file with gene information') # Provide the input file name
    parser.add_argument('-o', '--output_file', required=True, help='Set the output file name') # Provide the output file name
    parser.add_argument('-u', '--upstream', type=int, default=500, help='Set the length how long') # Provide the length we want to fetch
    parser.add_argument('-d', '--downstream', type=int, default=500, help='Set the length of downstream of gene') # Provide the length we want to fetch
    parser.add_argument('-s', '--species', type=str, required=True, help='Species name') # Provide species information

    # Set the args to use the argument with order (index)
    args = parser.parse_args()

    fetch_sequences(args.input_file, args.output_file, args.upstream, args.downstream, args.species)