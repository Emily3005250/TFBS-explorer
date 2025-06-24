
#### Scanning Exact Matches in Promoter Region ####

import requests
import csv
import argparse
import os

def fetch_sequences(input_file, output_file ,upstream, downstream, species):
    # Set the dict to store gene_id and gene_name
    gene_ids = {}

    with open(input_file, "r") as file:
        next(file)
        for line in file:
            data = line.split(",")
            gene_id = data[1]
            gene_name = data[2]
            gene_ids[gene_name] = gene_id

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
            # Set the endpoint to get 5' end expanded sequences
            endpoint = f"/sequence/id/{gene_id}?expand_5prime={upstream}&expand_3prime={downstream}"
            url = server + endpoint

            r = requests.get(url, headers=headers)

            if not r.ok:
                print(f"Error while retrieve {gene_id}: {r.status_code}") 
                continue
            
            data = r.json() # Store extract data in json format
            sequence = data.get('seq', '') # extract the sequence from data in json format
            
            promoter_seq = sequence[ :upstream]
            downstream_seq = sequence[-downstream:]
            gene_seq = sequence[upstream:-downstream]
            
            writer.writerow([species, gene_name, gene_id, promoter_seq, gene_seq ,downstream_seq])
    
    print(f'Save sequence to {output_file}')

## Initiate function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch the DNA sequence from Ensembl')

    parser.add_argument('-i', '--input_file', required=True, help='Input CSV file with gene information')
    parser.add_argument('-o', '--output_file', required=True, help='Set the output file name')
    parser.add_argument('-u', '--upstream', type=int, default=500, help='Set the length how long')
    parser.add_argument('-d', '--downstream', type=int, default=500, help='Set the length of downstream of gene')
    parser.add_argument('-s', '--species', type=str, required=True, help='Species name')

    args = parser.parse_args()

    fetch_sequences(args.input_file, args.output_file, args.upstream, args.downstream, args.species)