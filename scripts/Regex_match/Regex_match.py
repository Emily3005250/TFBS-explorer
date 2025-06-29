
import argparse
import re
import csv

def scan_with_regex(input_file, output_file):
    # Define the regex motif patterns
    forward_pattern = re.compile(r'GAAA[ATGC]{2}GAAA')
    reverse_pattern = re.complie(r'TTTC[ATGC]{2}TTTC')

    results = []

    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            gene_name = row[0]
            gene_id = row[1]
            sequence = row[2]

            # Scan the matches forward pattern
            for match in forward_pattern.finditer(sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '+', start, end])

            # Scan the matches reverse pattern
            for match in reverse_pattern.finditer(sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '-', start, end])

    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['Gene_Name', 'Gene_ID', 'Strand', 'Start', 'End'])
        writer.writerow(results)

def main():
    parser = argparse.ArgumentParser(description='Scan TF binding site using regex (forward + reverse pattern).')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file from Ensembl')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file name after scanning')

    args = parser.parse_args()

    scan_with_regex(args.input, args.output)

if __name__ == '__main__':
    main()