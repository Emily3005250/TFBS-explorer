
import argparse
import re
import csv
import sys 

csv.field_size_limit(sys.maxsize)

def scan_with_regex(input_file, output_file):
    # Define the regex motif patterns
    forward_pattern = re.compile(r'GAAA[ATGC]{2}GAAA')
    reverse_pattern = re.compile(r'TTTC[ATGC]{2}TTTC')

    results = []

    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)
        for data in reader:
            gene_name = data[0]
            gene_id = data[1]
            upstream = data[2]
            gene = data[3]
            downstream = data[4]
            utr5 = data[5]
            cds = data[6]
            utr3 = data[7]

            whole_sequence = upstream + gene + downstream # Contain Introns
            transcript_sequence = utr5 + cds + utr3 # With out Introns - mRNA

            # Scan the matches forward pattern
            for match in forward_pattern.finditer(whole_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Genomic'])

            # Scan the matches reverse pattern
            for match in reverse_pattern.finditer(whole_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Genomic'])

            for match in forward_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Transcript'])

            for match in reverse_pattern.finditer(transcript_sequence):
                start, end = match.start(), match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Transcript'])

    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['Gene_Name', 'Gene_ID', 'Strand', 'Start', 'End', 'Sequence_Type'])
        writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description='Scan TF binding site using regex (forward + reverse pattern).')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file from Ensembl')
    parser.add_argument('--output', '-o', required=True, help='Output CSV file name after scanning')

    args = parser.parse_args()

    scan_with_regex(args.input, args.output)

if __name__ == '__main__':
    main()