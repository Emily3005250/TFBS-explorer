
import argparse
import csv
import re
from Bio.Seq import Seq

def scan_exact_motif(csv_file, motif, output_file):

    motif = motif.upper()
    forward_pattern = re.compile(motif)
    reverse_motif = str(Seq(motif).reverse_complement())
    reverse_pattern = re.compile(reverse_motif)

    results = []

    with open(csv_file, 'r') as input_file:
        next(input_file)
        for line in input_file:
            data = line.strip().split(',')
            gene_name = data[0]
            gene_id = data[1]
            upstream = data[2]
            gene = data[3]
            downstream = data[4]
            utr5 = data[5]
            cds = data[6]
            utr3 = data[7]

            whole_sequence = upstream + gene + downstream # Contain Introns
            transcript_sequence = utr5 + cds + utr3 # With out Introns

            # Find matches in Genomic
            for match in forward_pattern.finditer(whole_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Genomic'])

            for match in reverse_pattern.finditer(whole_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Genomic'])

            # Find matches in Transcript(mRNA)
            for match in forward_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '+', start, end, 'Transcript'])

            for match in reverse_pattern.finditer(transcript_sequence):
                start = match.start()
                end = match.end()
                results.append([gene_name, gene_id, '-', start, end, 'Transcript'])
            
    
    with open(output_file, 'w', newline="") as output:
        writer = csv.writer(output)
        writer.writerow(['Gene_Name', 'Gene_ID', 'Strand', 'Start','End', 'Sequence_Type'])
        writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description='Scan exact match in DNA sequence')

    parser.add_argument('--input', '-i', required=True, help='Input CSV file - DNA sequence data for scanning')
    parser.add_argument('--motif', '-m', required=True, help='Input motif sequence')
    parser.add_argument('--output', '-o', required=True, help='Assign the output file name')

    args = parser.parse_args()

    scan_exact_motif(args.input, args.motif, args.output)

if __name__ == '__main__':
    main()