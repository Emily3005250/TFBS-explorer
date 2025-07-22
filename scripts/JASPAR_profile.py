import argparse
import csv
from pyjaspar import jaspardb
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
import sys

csv.field_size_limit(sys.maxsize)

# Fetch motifs from JASPAR
def fetch_motifs_from_jaspar(tf_list, release="JASPAR2024"):
    jdb = jaspardb(release=release)
    motif_dict = {}

    for tf in tf_list:
        motif_list = jdb.fetch_motifs_by_name(tf)
        if motif_list:
            motif_dict[tf] = motif_list[0]
            print(f"Fetched motif for {tf}: {motif_list[0].matrix_id}")
        else:
            print(f"No motif found for {tf}")
    return motif_dict

# Scan sequences
def scan_sequences(csv_file, motif_dict, threshold):
    results = []

    with open (csv_file, 'r') as csv_file:
        next(csv_file)
        for line in csv_file:
            data = line.split(',')
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

            for tf, motif in motif_dict.items():
                pssm_fwd = motif.pssm
                motif_rev = motif.reverse_complement()
                pssm_rev = motif_rev.pssm
                motif_len = len(motif)
 
                for pos, score in pssm_fwd.search(whole_sequence, threshold = threshold): # Iterate through matches
                    start = pos 
                    end = pos + motif_len
                    if start < 0 or end > len(whole_sequence):
                        continue # Skip when start value is negative or end value is over the promoter sequence length.
                    results.append([gene_name, gene_id, tf,"+", start, end, 'Genomic']) # Add matched data into empty list to write CSV file

                for pos, score in pssm_rev.search(whole_sequence, threshold = threshold):
                    start = pos
                    end = pos + motif_len
                    if start < 0 or end > len(whole_sequence):
                        continue
                    results.append([gene_name, gene_id, tf, '-', start, end, 'Genomic'])

                for pos, score in pssm_fwd.search(transcript_sequence, threshold = threshold):
                    start = pos
                    end = pos + motif_len
                    if start < 0 or end > len(transcript_sequence):
                        continue
                    results.append([gene_name, gene_id, '+', start, end, 'Transcript'])

                for pos, score in pssm_rev.search(transcript_sequence, threshold = threshold):
                    start = pos
                    end = pos + motif_len
                    if start < 0 or end > len(transcript_sequence):
                        continue
                    results.append([gene_name, gene_id, '-', start, end, 'Transscript'])

    return results

# Write CSV
def write_results_to_csv(results, output_file):
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_Name", "Gene_ID", "TF", "Strand", "Start", "End", "Sequence_Type"])
        writer.writerows(results)
    print(f"[INFO] Saved results to {output_file}")

# Main CLI logic
def main():
    parser = argparse.ArgumentParser(description="Scan DNA sequences for TF motifs using JASPAR profile and SeqIO PSSMs.")

    parser.add_argument('--input', '-i', required=True, help="Input CSV file with promoter sequences.")
    parser.add_argument('--transcription_factor', '-t', nargs='+', required=True, help="List of Transcription Factors to scan")
    parser.add_argument('--sensetivity', '-s', type=float, default=7.0, help="PSSM match threshold")
    parser.add_argument('--output', '-o', default="scan_results.csv", help="Output CSV filename")
    parser.add_argument('--jaspar_release', default="JASPAR2024", help="JASPAR version")

    args = parser.parse_args()

    print(f"Scanning {args.input} for motifs: {', '.join(args.tf)} at threshold {args.threshold}")
    motifs_dict = fetch_motifs_from_jaspar(args.tf, release=args.jaspar_release)
    scan_results = scan_sequences(args.input, motifs_dict, threshold=args.threshold)
    write_results_to_csv(scan_results, args.output)

if __name__ == "__main__":
    main()