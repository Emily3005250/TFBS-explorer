import argparse
import csv
from pyjaspar import jaspardb
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq

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
            gene_name = data[1]
            gene_id = data[2]
            sequence = data[3]

            for tf, motif in motif_dict.items():
                pssm_fwd = motif.pssm
                motif_rev = motif.reverse_complement()
                pssm_rev = motif_rev.pssm
                motif_len = len(motif)

                forward_matches = list(pssm_fwd.search(sequence, threshold=threshold)) # Scan Promoter sequences with forward motif patterns 
                for pos, score in forward_matches: # Iterate through matches
                    start = pos 
                    end = pos + motif_len
                    if start < 0 or end > len(sequence):
                        continue # Skip when start value is negative or end value is over the promoter sequence length.
                    results.append([gene_name, gene_id, tf,"+", start, end]) # Add matched data into empty list to write CSV file

                reverse_matches = list(pssm_rev.search(sequence, threshold=threshold))
                for pos, score in reverse_matches:
                    start = pos
                    end = pos + motif_len
                    if start < 0 or end > len(sequence):
                        continue
                    results.append([gene_name, gene_id, tf, '-', start, end])

    return results

# Write CSV
def write_results_to_csv(results, output_file):
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_Name", "Gene_ID", "TF", "Strand", "Start", "End"])
        writer.writerows(results)
    print(f"[INFO] Saved results to {output_file}")

# Main CLI logic
def main():
    parser = argparse.ArgumentParser(description="Scan DNA sequences for TF motifs using JASPAR profile and SeqIO PSSMs.")

    parser.add_argument('--input', '-i', required=True, help="Input CSV file with promoter sequences.")
    parser.add_argument('--tf', nargs='+', required=True, help="List of Transcription Factors to scan")
    parser.add_argument('--threshold', '-t', type=float, default=7.0, help="PSSM match threshold")
    parser.add_argument('--output', '-o', default="scan_results.csv", help="Output CSV filename")
    parser.add_argument('--jaspar_release', default="JASPAR2024", help="JASPAR version")

    args = parser.parse_args()

    print(f"Scanning {args.input} for motifs: {', '.join(args.tf)} at threshold {args.threshold}")
    motifs_dict = fetch_motifs_from_jaspar(args.tf, release=args.jaspar_release)
    scan_results = scan_sequences(args.input, motifs_dict, threshold=args.threshold)
    write_results_to_csv(scan_results, args.output)

if __name__ == "__main__":
    main()