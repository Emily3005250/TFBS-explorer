
import argparse
import csv
from pyjaspar import jaspardb
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq

# === Fetch motifs ===
def fetch_motifs_from_jaspar(tf_list, release="JASPAR2024"):
    jdb = jaspardb(release=release)
    motif_dict = {}
    for tf in tf_list:
        motif_list = jdb.fetch_motifs_by_name(tf)
        if motif_list:
            jaspar_motif = motif_list[0]
            motif_dict[tf] = jaspar_motif
            print(f"[INFO] Fetched motif for {tf}: {jaspar_motif.matrix_id}")
        else:
            print(f"[WARNING] No motif found for {tf}")
    return motif_dict

# Scan sequences
def scan_promoters_for_motifs(csv_file, motif_dict, threshold=3.0):
    results = []
    for record in SeqIO.parse(csv_file, "csv"):
        gene_info = record.id
        promoter_seq = str(record.seq)

        gene_name, gene_id = gene_info.split("_", 1) if "_" in gene_info else (gene_info, "NA")

        for tf, motif in motif_dict.items():
            motif_len = len(motif)
            pssm = motif.pssm
            forward_matches = list(pssm.search(promoter_seq, threshold=threshold))
            for pos, score in forward_matches:
                start = pos
                end = pos + motif_len
                if start < 0 or end > len(promoter_seq):
                    continue
                results.append([gene_name, gene_id, tf, "+", start, end])
    return results

# Write output
def write_results_to_csv(results, filename):
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Gene_Name", "Gene_ID", "TF", "Strand", "Start", "End"])
        writer.writerows(results)
    print(f"[INFO] Results written to: {filename}")

# === CLI entry ===
def main():
    parser = argparse.ArgumentParser(description="Scan promoter regions for TF motifs using JASPAR.")
    parser.add_argument("--csv", required=True, help="Input FASTA file with promoter sequences.")
    parser.add_argument("--tf", nargs="+", required=True, help="List of transcription factors to scan (e.g. IRF3 IRF7).")
    parser.add_argument("--threshold", type=float, default=3.0, help="PSSM threshold (default: 3.0)")
    parser.add_argument("--output", default="scan_results.csv", help="Output CSV filename (default: scan_results.csv)")
    parser.add_argument("--jaspar_release", default="JASPAR2024", help="JASPAR release version (default: JASPAR2024)")

    args = parser.parse_args()

    motif_dict = fetch_motifs_from_jaspar(args.tf, release=args.jaspar_release)
    results = scan_promoters_for_motifs(args.fasta, motif_dict, threshold=args.threshold)
    write_results_to_csv(results, args.output)

if __name__ == "__main__":
    main()