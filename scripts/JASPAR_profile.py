import argparse
import pandas as pd
from pyjaspar import jaspardb
from Bio import motifs

# Fetch motifs from JASPAR
def fetch_motifs_from_jaspar(tf_list, release="JASPAR2024"):
    jdb = jaspardb(release=release)
    motif_dict = {}
    # Fetch motifs for each transcription factor
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

    # Read the CSV file
    df = pd.read_csv(csv_file)
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
        # Check if genomic_sequence is empty
        
        # Scan genomic sequence for motifs
        # Iterate through each transcription factor motif
        # for each transcription factor in the motif dictionary
        # motif_dict is a dictionary where keys are transcription factor names and values are motif objects
        # motif is a Bio.motifs.Motif object
        # pssm_fwd is a Bio.motifs.PositionSpecificScoreMatrix object
        # pssm_rev is a Bio.motifs.PositionSpecificScoreMatrix object
        # motif_len is the length of the motif
        # motif_rev is the reverse complement of the motif
        # pssm_fwd.search() and pssm_rev.search() return an iterator
        # that yields positions and scores of matches in the sequence
        # threshold is the minimum score for a match to be considered valid
        for tf, motif in motif_dict.items():
            pssm_fwd = motif.pssm
            motif_rev = motif.reverse_complement()
            pssm_rev = motif_rev.pssm
            motif_len = len(motif)

            # Scan genomic sequence
            # Search for matches in the genomic sequence using the forward and reverse PSSMs
            # pssm_fwd.search() and pssm_rev.search() return an iterator
            # that yields positions and scores of matches in the sequence
            # threshold is the minimum score for a match to be considered valid
            # Iterate through matches in the genomic sequence
            for pos, score in pssm_fwd.search(genomic_sequence, threshold = threshold): # Iterate through matches
                start = pos 
                end = pos + motif_len
                if start < 0 or end > len(genomic_sequence):
                    continue # Skip when start value is negative or end value is over the promoter sequence length.
                results.append([gene_name, gene_id, tf,"+", start, end, 'Genomic']) # Add matched data into empty list to write CSV file

            for pos, score in pssm_rev.search(genomic_sequence, threshold = threshold):
                start = pos
                end = pos + motif_len
                if start < 0 or end > len(genomic_sequence):
                    continue
                results.append([gene_name, gene_id, tf, '-', start, end, 'Genomic'])

            # Scan transcript sequence

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
                results.append([gene_name, gene_id, '-', start, end, 'Transcript'])

    return results

# Write CSV
def write_results_to_csv(results, output_file):
    columns = ['Gene_Name', 'Gene_ID', 'Transcription_Factor', 'Strand', 'Start', 'End', 'Sequence_Type']
    results_df = pd.DataFrame(results, columns=columns)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Main CLI logic
def main():
    parser = argparse.ArgumentParser(description="Scan DNA sequences for TF motifs using JASPAR profile and SeqIO PSSMs.")

    parser.add_argument('--input', '-i', required=True, help="Input CSV file with promoter sequences.")
    parser.add_argument('--transcription_factor', '-t', nargs='+', required=True, help="List of Transcription Factors to scan")
    parser.add_argument('--sensitivity', '-s', type=float, default=7.0, help="PSSM match threshold")
    parser.add_argument('--output', '-o', default="scan_results.csv", help="Output CSV filename")
    parser.add_argument('--jaspar_release', default="JASPAR2024", help="JASPAR version")

    args = parser.parse_args()

    print(f"Scanning {args.input} for motifs: {', '.join(args.transcription_factor)} at threshold {args.sensitivity}")
    motifs_dict = fetch_motifs_from_jaspar(args.transcription_factor, release=args.jaspar_release)
    scan_results = scan_sequences(args.input, motifs_dict, threshold=args.sensitivity)
    write_results_to_csv(scan_results, args.output)

if __name__ == "__main__":
    main()