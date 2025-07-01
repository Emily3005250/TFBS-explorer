
# === plot_tfbs_hits_manual.py ===

import argparse
import os
import matplotlib.pyplot as plt

def count_threshold_hits(directory, prefix):
    """Count motif hits from JASPAR threshold results in a directory."""
    counts = []
    for threshold in range(1, 16):
        filename = os.path.join(directory, f"{prefix}_scan_jaspar_{threshold}.0.csv")
        if os.path.exists(filename):
            count = sum(1 for _ in open(filename)) - 1  # skip header
        else:
            count = 0
        counts.append(count)
    return counts

def main():
    parser = argparse.ArgumentParser(description="Plot TFBS match counts with threshold lines and manual exact/regex inputs.")

    parser.add_argument('--pos_dir', required=True, help="Directory of JASPAR positive match files (pos_scan_jaspar_X.0.csv)")
    parser.add_argument('--neg_dir', required=True, help="Directory of JASPAR negative match files (neg_scan_jaspar_X.0.csv)")

    parser.add_argument('--exact_pos', type=int, required=True, help="Exact match count for positive control")
    parser.add_argument('--exact_neg', type=int, required=True, help="Exact match count for negative control")
    parser.add_argument('--regex_pos', type=int, required=True, help="Regex match count for positive control")
    parser.add_argument('--regex_neg', type=int, required=True, help="Regex match count for negative control")

    parser.add_argument('-o', '--output', default='match_plot.png', help="Output filename for the plot")

    args = parser.parse_args()

    # Step 1: Count JASPAR threshold hits
    pos_hits = count_threshold_hits(args.pos_dir, prefix="pos")
    neg_hits = count_threshold_hits(args.neg_dir, prefix="neg")
    thresholds = list(range(1, 16))

    # Step 2: Plot
    plt.figure(figsize=(10, 6))

    # JASPAR threshold curves
    plt.plot(thresholds, pos_hits, label="JASPAR Positive", color="green", marker='o')
    plt.plot(thresholds, neg_hits, label="JASPAR Negative", color="red", marker='o')

    # Exact/Regex match overlays
    plt.plot([1, 15], [args.exact_pos] * 2, linestyle='--', color='green', label="Exact Positive")
    plt.plot([1, 15], [args.regex_pos] * 2, linestyle=':', color='green', label="Regex Positive")
    plt.plot([1, 15], [args.exact_neg] * 2, linestyle='--', color='red', label="Exact Negative")
    plt.plot([1, 15], [args.regex_neg] * 2, linestyle=':', color='red', label="Regex Negative")

    # Final styling
    plt.title("TFBS Match Distribution: JASPAR Thresholds + Exact & Regex (Manual Count)")
    plt.xlabel("PSSM Threshold")
    plt.ylabel("Motif Match Count")
    plt.xticks(thresholds)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output)
    print(f"[INFO] Plot saved to {args.output}")

if __name__ == "__main__":
    main()