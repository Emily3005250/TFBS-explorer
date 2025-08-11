
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

# This script generates an elbow plot of motif matches across thresholds for a given transcription factor (TF).
def count_hits_across_thresholds(tf_name, dir, thresholds):
    pos_hits = []
    neg_hits = []
    # Check if the directories exist
    for threshold in thresholds:
        pos_file = os.path.join(dir, f'pos_jaspar_{tf_name}_{threshold}.0.csv')
        neg_file = os.path.join(dir, f'neg_jaspar_{tf_name}_{threshold}.0.csv')
        # Count the number of hits for each threshold
        if os.path.exists(pos_file):
            pos_df = pd.read_csv(pos_file)
            pos_hits.append(len(pos_df))
        else:
            pos_hits.append(0)
        # Count the number of hits for negative scans
        if os.path.exists(neg_file):
            neg_df = pd.read_csv(neg_file)
            neg_hits.append(len(neg_df))
        else:
            neg_hits.append(0)
    # Return the counts of hits for positive and negative scans    
    return pos_hits, neg_hits

# Function to plot the elbow plot using matplotlib
def plot_elbow(thresholds, pos_hits, neg_hits, tf_name, output_file, pos_exact=None, pos_regex=None, neg_exact=None, neg_regex=None):
    # Create a plot for the elbow plot
    # Set the figure size and plot the positive and negative hits
    plt.figure(figsize=(10, 6))
    plt.plot(thresholds, pos_hits, label='JASPAR Positive Hits', color='green')
    plt.plot(thresholds, neg_hits, label='JASPAR Negative Hits', color='red')

    # Optional lines for exact and regex hits
    # If provided, these will be plotted as horizontal lines
    if pos_exact is not None:
        plt.axhline(y=pos_exact, color='green', linestyle='--', label='Exact Positive Hits')
    if pos_regex is not None:
        plt.axhline(y=pos_regex, color='green', linestyle=':', label='Regex Positive Hits')
    if neg_exact is not None:
        plt.axhline(y=neg_exact, color='red', linestyle='--', label='Exact Negative Hits')
    if neg_regex is not None:
        plt.axhline(y=neg_regex, color='red', linestyle=':', label='Regex Negative Hits')

    # Set the title and labels for the plot
    # Add title, labels, and legend
    plt.title(f'Elbow Plot for {tf_name}')
    plt.xlabel('PSSM Thresholds')
    plt.ylabel('Number of Hits')
    plt.xticks(thresholds)
    plt.legend()
    plt.grid()
    plt.savefig(output_file) # Save the plot to the specified output file
    plt.show() # Display the plot
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate an elbow plot of motif matches across thresholds for Transcription Factor.')
    
    # Required arguments
    parser.add_argument('--tf', '-t', required=True, help='Transcription Factor name')
    parser.add_argument('--dir', '-d', required=True, help='Directory containing JSAPAR positive scan results')
    parser.add_argument('--output', '-o',default='elbow_plot.png', help='Output file for the elbow plot')
    parser.add_argument('--min_thresh', type=int, default=1, help="Minimum PSSM threshold (default: 1)")
    parser.add_argument('--max_thresh', type=int, default=15, help="Maximum PSSM threshold (default: 15)")

    # Optional arguments for exact and regex hits
    parser.add_argument('--pos_exact', type=int, default = None, help='Exact positive hits for the transcription factor')
    parser.add_argument('--pos_regex', type=int, default = None, help='Regex positive hits for the transcription factor')
    parser.add_argument('--neg_exact', type=int, default = None, help='Exact negative hits for the transcription factor')
    parser.add_argument('--neg_regex', type=int, default = None, help='Regex negative hits for the transcription factor')

    # Parse the arguments
    # Parse the command line arguments
    args = parser.parse_args()
    thresholds = range(args.min_thresh, args.max_thresh + 1)
    pos_hits, neg_hits = count_hits_across_thresholds(args.tf, args.dir, thresholds)

    # Generate the elbow plot
    # Call the plot_elbow function to create the elbow plot
    plot_elbow(thresholds, pos_hits, neg_hits, args.tf, args.output,
               pos_exact=args.pos_exact, pos_regex=args.pos_regex,
               neg_exact=args.neg_exact, neg_regex=args.neg_regex)
    
if __name__ == '__main__':
    main()
