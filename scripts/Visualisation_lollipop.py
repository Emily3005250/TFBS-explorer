
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os

## Function to merge multiple CSV files into one
## This function reads multiple CSV files (TFBS scanning result), concatenates them into a single DataFrame,
def merge_csv_files(input_files, output_file):
    dataframes = []
    for file in input_files:
        df = pd.read_csv(file)
        dataframes.append(df)
    
    merged_df = pd.concat(dataframes, ignore_index=True)
    merged_df.to_csv(output_file, index=False)
    
    return merged_df
## Function to parse command line arguments
def generate_lollipop_plot(df, output_dir, upstream = 500, downstream = 500):
    os.makedirs(output_dir, exist_ok=True)

    unique_tfs = sorted(df['Transcription_Factor'].unique())
    cmap = plt.get_cmap('Set1')
    colors = [cmap(i % 9) for i in range(len(unique_tfs))]
    tf_color = dict(zip(unique_tfs, colors))
    heights = np.linspace(0.15, 0.15 * len(unique_tfs), len(unique_tfs))
    tf_height = dict(zip(unique_tfs, heights))

    df['tf_color'] = df['Transcription_Factor'].map(tf_color)
    df['line_height'] = df['Transcription_Factor'].map(tf_height)

    ## Calculate the length of each region
    for (species, gene_id, gene_name), gene_df in df.groupby(['Species', 'Gene_ID', 'Gene_Name']):
        row = gene_df.iloc[0]
        utr5_length = len(str(row['UTR5'])) if pd.notna(row['UTR5']) else 0
        utr3_length = len(str(row['UTR3'])) if pd.notna(row['UTR3']) else 0
        cds_length = len(str(row['CDS'])) if pd.notna(row['CDS']) else 0
        gene_len = len(str(row['Gene'])) if pd.notna(row['Gene']) else 0

        # length of each region
        region_length = {'upstream': upstream, 'gene': gene_len, 'UTR5': utr5_length, 'CDS': cds_length, 'UTR3': utr3_length, 'downstream': downstream}

        # calculate the start position of each region 
        region_start = {}
        cumulative = 0
        for region, length in region_length.items():
            region_start[region] = cumulative
            cumulative += length

        # Create a DataFrame for the gene with the calculated positions
        def compute_positions(row):
            seq_type = row['Sequence_Type']
            start = row['Start']
            if seq_type == 'Genomic': # Genomic sequence
                if start < upstream: # Upstream region
                    return start - upstream
                elif start < upstream + gene_len: # Gene region
                    return start - upstream
                else: # Downstream region
                    return start - upstream
            elif seq_type == 'Transcript':
                if start < utr5_length: # UTR5 region
                    return region_start['utr5'] + start
                elif start < utr5_length + cds_length: # CDS region
                    return region_start['cds'] + (start - utr5_length)
                else: # UTR3 region
                    return region_start['utr3'] + (start - utr5_length - cds_length)
            return
        
        # Create a DataFrame for the gene with the calculated positions
        gene_df['Position'] = gene_df.apply(compute_positions, axis=1)

        # Create the plot
        plt.figure(figsize=(12, 6))
        for _, row in gene_df.iterrows():
            plt.vlines(x=row['Position'], ymin=0, ymax=row['line_height'], color=row['tf_color'], linewidth=2) # vertical line for each TFBS
            plt.scatter(x=row['Position'], y=row['line_height'], color=row['tf_color'], s=50, edgecolor='black', zorder=3) # lollipop head for each TFBS
        
        # Set the x-ticks to the midpoints of each region
        region_labels = ['Upstream', 'Gene', 'Downstream', 'UTR5', 'CDS', 'UTR3']
        region_midpoints = {region : region_start[region] + region_length[region] / 2 for region in region_labels if region in region_start}

        # Set the x-ticks and labels
        plt.xticks(list(region_midpoints.values()), list(region_midpoints.keys()), fontsize=10)
        plt.xlabel('Gene Region (Genomic + Transcript)', fontsize=10)
        plt.yticks([])
        plt.ylim(0, max(heights) + 0.2)
        plt.title(f'Lollipop Plot for {gene_name}', fontsize=14)
        plt.grid(axis='x', linestyle=':', linewidth=0.5)
        
        # Add a legend for transcription factors
        legend_patches = [mpatches.Patch(color=color, label=tf) for tf, color in tf_color.items()]
        plt.legend(handles=legend_patches,title='Transcription Factors', loc='upper right', fontsize=10)
        plt.tight_layout()

        # Save the plot
        file_name = f"{species}_{gene_id}_{gene_name}.png"
        out_path = os.path.join(output_dir, file_name)
        plt.savefig(out_path, dpi = 300 ,bbox_inches='tight')
        plt.close()
        print(f"Plot saved to {out_path}")

# Main function to handle command line arguments and execute the script
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate lollipop plots from TFBS scanning results.')
    # Add command line arguments for input files, output file, and output directory
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='Input CSV files containing TFBS scanning results.')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for the merged results.')
    parser.add_argument('-d', '--output_dir', required=True, help='Output directory to save the lollipop plots.')
    parser.add_argument('--upstream', type=int, default=500, help='Length of upstream region (default: 500).')
    parser.add_argument('--downstream', type=int, default=500, help='Length of downstream region (default: 500).')

    # Parse the arguments
    args = parser.parse_args()

    # Merge the input CSV files and generate lollipop plots
    merged_df = merge_csv_files(args.input_files, 'merged_results.csv')
    generate_lollipop_plot(merged_df, args.output_dir, args.upstream, args.downstream)

# This script can be run from the command line with the required arguments.
if __name__ == '__main__':
    main()