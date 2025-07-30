
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os

# Read csv files from scanning and merge them
def merge_csv_files(input_files, output_file):
    dataframes = [pd.read_csv(file) for file in input_files] # Read each file into a DataFrame
    merged_df = pd.concat(dataframes, ignore_index=True) # Combine all dataframes into one
    merged_df.to_csv(output_file, index=False)
    return merged_df


def generate_lollipop_plot(df, output_dir, upstream=500, downstream=500):
    os.makedirs(output_dir, exist_ok=True)

    unique_tfs = sorted(df['Transcription_Factor'].dropna().unique())
    # Assign colors to each transcription factor
    # Using a colormap to ensure distinct colors for each TF
    cmap = plt.get_cmap('Set1')
    colors = [cmap(i % 9) for i in range(len(unique_tfs))]
    tf_color = dict(zip(unique_tfs, colors))

    legend_handles = [mpatches.Patch(color=tf_color[tf], label = tf) for tf in unique_tfs]

    for (species, gene_id, gene_name), gene_df in df.groupby(['Species', 'Gene_ID', 'Gene_Name']):
        utr5_length = len(str(gene_df.iloc[0]['5_UTR'])) if pd.notna(gene_df.iloc[0]['5_UTR']) else 0
        cds_length = len(str(gene_df.iloc[0]['CDS'])) if pd.notna(gene_df.iloc[0]['CDS']) else 0
        utr3_length = len(str(gene_df.iloc[0]['3_UTR'])) if pd.notna(gene_df.iloc[0]['3_UTR']) else 0
        gene_len = len(str(gene_df.iloc[0]['Gene_Sequence'])) if pd.notna(gene_df.iloc[0]['Gene_Sequence']) else 0
        
        #### Genomic lollipop plot ####
        # Filter for genomic sequence type
        # This assumes the DataFrame has a 'Sequence_Type' column indicating 'Genomic'
        gdf = gene_df[gene_df['Sequence_Type'] == 'Genomic'].copy()
        if not gdf.empty:
            
            gdf['abs_position'] = gdf['Start'] - upstream
            
            ## Divide the genomic positions into bins
            # Bin upstream region
            # This assumes the upstream region is divided into 5 equal parts, adjust as necessary
            upstream_bins = list(range(-upstream, 0, 100)) # e.g., -500, -400, ..., -100

            # Bin gene region with 20 bins
            # This assumes the gene is divided into 20 equal parts, adjust as necessary
            gene_bins = 20
            gene_edges = list(np.linspace(0, gene_len, gene_bins + 1))

            # Bin downstream region
            # This assumes the downstream region is divided into 5 equal parts, adjust as necessary
            downstream_bins = list(range(gene_len + 100, gene_len + downstream + 1, 100))

            # Combine all bins
            all_bins = upstream_bins + gene_edges + downstream_bins

            gdf['bin'] = pd.cut(gdf['abs_position'], bins= all_bins, include_lowest=True)
            gdf['bin_mid'] = gdf['bin'].apply(lambda x: x.mid) 
            # Assign colors to transcription factors
            gdf['tf_color'] = gdf['Transcription_Factor'].map(tf_color)

            # Sort the DataFrame by bin_mid and Transcription_Factor
            # This ensures that the lollipop plot is organized by bins and transcription factors
            gdf.sort_values(by=['bin_mid', 'Transcription_Factor'], inplace=True) # Sort by bin_mid and Transcription_Factor
            gdf['stack_level'] = gdf.groupby(['bin_mid', 'Transcription_Factor']).cumcount() # This creates a stacking effect for overlapping TFBS
            gdf['line_height'] = gdf['stack_level'] * 0.15 + 0.15 # Adjust the height of the lines for visibility
            
            # Create the plot
            plt.figure(figsize=(10, 3))
            
            # Plot each transcription factor binding site as a vertical line
            # This creates the lollipop effect with lines and points
            for _, row in gdf.iterrows():
                plt.vlines(x=row['bin_mid'], ymin=0, ymax=row['line_height'], color=row['tf_color'])
                plt.scatter(row['bin_mid'], row['line_height'], color=row['tf_color'],s = 60, zorder = 3)
            
            plt.title(f"{species} - {gene_name} ({gene_id}) - Genomic Lollipop Plot")
            plt.vlines(x=0, ymin=0, ymax= 0.3, color='black', linestyle='--') # Vertical line at TSS
            plt.xlabel('TFBS position relative to TSS') # X-axis label
            # Tick positions and labels for upstream (every 100 bp)
            upstream_ticks = list(range(-upstream, 0, 100)) # e.g., -500, -400, ..., -100
            upstream_labels = [f"{x}" for x in upstream_ticks]  # e.g., -500, -400, ..., -100

            # TSS tick
            gene_ticks = [0] # TSS position
            gene_labels = ['+1'] # e.g., +1 for TSS

            # Downstream ticks (every 100 bp after gene end)
            downstream_ticks = list(range(gene_len + 100, gene_len + downstream + 1, 100)) # e.g., 100, 200, ..., 500
            downstream_labels = [f"+{x - gene_len}" for x in downstream_ticks]  # e.g., +100, +200, ...

            # Combine all ticks and labels
            tick_locs = upstream_ticks + gene_ticks + downstream_ticks
            tick_labels = upstream_labels + gene_labels + downstream_labels
            # Set the ticks and labels on the x-axis
            plt.xticks(tick_locs, labels=tick_labels, rotation=45)
            plt.xlim(-upstream - 50, gene_len + downstream + 50) # Set x-axis limits
            # Set y-axis labels and limits
            plt.ylabel('TFBS Count')
            plt.ylim(0, gdf['line_height'].max() + 0.2) # Set y-axis limits
            
            # Add a legend for transcription factors
            plt.legend(handles=legend_handles, loc='upper right', title='Transcription Factors')
            
            # Save the plot
            plt.tight_layout()
            output_file = os.path.join(output_dir, f"{species}_{gene_id}_genomic_lollipop.png")
            plt.savefig(output_file, dpi=300)
            plt.close()

        print(f"Genomic lollipop plot saved to {output_file}")
        
        #### Transcript lollipop plot ####
        # Filter for transcript sequence type
        # This assumes the DataFrame has a 'Sequence_Type' column indicating 'Transcript'
        tdf = gene_df[gene_df['Sequence_Type'].str.strip().str.lower() == 'transcript'].copy()
        if not tdf.empty:
            tdf['abs_position'] = tdf['Start'] - utr5_length
            
            # Divide the transcript positions into bins
            # Bin 5' UTR region & 3' UTR region
            utr5_bins = list(range(0, utr5_length + 100, 100))
            utr3_bins = list(range(utr5_length + cds_length, utr5_length + cds_length + utr3_length + 100, 100))
            # Bin CDS region
            cds_bin = 20
            cds_bins = list(np.linspace(utr5_length, utr5_length + cds_length, cds_bin + 1))

            # Combine all bins
            all_bins = utr5_bins + cds_bins + utr3_bins

            # Create bin columns and falls hits into bins
            # Calculate the mid-point of each bin for plotting
            tdf['bin'] = pd.cut(tdf['abs_position'], bins=all_bins, include_lowest=True)
            tdf['bin_mid'] = tdf['bin'].apply(lambda x: x.mid)
            # Assign colors to transcription factors
            tdf['tf_color'] = tdf['Transcription_Factor'].map(tf_color)

            # Sort the DataFrame by bin_mid and Transcription_Factor
            tdf.sort_values(by=['bin_mid', 'Transcription_Factor'], inplace=True)
            tdf['stack_level'] = tdf.groupby(['bin_mid', 'Transcription_Factor']).cumcount()
            tdf['line_height'] = tdf['stack_level'] * 0.15 + 0.15

            # Create the plot
            plt.figure(figsize=(10, 3))

            # Plot each transcription factor binding site as a vertical line
            for _, row in tdf.iterrows():
                plt.vlines(x=row['bin_mid'], ymin=0, ymax=row['line_height'], color=row['tf_color'])
                plt.scatter(row['bin_mid'], row['line_height'], color=row['tf_color'], s=60, zorder=3)
            
            plt.title(f"{species} - {gene_name} ({gene_id}) - Transcript Lollipop Plot")
            plt.vlines(x= 0, ymin= 0, ymax= 0.3, color='black', linestyle='--') # Vertical line at CDS start
            plt.xlabel('TFBS position relative to CDS Start') # X-axis label
            
            # Tick positions and labels for 5' UTR (every 100 bp)
            utr5_ticks = list(range(0, utr5_length + 100, 100)) # e.g., 0, 100, ..., 400
            utr5_labels = [f"{x}" for x in utr5_ticks]  # e.g., 0, 100, ..., 400
            # CDS tick
            cds_ticks = [utr5_length] # CDS start position
            cds_labels = ['CDS Start'] # e.g., CDS Start
            # Downstream ticks (every 100 bp after CDS end)
            utr3_ticks = list(range(utr5_length + cds_length + 100, utr5_length + cds_length + utr3_length + 1, 100)) # e.g., 100, 200, ..., 500
            utr3_labels = [f"{x - (utr5_length + cds_length)}" for x in utr3_ticks]  # e.g., +100, +200, ...
            # Combine all ticks and labels
            tick_locs = utr5_ticks + cds_ticks + utr3_ticks
            tick_labels = utr5_labels + cds_labels + utr3_labels
            
            # Set the ticks and labels on the x-axis
            plt.xticks(tick_locs, labels=tick_labels, rotation=45)
            plt.xlim(-50, utr5_length + cds_length + utr3_length + 50) # Set x-axis limits
            # Set y-axis labels and limits
            plt.ylabel('TFBS Count')
            plt.ylim(0, tdf['line_height'].max() + 0.2) # Set y-axis limits
            
            # Add a legend for transcription factors
            plt.legend(handles=legend_handles, loc='upper right', title='Transcription Factors')
            
            # Save the plot
            plt.tight_layout()
            output_file = os.path.join(output_dir, f"{species}_{gene_id}_transcript_lollipop.png")
            plt.savefig(output_file, dpi=300)
            plt.close()
        print(f"Transcript lollipop plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate lollipop plots of TFBSs from motif scanning results.")

    parser.add_argument( "-i", "--input_files", nargs='+', required=True, help="List of input CSV files (scanning results) to merge and plot")
    parser.add_argument("-o", "--output_file", required=True, help="Filename for the merged CSV output (e.g. merged_results.csv)")
    parser.add_argument("-d", "--output_dir", required=True, help="Directory where plots will be saved")
    parser.add_argument("--upstream", type=int, default=500, help="Upstream region length in bp (default: 500)")
    parser.add_argument("--downstream", type=int, default=500, help="Downstream region length in bp (default: 500)")

    args = parser.parse_args()

    # Merge input files
    merged_df = merge_csv_files(args.input_files, args.output_file)

    # Generate lollipop plots
    generate_lollipop_plot(merged_df, output_dir=args.output_dir, upstream=args.upstream, downstream=args.downstream)

if __name__ == "__main__":
    main()