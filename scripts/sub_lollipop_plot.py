
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os
import numpy as np
import math

def generate_subplots(match_file, seq_file, output_dir, gene_name, chunk_size):
    # Load data
    match_df = pd.read_csv(match_file)
    seq_df = pd.read_csv(seq_file)

    # Filter to the selected gene and genomic only
    gdf = match_df[(match_df['Sequence_Type'] == 'Genomic') & (match_df['Gene_Name'] == gene_name)].copy()
    seq_row = seq_df[seq_df['gene_name'] == gene_name]
    if gdf.empty or seq_row.empty:
        print(f"No data found for {gene_name}")
    
    else:
        # Assign colours and heights
        unique_tfs = gdf['Transcription_Factor'].unique()
        cmap = plt.get_cmap('Set1')
        colors = [cmap(i) for i in range(len(unique_tfs))]
        tf_colors = dict(zip(unique_tfs, colors))
        heights = np.linspace(0.15, 0.15*len(unique_tfs), len(unique_tfs))
        tf_height = dict(zip(unique_tfs, heights))

        gdf['Color'] = gdf['Transcription_Factor'].map(tf_colors)
        gdf['line_height'] = gdf['Transcription_Factor'].map(tf_height)

        if 'Strand' in gdf.columns:
            gdf['Strand_Symbol'] = gdf['Strand']
        else:
            gdf['Strand_Symbol'] = ''

        # Sequence length
        upstream_len = len(str(seq_row.iloc[0]['Upstream']))
        gene_len = len(str(seq_row.iloc[0]['Gene_seq']))
        downstream_len = len(str(seq_row.iloc[0]['Downstream']))
        total_len = upstream_len + gene_len + downstream_len

        # Adjust match positions relative to TSS
        gdf['Rel_Start'] = gdf['Start'] - upstream_len

        domain_min = -upstream_len
        domain_max = gene_len + downstream_len
        sum = domain_max - domain_min

        # Divide into chunks
        n_chunks = (total_len // chunk_size) + 1
        os.makedirs(output_dir, exist_ok=True)

        # Loop through chunks
        # Filter for the current chunk and make the dataframe
        for i in range(n_chunks):
            start = domain_min + i * chunk_size
            end = min(domain_min + (i+1) * chunk_size, domain_max)
            df_chunk = gdf[(gdf['Rel_Start'] >= start) & (gdf['Rel_Start'] < end)].copy()

            # Make a figure for the chunk
            plt.figure(figsize=(9, 4))

            # Plotting the lollipop plot for the chunk
            for _, row in df_chunk.iterrows():
                plt.vlines(row['Rel_Start'], 0, row['line_height'], color=row['Color'])
            
            plt.scatter(df_chunk['Rel_Start'], df_chunk['line_height'],color=df_chunk['Color'], s=60, zorder=3, alpha=0.8)

            # Add strand symbols
            for _, row in df_chunk.iterrows():
                x = float(row['Rel_Start'])
                y = float(row['line_height'])

                if row['Strand_Symbol']:
                    plt.text(x, y, str(row['Strand_Symbol']), ha='center', va = 'center', color = 'black', fontsize=8, fontweight='bold', zorder = 4)

                plt.annotate(str(int(row['Rel_Start'])), xy=(x, y), xytext=(5, 0), textcoords='offset points', ha='left', va='center', fontsize=8, color='black', zorder=4)

            # Labels and title
            species = seq_row.iloc[0]['species']
            gene_id = seq_row.iloc[0]['ensembl_id']
            plt.title(f'{species} - {gene_name} ({gene_id})\nChunk {i+1}: {start} to {end} bp')
            plt.yticks([])
            plt.ylim(0, 1.5)
            plt.xlabel('Position (bp)')
            plt.xlim(start, end)

            # Legend
            legend_handles = [mpatches.Patch(color=color, label=tf) for tf, color in tf_colors.items()]
            plt.legend(handles=legend_handles, title='TFs', fontsize='small', loc='upper left')

            # Aesthetics
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)

            # Save
            output_path = os.path.join(output_dir, f'{gene_name}_chunk_{i+1}.svg')
            plt.tight_layout()
            plt.savefig(output_path, format='svg', bbox_inches='tight')
            plt.close()
            print(f"Saved {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate sub-lollipop plots by chunking sequence length.")
    parser.add_argument("-m", "--match_file", required=True, help="CSV file from Pipeline 2 TFBS scanning")
    parser.add_argument("-s", "--seq_file", required=True, help="CSV file with sequences from Ensembl (Pipeline 1)")
    parser.add_argument("-d", "--output_dir", required=True, help="Directory where plots will be saved")
    parser.add_argument("-g", "--gene", required=True, help="Gene name to plot")
    parser.add_argument("-c", "--chunk_size", type=int, default=3000, help="Chunk size in bp (default: 3000)")
    args = parser.parse_args()

    generate_subplots(args.match_file, args.seq_file, args.output_dir, args.gene, args.chunk_size)

if __name__ == "__main__":
    main()