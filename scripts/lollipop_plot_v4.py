
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os

def generate_lollipop_plots(match_file, seq_file, output_dir, upstream, downstream):
    # Load data
    match_df = pd.read_csv(match_file)
    seq_df = pd.read_csv(seq_file)

    # Filter to genomic only
    match_df = match_df[match_df['Sequence_Type'] == 'Genomic']

    # Assign colors and heights for transcription factors
    unique_tfs = match_df['Transcription_Factor'].unique()
    cmap = plt.get_cmap('Set1')
    colors = [cmap(i) for i in range(len(unique_tfs))]
    tf_colors = dict(zip(unique_tfs, colors))
    heights = np.linspace(0.15, 0.15*len(unique_tfs), len(unique_tfs))
    tf_height = dict(zip(unique_tfs, heights))

    match_df['Color'] = match_df['Transcription_Factor'].map(tf_colors)
    match_df['line_height'] = match_df['Transcription_Factor'].map(tf_height)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through genes
    for gene_name in match_df['Gene_Name'].unique():
        df = match_df[match_df['Gene_Name'] == gene_name].copy()
        seq_row = seq_df[seq_df['gene_name'] == gene_name]

        if seq_row.empty:
            continue

        # Get sequence lengths
        upstream_len = len(str(seq_row.iloc[0]['Upstream']))
        gene_len = len(str(seq_row.iloc[0]['Gene_seq']))
        downstream_len = len(str(seq_row.iloc[0]['Downstream']))
        total_len = upstream_len + gene_len + downstream_len
        
        # Relative position to TSS
        df['Rel_Start'] = df['Start'] - upstream_len

        plt.figure(figsize=(9, 4))

        # Plotting the lollipop plot
        for _, row in df.iterrows():
            plt.vlines(row['Rel_Start'], 0, row['line_height'], color=row['Color'])

        plt.scatter(df['Rel_Start'], df['line_height'], color=df['Color'], s=60, zorder=3, alpha=0.8)

        # Labels and title
        species = seq_row.iloc[0]['species']
        gene_id = seq_row.iloc[0]['ensembl_id']
        plt.title(f'Lollipop Plot {species} - {gene_name} ({gene_id})')
        plt.yticks([])
        plt.ylim(0, 1.5)
        plt.xlabel('Relative Position to TSS (bp)')
        plt.xlim(-upstream_len, gene_len + downstream_len)
        tick_locs = [-upstream_len, 0, upstream_len + gene_len + downstream_len]
        tick_labels = [f'-{upstream_len}','+1', str(total_len - upstream_len)]
        plt.xticks(tick_locs, tick_labels, rotation=45)
        #plt.text(0, 0.02, '+1', transform=plt.gca().get_xaxis_transform(), ha='center', va='bottom', fontsize=10)

        # Vertical lines for gene start and gene end
        plt.axvline(x=0, color='black', linestyle='--', linewidth=1, ymin=0, ymax=0.5/1.5)
        plt.axvline(x=gene_len, color='black', linestyle='--', linewidth=1, ymin=0, ymax=0.5/1.5)

        # Legend
        legend_handles = [mpatches.Patch(color=color, label=tf) for tf, color in tf_colors.items()]
        plt.legend(handles=legend_handles, title='Transcription Factors', fontsize='small', loc='upper left')

        # Aesthetics
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

        # Save and close
        output_path = os.path.join(output_dir, f'{gene_name}_tfbs_positions_lollipop_plot.svg')
        plt.tight_layout()
        plt.savefig(output_path, format='svg', bbox_inches='tight')
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate lollipop plots of TFBSs from motif scanning results.")
    parser.add_argument("-m", "--match_file", required=True, help="CSV file from Pipeline 2 TFBS scanning (e.g. pos_jaspar.csv)")
    parser.add_argument("-s", "--seq_file", required=True, help="CSV file with sequences from Ensembl (Pipeline 1)")
    parser.add_argument("-d", "--output_dir", required=True, help="Directory where plots will be saved")
    parser.add_argument("--upstream", type=int, default=3000, help="Upstream region length in bp (default: 3000)")
    parser.add_argument("--downstream", type=int, default=3000, help="Downstream region length in bp (default: 3000)")
    args = parser.parse_args()

    generate_lollipop_plots(match_file=args.match_file, 
                            seq_file=args.seq_file, 
                            output_dir=args.output_dir, 
                            upstream=args.upstream, 
                            downstream=args.downstream)

if __name__ == "__main__":
    main()