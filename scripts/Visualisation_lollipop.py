import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os

def merge_csv_files(input_files, output_file):
    dataframes = [pd.read_csv(file) for file in input_files]
    merged_df = pd.concat(dataframes, ignore_index=True)
    merged_df.to_csv(output_file, index=False)
    return merged_df

def generate_lollipop_plot(df, output_dir, upstream=500, downstream=500):
    os.makedirs(output_dir, exist_ok=True)

    unique_tfs = sorted(df['Transcription_Factor'].unique())
    cmap = plt.get_cmap('Set1')
    colors = [cmap(i % 9) for i in range(len(unique_tfs))]
    tf_color = dict(zip(unique_tfs, colors))
    tf_height = dict(zip(unique_tfs, np.linspace(0.2, 1.5, len(unique_tfs))))

    df['tf_color'] = df['Transcription_Factor'].map(tf_color)
    df['line_height'] = df['Transcription_Factor'].map(tf_height)

    for (species, gene_id, gene_name), gene_df in df.groupby(['Species', 'Gene_ID', 'Gene_Name']):
        utr5_length = len(str(gene_df.iloc[0]['5_UTR'])) if pd.notna(gene_df.iloc[0]['5_UTR']) else 0
        cds_length = len(str(gene_df.iloc[0]['CDS'])) if pd.notna(gene_df.iloc[0]['CDS']) else 0
        utr3_length = len(str(gene_df.iloc[0]['3_UTR'])) if pd.notna(gene_df.iloc[0]['3_UTR']) else 0
        gene_len = len(str(gene_df.iloc[0]['Gene_Sequence'])) if pd.notna(gene_df.iloc[0]['Gene_Sequence']) else 0

        region_lengths_genomic = {'Upstream': upstream, 'Gene': gene_len, 'Downstream': downstream}
        region_lengths_transcript = {'5_UTR': utr5_length, 'CDS': cds_length, '3_UTR': utr3_length}

        def compute_genomic_pos(row):
            start = row['Start']
            if start < upstream:
                return start - upstream  # Upstream → negative
            elif start < upstream + gene_len:
                return start - upstream  # Gene → starts at 0
            else:
                return start - upstream  # Downstream → positive

        def compute_transcript_pos(row):
            start = row['Start']
            if start < utr5_length:
                return start
            elif start < utr5_length + cds_length:
                return utr5_length + (start - utr5_length)
            else:
                return utr5_length + cds_length + (start - utr5_length - cds_length)

        # Genomic lollipop plot
        gdf = gene_df[gene_df['Sequence_Type'].str.lower() == 'genomic'].copy()
        if not gdf.empty:
            gdf['Position'] = gdf.apply(compute_genomic_pos, axis=1)
            fig, ax = plt.subplots(figsize=(10, 3))
            for _, row in gdf.iterrows():
                ax.vlines(x=row['Position'], ymin=0, ymax=row['line_height'], color=row['tf_color'])
                ax.scatter(row['Position'], row['line_height'], color=row['tf_color'], s=60, edgecolor='black', zorder=3)
            ax.axvline(x=0, color='black', linestyle='--', label='Gene Start')
            ax.set_title(f"Lollipop Plot: {gene_name} - Genomic")
            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_xlabel("Genomic Region (TSS = 0)")
            ax.set_xlim(min(gdf['Position']) - 50, max(gdf['Position']) + 50)
            ax.grid(axis='x', linestyle=':')
            patches = [mpatches.Patch(color=color, label=tf) for tf, color in tf_color.items()]
            ax.legend(handles=patches, title="TF")
            plt.tight_layout()
            out_path = os.path.join(output_dir, f"{gene_id}_{gene_name}_genomic.png")
            plt.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved plot: {out_path}")

        # Transcript lollipop plot #
        tdf = gene_df[gene_df['Sequence_Type'].str.lower() == 'transcript'].copy()
        if not tdf.empty:
            tdf['Position'] = tdf.apply(compute_transcript_pos, axis=1)
            fig, ax = plt.subplots(figsize=(10, 3))
            for _, row in tdf.iterrows():
                ax.vlines(x=row['Position'], ymin=0, ymax=row['line_height'], color=row['tf_color'])
                ax.scatter(row['Position'], row['line_height'], color=row['tf_color'], s=60, edgecolor='black', zorder=3)
            # Annotate the regions
            region_ends = {
                '5\'UTR': utr5_length,
                'CDS': utr5_length + cds_length,
                '3\'UTR': utr5_length + cds_length + utr3_length
            }
            for label, x in region_ends.items():
                ax.axvline(x=x, linestyle=':', color='gray')
                ax.text(x - 20, max(tf_height.values()) + 0.1, label, rotation=0, fontsize=9)

            ax.set_title(f"Lollipop Plot: {gene_name} - Transcript")
            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_xlabel("Transcript Region")
            ax.set_xlim(-20, tdf['Position'].max() + 50)
            ax.grid(axis='x', linestyle=':')
            patches = [mpatches.Patch(color=color, label=tf) for tf, color in tf_color.items()]
            ax.legend(handles=patches, title="TF")
            plt.tight_layout()
            out_path = os.path.join(output_dir, f"{gene_id}_{gene_name}_transcript.png")
            plt.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved plot: {out_path}")

def main():
    parser = argparse.ArgumentParser(description='Generate lollipop plots from TFBS scanning results.')
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='Input CSV files with TFBS scanning results.')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for merged data.')
    parser.add_argument('-d', '--output_dir', required=True, help='Directory to save plots.')
    parser.add_argument('--upstream', type=int, default=500, help='Upstream length')
    parser.add_argument('--downstream', type=int, default=500, help='Downstream length')

    args = parser.parse_args()
    merged_df = merge_csv_files(args.input_files, args.output_file)
    generate_lollipop_plot(merged_df, args.output_dir, upstream=args.upstream, downstream=args.downstream)

if __name__ == '__main__':
    main()