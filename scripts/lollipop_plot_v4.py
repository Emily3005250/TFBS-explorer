import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# === Input files ===
match_file = '/Users/sang/bioinformatics_msc_project_2025/TFBS-explorer-19/scripts/pos_exact.csv'
seq_file = '/Users/sang/bioinformatics_msc_project_2025/TFBS-explorer-19/scripts/human_IFN_up_regulation_extract.csv'

# === Load data ===
match_df = pd.read_csv(match_file)
seq_df = pd.read_csv(seq_file)

# === Filter to genomic only ===
match_df = match_df[match_df['Sequence_Type'] == 'Genomic']

# === Assign colors ===
unique_tfs = match_df['Transcription_Factor'].unique()
cmap = plt.get_cmap('Set1')
colors = [cmap(i) for i in range(len(unique_tfs))]
tf_colors = dict(zip(unique_tfs, colors))
heights = np.linspace(0.15, 0.15*len(unique_tfs), len(unique_tfs))
tf_height = dict(zip(unique_tfs, heights))

match_df['Color'] = match_df['Transcription_Factor'].map(tf_colors)
match_df['line_height'] = match_df['Transcription_Factor'].map(tf_height)

# === Loop through genes ===
for gene_name in match_df['Gene_Name'].unique():
    df = match_df[match_df['Gene_Name'] == gene_name].copy()
    seq_row = seq_df[seq_df['gene_name'] == gene_name]

    if seq_row.empty:
        continue

    # === Get sequence lengths ===
    upstream_len = len(str(seq_row.iloc[0]['Upstream']))
    gene_len = len(str(seq_row.iloc[0]['Gene_seq']))
    downstream_len = len(str(seq_row.iloc[0]['Downstream']))
    total_len = upstream_len + gene_len + downstream_len
    
    # === Relative position to TSS ===
    df['Rel_Start'] = df['Start'] - upstream_len

    plt.figure(figsize=(9, 4))

    for -, row in df.iterrows():
        plt.vlines(row['Rel_Start'], 0, row['line_height'], color=row['Color'], linewidth=2)


