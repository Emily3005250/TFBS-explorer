# import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# read input file
input_file = "input.csv"
input_df = pd.read_csv(input_file)

# set colour and length for each TF
unique_tfs = input_df['TF'].unique()

cmap = plt.get_cmap('Set1')
colors = [cmap(i) for i in range(len(unique_tfs))]

tf_colour = dict(zip(unique_tfs, colors))
heights = np.linspace(0.15, 0.15*len(unique_tfs), len(unique_tfs))
tf_height = dict(zip(unique_tfs, heights))

input_df['tf_colour'] = input_df['TF'].map(tf_colour)
input_df['line_height'] = input_df['TF'].map(tf_height)

# make a dataframe subset for a gene of interest
gene_name = 'RTP4'
df = input_df[input_df['Gene_Name'] == gene_name]

# make a plot
plt.figure(figsize=(6, 2))

for _, row in df.iterrows():
    plt.vlines(row['Start'], 0, row['line_height'], color=row['tf_colour'])

plt.scatter(df['Start'], df['line_height'], color=df['tf_colour'], s=60, zorder=3, alpha=0.8)

plt.yticks([])
plt.ylim(0, 1.5)

plt.xlabel('TFBS position relative to TSS')
plt.xlim(0, 500)
plt.xticks(
    ticks=[0, 100, 200, 300, 400, 500],
    labels=['-500', '-400', '-300', '-200', '-100', '+1']
)

legend_handles = [mpatches.Patch(color=color, label=tf) for tf, color in tf_colour.items()]
plt.legend(handles=legend_handles, title="Transcription factor", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(f'{gene_name}_tfbs_positions_lollipop_plot.svg', format='svg', bbox_inches='tight')
plt.show()