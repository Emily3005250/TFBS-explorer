# TFBS-explorer

Explores the presence of transcription factor binding sites (TFBSs) in the promoter regions of human interferon-stimulated genes (ISGs). It identifies and compares motif occurrences of IRF3, IRF7, and IRF9 using JASPAR motif profiles and Biopython-based scanning.

## Clone repository
```
git clone https://github.com/Sanghee-L/TFBS-explorer.git
```

## Installation
To perfrom the analysis, install the required Python libraries:

```bash
pip install requests biopython pyjaspar pandas matplotlib numpy 
```

## Pipeline Module

### Module 1: Sequence Retrieval
'Fetch_DNA_Sequence.py'
- Fetches 3,000 bp upstream, full gene sequence, and 3,000 bp downstream sequence using the Ensemb REST API

### Module 2: TFBS Scanning
- 'Exact_match.py' - Identifies perfect matches to core motifs
- 'Regex_match.py' - Uses flexiblel pattern matching
- 'JASPAR_profile.py' - Scans using PSSMs from the JASPAR 2024 database
- 'Elbow_plot.py' - Helps select the optimal threshold for PSSM matching

### Module 3: Visualisation
- 'Lollipop_plot.py' - Creates lollipop plots of TFBS locations per gene
- 'sub_lollipop_plot.py' - Generates zoomed-in, windowed TFBS plots for detailed inspection



## Inputs

human_IFN_up-regulation.csv : List of 100 genes that are upregulated by IFN (positive control)

human_no_IFN-regulation.csv : List of 100 genes not regulated by IFN (negative control)

## Input file structure for module 1 (Fetch_DNA_Sequence.py) - positive control

| Species       | Ensembl ID       | Gene       | Expression  | Orthologous Cluster ID |
|----------------|----------------|----------------|----------|-------------|
| Homo sapiens   | ENSG00000157601   | MX1   | up_regulated | HS6198 |
| Homo sapiens   | ENSG00000135114   | OASL   | up_regulated |HS1036 | 

Visit the [Orthologous Clusters of Interferon-Stimulated Genes (ISGs)](https://isg.data.cvr.ac.uk/) database for more information.

## Outputs
Module 1 (CSV file)
- tabular format csv file which contains sequence data (upstream, gene, downstream, transcript)

| Species | gene_Name | ensembl_ID | Upstream | Gene_seq | Downstream | Transcript_ID | 5_UTR | CDS | 3_UTR |
|---------|-----------|------------|----------|----------|------------|---------------|-------|-----|-------|
| Homo sapience | MFSD10 | ENSG00000109736 | 3kb seq | gene seq | 3kb seq | ENST00000355443 | 5' UTR seq | CDS seq | 3' UTR seq |


Module 2 (CSV file)
- tabular format csv file which contains matches coordinates and metadata

| Species | Gene_Name | ID | Transcription_Factor | Strand | Start | End | Sequence_Type |
|---------|-----------|----|----------------------|--------|-------|-----|---------------|
| Homo sapience | COL23A1 | ENSG00000050767 | IRF3 | + | 194685 | 194695 | Genomic |
| Homo sapience | COL23A1 | ENSG00000050767 | IRF9 | + | 194685 | 194695 | Genomic |
| Homo sapience | STK4 | NSG00000101109 | IRF7 | - | 61826 | 61836 | Genomic |


Module 3 (SVG)

- Lollipop plots to show the TFBS distribution
![Example Lollipop Plot](/TFBS-explorer/outputs/lollipop_plot/TRIM5_tfbs_positions_lollipop_plot.svg)


- Zoomed (windowed) plots for detailed TFBS distribution and features
![Zoomed TFBS plot](/TFBS-explorer/outputs/sub_lollipop_plot/TRIM5_chunk_6.svg)



## Example Command-Line Usage

#### Module 1 (Fetch_DNA_Sequence.py)
```bash
python3 Fetch_DNA_Sequence.py -i human_IFN_up-regulation.csv -o human_IFN_up_regulation_seq.csv -u 3000 -d 3000
```

#### Module 2 (Scanning)

- Exact_match.py
```bash
python3 Exact_match.py -i ../output/human_IFN_up_regulation_seq.csv -m motif.seq -o pos_exact.csv
```

- Regex_match.py
```bash
python3 Regex_match.py -i ../output/human_IFN_up_regulation_seq.csv -m motif.seq -o pos_regex.csv
```

- JASPAR_profile.py
```bash
python3 JASPAR_profile.py -i ../output/human_IFN_up_regulation_seq.csv -t IRF3 IRF7 IRF9 -s 14.0 -o pos_jaspar.csv
```

- Elbow_plot.py
```bash
python3 Elbow_plot.py -t IRF3 -d ./outputs --min_thresh 1.0 --max_thresh 15.0 --pos_exact 20 --pos_regex 150 --neg_exact 5 --neg_regex 90
```

#### Module 3 (Visualisation)
- Lollipop_plot.py
```bash
python3 Lollipop_plot.py -m pos_jaspar.csv -s human_IFN_up_regulation_seq.csv -d lollipop_plots --upstream 3000 --downstream 3000
```

- sub_lollipop_plot.py
```bash
python3 sub_lollipop_plot.py -m pos_jaspar.csv -s human_IFN_up_regulation_seq.csv -d sub_lollipop_plot -g TRIM5 -c 1000
```