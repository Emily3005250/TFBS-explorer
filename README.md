# TFBS-explorer

Explores the presence of transcription factor binding sites (TFBSs) in the promoter regions of human interferon-stimulated genes (ISGs). It identifies and compares motif occurrences of IRF3, IRF7, and IRF9 using JASPAR motif profiles and Biopython-based scanning.

## Clone repository
```
git clone https://github.com/Sanghee-L/TFBS-explorer.git
```

## Scripts

[TFBS_explorer](scripts/TFBS_explorer_human.ipynb)

This notebook performs the following steps:

1. Retrieves the promoter sequences from ensembl using REST API
2. Fetchs the Trasncription Factor Binding Site profile from JASPAR database by pyJASPAR
3. Scans for TFBS motifs (IRF3, IRF7, IRF9) in the promoter sequences
4. Output CSV files containing match positions
5. Visualizes match counts across thresholds through elbow plot

## Installation
To perfrom the analysis, install the required Python libraries:

```bash
pip install requests biopython pyjaspar pandas matplotlib
```

## Inputs

human_IFN_up-regulation.csv : List of genes that are upregulated by IFN (positive control)

human_no_IFN-regulation.csv : List of genes not regulated by IFN (negative control)

## Outputs
Promoter_sequence.fasta : 3000bp upstream DNA for each gene

motif_scan_positive_threshold_X.csv : scan results at PSSM thresholds from 1 to 10 : Hit number for motif

motif_scan_negative_threshold_X.csv : scan results at PSSM thresholds from 1 to 10 : Hit number for motif

Elbow plot : Visualise the relationship between pssm threshold and number of motif matches