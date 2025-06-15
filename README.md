# TFBS-explorer

Exploration of Transcription Factor Binding Sites in promoter regions of IFN-stimulated genes

## Scripts

[TFBS_explorer](scripts/TFBS_explorer_human.ipynb)

Work-Flow

1. Retrieve the promoter sequence from ensembl using rest api
2. Fetch the Trasncription Factor Binding Site sequence from JASPAR database using pyJASPAR
3. Scan TFBS in the promoter sequences
4. Output CSV with match positions
5. Visualizes match counts (elbow plot)

## Installation
To perfrom the script you need to install python libraries.

```bash
pip install requests biopython pyjaspar pandas matplotlib
```

## Inputs
Genes that are regulated by IFN (Positive Control) and not regulated by IFN (Negative Control)

human_IFN_up-regulation.csv 

human_no_IFN-regulation.csv

## Outputs
Promoter_sequence.fasta : 3000bp upstream DNA for each gene

motif_scan_positive_threshold_X.csv : scan results at PSSM thresholds from 1 to 10 : Hit number for motif

motif_scan_negative_threshold_X.csv : scan results at PSSM thresholds from 1 to 10 : Hit number for motif

Elbow plots of Transscription Factor Binding Site : To check the sensitivity