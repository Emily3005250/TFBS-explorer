# TFBS-explorer

Exploration of Transcription Factor Binding Sites in promoter regions of IFN-stimulated genes

## Scripts

[TFBS_explorer](scripts/TFBS_explorer_human.ipynb)

1. Retrieve the promoter sequence from ensembl using rest api
2. Fetch the Trasncription Factor Binding Site sequence from JASPAR database using pyJASPAR
3. Scan TFBS in the promoter sequences
4. Output CSV with match positions
5. Visualizes match counts (elbow plot)

## Installation
To perfrom the script you need to install python libraries.

<pre>
'''bash
pip install request biopython pyjaspar pandas matplotlib
'''
<pre>

## Outputs
