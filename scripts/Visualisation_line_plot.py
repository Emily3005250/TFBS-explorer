
#### Visualisation with matches ####

## Import the libraries
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import sys

csv.field_size_limit(sys.maxsize)
sns.set_theme(style="whitegrid")



## Set the Commnadline Interfacae work
def main():
    parser = argparse.ArgumentParser(description='Visualise the location and number of matches')

    parser.add_argument('--dna', '-d', required=True, help='file for DNAseq')
    parser.add_argument('--scan', '-s', required=True, help='file for scanned with motif')

    args = parser.parse_args()

    line_plot(args.dna, args.scan)

if __name__ == '__main__':
    main()