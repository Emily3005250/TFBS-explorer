
import pandas as pd

header = ['Transcription_Factor', 'exact_forward', 'exact_reverse','regex_forward', 'regex_reverse']

sequence = [
            ['IRF3', 'GAAACCGAAA','TTTCGGTTTC', 'GAAA[ATGC]{2}GAAA', 'TTTC[ATGC]{2}TTTC'], 
            ['IRF7', 'GAAAGCGAAA','TTTCGCTTTC' , 'GAAA[ATGC]{2}GAAA', 'TTTC[ATGC]{2}TTTC'], 
            ['IRF9', 'GAAACCGAAA', 'TTTCGGTTTC', 'GAAA[ATGC]{2}GAAA', 'TTTC[ATGC]{2}TTTC']
            ]

df = pd.DataFrame(sequence, columns=header)

df.to_csv("motif_seq.csv", index=False)