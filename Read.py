import numpy as np
import pandas as pd
from Bio import SeqIO



# df = pd.read_csv('SRR24133807.fastq', delimiter='\t')

df = pd.DataFrame

for seqrecord in SeqIO.parse("SRR24133807.fastq", "fastq"):
    df = pd.DataFrame(seqrecord.seq)
    # print (seqrecord.seq)
# print(df)