import numpy as np
import pandas as pd
#from Bio import SeqIO



# df = pd.read_csv('SRR24133807.fastq', delimiter='\t')

df = pd.DataFrame(pd.read_csv('SRR24133807.fastq', sep='\t', header=None).values.reshape(-1, 4), columns=['read_id', 'seq', '+', 'qual'])
print(df.qual.to_string(index=False))

'''
df = pd.DataFrame

for seqrecord in SeqIO.parse("SRR24133807.fastq", "fastq"):
    df = pd.DataFrame(seqrecord.seq)
  ''' 

    # print (seqrecord.seq)
# print(df)