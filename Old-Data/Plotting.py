import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.cluster import KMeans
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

# Read sequences from FASTQ file
rawData = []
for record in SeqIO.parse('SRR24133807.fastq', 'fastq'):
    rawData.append(str(record.seq))

# Create k-mer count matrix from list of sequences
def kmer_count_matrix(sequences, k):
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(k, k))
    return vectorizer.fit_transform(sequences)

# Perform clustering on sequences using k-mer count matrix
def cluster_sequences(sequences, k, n_clusters):
    kmer_matrix = kmer_count_matrix(sequences, k)
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=0)
    kmeans.fit(kmer_matrix)
    return kmeans.labels_

# Parameters for clustering and dimensionality reduction
k = 4
n_clusters = 5

# Performing clustering on sequences
labels = cluster_sequences(rawData, k, n_clusters)

# Perform t-SNE dimensionality reduction on k-mer count matrix
kmer_matrix = kmer_count_matrix(rawData, k)
tsne = TSNE(n_components=2, random_state=0)
reduced_data = tsne.fit_transform(kmer_matrix.toarray())

# Plot t-SNE results with cluster labels
plt.figure(figsize=(10, 8))
for i in range(n_clusters):
    cluster_data = reduced_data[labels == i]
    plt.scatter(cluster_data[:, 0], cluster_data[:, 1], label=f'Cluster {i}')
plt.legend()
plt.title('t-SNE visualization of clustered sequences')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.show()