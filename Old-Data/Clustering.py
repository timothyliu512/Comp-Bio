import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.cluster import KMeans

# Function to generate k-mers
def build_kmers(sequence, k_size):
    kmers = []
    n_kmers = len(sequence) - k_size + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + k_size]
        kmers.append(kmer)
    
    return ' '.join(kmers)

# Load sequences
sequences = [str(record.seq) for record in SeqIO.parse("SRR24133807.fastq", "fastq")]

# Generate 4-mers
k_size = 4
sequences_kmers = [build_kmers(seq, k_size) for seq in sequences]

# Convert k-mers to matrix format
vectorizer = CountVectorizer(analyzer='word', ngram_range=(4, 4), max_features=5000)
seq_vectors = vectorizer.fit_transform(sequences_kmers)

# Perform k-means clustering
kmeans = KMeans(n_clusters=10, random_state=0)
kmeans.fit(seq_vectors)

# Get cluster assignments
clusters = kmeans.predict(seq_vectors)

# Print cluster assignments
for i, cluster in enumerate(clusters):
    print(f"{i} : Cluster {cluster}")
