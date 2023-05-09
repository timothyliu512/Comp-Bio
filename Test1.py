import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

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

# k-means clustering
kmeans = KMeans(n_clusters=10, random_state=0)
kmeans.fit(seq_vectors)

# Get cluster assignments
clusters = kmeans.predict(seq_vectors)

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=0)
seq_vectors_2d = tsne.fit_transform(seq_vectors.toarray())

# Create a scatter plot
plt.figure(figsize=(8, 8))
scatter = plt.scatter(seq_vectors_2d[:, 0], seq_vectors_2d[:, 1], c=clusters, cmap='viridis')

# Add a color bar
plt.colorbar(scatter)

plt.show()
