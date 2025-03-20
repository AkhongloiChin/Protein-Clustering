import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from itertools import product
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import base64
from io import BytesIO
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.datasets import make_blobs

import numpy as np
from itertools import product

def get_kmer_counts(sequence, k=2):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    possible_kmers = [''.join(p) for p in product(amino_acids, repeat=k)]
    kmer_counts = {kmer: 0 for kmer in possible_kmers}
    total = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
            total += 1
    if total == 0:
        return [0.0] * len(possible_kmers)
    return [kmer_counts[kmer] / total for kmer in possible_kmers]

def get_optimal_k(X, is_Kmeans : bool):
    '''
    is_Kmeans == True : K-means
    is_Kmeans == False : hierarchical
    ''' 
    k_values = range(2, 11)
    # Store WCSS and Silhouette Scores
    wcss = []
    silhouette_scores = []

    # Compute WCSS and Silhouette Score for each k
    for k in k_values:
        # Fit k-means model
        if is_Kmeans :
            model = KMeans(n_clusters=k, random_state=42)
        else :
            model = AgglomerativeClustering(n_clusters = k)   
        model.fit_predict(X)
        # Compute Silhouette Score
        labels = model.labels_
        silhouette_avg = silhouette_score(X, labels)
        silhouette_scores.append(silhouette_avg)
        # Tránh lỗi silhouette_score khi có cụm chỉ có 1 điểm
        #if len(set(labels)) > 1:
        #    silhouette_avg = silhouette_score(X, labels)
        #    silhouette_scores.append(silhouette_avg)
        #else:
        #    silhouette_scores.append(-1)
    optimal_k = k_values[np.argmax(silhouette_scores)]
    return optimal_k
