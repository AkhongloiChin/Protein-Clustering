import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from itertools import product
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import base64
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as sch

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

def get_optimal_k(X, is_Kmeans : bool, linkage = None):
    '''
    is_Kmeans == True : K-means
    is_Kmeans == False : hierarchical
    ''' 
    k_values = range(2, 11)
    silhouette_scores = []

    # Compute WCSS and Silhouette Score for each k
    for k in k_values:
        # Fit k-means model
        if linkage == "centroid":
            linkage_matrix = sch.linkage(X, method="centroid")
            labels = sch.fcluster(linkage_matrix, t= k, criterion="maxclust")
        elif is_Kmeans :
            model = KMeans(n_clusters=k, random_state=42)
            model.fit_predict(X)
            labels = model.labels_
        else:
            model = AgglomerativeClustering(n_clusters= k, linkage=linkage)
            model.fit_predict(X)
            labels = model.labels_
        # Compute Silhouette Score
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

def kmeans_clustering(X_pca, labels, sequences, n_clusters=3):
    clusters = {}
    for i, (header, _) in enumerate(sequences):
        cluster_id = labels[i]
        clusters.setdefault(cluster_id, []).append(header)
    sorted_clusters = dict(sorted(clusters.items()))
    # Scatter plot
    plt.figure(figsize=(8, 6))
    for cluster_id in range(n_clusters):
        plt.scatter(X_pca[labels == cluster_id, 0], X_pca[labels == cluster_id, 1], label=f'Cluster {cluster_id}')
    plt.title('K-Means Protein Clusters')
    plt.legend()

    # Save and encode the plot
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    cluster_plot = base64.b64encode(buffer.getvalue()).decode('utf-8')
    plt.close()

    return sorted_clusters, cluster_plot

def hierarchical_clustering(X_pca, labels, sequences, linkage_method='ward'):
    """Performs Hierarchical clustering and generates a dendrogram."""
    clusters = {}
    for i, (header, _) in enumerate(sequences):
        cluster_id = labels[i]
        clusters.setdefault(cluster_id, []).append(header)
    sorted_clusters = dict(sorted(clusters.items()))
    # Generate dendrogram
    plt.figure(figsize=(10, 6))
    Z = linkage(X_pca, method=linkage_method)
    dendrogram(Z, labels=[header for header, _ in sequences], leaf_rotation=90, leaf_font_size=10)
    plt.title(f'Hierarchical Clustering Dendrogram ({linkage_method} linkage)')

    # Save and encode the dendrogram
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    dendrogram_plot = base64.b64encode(buffer.getvalue()).decode('utf-8')
    plt.close()

    return sorted_clusters, dendrogram_plot
