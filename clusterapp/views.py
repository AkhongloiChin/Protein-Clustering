import matplotlib
matplotlib.use('Agg')
from django.shortcuts import render, redirect
from .forms import FastaUploadForm
from .models import FastaFile
from django.views.generic import DetailView
import numpy as np
from sklearn.cluster import KMeans
from itertools import product
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import base64
from io import BytesIO


def parse_fasta(content):
    sequences = []
    headers = []
    current_header = None
    current_seq = []
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header is not None:
                headers.append(current_header)
                sequences.append(''.join(current_seq))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_header is not None:
        headers.append(current_header)
        sequences.append(''.join(current_seq))
    return list(zip(headers, sequences))

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

def upload_fasta(request):
    if request.method == 'POST':
        form = FastaUploadForm(request.POST, request.FILES)
        if form.is_valid():
            fasta = form.save()
            return redirect('cluster_results', pk=fasta.pk)
    else:
        form = FastaUploadForm()
    return render(request, 'upload.html', {'form': form})

class ClusterResultsView(DetailView):
    model = FastaFile
    template_name = 'results.html'
    context_object_name = 'fasta'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        fasta = self.object
        
        # Read and parse FASTA file
        content = fasta.file.read().decode('utf-8')
        sequences = parse_fasta(content)
        
        # Feature extraction
        k = 2  # Using 2-mers
        features = [get_kmer_counts(seq, k) for header, seq in sequences]
        X = np.array(features)
        
        n_clusters = fasta.num_clusters  
        kmeans = KMeans(n_clusters=n_clusters)
        labels = kmeans.fit_predict(X)
        
        clusters = {}
        for i, (header, _) in enumerate(sequences):
            cluster_id = labels[i]
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(header)
        
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X)
        
        plt.figure(figsize=(8, 6))
        for cluster_id in range(n_clusters):
            plt.scatter(X_pca[labels == cluster_id, 0], X_pca[labels == cluster_id, 1], label=f'Cluster {cluster_id}')
        plt.title('Protein Clusters (2D PCA)')
        plt.xlabel('PCA Component 1')
        plt.ylabel('PCA Component 2')
        plt.legend()

        buffer = BytesIO()
        plt.savefig(buffer, format='png')
        buffer.seek(0)
        cluster_plot = base64.b64encode(buffer.getvalue()).decode('utf-8')
        plt.close()
        
        context['clusters'] = clusters
        context['n_clusters'] = n_clusters
        context['cluster_plot'] = cluster_plot
        return context
        