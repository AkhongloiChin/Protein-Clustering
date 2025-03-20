from io import BytesIO
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from django.shortcuts import render, redirect
from django.views.generic import DetailView
from .forms import FastaUploadForm
from .models import FastaFile

from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats.mstats import winsorize
from .utils import get_optimal_k, get_kmer_counts
import numpy as np

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

def upload_fasta(request):
    if request.method == 'POST':
        print("Received POST request")  
        form = FastaUploadForm(request.POST, request.FILES)
        if form.is_valid():
            print("Form is valid, saving...")
            fasta = form.save()
            print(f"Saved file: {fasta.file.name}")
            return redirect('cluster_results', pk=fasta.pk)
        else:
            print("Form is invalid:", form.errors) 
    else:
        form = FastaUploadForm()
    return render(request, 'upload.html', {'form': form})


class ClusterResultsView(DetailView):
    model = FastaFile
    template_name = 'results.html'
    context_object_name = 'fasta'

    def get_context_data(self, **kwargs):
        print("got it")
        context = super().get_context_data(**kwargs)
        fasta = self.object

        # Read and parse FASTA file
        with fasta.file.open() as f:
            content = f.read().decode('utf-8')

        sequences = parse_fasta(content)
        
        ########################
        k = 2  # Using 2-mers
        features = [get_kmer_counts(seq) for _, seq in sequences]
        X = np.array(features)
        if X.shape[0] == 0:
            context['error'] = "No valid sequences found."
            return context
        winsorized_X = winsorize(X, limits=[0.05, 0.05])
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(winsorized_X)  

        #umap_model = umap.UMAP(n_components=50, random_state=42)
        #X_umap = umap_model.fit_transform(X_scaled)
        pca = PCA(n_components=50)
        X_pca = pca.fit_transform(X_scaled)

        if fasta.model_choice == 'kmeans':
            n_clusters = get_optimal_k(X_pca, is_Kmeans = True)
            model = KMeans(n_clusters=n_clusters)

        elif fasta.model_choice == 'hierarchical':
            n_clusters = get_optimal_k(X_pca, is_Kmeans= False)
            model = AgglomerativeClustering(n_clusters=3)

        else:
            context['error'] = "Invalid clustering method."
        labels = model.fit_predict(X_pca)
        clusters = {}
        for i, (header, _) in enumerate(sequences):
            cluster_id = labels[i]
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(header)
        
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
        context['cluster_plot'] = cluster_plot
        return context
