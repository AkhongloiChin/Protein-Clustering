from io import BytesIO
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from django.shortcuts import render, redirect
from django.views.generic import DetailView
from .forms import FastaUploadForm
from .models import FastaFile

import scipy.cluster.hierarchy as sch
from scipy.stats.mstats import winsorize

from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

from django.apps import apps
from django.http import JsonResponse
from .utils import get_optimal_k, get_kmer_counts, hierarchical_clustering, kmeans_clustering
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

        k = 2  # Using 2-mers
        features = [get_kmer_counts(seq) for _, seq in sequences]
        X = np.array(features)
        if X.shape[0] == 0:
            context['error'] = "No valid sequences found."
            return context
        winsorized_X = winsorize(X, limits=[0.05, 0.05])
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(winsorized_X)  

        pca = PCA(n_components=3)
        X_pca = pca.fit_transform(X_scaled)

        cluster_app = apps.get_app_config('clusterapp')

        if fasta.model_choice == 'kmeans':
            n_clusters = get_optimal_k(X_pca, is_Kmeans = True)
            model_name = 'kmeans_'+str(n_clusters)
            #model = KMeans(n_clusters=n_clusters)
            model = cluster_app.cluster_models[model_name]
            labels = model.fit_predict(X_pca)
            context['clusters'], context['cluster_plot'] = kmeans_clustering(X_pca, labels, sequences, n_clusters)
            context['dendrogram_plot'] = None

        elif fasta.model_choice == 'hierarchical':
            n_clusters = get_optimal_k(X_pca, is_Kmeans=False, linkage=fasta.linkage)
            if fasta.linkage == "centroid":
                linkage_matrix = sch.linkage(X_pca, method="centroid")
                labels = sch.fcluster(linkage_matrix, t=n_clusters, criterion="maxclust")
                labels -= labels.min() 
                _ , context['dendrogram_plot'] = hierarchical_clustering(X_pca, labels, sequences, linkage_method= "centroid", n_clusters = n_clusters)
                context['clusters'], context['cluster_plot'] = kmeans_clustering(X_pca, labels,sequences,n_clusters)
            else:
                linkage_name = "hier_" + fasta.linkage + '_' + str(n_clusters)
                #model = AgglomerativeClustering(n_clusters=n_clusters, linkage=fasta.linkage)
                model = cluster_app.cluster_models[linkage_name]
                labels = model.fit_predict(X_pca)
                context['clusters'], context['cluster_plot'] = kmeans_clustering(X_pca, labels, sequences, n_clusters = n_clusters)
                _ , context['dendrogram_plot'] = hierarchical_clustering(X_pca, labels, sequences, linkage_method= fasta.linkage, n_clusters = n_clusters)
        else:
            context['error'] = "Invalid clustering method."
        if len(set(labels)) > 1:  # Kiểm tra nếu có hơn 1 cluster
            silhouette_avg = silhouette_score(X_pca, labels)
        else:
            silhouette_avg = None  # Tránh lỗi khi chỉ có 1 cụm
        context['silhouette_score'] = silhouette_avg
        return context
