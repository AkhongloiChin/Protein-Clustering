from django.db import models

class FastaFile(models.Model):
    file = models.FileField(upload_to='fasta_files/')
    model_choice = models.CharField(
        max_length=20,
        choices=[
            ('kmeans', 'K-Means Clustering'),
            ('hierarchical', 'Hierarchical Clustering'),
        ],
        default='kmeans'
    )
