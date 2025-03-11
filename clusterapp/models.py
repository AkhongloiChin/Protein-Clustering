from django.db import models

class FastaFile(models.Model):
    file = models.FileField(upload_to='fasta_files/')
    num_clusters = models.IntegerField()
    uploaded_at = models.DateTimeField(auto_now_add=True)