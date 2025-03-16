from django import forms
from .models import FastaFile

class FastaUploadForm(forms.ModelForm):
    class Meta:
        model = FastaFile
        fields = ['file', 'num_clusters']   