from django import forms
from .models import FastaFile
import os

class FastaUploadForm(forms.ModelForm):
    class Meta:
        model = FastaFile
        fields = ['file', 'model_choice']

    model_choice = forms.ChoiceField(
        choices=[('kmeans', 'K-Means'), ('hierarchical', 'Hierarchical')],
        required=True
    )

    def clean_file(self):
        file = self.cleaned_data.get('file')

        if file:
            ext = os.path.splitext(file.name)[1].lower()  # Get file extension
            allowed_extensions = ['.fasta', '.fa']

            if ext not in allowed_extensions:
                raise forms.ValidationError("Only FASTA files (.fasta, .fa) are allowed.")

        return file
