from django import forms
from .models import FastaFile
import os

class FastaUploadForm(forms.ModelForm):
    class Meta:
        model = FastaFile
        fields = ['file', 'model_choice', 'linkage'] 

    model_choice = forms.ChoiceField(
        choices=[('kmeans', 'K-Means'), ('hierarchical', 'Hierarchical')],
        required=True
    )

    linkage = forms.ChoiceField(
        choices=[
            ('single', 'Single Linkage'),
            ('complete', 'Complete Linkage'),
            ('average', 'Average Linkage'),
            ('centroid', 'Centroid Linkage')
        ],
        required=False  
    )

    def clean_file(self):
        """ Validate that the uploaded file is a FASTA file """
        file = self.cleaned_data.get('file')

        if file:
            ext = os.path.splitext(file.name)[1].lower()
            allowed_extensions = ['.fasta', '.fa']

            if ext not in allowed_extensions:
                raise forms.ValidationError("Only FASTA files (.fasta, .fa) are allowed.")

        return file

    def clean(self):
        """ Ensure linkage is selected only for Hierarchical clustering """
        cleaned_data = super().clean()
        model_choice = cleaned_data.get('model_choice')
        linkage = cleaned_data.get('linkage')

        if model_choice == 'hierarchical' and not linkage:
            raise forms.ValidationError("Please select a linkage method for hierarchical clustering.")

        if model_choice == 'kmeans':
            cleaned_data['linkage'] = None

        return cleaned_data
