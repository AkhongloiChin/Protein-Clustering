from django.apps import AppConfig
import joblib
import os

class ClusterAppConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'clusterapp'

    def ready(self):
        """Load clustering models once when Django starts"""
        from django.conf import settings
        models_dir = os.path.join(settings.BASE_DIR, 'models')  # Path to models folder

        # Load multiple models (adjust based on what you have)
        self.cluster_models = {
            'kmeans': joblib.load(os.path.join(models_dir, 'k_means26.joblib')),
            'hier_average': joblib.load(os.path.join(models_dir, 'average27.joblib')),
            'hier_complete': joblib.load(os.path.join(models_dir, 'complete26.joblib')),
            'hier_single': joblib.load(os.path.join(models_dir, 'single22.joblib')),
        }
        if hasattr(self.cluster_models['kmeans'], 'n_init'):
            self.cluster_models['kmeans'].n_init = 'auto'
        for key in ['hier_average', 'hier_complete', 'hier_single']:
            if hasattr(self.cluster_models[key], 'metric') and self.cluster_models[key].metric is None:
                self.cluster_models[key].metric = 'euclidean'
        print("✅ Clustering models loaded successfully!")
