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
            'kmeans_2': joblib.load(os.path.join(models_dir, 'k_means2.joblib')),
            'kmeans_3': joblib.load(os.path.join(models_dir, 'k_means3.joblib')),
            'kmeans_4': joblib.load(os.path.join(models_dir, 'k_means4.joblib')),
            'kmeans_5': joblib.load(os.path.join(models_dir, 'k_means5.joblib')),
            'kmeans_6': joblib.load(os.path.join(models_dir, 'k_means6.joblib')),
            'kmeans_7': joblib.load(os.path.join(models_dir, 'k_means7.joblib')),
            'kmeans_8': joblib.load(os.path.join(models_dir, 'k_means8.joblib')),
            'kmeans_9': joblib.load(os.path.join(models_dir, 'k_means9.joblib')),
            'kmeans_10': joblib.load(os.path.join(models_dir, 'k_means10.joblib')),

            'hier_average_2': joblib.load(os.path.join(models_dir, 'hier_average2.joblib')),
            'hier_average_3': joblib.load(os.path.join(models_dir, 'hier_average3.joblib')),
            'hier_average_4': joblib.load(os.path.join(models_dir, 'hier_average4.joblib')),
            'hier_average_5': joblib.load(os.path.join(models_dir, 'hier_average5.joblib')),
            'hier_average_6': joblib.load(os.path.join(models_dir, 'hier_average6.joblib')),
            'hier_average_7': joblib.load(os.path.join(models_dir, 'hier_average7.joblib')),
            'hier_average_8': joblib.load(os.path.join(models_dir, 'hier_average8.joblib')),
            'hier_average_9': joblib.load(os.path.join(models_dir, 'hier_average9.joblib')),
            'hier_average_10': joblib.load(os.path.join(models_dir, 'hier_average10.joblib')),

            'hier_complete_2': joblib.load(os.path.join(models_dir, 'hier_complete2.joblib')),
            'hier_complete_3': joblib.load(os.path.join(models_dir, 'hier_complete3.joblib')),
            'hier_complete_4': joblib.load(os.path.join(models_dir, 'hier_complete4.joblib')),
            'hier_complete_5': joblib.load(os.path.join(models_dir, 'hier_complete5.joblib')),
            'hier_complete_6': joblib.load(os.path.join(models_dir, 'hier_complete6.joblib')),
            'hier_complete_7': joblib.load(os.path.join(models_dir, 'hier_complete7.joblib')),
            'hier_complete_8': joblib.load(os.path.join(models_dir, 'hier_complete8.joblib')),
            'hier_complete_9': joblib.load(os.path.join(models_dir, 'hier_complete9.joblib')),
            'hier_complete_10': joblib.load(os.path.join(models_dir, 'hier_complete10.joblib')),

            'hier_single_2': joblib.load(os.path.join(models_dir, 'hier_single2.joblib')),
            'hier_single_3': joblib.load(os.path.join(models_dir, 'hier_single3.joblib')),
            'hier_single_4': joblib.load(os.path.join(models_dir, 'hier_single4.joblib')),
            'hier_single_5': joblib.load(os.path.join(models_dir, 'hier_single5.joblib')),
            'hier_single_6': joblib.load(os.path.join(models_dir, 'hier_single6.joblib')),
            'hier_single_7': joblib.load(os.path.join(models_dir, 'hier_single7.joblib')),
            'hier_single_8': joblib.load(os.path.join(models_dir, 'hier_single8.joblib')),
            'hier_single_9': joblib.load(os.path.join(models_dir, 'hier_single9.joblib')),
            'hier_single_10': joblib.load(os.path.join(models_dir, 'hier_single10.joblib')),
        }

        for key in self.cluster_models:
            if key.startswith('kmeans_') and hasattr(self.cluster_models[key], 'n_init'):
                self.cluster_models[key].n_init = 'auto'

        for key in self.cluster_models:
            if key.startswith(('hier_average_', 'hier_complete_', 'hier_single_')):
                if hasattr(self.cluster_models[key], 'metric') and self.cluster_models[key].metric is None:
                    self.cluster_models[key].metric = 'euclidean'

        print("✅ Clustering models loaded successfully!")
