# Generated by Django 5.1.7 on 2025-03-20 07:57

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("clusterapp", "0001_initial"),
    ]

    operations = [
        migrations.AddField(
            model_name="fastafile",
            name="linkage",
            field=models.CharField(
                blank=True,
                choices=[
                    ("single", "Single Linkage"),
                    ("complete", "Complete Linkage"),
                    ("average", "Average Linkage"),
                    ("centroid", "Centroid Linkage"),
                ],
                max_length=20,
                null=True,
            ),
        ),
    ]
