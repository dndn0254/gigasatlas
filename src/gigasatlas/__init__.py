"""gigasatlas: name the cells of new Pacific oyster single-cell data with the harmonised atlas.

Quick start
-----------
>>> import anndata as ad, gigasatlas as ga
>>> A = ga.read("my_sample.h5ad")               # .h5ad, 10x folder or count table; raw counts, xbMagGiga1.1 gene IDs
>>> ga.gene_coverage(A, "hemocyte")             # should be > 0.8
>>> ga.predict(A, "hemocyte", groupby="seurat_clusters")  # adds gigas_label, gigas_prob, gigas_status
"""
from .annotate import gene_coverage, log_normalize, predict
from .datasets import example_data, example_gene_list, example_old_ids, example_path
from .io import read, read_10x_mtx, read_table
from .enrich import cell_type_enrichment
from .ids import convert, rename_var
from .markers import score_clusters
from .model import AtlasModel, available_models, load_model, train_model

__version__ = "0.1.0"
__all__ = ["read", "read_10x_mtx", "read_table", "predict", "gene_coverage", "log_normalize", "cell_type_enrichment", "convert", "rename_var",
           "example_data", "example_gene_list", "example_old_ids", "example_path",
           "score_clusters", "AtlasModel", "available_models", "load_model", "train_model"]
