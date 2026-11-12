"""Small example data shipped with the package, so every command can be tried at once.

example_data()        400 gastrula cells (raw counts) from Gavery et al. 2025 (PRJNA909306, library E4),
                      with the authors' cell-type name in obs['author_label'].
example_gene_list()   200 genes higher in sorted granulocytes than hyalinocytes
                      (bulk RNA-seq of Mao et al. 2020, re-analysed with DESeq2 on xbMagGiga1.1).
example_old_ids()     Roslin v1 gene IDs (Ensembl G##### and RefSeq LOC) as used in older papers.

Note: the bundled 'development' model was trained on all four gastrula libraries, so names on the
example cells look better than on truly new data (held-out accuracy 0.91-0.92, see README).
"""
from __future__ import annotations

import os

from .model import DATA_DIR

EX_DIR = os.path.join(DATA_DIR, "examples")


def example_path(name: str) -> str:
    """Path of a bundled example file: 'gastrula_400.h5ad', 'granulocyte_up_genes.txt', 'roslin_ids.txt'."""
    p = os.path.join(EX_DIR, name)
    if not os.path.exists(p):
        raise FileNotFoundError(f"{name} not in {EX_DIR}: {sorted(os.listdir(EX_DIR))}")
    return p


def example_data():
    import anndata as ad
    from .model import UNNAMED
    A = ad.read_h5ad(example_path("gastrula_400.h5ad"))
    A.obs["author_label"] = A.obs["author_label"].astype(str).replace("nan", UNNAMED).astype("category")
    return A


def _lines(name):
    return [x.strip() for x in open(example_path(name)) if x.strip() and not x.startswith("#")]


def example_gene_list() -> list:
    return _lines("granulocyte_up_genes.txt")


def example_old_ids() -> list:
    return _lines("roslin_ids.txt")
