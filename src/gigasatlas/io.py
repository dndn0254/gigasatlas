"""Read single-cell count data in the common formats, so a user does not need scanpy to start.

Supported:
  .h5ad                         AnnData (raw counts in .X, or in a layer)
  10x folder                    matrix.mtx(.gz) + features.tsv(.gz) or genes.tsv + barcodes.tsv(.gz)
                                (Cell Ranger / STARsolo 'filtered' output)
  .csv / .tsv (.gz)             a count table; genes as rows and cells as columns (the usual export),
                                or cells as rows (detected from which side matches the model genes)
"""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .model import load_model


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _first(folder, names):
    for n in names:
        p = os.path.join(folder, n)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"none of {names} in {folder}")


def _model_genes():
    return set(load_model("lineage").genes) | set(load_model("development").genes)


def read_10x_mtx(folder: str):
    """Cell Ranger / STARsolo folder -> AnnData (cells x genes, raw counts).

    features.tsv has a gene ID column and a gene name column; the one that matches more
    xbMagGiga1.1 IDs (LOC... / symbols) is used as var_names.
    """
    import anndata as ad
    mtx = _first(folder, ["matrix.mtx.gz", "matrix.mtx"])
    feat = _first(folder, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
    bc = _first(folder, ["barcodes.tsv.gz", "barcodes.tsv"])
    with _open(mtx) as fh:
        X = sp.csr_matrix(scipy.io.mmread(fh).T, dtype=np.float32)       # file is genes x cells
    F = pd.read_csv(feat, sep="\t", header=None, dtype=str)
    B = pd.read_csv(bc, sep="\t", header=None, dtype=str)[0].values
    known = _model_genes()
    cols = [c for c in F.columns[:2]]
    best = max(cols, key=lambda c: F[c].isin(known).sum())
    A = ad.AnnData(X=X, obs=pd.DataFrame(index=B), var=pd.DataFrame({"gene_ids": F[0].values}, index=F[best].values))
    A.var_names_make_unique()
    return A


def read_table(path: str):
    """Count table (.csv/.tsv, optionally .gz) -> AnnData. Orientation is detected from gene names."""
    import anndata as ad
    sep = "," if ".csv" in os.path.basename(path) else "\t"
    D = pd.read_csv(path, sep=sep, index_col=0)
    known = _model_genes()
    if D.columns.isin(known).sum() > D.index.isin(known).sum():
        D = D.T                                   # was cells x genes
    D = D.apply(pd.to_numeric, errors="coerce").fillna(0)
    return ad.AnnData(X=sp.csr_matrix(D.T.values.astype(np.float32)), obs=pd.DataFrame(index=D.columns.astype(str)),
                      var=pd.DataFrame(index=D.index.astype(str)))


def read(path: str):
    """Read .h5ad, a 10x folder, or a .csv/.tsv count table into AnnData (cells x genes)."""
    import anndata as ad
    if os.path.isdir(path):
        return read_10x_mtx(path)
    low = path.lower()
    if low.endswith(".h5ad"):
        return ad.read_h5ad(path)
    if any(low.endswith(e) for e in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")):
        return read_table(path)
    raise ValueError(f"unknown input format: {path} (use .h5ad, a 10x folder, or a .csv/.tsv count table)")
