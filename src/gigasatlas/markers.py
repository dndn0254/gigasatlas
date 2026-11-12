"""Marker-set scoring for clusters (use when the tissue is not in the reference, e.g. gill)."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .annotate import _looks_like_counts, log_normalize
from .model import DATA_DIR


def load_markers(path: str | None = None) -> pd.DataFrame:
    """Marker table: columns `set` (cell type), `gene` (xbMagGiga1.1 ID). Default = atlas top-specific genes."""
    path = path or os.path.join(DATA_DIR, "label_markers.tsv")
    return pd.read_csv(path, sep="\t")


def score_clusters(adata, groupby: str, markers: pd.DataFrame | None = None, top_n: int = 50,
                   min_genes: int = 5) -> pd.DataFrame:
    """Mean log expression of each marker set minus the panel mean, per cluster, z-scored across clusters.

    Returns clusters x sets z-scores; the best set per cluster is `df.idxmax(1)`.
    """
    mk = load_markers() if markers is None else markers
    if "rank" in mk:
        mk = mk[mk["rank"] <= top_n]
    X = adata.X
    if _looks_like_counts(X):
        X = log_normalize(X)
    g = pd.Series(np.arange(adata.n_vars), index=adata.var_names)
    groups = adata.obs[groupby].astype(str).values
    ug = np.unique(groups)
    G = sp.csr_matrix((np.ones(len(groups)), ([np.where(ug == x)[0][0] for x in groups], np.arange(len(groups)))),
                      shape=(len(ug), len(groups)))
    G = sp.diags(1 / np.asarray(G.sum(1)).ravel()) @ G
    M = G @ X                                   # clusters x genes mean
    M = M.toarray() if sp.issparse(M) else np.asarray(M)
    bg = M.mean(1, keepdims=True)
    out = {}
    for s, d in mk.groupby("set"):
        ii = g.reindex(d["gene"]).dropna().astype(int).values
        if len(ii) >= min_genes:
            out[s] = M[:, ii].mean(1) - bg[:, 0]
    S = pd.DataFrame(out, index=ug)
    return (S - S.mean()) / S.std(ddof=0).replace(0, 1)
