"""Annotate new oyster single-cell data with the atlas reference models."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .model import AtlasModel, load_model


def log_normalize(X, target_sum: float = 1e4):
    """Raw counts -> log1p(counts per 10,000). Works on sparse or dense matrices."""
    if sp.issparse(X):
        X = sp.csr_matrix(X, dtype=np.float32)
        tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
        X = sp.diags(target_sum / tot) @ X
        X.data = np.log1p(X.data)
        return X
    X = np.asarray(X, dtype=np.float32)
    tot = X.sum(1, keepdims=True); tot[tot == 0] = 1
    return np.log1p(X / tot * target_sum)


def _looks_like_counts(X) -> bool:
    d = X.data if sp.issparse(X) else np.asarray(X).ravel()
    d = d[: 100000]
    return d.size > 0 and np.allclose(d, np.round(d)) and d.max() > 20


def predict(adata, model: "str | AtlasModel" = "lineage", layer: str | None = None, min_prob: float = 0.5,
            groupby: str | None = None, key: str = "gigas") -> pd.DataFrame:
    """Predict a cell-type name for every cell of an AnnData object.

    adata.X (or `layer`) may be raw counts (recommended) or log1p(CP10k); raw counts are detected
    and normalised automatically. Gene names must be xbMagGiga1.1 RefSeq symbols / LOC IDs; use
    `gigasatlas.ids.convert` first if your data use Roslin v1 IDs.

    Adds to adata.obs:
      {key}_label      best name
      {key}_prob       its probability
      {key}_lineage    lineage of the name (if known)
      {key}_status     'confident' (prob >= min_prob) or 'uncertain'
      {key}_label_group  majority name of each `groupby` group (if groupby is given)
    Returns the full cells x names probability table.
    """
    m = load_model(model) if isinstance(model, str) else model
    X = adata.layers[layer] if layer else adata.X
    if _looks_like_counts(X):
        X = log_normalize(X)
    P = m.decision(X, adata.var_names)
    best = P.argmax(1)
    lab = m.labels[best]; pr = P[np.arange(len(best)), best]
    adata.obs[f"{key}_label"] = pd.Categorical(lab)
    adata.obs[f"{key}_prob"] = pr
    adata.obs[f"{key}_lineage"] = [m.lineage_of(x) for x in lab]
    adata.obs[f"{key}_status"] = np.where(pr >= min_prob, "confident", "uncertain")
    if groupby is not None:
        maj = adata.obs.groupby(groupby, observed=True)[f"{key}_label"].agg(lambda s: s.value_counts().index[0])
        adata.obs[f"{key}_label_group"] = adata.obs[groupby].map(maj).astype("category")
    frac = float((pr >= min_prob).mean())
    adata.uns[f"{key}_model"] = {"model": m.name, "n_genes_model": int(len(m.genes)),
                                 "n_genes_found": int(np.isin(m.genes, adata.var_names).sum()), "frac_confident": frac}
    return pd.DataFrame(P, index=adata.obs_names, columns=m.labels)


def gene_coverage(adata, model: "str | AtlasModel" = "lineage") -> float:
    """Fraction of the model's genes present in adata.var_names (check this first: < 0.5 means wrong gene IDs)."""
    m = load_model(model) if isinstance(model, str) else model
    return float(np.isin(m.genes, adata.var_names).mean())
