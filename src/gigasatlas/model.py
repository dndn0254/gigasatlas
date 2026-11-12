"""Reference models: a multinomial logistic regression per atlas, stored as one small .npz file.

The model works on log1p(CP10k) expression of a fixed gene panel, standardised with the
training mean and standard deviation (the same idea as CellTypist, Dominguez Conde et al. 2022).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MODEL_DIR = os.path.join(DATA_DIR, "models")
UNNAMED = "Unnamed by authors"


@dataclass
class AtlasModel:
    name: str
    genes: np.ndarray            # (G,) gene IDs (xbMagGiga1.1 RefSeq symbols / LOC IDs)
    labels: np.ndarray           # (L,) cell-type names
    coef: np.ndarray             # (L, G) float32
    intercept: np.ndarray        # (L,)
    mean: np.ndarray             # (G,) training mean of log1p(CP10k)
    scale: np.ndarray            # (G,) training standard deviation
    meta: dict = field(default_factory=dict)   # lineage of each label, version, training notes

    # ------------------------------------------------------------------ io
    def save(self, path: str) -> None:
        np.savez_compressed(path, genes=self.genes.astype(str), labels=self.labels.astype(str),
                            coef=self.coef.astype(np.float32), intercept=self.intercept.astype(np.float32),
                            mean=self.mean.astype(np.float32), scale=self.scale.astype(np.float32),
                            meta=np.array(json.dumps({"name": self.name, **self.meta})))

    @classmethod
    def load(cls, path: str) -> "AtlasModel":
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        # cells the authors left unnamed were stored as the text 'nan', which spreadsheets and pandas read as missing
        labels = np.array([UNNAMED if str(x) == "nan" else str(x) for x in z["labels"]])
        if "nan" in meta.get("lineage", {}):
            meta["lineage"][UNNAMED] = meta["lineage"].pop("nan")
        return cls(name=meta.pop("name"), genes=z["genes"], labels=labels, coef=z["coef"],
                   intercept=z["intercept"], mean=z["mean"], scale=z["scale"], meta=meta)

    # ------------------------------------------------------------------ prediction
    def decision(self, X_log: "sp.spmatrix | np.ndarray", genes) -> np.ndarray:
        """Class probabilities for cells x genes log1p(CP10k) matrix with gene names `genes`."""
        X = _align(X_log, list(genes), list(self.genes))
        X = (X - self.mean) / self.scale
        np.clip(X, -10, 10, out=X)
        z = X @ self.coef.T + self.intercept
        z -= z.max(1, keepdims=True)
        p = np.exp(z)
        p /= p.sum(1, keepdims=True)
        return p

    def lineage_of(self, label: str) -> str:
        return self.meta.get("lineage", {}).get(label, "")


def _align(X, have: list, want: list) -> np.ndarray:
    """Columns of X (genes `have`) reordered to `want`; genes not measured become 0."""
    idx = {g: i for i, g in enumerate(have)}
    cols = np.array([idx.get(g, -1) for g in want])
    ok = cols >= 0
    out = np.zeros((X.shape[0], len(want)), dtype=np.float32)
    sub = X[:, cols[ok]]
    out[:, ok] = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    return out


def available_models() -> list[str]:
    if not os.path.isdir(MODEL_DIR):
        return []
    return sorted(f[:-4] for f in os.listdir(MODEL_DIR) if f.endswith(".npz"))


def load_model(name_or_path: str) -> AtlasModel:
    """Load a bundled model by name (e.g. 'hemocyte', 'lineage') or any .npz path."""
    path = name_or_path if name_or_path.endswith(".npz") else os.path.join(MODEL_DIR, f"{name_or_path}.npz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"model '{name_or_path}' not found; bundled models: {available_models()}")
    return AtlasModel.load(path)


def train_model(X_log, genes, labels, name: str, n_genes: int = 3000, max_per_label: int = 600,
                C: float = 0.1, seed: int = 0, lineage: dict | None = None, gene_panel=None) -> AtlasModel:
    """Train a reference model from a log1p(CP10k) matrix and one label per cell.

    Gene panel: `gene_panel` if given, else the `n_genes` most variable genes (dispersion of
    log expression). Each label is capped at `max_per_label` cells so large labels do not dominate.
    Requires scikit-learn.
    """
    from sklearn.linear_model import LogisticRegression

    labels = np.asarray(labels).astype(str)
    genes = np.asarray(genes).astype(str)
    rng = np.random.default_rng(seed)
    keep = []
    for lab in np.unique(labels):
        ii = np.where(labels == lab)[0]
        keep.append(ii if len(ii) <= max_per_label else rng.choice(ii, max_per_label, replace=False))
    keep = np.sort(np.concatenate(keep))
    Xs = sp.csr_matrix(X_log[keep]) if sp.issparse(X_log) else np.asarray(X_log[keep])
    if gene_panel is None:          # gene choice on the sparse matrix, so only the panel is ever made dense
        m = np.asarray(Xs.mean(0)).ravel()
        sq = Xs.multiply(Xs).mean(0) if sp.issparse(Xs) else (Xs ** 2).mean(0)
        v = np.asarray(sq).ravel() - m ** 2
        disp = np.where(m > 0.01, v / np.maximum(m, 1e-9), 0)
        top = np.argsort(-disp)[:n_genes]
    else:
        gi = {g: i for i, g in enumerate(genes)}
        top = np.array([gi[g] for g in gene_panel if g in gi])
    Xd = Xs[:, top]
    Xd = (Xd.toarray() if sp.issparse(Xd) else np.asarray(Xd)).astype(np.float32)
    mean = Xd.mean(0); scale = Xd.std(0); scale[scale == 0] = 1
    Z = np.clip((Xd - mean) / scale, -10, 10)
    clf = LogisticRegression(C=C, max_iter=2000, n_jobs=None)
    clf.fit(Z, labels[keep])
    return AtlasModel(name=name, genes=genes[top], labels=clf.classes_.astype(str), coef=clf.coef_.astype(np.float32),
                      intercept=clf.intercept_.astype(np.float32), mean=mean, scale=scale,
                      meta={"lineage": lineage or {}, "n_train_cells": int(len(keep)), "C": C,
                            "max_per_label": max_per_label, "seed": seed})
