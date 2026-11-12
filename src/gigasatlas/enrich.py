"""Which atlas cell types does a gene list (for example bulk DE genes) point to?"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from .model import DATA_DIR


def cell_type_enrichment(genes, top_n: int = 100, universe: int | None = None) -> pd.DataFrame:
    """Hypergeometric test of `genes` against each cell type's top-`top_n` specific genes.

    Universe = all genes in the marker table (or `universe`). Returns one row per cell type with
    overlap, fold enrichment, p and Benjamini-Hochberg q, sorted by p.
    """
    mk = pd.read_csv(os.path.join(DATA_DIR, "label_markers.tsv"), sep="\t")
    mk = mk[mk["rank"] <= top_n]
    gset = set(map(str, genes))
    U = universe or pd.read_csv(os.path.join(DATA_DIR, "gene_universe.tsv"), sep="\t").shape[0]
    q = len(gset)
    rows = []
    for s, d in mk.groupby("set"):
        ms = set(d["gene"]); k = len(ms & gset)
        p = hypergeom.sf(k - 1, U, len(ms), q) if k else 1.0
        rows.append((s, d["dataset"].iat[0], d["lineage"].iat[0], len(ms), k, (k / q) / (len(ms) / U) if q else 0, p,
                     ", ".join(sorted(ms & gset)[:15])))
    r = pd.DataFrame(rows, columns=["cell_type", "dataset", "lineage", "set_size", "overlap", "fold", "p", "genes"])
    r = r.sort_values("p").reset_index(drop=True)
    n = len(r); rk = np.arange(1, n + 1)
    r["q"] = np.minimum.accumulate((r["p"] * n / rk)[::-1])[::-1].clip(upper=1)
    return r
