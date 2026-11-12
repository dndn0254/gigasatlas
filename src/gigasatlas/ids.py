"""Gene ID conversion: Roslin v1 (Ensembl Metazoa G##### or RefSeq LOC) -> xbMagGiga1.1 RefSeq."""
from __future__ import annotations

import os

import pandas as pd

from .model import DATA_DIR


def crosswalk() -> pd.DataFrame:
    """Columns: old_id, source ('ensembl_rbh' | 'ncbi_compare_prev'), new_id, category, note."""
    return pd.read_csv(os.path.join(DATA_DIR, "id_crosswalk.tsv.gz"), sep="\t", dtype=str)


def convert(ids, keep_unmapped: bool = False) -> pd.DataFrame:
    """Convert a list of old gene IDs. One row per (old, new) pair; 1-to-many mappings give several rows."""
    cw = crosswalk()
    q = pd.DataFrame({"old_id": [str(x).strip() for x in ids]})
    out = q.merge(cw, on="old_id", how="left")
    if not keep_unmapped:
        out = out[out["new_id"].notna() & (out["new_id"] != "")]
    return out.reset_index(drop=True)


def rename_var(adata, drop_unmapped: bool = True, how: str = "first"):
    """Return a copy of adata whose var_names are xbMagGiga1.1 IDs (old IDs mapping to one gene only)."""
    cw = convert(list(adata.var_names))
    cw = cw.groupby("old_id").filter(lambda d: d["new_id"].nunique() == 1).drop_duplicates("old_id")
    m = dict(zip(cw["old_id"], cw["new_id"]))
    keep = [v in m for v in adata.var_names] if drop_unmapped else [True] * adata.n_vars
    a = adata[:, keep].copy()
    a.var["old_id"] = a.var_names
    a.var_names = [m.get(v, v) for v in a.var_names]
    a.var_names_make_unique()
    return a
