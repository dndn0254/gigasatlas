import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

import pytest

import gigasatlas as ga
from gigasatlas.model import train_model


def _toy(n=300, g=60, seed=0):
    rng = np.random.default_rng(seed)
    lab = np.repeat(["A", "B", "C"], n // 3)
    mu = np.ones((3, g)); mu[0, :10] = 20; mu[1, 10:20] = 20; mu[2, 20:30] = 20
    X = rng.poisson(mu[np.searchsorted(["A", "B", "C"], lab)])
    return ad.AnnData(X=sp.csr_matrix(X.astype(np.float32)), obs=pd.DataFrame({"truth": lab}, index=[f"c{i}" for i in range(n)]),
                      var=pd.DataFrame(index=[f"g{i}" for i in range(g)])), lab


def test_train_predict_roundtrip(tmp_path):
    pytest.importorskip("sklearn")          # training needs the optional [train] extra
    A, lab = _toy()
    X = ga.log_normalize(A.X)
    m = train_model(X, A.var_names, lab, name="toy", n_genes=40)
    f = tmp_path / "toy.npz"; m.save(str(f))
    m2 = ga.load_model(str(f))
    B, _ = _toy(seed=1)
    ga.predict(B, m2)
    assert (B.obs["gigas_label"].astype(str) == B.obs["truth"]).mean() > 0.95


def test_missing_genes_are_zero():
    pytest.importorskip("sklearn")          # training needs the optional [train] extra
    A, lab = _toy()
    m = train_model(ga.log_normalize(A.X), A.var_names, lab, name="toy", n_genes=40)
    B, _ = _toy(seed=2)
    B = B[:, 5:].copy()                      # drop some model genes
    assert ga.gene_coverage(B, m) < 1
    ga.predict(B, m)
    assert "gigas_prob" in B.obs


def test_bundled_models_load():
    for name in ga.available_models():
        m = ga.load_model(name)
        assert m.coef.shape == (len(m.labels), len(m.genes))


def test_read_10x_and_table(tmp_path):
    import gzip, scipy.io
    A, _ = _toy(n=30, g=12)
    d = tmp_path / "tenx"; d.mkdir()
    with gzip.open(d / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, sp.csr_matrix(A.X).T.tocoo())
    pd.DataFrame({0: A.var_names, 1: A.var_names, 2: "Gene Expression"}).to_csv(d / "features.tsv.gz", sep="\t", header=False, index=False)
    pd.Series(A.obs_names).to_csv(d / "barcodes.tsv.gz", sep="\t", header=False, index=False)
    B = ga.read(str(d))
    assert B.shape == A.shape and np.allclose(B.X.toarray(), A.X.toarray())
    t = tmp_path / "counts.csv"
    pd.DataFrame(A.X.toarray().T, index=A.var_names, columns=A.obs_names).to_csv(t)   # genes x cells
    C = ga.read(str(t))
    assert C.shape == A.shape


def test_examples_and_cli(tmp_path, capsys):
    from gigasatlas.cli import main
    E = ga.example_data()
    assert E.n_obs == 400 and "author_label" in E.obs
    assert ga.gene_coverage(E, "development") > 0.9
    assert len(ga.example_gene_list()) == 200 and len(ga.example_old_ids()) > 3
    main(["demo"])
    out = capsys.readouterr().out
    assert "Name cells" in out and "vesicular" in out.lower()
    main(["annotate", ga.example_path("gastrula_400.h5ad"), "--model", "development", "--table", str(tmp_path / "n.tsv")])
    t = pd.read_csv(tmp_path / "n.tsv", sep="\t")
    assert len(t) == 400 and {"gigas_label", "gigas_prob", "gigas_status"} <= set(t.columns)
