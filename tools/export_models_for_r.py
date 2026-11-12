"""export_models_for_r.py — gigasatlas 모델 (.npz) 을 R 에서 읽는 표로 내보낸다 (2026-11-12)

왜: Seurat (R) 사용자가 파이썬 없이 R 안에서 바로 이름을 붙이게 (R/gigasatlas.R). 모델은 선형 (multinomial logistic regression) 이라
    계수, 절편, 평균, 표준편차만 있으면 R 에서 같은 확률을 계산할 수 있다. 값은 .npz 와 같다 (float32 → 텍스트 9 자리).
사용: python tools/export_models_for_r.py   (패키지 폴더에서)
출력: R/models/<model>_genes.tsv.gz (gene, mean, scale, 이름별 계수 열), R/models/<model>_labels.tsv (label, intercept, lineage)
"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from gigasatlas.model import load_model, available_models
OUT = os.path.join(HERE, "R", "models"); os.makedirs(OUT, exist_ok=True)
for name in available_models():
    m = load_model(name)
    G = pd.DataFrame(m.coef.T.astype(np.float64), columns=[f"coef_{i}" for i in range(len(m.labels))])
    G.insert(0, "scale", m.scale.astype(np.float64)); G.insert(0, "mean", m.mean.astype(np.float64)); G.insert(0, "gene", m.genes.astype(str))
    G.to_csv(os.path.join(OUT, f"{name}_genes.tsv.gz"), sep="\t", index=False, float_format="%.9g", compression="gzip")
    pd.DataFrame({"label": m.labels, "intercept": m.intercept.astype(np.float64), "lineage": [m.lineage_of(x) for x in m.labels]}).to_csv(
        os.path.join(OUT, f"{name}_labels.tsv"), sep="\t", index=False, float_format="%.9g")
    print(name, len(m.genes), "genes", len(m.labels), "labels")
