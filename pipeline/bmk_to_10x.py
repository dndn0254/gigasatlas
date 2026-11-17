"""bmk_to_10x.py: let Cell Ranger read a non-10x droplet library whose R1 is 16 nt barcode + 12 nt UMI (10x v3 layout)
but whose barcodes are not on the 10x list (for example BMKMANU DG1000). Same idea as BMK's own converter and UniverSC.

  map : take the barcodes seen in the library (STARsolo run without a whitelist, raw output), keep those with at
        least MIN_UMI UMIs, and give each one its own 10x v3 barcode (one to one, most UMIs first) -> map.tsv
  conv: rewrite R1: barcode -> its 10x partner (exact, or the single candidate one base away), otherwise N x 16
        (Cell Ranger drops it). UMI and R2 are not changed.
  back: after Cell Ranger, write barcodes.tsv.gz with the original barcodes again (to match the authors' cell names)

usage: python3 bmk_to_10x.py map  <Solo.out/GeneFull/raw> <3M-february-2018.txt> <map.tsv> [MIN_UMI=10]
       python3 bmk_to_10x.py conv <map.tsv> <in_R1.fastq.gz> <out_R1.fastq.gz>
       python3 bmk_to_10x.py back <map.tsv> <outs/filtered_feature_bc_matrix/barcodes.tsv.gz> <out_barcodes.tsv.gz>
needs: numpy, pandas, scipy
"""
import sys, gzip, json
import numpy as np, pandas as pd, scipy.io

def make_map(raw, wl, out, min_umi=10):
    bc = pd.read_csv(f"{raw}/barcodes.tsv", header=None)[0].values
    u = np.asarray(scipy.io.mmread(f"{raw}/matrix.mtx").tocsc().sum(0)).ravel()
    keep = np.where(u >= min_umi)[0]; keep = keep[np.argsort(-u[keep])]
    tenx = []
    with open(wl) as f:
        for line in f:
            tenx.append(line.strip())
            if len(tenx) >= len(keep): break
    pd.DataFrame({"bmk": bc[keep], "tenx": tenx, "umi": u[keep].astype(int)}).to_csv(out, sep="\t", index=False)
    print(json.dumps({"barcodes_seen": int(len(bc)), "min_umi": min_umi, "mapped": int(len(keep))}))

def conv(mp, fin, fout):
    m = pd.read_csv(mp, sep="\t"); D = dict(zip(m.bmk, m.tenx)); memo = {}; N16 = "N" * 16
    st = dict(reads=0, exact=0, one_mismatch=0, no_match=0)
    def fix(b):
        if b in memo: return memo[b]
        hits = {D[x] for i in range(16) for c in "ACGT" if c != b[i] for x in [b[:i] + c + b[i + 1:]] if x in D}
        r = next(iter(hits)) if len(hits) == 1 else N16
        if len(memo) < 20_000_000: memo[b] = r
        return r
    with gzip.open(fin, "rt") as fi, gzip.open(fout, "wt", compresslevel=3) as fo:
        for h in fi:
            s = next(fi); p = next(fi); q = next(fi); st["reads"] += 1; b = s[:16]
            t = D.get(b)
            if t is not None: st["exact"] += 1
            else:
                t = fix(b); st["one_mismatch" if t != N16 else "no_match"] += 1
            fo.write(h + t + s[16:28] + "\n" + p + q[:28] + "\n")
    print(json.dumps(st))

def back(mp, fin, fout):
    m = pd.read_csv(mp, sep="\t"); B = dict(zip(m.tenx, m.bmk))
    bc = pd.read_csv(fin, header=None)[0].str.replace(r"-1$", "", regex=True).map(B)
    assert bc.notna().all()
    with gzip.open(fout, "wt") as g: g.write("\n".join(bc + "-1") + "\n")

if __name__ == "__main__":
    k = sys.argv[1]
    if k == "map": make_map(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 10)
    elif k == "conv": conv(*sys.argv[2:5])
    else: back(*sys.argv[2:5])
