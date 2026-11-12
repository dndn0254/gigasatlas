"""Command line: gigasatlas demo | example | annotate | enrich | convert | models

    gigasatlas demo                                   run every tool on the bundled examples (about 10 s)
    gigasatlas example --out my_folder                copy the example files to try the commands yourself
    gigasatlas annotate INPUT --model hemocyte --table names.tsv [--out named.h5ad] [--groupby seurat_clusters]
    gigasatlas enrich genes.txt
    gigasatlas convert old_ids.txt
    gigasatlas models
INPUT may be .h5ad, a 10x folder (matrix.mtx + features.tsv + barcodes.tsv) or a .csv/.tsv count table.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

MODEL_HELP = {
    "development": "cleavage, blastula, gastrula embryos and trochophore larvae (authors' names), 73 types",
    "hemocyte": "circulating blood cells (hemolymph), 15 types",
    "gonad": "female, male, diploid and triploid gonads (nuclei), 28 types",
    "ganglia": "cerebral and visceral ganglia (nuclei), 49 types",
    "mantle": "mantle, 21 types",
    "lineage": "7 broad lineages only; for a tissue not in the atlas use score_clusters instead",
}


def _lines(path):
    return [x.strip() for x in open(path) if x.strip() and not x.startswith("#")]


def _annotate(inp, model, groupby=None, min_prob=0.5, out=None, table=None, quiet=False):
    from .annotate import gene_coverage, predict
    from .io import read
    A = read(inp) if isinstance(inp, str) else inp
    cov = gene_coverage(A, model)
    print(f"{A.n_obs:,} cells, {A.n_vars:,} genes; {cov:.0%} of the '{model}' model genes found", file=sys.stderr)
    if cov < 0.5:
        print("warning: fewer than half of the model genes found. Are the genes xbMagGiga1.1 IDs (LOC... or symbols)?\n"
              "         Roslin v1 IDs: convert first (python: A = gigasatlas.rename_var(A)).", file=sys.stderr)
    predict(A, model, groupby=groupby, min_prob=min_prob)
    cols = [c for c in A.obs.columns if c.startswith("gigas_")]
    if table:
        A.obs[cols].to_csv(table, sep="\t", index_label="cell"); print(f"wrote {table}", file=sys.stderr)
    if out:
        A.write_h5ad(out); print(f"wrote {out}", file=sys.stderr)
    if not quiet:
        vc = A.obs["gigas_label"].value_counts()
        print(f"\n{(A.obs['gigas_status'] == 'confident').mean():.0%} of cells confident (probability >= {min_prob})")
        print("most common names:"); print(vc.head(15).to_string())
    return A


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gigasatlas", description="Pacific oyster single-cell atlas tools. Start with: gigasatlas demo",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("\n\n", 1)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("demo", help="run every tool on the bundled examples")
    e0 = sub.add_parser("example", help="copy the example files into a folder"); e0.add_argument("--out", default="gigasatlas_example")
    a = sub.add_parser("annotate", help="name the cells of .h5ad / 10x folder / count table")
    a.add_argument("input"); a.add_argument("--model", default="lineage", choices=sorted(MODEL_HELP))
    a.add_argument("--table", help="write cell names as a .tsv table"); a.add_argument("--out", help="write an .h5ad with the names added")
    a.add_argument("--groupby", help="obs column with clusters: also gives one name per cluster"); a.add_argument("--min-prob", type=float, default=0.5)
    e = sub.add_parser("enrich", help="cell types a gene list points to (one gene per line)")
    e.add_argument("genes"); e.add_argument("--top", type=int, default=20); e.add_argument("--out", help="write the full table (.tsv)")
    c = sub.add_parser("convert", help="Roslin v1 IDs -> xbMagGiga1.1 (one ID per line)")
    c.add_argument("ids"); c.add_argument("--out", help="write the table (.tsv)")
    sub.add_parser("models", help="list bundled reference models")
    args = ap.parse_args(argv)

    if args.cmd == "models":
        from .model import available_models
        for m in available_models():
            print(f"{m:12s} {MODEL_HELP.get(m, '')}")
    elif args.cmd == "example":
        from .datasets import EX_DIR
        os.makedirs(args.out, exist_ok=True)
        for f in sorted(os.listdir(EX_DIR)):
            shutil.copy(os.path.join(EX_DIR, f), os.path.join(args.out, f)); print(os.path.join(args.out, f))
        print(f"\nnext:\n  gigasatlas annotate {args.out}/gastrula_400.h5ad --model development --table names.tsv\n"
              f"  gigasatlas enrich {args.out}/granulocyte_up_genes.txt\n  gigasatlas convert {args.out}/roslin_ids.txt")
    elif args.cmd == "annotate":
        _annotate(args.input, args.model, args.groupby, args.min_prob, args.out, args.table)
    elif args.cmd == "enrich":
        from .enrich import cell_type_enrichment
        r = cell_type_enrichment(_lines(args.genes))
        if args.out: r.to_csv(args.out, sep="\t", index=False)
        print(r.head(args.top)[["cell_type", "overlap", "set_size", "fold", "q"]].to_string(index=False))
    elif args.cmd == "convert":
        from .ids import convert
        r = convert(_lines(args.ids), keep_unmapped=True)
        if args.out: r.to_csv(args.out, sep="\t", index=False)
        print(r.to_csv(sep="\t", index=False), end="")
    elif args.cmd == "demo":
        from .datasets import example_data, example_gene_list, example_old_ids
        from .enrich import cell_type_enrichment
        from .ids import convert
        print("== 1. Name cells: 400 gastrula cells, model 'development'")
        A = _annotate(example_data(), "development", quiet=True)
        agree = (A.obs["gigas_label"].astype(str) == A.obs["author_label"].astype(str)).mean()
        print(f"   {(A.obs['gigas_status'] == 'confident').mean():.0%} confident; same name as the authors for {agree:.0%} of cells")
        print("   (optimistic: these cells were in the training data; held-out accuracy is 0.91-0.92)")
        print(A.obs[["author_label", "gigas_label", "gigas_prob"]].head(5).to_string())
        print("\n== 2. Which cell types does a gene list point to? 200 genes up in sorted granulocytes")
        print(cell_type_enrichment(example_gene_list()).head(5)[["cell_type", "overlap", "set_size", "q"]].to_string(index=False))
        print("\n== 3. Convert old Roslin v1 gene IDs")
        print(convert(example_old_ids(), keep_unmapped=True)[["old_id", "new_id", "category"]].to_string(index=False))
        print("\nNow try your own data: gigasatlas annotate my_data.h5ad --model hemocyte --table names.tsv")


if __name__ == "__main__":
    main()
