# gigasatlas

Name the cells of new Pacific oyster (*Magallana gigas*, syn. *Crassostrea gigas*) single-cell data in one command,
using the Gigas Cell Atlas: 186,414 cells from 8 public studies, re-processed from raw reads on the current genome
(xbMagGiga1.1, RefSeq annotation RS_2024_06). Website: https://dndn0254.github.io/gigas-atlas/

| What you have | Command | What you get |
|---|---|---|
| a Seurat object or a count matrix from your single-cell experiment | `gigasatlas_annotate()` in R, `gigasatlas annotate` in Python | a cell-type name, probability and lineage for every cell |
| a gene list, for example DE genes from bulk RNA-seq | `gigasatlas enrich` | which atlas cell types the list points to |
| gene IDs from an older paper (Roslin v1) | `gigasatlas convert` | current xbMagGiga1.1 IDs |

## 0. Only have FASTQ files? Follow the pipeline

[`pipeline/`](pipeline/) takes one 10x library from raw reads to named cells, with the same settings as every
library of the Gigas Cell Atlas:

| Step | Program | Script |
|---|---|---|
| Reference genome + gene annotation (GTF), STAR index | NCBI xbMagGiga1.1 (RS_2024_06), STAR 2.7.11b | `pipeline/01_get_reference.sh` |
| Reads from SRA (or your own FASTQ) | SRA Toolkit | `pipeline/02_get_reads.sh` |
| Align, count, call cells (Cell Ranger-matched) | STARsolo 2.7.11b | `pipeline/03_starsolo.sh` |
| QC (median - 3 MAD), doublets, Seurat clusters, UMAP, markers | R, Seurat 5, scDblFinder | `pipeline/04_seurat.R` |
| Cell names from the atlas | gigasatlas (R) | `pipeline/04_seurat.R` |

```bash
bash pipeline/01_get_reference.sh ref/                                   # once
bash pipeline/03_starsolo.sh R1.fastq.gz R2.fastq.gz ref/ out/sample/ v3
Rscript pipeline/04_seurat.R out/sample/Solo.out/Gene/filtered gonad sample_named.rds
```

## 1. Seurat (R) users

No installation and no Python: `R/gigasatlas.R` runs the models directly on a Seurat object
(Seurat >= 5; the input is Seurat's `NormalizeData` LogNormalize layer, which is exactly what the models expect).
Gene names must be xbMagGiga1.1 RefSeq IDs (LOC... or symbols); for Roslin v1 IDs run `gigasatlas_convert_ids()` first.

```r
source("https://raw.githubusercontent.com/dndn0254/gigasatlas/main/R/gigasatlas.R")   # or source("R/gigasatlas.R") in a clone
seu <- gigasatlas_annotate(seu, model = "gonad", group.by = "seurat_clusters")
head(seu[[c("gigas_label", "gigas_prob", "gigas_lineage", "gigas_status", "gigas_label_group")]])
DimPlot(seu, group.by = "gigas_label", label = TRUE)

# a tissue that is not in the atlas (gill, digestive gland): score every atlas marker set per cluster
z <- gigasatlas_score_clusters(seu, group.by = "seurat_clusters"); colnames(z)[max.col(z)]
# old Roslin v1 gene IDs (G12345) -> xbMagGiga1.1
seu <- gigasatlas_convert_ids(seu)
```

The R and Python versions give the same names (tested on 3,815 hemocytes: 100% identical names, probabilities
within 3e-7). Model coefficients for R are plain tables in `R/models/` (made by `tools/export_models_for_r.py`).

| R function | Purpose |
|---|---|
| `gigasatlas_annotate(seu, model, group.by = NULL, min_prob = 0.5)` | name every cell; adds `gigas_label`, `gigas_prob`, `gigas_lineage`, `gigas_status` (and `gigas_label_group` per cluster) |
| `gigasatlas_score_clusters(seu, group.by)` | marker-set score per cluster, for a tissue not in the atlas |
| `gigasatlas_convert_ids(seu)` | Roslin v1 -> xbMagGiga1.1 gene IDs |
| `gigasatlas_markers(top_n = 50)` | marker table of every atlas cell type |
| `gigasatlas_models()`, `gigasatlas_load(model)` | bundled models |

## 2. Python: install


You need Python 3.9 or newer. In a terminal (Linux, macOS, or Windows PowerShell):

```bash
python -m venv gigas                 # a clean environment (optional but recommended)
source gigas/bin/activate            # Windows: gigas\Scripts\activate
pip install git+https://github.com/dndn0254/gigasatlas.git
```

This also installs numpy, pandas, scipy and anndata. In a conda environment the same `pip install` line works.
Check it: `gigasatlas models` should list six models.

## 3. Python: try it in 10 seconds

The package ships a small example (400 gastrula cells, a 200-gene list, some old gene IDs), so you can see
every tool work before using your own data:

```bash
gigasatlas demo
```

Real output (Linux, Python 3.11; takes about 10 s):

```text
400 cells, 19,395 genes; 92% of the 'development' model genes found
== 1. Name cells: 400 gastrula cells, model 'development'
   100% confident; same name as the authors for 92% of cells
   (optimistic: these cells were in the training data; held-out accuracy is 0.91-0.92)
                           author_label gigas_label  gigas_prob
E4_AGACCATAGCAAGTCG                Un-M        Un-M    0.999812
E4_TGCGGCAAGATTGCGG                Un-M        Un-M    0.999806
E4_ACCAACATCCACCTGT              Ec-SF1      Ec-SF1    0.999965
E4_TTACTGTTCACACGAT  Unnamed by authors    EnM-Mus2    0.974146
E4_CTGCCATCACATAACC             EnM-Un1     EnM-Un1    0.998173

== 2. Which cell types does a gene list point to? 200 genes up in sorted granulocytes
                                                                      cell_type  overlap  set_size            q
                                     hemocyte | Hemocyte — vesicular cells (VC)       19       100 5.345133e-21
                      mantle | Hemocyte-like (VC-like: Rab-20/NPFF receptor 2+)       15       100 5.254497e-15
                                      priming | Hemocyte — vesicular cells (VC)       14       100 1.079704e-13
ganglia | Hemocyte-like (VC-like: HAVCR1+/NPFF receptor 2+) [author: Hemocytes]       12       100 5.922497e-11
                      gonad | Hemocyte-like (VC-like: NPFF receptor 2+/HAVCR1+)       11       100 1.112499e-09

== 3. Convert old Roslin v1 gene IDs
      old_id       new_id          category
      G11998 LOC105343850           RBH 1:1
       G5032 LOC105337090           RBH 1:1
      G29424 LOC105342567           RBH 1:1
LOC105348812 LOC105348812         Identical
LOC117688139 LOC105333716            Merged
LOC105348378 LOC136269620             Split
LOC105348378 LOC105348378             Other
LOC117690076          NaN Previous-unmapped

Now try your own data: gigasatlas annotate my_data.h5ad --model hemocyte --table names.tsv
```

To run the commands on the example files yourself:

```bash
gigasatlas example --out ex          # copies gastrula_400.h5ad, granulocyte_up_genes.txt, roslin_ids.txt into ./ex
gigasatlas annotate ex/gastrula_400.h5ad --model development --table names.tsv
gigasatlas enrich ex/granulocyte_up_genes.txt
gigasatlas convert ex/roslin_ids.txt
```

## 4. Python: name the cells of your own data

```bash
gigasatlas annotate MY_INPUT --model hemocyte --table names.tsv --out named.h5ad
```

**Input (`MY_INPUT`)**, any of:

| Format | Example |
|---|---|
| AnnData `.h5ad` (Python; Seurat users: use R/gigasatlas.R on the Seurat object) | `sample.h5ad` |
| a 10x folder: `matrix.mtx(.gz)`, `features.tsv(.gz)` or `genes.tsv`, `barcodes.tsv(.gz)` (Cell Ranger, STARsolo) | `outs/filtered_feature_bc_matrix` |
| a count table `.csv` / `.tsv` (genes as rows, cells as columns, or the other way round) | `counts.csv` |

Two things must be true:

1. **Raw counts** (integers). Log-normalised data also work, but counts are safer.
2. **Gene names are xbMagGiga1.1 IDs** (`LOC105347058`, or gene symbols as in the NCBI RefSeq annotation),
   which is what you get when you align to the current genome GCF_963853765.1. The command prints how many
   of the model's genes it found; above 80% is good. Below 50% usually means old Roslin v1 IDs (see section 5).

**Choose the model** closest to your sample (`gigasatlas models`):

| Model | Cell types | For |
|---|---:|---|
| `development` | 73 | cleavage, blastula and gastrula embryos, trochophore larvae (authors' names) |
| `hemocyte` | 15 | circulating blood cells (hemolymph) |
| `gonad` | 28 | female, male, diploid and triploid gonads (nuclei) |
| `ganglia` | 49 | cerebral and visceral ganglia (nuclei) |
| `mantle` | 21 | mantle |
| `lineage` | 7 | broad lineage only (hemocyte, neural, epithelial, muscle, germline, ...) |

**Options**: `--table names.tsv` writes one row per cell; `--out named.h5ad` writes your data with the names added;
`--groupby cluster` (a cluster column in your `.h5ad`) also gives the majority name of each cluster; `--min-prob 0.5` sets the
probability below which a cell is called `uncertain`.

**Output columns** (in the table and in `.obs`):

| Column | Meaning |
|---|---|
| `gigas_label` | the most likely atlas cell type |
| `gigas_prob` | its probability (0-1) |
| `gigas_lineage` | lineage of that cell type |
| `gigas_status` | `confident` (prob >= 0.5) or `uncertain` |
| `gigas_label_group` | majority name of the cell's cluster (only with `--groupby`) |

Cells of a type the reference does not contain get some name with a low probability: look at `uncertain` cells
and at clusters with mixed names before trusting a label.

## 5. Read a gene list (bulk RNA-seq, proteomics, a GWAS gene set)

Put one gene per line in a text file and run:

```bash
gigasatlas enrich my_de_genes.txt --out enrich.tsv
```

Each atlas cell type's 100 most specific genes are tested for overlap with your list (hypergeometric test,
Benjamini-Hochberg q). Columns: `cell_type` (dataset | name), `overlap`, `set_size`, `fold`, `p`, `q`, `genes`.
Example: the 200 genes higher in sorted granulocytes point to vesicular cells (VC) in five datasets (demo above).
A significant cell type means your genes are typical of that cell type; it does not tell whether the cells became
more numerous or changed their expression.

## 6. Convert old gene IDs

```bash
gigasatlas convert old_ids.txt --out converted.tsv
```

Roslin v1 Ensembl IDs (`G11998`) are matched by reciprocal best protein hit (1:1 only). Roslin v1 RefSeq IDs
(`LOC...`) use NCBI's own comparison of the two annotations, which also says whether a gene was `Identical`,
`Merged`, `Split`, `Changed` or is absent from the new genome (`Previous-unmapped`).

For a whole data set whose genes are Roslin v1 IDs, convert inside Python before annotating:

```python
import gigasatlas as ga
A = ga.read("old_ids_data.h5ad")
A = ga.rename_var(A)              # keeps genes with one unambiguous new ID
ga.predict(A, "hemocyte")
```

## 7. Python functions

```python
import gigasatlas as ga

A = ga.read("sample.h5ad")                 # or a 10x folder, or counts.csv
ga.gene_coverage(A, "hemocyte")            # share of model genes found; expect > 0.8
probs = ga.predict(A, "hemocyte", groupby="cluster")   # adds gigas_* columns to A.obs; returns cells x types probabilities
A.obs[["gigas_label", "gigas_prob", "gigas_status"]].head()

ga.cell_type_enrichment(["LOC105347058", "LOC105332995"]).head()   # gene list -> cell types
ga.convert(["G11998", "LOC117688139"])                             # old -> new IDs
```

| Function | Purpose |
|---|---|
| `read(path)` | `.h5ad`, 10x folder or count table -> AnnData |
| `predict(adata, model, groupby=None, min_prob=0.5)` | name every cell |
| `gene_coverage(adata, model)` | share of the model's genes present |
| `cell_type_enrichment(genes, top_n=100)` | gene list -> cell types |
| `convert(ids)`, `rename_var(adata)` | Roslin v1 -> xbMagGiga1.1 |
| `score_clusters(adata, groupby)` | marker-set score per cluster, for a tissue not in the atlas |
| `available_models()`, `load_model(name)` | bundled models |
| `train_model(X_log, genes, labels, name)` | build your own reference (needs `pip install "gigasatlas[train] @ git+https://github.com/dndn0254/gigasatlas.git"`) |
| `example_data()`, `example_gene_list()`, `example_old_ids()` | the bundled examples |

## 8. How accurate is it? (cells never used for training)

| Test | Truth | Cells | Accuracy | Guessing the most common name |
|---|---|---:|---:|---:|
| Gastrula embryo, hold out 1 of 4 libraries (4 runs) | authors' names | 4,499-4,752 | 0.91-0.92 | 0.11-0.12 |
| Ganglia, cerebral -> visceral library | atlas names | 19,672 | 0.74 | 0.20 |
| Ganglia, visceral -> cerebral library | atlas names | 13,457 | 0.73 | 0.25 |
| Gonad, 2024 study -> 2025 study (diploid and triploid) | atlas names | 26,369 | 0.74 | 0.21 |
| Blood, healthy oyster -> primed hemocytes (another lab) | atlas names | 41,453 | 0.83 | 0.18 |

Cells called `confident` are 83-99% of cells and are right more often (0.77-0.93).
The example cells in `gigasatlas demo` were part of the training data, so their agreement with the authors is
higher than you should expect on new data.

Each model is a multinomial logistic regression on 3,000 variable genes (log1p counts per 10,000, standardised),
trained with at most 600 cells per name, the same idea as CellTypist (Dominguez Conde et al. 2022, Science).

## 9. Limits

- **A tissue that is not in the atlas** (gill, digestive gland, adductor muscle): the `lineage` model does not transfer well
  between tissues (held-out tissue accuracy 0.35-0.76), because every tissue in the atlas comes from a different lab and
  protocol. For a new tissue, cluster first, then use `gigasatlas_score_clusters` (R) or `score_clusters` (Python) and check the result by hand.
- Names that are rare in the reference (for example triploid-specific cells absent from the 2024 gonad study) are hard to predict.
- Whole cells and nuclei give different absolute levels; use the model built from the same kind of data when you can.
- `Unnamed by authors` is a group of embryo cells the original study did not name.

## Data and citation

Reference data: Gavery et al. 2025; Piovani et al. 2023; de La Forest Divonne et al. 2025; PRJNA1236035; Wang et al. 2024;
Wang et al. 2025; Yang et al. 2026; Bai et al. 2026 (preprint). Gene ID crosswalk: NCBI annotation comparison
(GCF_963853765.1-RS_2024_06) and reciprocal best protein hits to Roslin v1 Ensembl.
Example cells: Gavery et al. 2025 (PRJNA909306, gastrula library E4). Example gene list: sorted granulocytes vs hyalinocytes
(Mao et al. 2020), re-analysed on xbMagGiga1.1.

Manuscript in preparation; until then cite this repository. Web atlas, cell annotation tables and cell-type families:
https://dndn0254.github.io/gigas-atlas/ (source: https://github.com/dndn0254/gigas-atlas).

License: MIT.
