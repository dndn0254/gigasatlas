# From FASTQ to named oyster cells: the Gigas Cell Atlas pipeline

This folder shows, step by step, how to process a new Pacific oyster (*Crassostrea gigas*) single-cell or
single-nucleus RNA-seq library the same way as every library in the Gigas Cell Atlas.
**Cell Ranger is the default counter. STARsolo is an option**, set up to match Cell Ranger, for data that Cell Ranger
cannot read, for computers without Cell Ranger, or as a second count to compare with.
You need: the FASTQ files of one library, a Linux computer or cluster, and R.

| Step | What it does | Program (version used in the atlas) | Script |
|---|---|---|---|
| 1-2 | Download the current genome and gene annotation (GTF), build the reference (once) | NCBI RefSeq xbMagGiga1.1, annotation RS_2024_06; `cellranger mkref` (STAR index optional) | `01_get_reference.sh` |
| 3 | Get the reads (skip if you have your own FASTQ) | SRA Toolkit | `02_get_reads.sh` |
| 4 | Align, correct barcodes and UMIs, count genes, call cells | **Cell Ranger 7.0.1** | `03_cellranger.sh` |
| 4, option | the same with STARsolo, Cell Ranger-matched settings | STARsolo 2.7.11b | `03b_starsolo.sh` |
| 4, non-10x | droplet data with non-10x barcodes (for example BMKMANU) | Cell Ranger after barcode conversion, and STARsolo | `03c_non10x.sh`, `bmk_to_10x.py` |
| 5-6 | Quality control, doublets, normalisation, clusters, map, markers, atlas names | R 4.5, Seurat 5.5, scDblFinder, gigasatlas | `04_seurat.R` |
| 7 | Compare two counts of the same library | R, Seurat | `05_compare_counters.R` |

## Step by step

**1-2. Reference (once).** `bash 01_get_reference.sh ref/ 16` downloads the genome FASTA and GTF from NCBI
(GCF_963853765.1) and builds the Cell Ranger reference (33,068 genes: every exon that belongs to a transcript;
no gene-type filter). Add `star` (`bash 01_get_reference.sh ref/ 16 star`) for the STARsolo option and for non-10x data.
The mitochondrial genome is not in the RefSeq FASTA, so mitochondrial reads are not counted.

**3. Reads.** `bash 02_get_reads.sh SRR22554543 fastq/` (public data), or use your own FASTQ.
Cell Ranger wants the 10x file names, so link (do not copy) the reads:
`fq/<sample>_S1_L001_R1_001.fastq.gz` (barcode + UMI) and `fq/<sample>_S1_L001_R2_001.fastq.gz` (cDNA);
more runs or lanes of the same library become L002, L003, ...

**4. Counting with Cell Ranger (default).** `bash 03_cellranger.sh fq/ E1 ref/ out/ cells`
- `cells`: whole cells, exons only (`--include-introns=false`), as for the atlas embryo and larva libraries;
- `nuclei`: single-nucleus data and tissue cells with many unspliced reads (exons + introns), as for the atlas
  hemocyte, gonad and ganglia libraries;
- `--chemistry=auto` detects 10x v2, v3 and v3.1; `--no-bam` saves disk.

**4, option. STARsolo.** `bash 03b_starsolo.sh R1.fastq.gz R2.fastq.gz ref/ out/sample/ v3` runs STARsolo with
settings that reproduce Cell Ranger (10x whitelist with one mismatch, Cell Ranger UMI deduplication and
multi-gene filtering, EmptyDrops cell calling). Use the `Gene` output for whole cells and `GeneFull` for nuclei.
In the atlas (20 libraries) Cell Ranger kept 203,071 cells and STARsolo 186,414 after the same quality control;
176,851 cells were kept by both. On these shared cells the Seurat clusters of the two counts agreed with an adjusted Rand
index of 0.76 to 0.95 within tissues (0.61 for all tissues together), 86 of the 97 cross-dataset cell-type matches
(MetaNeighbor) were found again, and the vesicular-cell score of 64 shared cell types correlated at r = 0.995.

**4, non-10x droplet data.** `bash 03c_non10x.sh R1.fastq.gz R2.fastq.gz ref/ out/mantle/` runs both routes and
compares them. Cell Ranger alone stops with "an extremely low rate of correct barcodes" (0.9% for the atlas mantle
library, BMKMANU DG1000). Route A gives each barcode seen in the library its own 10x barcode, so Cell Ranger can run
(the same idea as the platform's own converter and as UniverSC); route C gives STARsolo the library's own barcode list.
In the atlas (mantle, 290.8 M reads), 99.2% of reads kept a valid barcode after conversion, Cell Ranger called all
28,913 mantle cells of the STARsolo atlas plus about 18,600 shallower barcodes (after quality control and doublet
removal 40,402 cells remain, 25,299 of them atlas cells; check doublets carefully, 15% here), and
the two routes agreed on the same cells (UMI rank correlation 0.999, pseudobulk Pearson 0.983).

**5-6. In R.** `Rscript 04_seurat.R out/E1/outs/filtered_feature_bc_matrix development E1_named.rds`

- quality control: remove cells whose log UMI count or log gene count is below the library median minus 3 MAD
  (no fixed cut-offs, because oyster tissues and protocols differ a lot);
- doublets: scDblFinder (the atlas removed a cell when 2 of 3 tools agreed: Scrublet, scDblFinder, SOLO);
- Seurat: NormalizeData (LogNormalize), FindVariableFeatures (vst, 3,000), ScaleData, RunPCA (50),
  FindNeighbors (k 20), FindClusters (Louvain, resolution 1), RunUMAP, FindAllMarkers;
- names: `gigasatlas_annotate(seu, model = "development", group.by = "seurat_clusters")` adds `gigas_label`,
  `gigas_prob`, `gigas_lineage`, `gigas_status` and the majority name of each cluster. Models: `lineage`,
  `development`, `hemocyte`, `gonad`, `ganglia`, `mantle`. For a tissue that is not in the atlas use
  `gigasatlas_score_clusters(seu)` and check the best marker sets by hand.

**7. Two counters, one library.** `Rscript 05_compare_counters.R out/E1/outs/filtered_feature_bc_matrix out/E1_star/Solo.out/Gene/filtered E1_cmp CellRanger STARsolo`
reports shared cells (Jaccard), UMI ratio and rank correlation per shared cell, pseudobulk correlation, the genes
that differ most, and cluster agreement (ARI) after the same Seurat steps. Use it when you had to count a library
with STARsolo (non-10x data, no Cell Ranger) and want to know how close it is to the Cell Ranger result.

**Several libraries of the same tissue**: merge them, then integrate with Seurat RPCA
(see the commented block in `04_seurat.R`). Do not integrate different tissues: in public oyster data each tissue
comes from a different lab, so correction would also remove real tissue differences.

**Old gene IDs** (Roslin v1, `G12345`): `seu <- gigasatlas_convert_ids(seu)` before naming, or paste them into the
ID converter of the website, which also shows what each gene is and where it is expressed.

## Other platforms

| Platform | Default | Option |
|---|---|---|
| 10x Chromium 3' v2, v3, v3.1 | Cell Ranger (`03_cellranger.sh`) | STARsolo, Cell Ranger-matched (`03b_starsolo.sh`) |
| droplet platforms with a 10x-like R1 but other barcodes (for example BMKMANU DG1000) | Cell Ranger after barcode conversion (`03c_non10x.sh`, route A) | STARsolo with the library's own barcode list (route C) |
| Parse Evercode / SPLiT-seq (combinatorial barcodes) | the company pipeline (split-pipe) | STARsolo `--soloType CB_UMI_Complex` with one whitelist per barcode round |
| BD Rhapsody | the BD Rhapsody pipeline | STARsolo `CB_UMI_Complex` with the three BD cell-label lists |
| Smart-seq2 / plate data (one cell per FASTQ) | STARsolo `--soloType SmartSeq` | featureCounts per cell |
| many platforms through Cell Ranger | UniverSC (wraps Cell Ranger 3.0.2) | |

The atlas tested only the first two rows. Whatever counter you use, steps 5-7 are the same, and the atlas names
work on any count matrix with xbMagGiga1.1 gene IDs.

## Notes
- Counting is the heavy step (Cell Ranger: about 1-3 h per library with 16 threads and 54 GB; STARsolo is faster);
  everything in R runs on a laptop for one library.
- The settings above are the ones used for all 20 public libraries of the atlas (8 studies). The atlas scripts
  ran on a PBS cluster; these files are simplified, single-library versions of the same steps.
