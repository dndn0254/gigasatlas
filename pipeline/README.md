# From FASTQ to named oyster cells: the Gigas Cell Atlas pipeline

This folder shows, step by step, how to process a new Pacific oyster (*Crassostrea gigas*) single-cell or
single-nucleus RNA-seq library the same way as every library in the Gigas Cell Atlas.
You need: the FASTQ files of one 10x library, a Linux computer or cluster, and R.

| Step | What it does | Program (version used in the atlas) | Script |
|---|---|---|---|
| 1 | Download the current oyster genome and gene annotation (GTF) | NCBI RefSeq xbMagGiga1.1, annotation RS_2024_06 | `01_get_reference.sh` |
| 2 | Build the alignment index (once) | STAR 2.7.11b | `01_get_reference.sh` |
| 3 | Get the reads (skip if you have your own FASTQ) | SRA Toolkit | `02_get_reads.sh` |
| 4 | Align, correct barcodes and UMIs, count genes, call cells | STARsolo 2.7.11b, Cell Ranger-matched settings | `03_starsolo.sh` |
| 5 | Quality control, doublets, normalisation, clusters, map, markers | R 4.5, Seurat 5.5, scDblFinder | `04_seurat.R` |
| 6 | Name the cells with the atlas | gigasatlas (R version, this repository) | `04_seurat.R` |

## Step by step

**1-2. Reference (once).** `bash 01_get_reference.sh ref/` downloads the genome FASTA and GTF from NCBI
(GCF_963853765.1) and builds a STAR index. The mitochondrial genome is not in the RefSeq FASTA; mitochondrial
reads are therefore not counted (the atlas measured them separately and did not filter on them).

**3. Reads.** `bash 02_get_reads.sh SRR22554543 fastq/` (public data), or use your own `*_R1_001.fastq.gz`
(barcode + UMI) and `*_R2_001.fastq.gz` (cDNA).

**4. Counting.** `bash 03_starsolo.sh fastq/SRR_R1.fastq.gz fastq/SRR_R2.fastq.gz ref/ out/sample/ v3` runs STARsolo with the settings that
reproduce Cell Ranger (barcode whitelist with one mismatch, UMI deduplication and multi-gene filtering as in
Cell Ranger, EmptyDrops cell calling). Use `v2` for 10x v2 chemistry (10 nt UMI) and `v3` for v3 or v3.1 (12 nt UMI).
Single-nucleus data: use the `GeneFull` output (exons + introns); whole cells: `Gene`.

**5-6. In R.** `Rscript 04_seurat.R out/sample/Solo.out/Gene/filtered gonad sample_named.rds`

- quality control: remove cells whose log UMI count or log gene count is below the library median minus 3 MAD
  (no fixed cut-offs, because oyster tissues and protocols differ a lot);
- doublets: scDblFinder (the atlas removed a cell when 2 of 3 tools agreed: Scrublet, scDblFinder, SOLO);
- Seurat: NormalizeData (LogNormalize), FindVariableFeatures (vst, 3,000), ScaleData, RunPCA (50),
  FindNeighbors (k 20), FindClusters (Louvain, resolution 1), RunUMAP, FindAllMarkers;
- names: `gigasatlas_annotate(seu, model = "gonad", group.by = "seurat_clusters")` adds `gigas_label`,
  `gigas_prob`, `gigas_lineage`, `gigas_status` and the majority name of each cluster. Models: `lineage`,
  `development`, `hemocyte`, `gonad`, `ganglia`, `mantle`. For a tissue that is not in the atlas use
  `gigasatlas_score_clusters(seu)` and check the best marker sets by hand.

**Several libraries of the same tissue**: merge them, then integrate with Seurat RPCA
(see the commented block in `04_seurat.R`). Do not integrate different tissues: in public oyster data each tissue
comes from a different lab, so correction would also remove real tissue differences.

**Old gene IDs** (Roslin v1, `G12345`): `seu <- gigasatlas_convert_ids(seu)` before naming.

## Notes
- Counting is the heavy step (about 20-40 min per library with 8-16 threads); everything in R runs on a laptop
  for one library.
- The settings above are the ones used for all 20 public libraries of the atlas (8 studies). The atlas scripts
  ran on a PBS cluster; these files are simplified, single-library versions of the same steps.
