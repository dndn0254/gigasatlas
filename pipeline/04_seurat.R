# Step 5-6: quality control, doublets, Seurat analysis and atlas names for one oyster library.
# usage: Rscript 04_seurat.R <STARsolo filtered folder> <gigasatlas model> <out.rds> [path to R/gigasatlas.R]
#   folder: out/sample/Solo.out/Gene/filtered (cells) or Solo.out/GeneFull/filtered (nuclei)
#   model : lineage, development, hemocyte, gonad, ganglia or mantle
# needs: R >= 4.3, Seurat >= 5; optional scDblFinder (Bioconductor) for doublets.
suppressPackageStartupMessages(library(Seurat))
a <- commandArgs(trailingOnly = TRUE)
IN <- a[1]; MODEL <- if (length(a) > 1) a[2] else "lineage"; OUT <- if (length(a) > 2) a[3] else "named.rds"
GA <- if (length(a) > 3) a[4] else "https://raw.githubusercontent.com/dndn0254/gigasatlas/main/R/gigasatlas.R"
source(GA)
set.seed(0)

# --- counts (STARsolo writes the same three files as Cell Ranger, gzipped or not)
read_solo <- function(d) {
  f <- function(n) { p <- file.path(d, n); if (file.exists(p)) p else paste0(p, ".gz") }
  m <- Matrix::readMM(f("matrix.mtx"))
  rownames(m) <- read.delim(f("features.tsv"), header = FALSE)[, 1]    # column 1 = gene ID of the GTF (LOC numbers)
  colnames(m) <- read.delim(f("barcodes.tsv"), header = FALSE)[, 1]
  methods::as(m, "CsparseMatrix")
}
m <- read_solo(IN)
seu <- CreateSeuratObject(m, min.cells = 0, min.features = 0)
cat("cells called by STARsolo:", ncol(seu), "\n")

# --- quality control: drop cells below the library median - 3 MAD of log UMIs or log genes (lower tail only)
low <- function(x) { lx <- log1p(x); median(lx) - 3 * mad(lx) }
keep <- log1p(seu$nCount_RNA) >= low(seu$nCount_RNA) & log1p(seu$nFeature_RNA) >= low(seu$nFeature_RNA)
seu <- subset(seu, cells = colnames(seu)[keep])
cat("after QC:", ncol(seu), "\n")

# --- doublets (the atlas removed a cell when 2 of 3 tools agreed: Scrublet, scDblFinder, SOLO; here scDblFinder only)
if (requireNamespace("scDblFinder", quietly = TRUE)) {
  sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(seu))
  seu$doublet <- sce$scDblFinder.class
  seu <- subset(seu, doublet == "singlet"); cat("after doublets:", ncol(seu), "\n")
} else message("scDblFinder not installed: BiocManager::install('scDblFinder') to remove doublets")

# --- Seurat
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:50, k.param = 20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:50, verbose = FALSE)

# Several libraries of the SAME tissue: merge first, then integrate with RPCA instead of the block above
#   seu <- merge(seu1, list(seu2, seu3)); seu <- NormalizeData(seu); seu <- FindVariableFeatures(seu, nfeatures = 3000)
#   seu <- ScaleData(seu); seu <- RunPCA(seu, npcs = 50)        # layers are already split by library after merge
#   seu <- IntegrateLayers(seu, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca")
#   seu <- JoinLayers(seu); seu <- FindNeighbors(seu, reduction = "integrated.rpca", dims = 1:50)
#   seu <- FindClusters(seu, resolution = 1); seu <- RunUMAP(seu, reduction = "integrated.rpca", dims = 1:50)

# --- names from the atlas
seu <- gigasatlas_annotate(seu, model = MODEL, group.by = "seurat_clusters")
print(head(sort(table(seu$gigas_label_group), decreasing = TRUE), 15))

# --- markers of each cluster (check the names by eye: read the product names of the top genes)
mk <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, max.cells.per.ident = 1500, verbose = FALSE)
write.table(mk, sub("\\.rds$", "_markers.tsv", OUT), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(seu, OUT)
cat("saved", OUT, "\n")
# plots: DimPlot(seu, group.by = "gigas_label_group", label = TRUE); FeaturePlot(seu, "gigas_prob")
