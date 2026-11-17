# Step 7 (optional): compare two counts of the SAME library, for example Cell Ranger (03_cellranger.sh) and STARsolo
# (03b_starsolo.sh), or Cell Ranger after barcode conversion and STARsolo for a non-10x library (03c_non10x.sh).
# usage: Rscript 05_compare_counters.R <filtered folder A> <filtered folder B> <out_prefix> [name A] [name B]
# writes <out_prefix>_summary.tsv and <out_prefix>_genes.tsv, and prints the summary.
# What it reports, and why:
#   cells A, cells B, shared, Jaccard   do both counters call the same barcodes as cells?
#   UMI per shared cell: median ratio B/A and Spearman   same depth and same ranking of cells?
#   pseudobulk Pearson (log CPM of the summed shared cells)   same gene expression profile?
#   per-gene ratio (genes.tsv)   which genes one counter counts much more (intron reads, multi-mapping)?
#   cluster agreement (ARI) of shared cells after the same Seurat steps on both   same cell types?
suppressPackageStartupMessages(library(Seurat))
a <- commandArgs(trailingOnly = TRUE)
A <- a[1]; B <- a[2]; OUT <- a[3]; NA_ <- if (length(a) > 3) a[4] else "A"; NB_ <- if (length(a) > 4) a[5] else "B"
rd <- function(d) {
  f <- function(n) { p <- file.path(d, n); if (file.exists(p)) p else paste0(p, ".gz") }
  m <- Matrix::readMM(f("matrix.mtx"))
  rownames(m) <- read.delim(f("features.tsv"), header = FALSE)[, 1]
  colnames(m) <- sub("-1$", "", read.delim(f("barcodes.tsv"), header = FALSE)[, 1])   # Cell Ranger adds -1
  methods::as(m, "CsparseMatrix")
}
ma <- rd(A); mb <- rd(B)
g <- intersect(rownames(ma), rownames(mb)); sh <- intersect(colnames(ma), colnames(mb))
ua <- Matrix::colSums(ma[, sh]); ub <- Matrix::colSums(mb[, sh])
pa <- Matrix::rowSums(ma[g, sh]); pb <- Matrix::rowSums(mb[g, sh])
lc <- function(x) log1p(x / sum(x) * 1e6)
clus <- function(m) {
  s <- CreateSeuratObject(m); s <- NormalizeData(s, verbose = FALSE); s <- FindVariableFeatures(s, nfeatures = 3000, verbose = FALSE)
  s <- ScaleData(s, verbose = FALSE); s <- RunPCA(s, npcs = 50, verbose = FALSE)
  s <- FindNeighbors(s, dims = 1:50, verbose = FALSE); s <- FindClusters(s, resolution = 1, verbose = FALSE); s$seurat_clusters
}
ari <- function(x, y) {   # adjusted Rand index, no extra package
  t <- table(x, y); n <- sum(t); s2 <- function(v) sum(choose(v, 2))
  e <- s2(rowSums(t)) * s2(colSums(t)) / choose(n, 2); (s2(t) - e) / ((s2(rowSums(t)) + s2(colSums(t))) / 2 - e)
}
set.seed(0); ca <- clus(ma[g, sh]); cb <- clus(mb[g, sh])
res <- data.frame(item = c(paste("cells", NA_), paste("cells", NB_), "shared cells", "Jaccard of called cells",
                           "median UMI ratio B/A (shared)", "Spearman UMI (shared)", "pseudobulk Pearson (log CPM)",
                           "cluster ARI (shared cells, same Seurat steps)", "genes compared"),
                  value = c(ncol(ma), ncol(mb), length(sh), length(sh) / length(union(colnames(ma), colnames(mb))),
                            median((ub + 1) / (ua + 1)), cor(ua, ub, method = "spearman"), cor(lc(pa), lc(pb)),
                            ari(ca, cb), length(g)))
write.table(res, paste0(OUT, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
gt <- data.frame(gene = g, umi_A = pa, umi_B = pb, log2_ratio_B_A = log2((pb + 1) / (pa + 1)))
write.table(gt[order(-abs(gt$log2_ratio_B_A)), ], paste0(OUT, "_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
