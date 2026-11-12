# gigasatlas.R: name the cells of a Seurat object with the Gigas Cell Atlas, in R only (no Python needed)
#
#   source("https://raw.githubusercontent.com/dndn0254/gigasatlas/main/R/gigasatlas.R")   # or source("R/gigasatlas.R") in a clone
#   seu <- gigasatlas_annotate(seu, model = "gonad", group.by = "seurat_clusters")
#   head(seu[[c("gigas_label", "gigas_prob", "gigas_status", "gigas_label_group")]])
#
# The models are the same as in the Python package (multinomial logistic regression on log1p(counts per 10,000),
# standardised with the training mean and SD, values clipped to +-10). Seurat's NormalizeData (LogNormalize, 1e4)
# gives exactly that input, so R and Python return the same names and probabilities.
# Needs: Seurat >= 5 (SeuratObject), Matrix. Gene names must be xbMagGiga1.1 RefSeq IDs (LOC... / symbols);
# use gigasatlas_convert_ids() first for Roslin v1 IDs (G12345).

.gigas_here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
.gigas_url <- "https://raw.githubusercontent.com/dndn0254/gigasatlas/main"

.gigas_read <- function(rel, ...) {
  loc <- if (!is.null(.gigas_here)) file.path(.gigas_here, "..", rel) else NULL
  src <- if (!is.null(loc) && file.exists(loc)) loc else paste0(.gigas_url, "/", rel)
  if (grepl("^https?://", src) && grepl("\\.gz$", src)) src <- gzcon(url(src))   # read.delim opens and closes it
  utils::read.delim(src, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

gigasatlas_models <- function() c("lineage", "development", "hemocyte", "gonad", "ganglia", "mantle")

gigasatlas_load <- function(model = "lineage") {
  stopifnot(model %in% gigasatlas_models())
  G <- .gigas_read(sprintf("R/models/%s_genes.tsv.gz", model)); L <- .gigas_read(sprintf("R/models/%s_labels.tsv", model), na.strings = character(0))
  L$label[L$label == "nan"] <- "Unnamed by authors"
  W <- as.matrix(G[, grep("^coef_", names(G))]); colnames(W) <- L$label
  list(name = model, genes = G$gene, mean = G$mean, scale = G$scale, W = W, b = L$intercept, labels = L$label, lineage = setNames(L$lineage, L$label))
}

.gigas_lognorm <- function(object, assay) {
  X <- tryCatch(SeuratObject::LayerData(object, assay = assay, layer = "data"), error = function(e) NULL)
  C <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  if (is.null(X) || length(X@x) == 0 || identical(X, C)) {     # not normalised yet: log1p(counts per 10,000)
    message("gigasatlas: normalising counts (LogNormalize, 1e4)")
    X <- Seurat::NormalizeData(C, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  }
  X
}

gigasatlas_annotate <- function(object, model = "lineage", assay = NULL, min_prob = 0.5, group.by = NULL, prefix = "gigas",
                                chunk = 5000, return_probs = FALSE) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  m <- if (is.list(model)) model else gigasatlas_load(model)
  X <- .gigas_lognorm(object, assay)
  feat <- rownames(X); idx <- match(m$genes, feat)
  miss <- is.na(idx); idx[miss] <- match(gsub("_", "-", m$genes[miss]), feat)      # Seurat turns '_' into '-'
  found <- !is.na(idx); cov <- mean(found)
  message(sprintf("gigasatlas: model %s, %d of %d model genes found (%.0f%%)", m$name, sum(found), length(found), 100 * cov))
  if (cov < 0.5) warning("fewer than half of the model genes were found: are the gene IDs xbMagGiga1.1? see gigasatlas_convert_ids()")
  n <- ncol(X); P <- matrix(0, n, length(m$labels), dimnames = list(colnames(X), m$labels))
  for (s in seq(1, n, by = chunk)) {
    j <- s:min(n, s + chunk - 1)
    Z <- matrix(0, length(j), length(m$genes))                                    # genes not measured count as 0 (as in Python)
    Z[, found] <- as.matrix(Matrix::t(X[idx[found], j, drop = FALSE]))
    Z <- sweep(sweep(Z, 2, m$mean, "-"), 2, m$scale, "/"); Z[Z > 10] <- 10; Z[Z < -10] <- -10
    z <- Z %*% m$W; z <- sweep(z, 2, m$b, "+"); z <- z - apply(z, 1, max); e <- exp(z)
    P[j, ] <- e / rowSums(e)
  }
  best <- max.col(P, ties.method = "first"); lab <- m$labels[best]; pr <- P[cbind(seq_len(n), best)]
  md <- data.frame(lab, pr, unname(m$lineage[lab]), ifelse(pr >= min_prob, "confident", "uncertain"), row.names = colnames(X))
  names(md) <- paste0(prefix, c("_label", "_prob", "_lineage", "_status"))
  if (!is.null(group.by)) {
    g <- as.character(object[[group.by]][colnames(X), 1])
    maj <- tapply(lab, g, function(x) names(sort(table(x), decreasing = TRUE))[1])
    md[[paste0(prefix, "_label_group")]] <- unname(maj[g])
  }
  object <- SeuratObject::AddMetaData(object, md)
  object@misc[[paste0(prefix, "_model")]] <- list(model = m$name, n_genes_model = length(m$genes), n_genes_found = sum(found), frac_confident = mean(pr >= min_prob))
  if (return_probs) object@misc[[paste0(prefix, "_probs")]] <- P
  message(sprintf("gigasatlas: %.0f%% of cells confident (prob >= %.2f)", 100 * mean(pr >= min_prob), min_prob))
  object
}

# Marker sets of every atlas cell type (top specific genes); also the input of gigasatlas_score_clusters
gigasatlas_markers <- function(top_n = 50) { mk <- .gigas_read("src/gigasatlas/data/label_markers.tsv"); mk[mk$rank <= top_n, ] }

# For a tissue that is not in the atlas (gill, digestive gland): score every atlas marker set in every cluster.
# Same as the Python score_clusters: mean log expression of the set minus the cluster's panel mean, z-scored across clusters.
gigasatlas_score_clusters <- function(object, group.by = "seurat_clusters", assay = NULL, top_n = 50, min_genes = 5) {
  assay <- assay %||% SeuratObject::DefaultAssay(object); X <- .gigas_lognorm(object, assay)
  g <- as.character(object[[group.by]][colnames(X), 1]); ug <- sort(unique(g))
  G <- Matrix::sparseMatrix(i = match(g, ug), j = seq_along(g), x = 1, dims = c(length(ug), length(g)))
  G <- Matrix::Diagonal(x = 1 / Matrix::rowSums(G)) %*% G
  M <- as.matrix(G %*% Matrix::t(X)); rownames(M) <- ug                              # clusters x genes mean
  bg <- rowMeans(M); mk <- gigasatlas_markers(top_n); feat <- colnames(M)
  S <- sapply(split(mk$gene, mk$set), function(gs) {
    ii <- match(gs, feat); ii[is.na(ii)] <- match(gsub("_", "-", gs[is.na(ii)]), feat); ii <- ii[!is.na(ii)]
    if (length(ii) >= min_genes) rowMeans(M[, ii, drop = FALSE]) - bg else rep(NA_real_, nrow(M))
  })
  S <- S[, colSums(is.na(S)) == 0, drop = FALSE]
  sdv <- apply(S, 2, function(x) sqrt(mean((x - mean(x))^2))); sdv[sdv == 0] <- 1
  sweep(sweep(S, 2, colMeans(S), "-"), 2, sdv, "/")
}

# Roslin v1 gene IDs (G12345 or old LOC) -> xbMagGiga1.1; keeps genes that map to exactly one new gene
gigasatlas_convert_ids <- function(object, assay = NULL) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  cw <- .gigas_read("src/gigasatlas/data/id_crosswalk.tsv.gz", colClasses = "character"); cw <- cw[cw$new_id != "", ]
  one <- tapply(cw$new_id, cw$old_id, function(x) length(unique(x)) == 1); cw <- cw[cw$old_id %in% names(one)[one] & !duplicated(cw$old_id), ]
  C <- SeuratObject::LayerData(object, assay = assay, layer = "counts"); keep <- rownames(C) %in% cw$old_id
  C <- C[keep, ]; rownames(C) <- make.unique(cw$new_id[match(rownames(C), cw$old_id)])
  message(sprintf("gigasatlas: %d of %d genes converted", sum(keep), length(keep)))
  Seurat::CreateSeuratObject(C, meta.data = object[[]], assay = assay)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
