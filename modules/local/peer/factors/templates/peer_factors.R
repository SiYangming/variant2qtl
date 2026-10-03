#!/usr/bin/env Rscript
# PEER hidden factors from FastQTL BED (samples in columns 5+).
# Optional known covariates: sample×cov or cov×sample.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1]]
}

pheno_path <- get_arg("--phenotype")
cov_path <- get_arg("--covariates", default = "")
out_prefix <- get_arg("--out", default = "peer")
nk <- as.integer(get_arg("--nk", default = "10"))

if (is.null(pheno_path) || !file.exists(pheno_path)) {
  stop("--phenotype required and must exist")
}
if (!requireNamespace("peer", quietly = TRUE)) {
  stop("R package peer is required (bioconda::r-peer)")
}
library(peer)

open_tab <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "rt")
}

con <- open_tab(pheno_path)
header <- readLines(con, n = 1)
close(con)
hdr <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1]]
if (length(hdr) < 5) {
  stop("phenotype must be FastQTL BED (chr start end id + samples)")
}
samples <- hdr[seq(5, length(hdr))]

ph <- utils::read.delim(pheno_path, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
expr_cols <- colnames(ph)[seq(5, ncol(ph))]
expr <- t(as.matrix(ph[, expr_cols, drop = FALSE]))
storage.mode(expr) <- "double"
expr[!is.finite(expr)] <- 0

cov_mat <- NULL
cov_names <- character(0)
if (!is.null(cov_path) && nzchar(cov_path) && file.exists(cov_path)) {
  cv <- utils::read.delim(cov_path, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
  cn <- colnames(cv)
  sample_set <- samples
  if (all(cn[-1] %in% sample_set) || sum(cn[-1] %in% sample_set) > length(sample_set) / 2) {
    cov_names <- as.character(cv[[1]])
    cov_mat <- t(as.matrix(cv[, cn[-1], drop = FALSE]))
    cov_mat <- cov_mat[samples, , drop = FALSE]
  } else {
    row_ids <- as.character(cv[[1]])
    cov_names <- cn[-1]
    cov_mat <- as.matrix(cv[, cn[-1], drop = FALSE])
    rownames(cov_mat) <- row_ids
    cov_mat <- cov_mat[samples, , drop = FALSE]
  }
  storage.mode(cov_mat) <- "double"
}

nk <- max(1L, nk)
model <- PEER()
PEER_setPhenoMean(model, expr)
PEER_setNk(model, nk)
if (!is.null(cov_mat)) {
  PEER_setCovariates(model, cov_mat)
}
PEER_update(model)
X <- PEER_getX(model)

n_cov <- if (is.null(cov_mat)) 0L else ncol(cov_mat)
factor_names <- c(cov_names, paste0("PEER", seq_len(max(1L, ncol(X) - n_cov))))
if (length(factor_names) < ncol(X)) {
  factor_names <- c(factor_names, paste0("F", seq_len(ncol(X) - length(factor_names))))
}
factor_names <- factor_names[seq_len(ncol(X))]
colnames(X) <- factor_names
rownames(X) <- samples

tensor <- data.frame(id = samples, X, check.names = FALSE)
omiga <- data.frame(id = factor_names, t(X), check.names = FALSE)

utils::write.table(tensor, file = paste0(out_prefix, ".peer.tensor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(omiga, file = paste0(out_prefix, ".peer.omiga.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
