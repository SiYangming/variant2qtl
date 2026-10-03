#!/usr/bin/env Rscript
# Build SuSiE sumstats + LD from the official susieR N3finemapping dataset
# (same source as susieR vignettes). Falls back to the vignette simulation recipe
# if N3finemapping is unavailable.
#
# Usage:
#   Rscript scripts/make_finemap_susie_official_testdata.R
#   conda run -n susier Rscript scripts/make_finemap_susie_official_testdata.R

suppressPackageStartupMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "assets/testdata/finemap_susie_official"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

X <- NULL
y <- NULL
source_label <- NULL

tryCatch({
  utils::data("N3finemapping", package = "susieR", envir = environment())
  if (exists("N3finemapping", inherits = FALSE)) {
    dat <- get("N3finemapping")
    if (is.list(dat) && !is.null(dat$X) && !is.null(dat$y)) {
      X <- as.matrix(dat$X)
      y <- as.numeric(dat$y)
      source_label <- "susieR::N3finemapping"
    }
  }
}, error = function(e) invisible(NULL))

if (is.null(X)) {
  # Official vignette recipe (susie_rss.html): n=200; use p=100 for a compact asset
  n <- 200
  p <- 100
  beta <- rep(0, p)
  beta[1:4] <- 1
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- drop(X %*% beta + rnorm(n))
  source_label <- "susieR vignette simulation (n=200,p=100,causal=1:4)"
}

n <- nrow(X)
p <- ncol(X)
X <- scale(X, center = TRUE, scale = TRUE)
y <- as.numeric(scale(y, center = TRUE, scale = FALSE))

beta_hat <- numeric(p)
se_hat <- numeric(p)
z <- numeric(p)
for (j in seq_len(p)) {
  fit <- stats::lm(y ~ X[, j])
  sm <- summary(fit)$coefficients
  beta_hat[[j]] <- sm[2, 1]
  se_hat[[j]] <- sm[2, 2]
  z[[j]] <- sm[2, 3]
}

R <- stats::cor(X)
R[R > 1] <- 1
R[R < -1] <- -1
diag(R) <- 1

ids <- colnames(X)
if (is.null(ids) || any(!nzchar(ids))) {
  ids <- paste0("var", seq_len(p))
}

sumstats <- data.frame(
  variant_id = ids,
  beta = beta_hat,
  se = se_hat,
  z = z,
  n = n,
  stringsAsFactors = FALSE
)

utils::write.table(
  sumstats,
  file = file.path(out_dir, "sumstats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  R,
  file = file.path(out_dir, "ld.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

writeLines(
  c(
    paste0("source: ", source_label),
    paste0("n_samples: ", n),
    paste0("n_variants: ", p),
    paste0("generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  con = file.path(out_dir, "SOURCE.txt")
)

message("Wrote ", file.path(out_dir, "sumstats.tsv"), " and ld.txt (", p, " variants; ", source_label, ")")
