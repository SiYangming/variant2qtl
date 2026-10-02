#!/usr/bin/env Rscript
# Minimal SuSiE fine-mapping from summary statistics (susie_rss).
# Expected sumstats columns (case-insensitive): variant_id|snp|rsid, beta, se
# Optional: z, n. Optional LD matrix (--ld) as square numeric matrix without header.

suppressPackageStartupMessages({
  library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}

has_flag <- function(flag) !is.na(match(flag, args))

sumstats_path <- get_arg("--sumstats")
ld_path       <- get_arg("--ld", default = NULL)
out_prefix    <- get_arg("--out", default = "susie")
L             <- as.integer(get_arg("--L", default = "10"))
coverage      <- as.numeric(get_arg("--coverage", default = "0.95"))

if (is.null(sumstats_path) || !file.exists(sumstats_path)) {
  stop("--sumstats required and must exist")
}

sep <- if (grepl("\\.csv$", sumstats_path, ignore.case = TRUE)) "," else "\t"
df  <- utils::read.table(sumstats_path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
colnames(df) <- tolower(colnames(df))

id_col <- intersect(c("variant_id", "snp", "rsid", "id", "snpid"), colnames(df))
if (length(id_col) == 0) {
  df$variant_id <- paste0("var", seq_len(nrow(df)))
  id_col <- "variant_id"
} else {
  id_col <- id_col[[1]]
}

if ("z" %in% colnames(df)) {
  z <- as.numeric(df[["z"]])
} else if (all(c("beta", "se") %in% colnames(df))) {
  z <- as.numeric(df[["beta"]]) / as.numeric(df[["se"]])
} else {
  stop("sumstats need columns z, or beta+se")
}

n_var <- length(z)
if (!is.null(ld_path) && nzchar(ld_path) && file.exists(ld_path)) {
  R <- as.matrix(utils::read.table(ld_path, header = FALSE, sep = "", stringsAsFactors = FALSE))
  if (nrow(R) != n_var || ncol(R) != n_var) {
    stop(sprintf("LD matrix dim %dx%d != n_variants %d", nrow(R), ncol(R), n_var))
  }
} else {
  # Identity LD: valid for smoke / when LD unavailable; document limitation.
  R <- diag(n_var)
}

fit <- susie_rss(z = z, R = R, L = L, coverage = coverage)

pip <- susie_get_pip(fit)
cs  <- susie_get_cs(fit, Xcorr = R, coverage = coverage)

pip_df <- data.frame(
  variant_id = df[[id_col]],
  z = z,
  pip = as.numeric(pip),
  stringsAsFactors = FALSE
)
utils::write.table(pip_df, file = paste0(out_prefix, ".pip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cs_rows <- list()
if (!is.null(cs$cs) && length(cs$cs) > 0) {
  for (nm in names(cs$cs)) {
    idx <- cs$cs[[nm]]
    cs_rows[[length(cs_rows) + 1]] <- data.frame(
      credible_set = nm,
      variant_id = df[[id_col]][idx],
      pip = as.numeric(pip[idx]),
      stringsAsFactors = FALSE
    )
  }
  cs_df <- do.call(rbind, cs_rows)
} else {
  cs_df <- data.frame(credible_set = character(), variant_id = character(), pip = numeric(), stringsAsFactors = FALSE)
}
utils::write.table(cs_df, file = paste0(out_prefix, ".cs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(paste0(out_prefix, ".susie.log"))
cat("susieR version:", as.character(packageVersion("susieR")), "\n")
cat("n_variants:", n_var, "\n")
cat("L:", L, "\n")
cat("coverage:", coverage, "\n")
cat("ld:", if (!is.null(ld_path) && file.exists(ld_path)) ld_path else "identity", "\n")
cat("n_cs:", if (is.null(cs$cs)) 0 else length(cs$cs), "\n")
sink()
