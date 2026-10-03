#!/usr/bin/env Rscript
# mashr multi-condition shrinkage (Urbut / Stephens mash).
# Long-format: snp, condition|tissue|study, beta, se.

suppressPackageStartupMessages({
  library(mashr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1]]
}

in_path <- get_arg("--effects")
out_prefix <- get_arg("--out", default = "mashr")

if (is.null(in_path) || !file.exists(in_path)) {
  stop("--effects required and must exist")
}

sep <- if (grepl("\\.csv$", in_path, ignore.case = TRUE)) "," else "\t"
df <- utils::read.table(in_path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
colnames(df) <- tolower(gsub("[^a-zA-Z0-9_]", "_", colnames(df)))

id_col <- intersect(c("snp", "variant_id", "rsid", "id"), colnames(df))
cond_col <- intersect(c("condition", "tissue", "study", "trait", "cohort"), colnames(df))
if (length(id_col) == 0 || length(cond_col) == 0) {
  stop("effects need snp and condition/tissue/study")
}
if (!all(c("beta", "se") %in% colnames(df))) {
  stop("effects need beta and se")
}

df$snp <- as.character(df[[id_col[[1]]]])
df$condition <- as.character(df[[cond_col[[1]]]])
df$beta <- as.numeric(df$beta)
df$se <- as.numeric(df$se)
df <- df[!is.na(df$snp) & !is.na(df$condition) & !is.na(df$beta) & !is.na(df$se) & df$se > 0, ]
if (nrow(df) < 2) {
  stop("need at least two effect rows")
}

snps <- unique(df$snp)
conds <- unique(df$condition)
Bhat <- matrix(NA_real_, nrow = length(snps), ncol = length(conds), dimnames = list(snps, conds))
Shat <- matrix(NA_real_, nrow = length(snps), ncol = length(conds), dimnames = list(snps, conds))
for (i in seq_len(nrow(df))) {
  Bhat[df$snp[[i]], df$condition[[i]]] <- df$beta[[i]]
  Shat[df$snp[[i]], df$condition[[i]]] <- df$se[[i]]
}

ok <- rowSums(!is.na(Bhat) & !is.na(Shat)) == length(conds)
Bhat <- Bhat[ok, , drop = FALSE]
Shat <- Shat[ok, , drop = FALSE]
if (nrow(Bhat) < 2) {
  stop("need >=2 SNPs observed in every condition")
}

data <- mash_set_data(Bhat, Shat)
U.c <- cov_canonical(data)
m <- mash(data, U.c)

pm <- get_pm(m)
psd <- get_psd(m)
lfsr <- get_lfsr(m)

long <- data.frame(
  snp = rep(rownames(pm), times = ncol(pm)),
  condition = rep(colnames(pm), each = nrow(pm)),
  posterior_mean = as.vector(pm),
  posterior_sd = as.vector(psd),
  lfsr = as.vector(lfsr),
  stringsAsFactors = FALSE
)

logf <- file(paste0(out_prefix, ".mashr.log"), open = "wt")
writeLines(c(
  paste("nsnps", nrow(Bhat)),
  paste("nconditions", ncol(Bhat)),
  paste("conditions", paste(colnames(Bhat), collapse = ","))
), logf)
close(logf)

utils::write.table(long, paste0(out_prefix, ".mash.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(lfsr, paste0(out_prefix, ".lfsr.tsv"), sep = "\t", quote = FALSE, col.names = NA)
