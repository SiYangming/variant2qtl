#!/usr/bin/env Rscript
# Coloc ABF for one QTL locus vs one GWAS trait (Giambartolomei et al.).
# Expected columns (case-insensitive): snp|variant_id|rsid, beta, se
# Optional: maf|af, n|n_samples. type defaults to quant.

suppressPackageStartupMessages({
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    return(default)
  }
  args[[i + 1]]
}

qtl_path  <- get_arg("--qtl")
gwas_path <- get_arg("--gwas")
out_prefix <- get_arg("--out", default = "coloc")
type1 <- get_arg("--type1", default = "quant")
type2 <- get_arg("--type2", default = "quant")

if (is.null(qtl_path) || !file.exists(qtl_path)) {
  stop("--qtl required and must exist")
}
if (is.null(gwas_path) || !file.exists(gwas_path)) {
  stop("--gwas required and must exist")
}

read_ss <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  df <- utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- tolower(gsub("[^a-zA-Z0-9_]", "_", colnames(df)))
  id_col <- intersect(c("snp", "variant_id", "rsid", "id", "snpid", "rs"), colnames(df))
  if (length(id_col) == 0) {
    stop("sumstats need snp/variant_id/rsid")
  }
  df$snp <- as.character(df[[id_col[[1]]]])
  if (!all(c("beta", "se") %in% colnames(df))) {
    stop("sumstats need beta and se")
  }
  df$beta <- as.numeric(df$beta)
  df$se <- as.numeric(df$se)
  df$varbeta <- df$se * df$se
  maf_col <- intersect(c("maf", "af", "aaf", "eaf"), colnames(df))
  df$maf <- if (length(maf_col) > 0) as.numeric(df[[maf_col[[1]]]]) else 0.2
  n_col <- intersect(c("n", "n_samples", "sample_size"), colnames(df))
  df$n <- if (length(n_col) > 0) as.numeric(df[[n_col[[1]]]]) else 200
  df[!is.na(df$snp) & !is.na(df$beta) & !is.na(df$se) & df$se > 0, ]
}

qtl <- read_ss(qtl_path)
gwas <- read_ss(gwas_path)
merged <- merge(qtl, gwas, by = "snp", suffixes = c("_qtl", "_gwas"))
if (nrow(merged) < 2) {
  stop("need >=2 overlapping SNPs between QTL and GWAS sumstats")
}

ds1 <- list(
  snp = merged$snp,
  beta = merged$beta_qtl,
  varbeta = merged$varbeta_qtl,
  maf = merged$maf_qtl,
  N = median(merged$n_qtl, na.rm = TRUE),
  type = type1,
  sdY = 1
)
ds2 <- list(
  snp = merged$snp,
  beta = merged$beta_gwas,
  varbeta = merged$varbeta_gwas,
  maf = merged$maf_gwas,
  N = median(merged$n_gwas, na.rm = TRUE),
  type = type2,
  sdY = 1
)

res <- coloc.abf(dataset1 = ds1, dataset2 = ds2)
sm <- as.list(res$summary)
summary_df <- data.frame(
  nsnps = sm$nsnps,
  PP.H0.abf = sm[["PP.H0.abf"]],
  PP.H1.abf = sm[["PP.H1.abf"]],
  PP.H2.abf = sm[["PP.H2.abf"]],
  PP.H3.abf = sm[["PP.H3.abf"]],
  PP.H4.abf = sm[["PP.H4.abf"]],
  stringsAsFactors = FALSE
)
utils::write.table(summary_df, file = paste0(out_prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

results <- res$results
utils::write.table(results, file = paste0(out_prefix, ".abf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  sprintf("coloc.abf nsnps=%s PP.H4=%s", sm$nsnps, sm[["PP.H4.abf"]]),
  con = paste0(out_prefix, ".coloc.log")
)
