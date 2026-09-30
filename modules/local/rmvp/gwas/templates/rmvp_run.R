#!/usr/bin/env Rscript
# Thin wrapper around rMVP::MVP for Nextflow RMVP_GWAS
suppressPackageStartupMessages({
  library(rMVP)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list(
    vcf = NULL,
    phenotype = NULL,
    out = "rmvp",
    method = "FarmCPU",
    ncpus = max(1L, parallel::detectCores() - 1L)
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--vcf", "-g")) {
      out$vcf <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--phenotype", "-p")) {
      out$phenotype <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--out", "-o")) {
      out$out <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--method", "-m")) {
      out$method <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--ncpus", "-n")) {
      out$ncpus <- as.integer(args[[i + 1L]]); i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

opt <- parse_args(args)
if (is.null(opt$vcf) || is.null(opt$phenotype)) {
  stop("Usage: rmvp_run.R --vcf file.vcf --phenotype pheno.txt --out prefix [--method FarmCPU|MLM|GLM] [--ncpus N]")
}

methods <- strsplit(opt$method, ",", fixed = TRUE)[[1]]
methods <- trimws(methods)

# Prepare MVP binary files from VCF (creates .geno.map / .geno.desc etc. under out prefix dir)
mvp_prefix <- paste0(opt$out, ".mvp")
MVP.Data(
  fileVCF = opt$vcf,
  filePhe = opt$phenotype,
  fileKin = TRUE,
  filePC = TRUE,
  out = mvp_prefix,
  ncpus = opt$ncpus,
  maxLine = 10000
)

phe <- read.table(paste0(mvp_prefix, ".phe"), header = TRUE, stringsAsFactors = FALSE)
geno <- attach.big.matrix(paste0(mvp_prefix, ".geno.desc"))
map <- read.table(paste0(mvp_prefix, ".geno.map"), header = TRUE, stringsAsFactors = FALSE)
kinship <- attach.big.matrix(paste0(mvp_prefix, ".kin.desc"))
pcs <- read.table(paste0(mvp_prefix, ".pc.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

MVP(
  phe = phe,
  geno = geno,
  map = map,
  K = kinship,
  CV.GLM = pcs,
  CV.MLM = pcs,
  CV.FarmCPU = pcs,
  nPC.GLM = min(3L, ncol(pcs)),
  nPC.MLM = min(3L, ncol(pcs)),
  nPC.FarmCPU = min(3L, ncol(pcs)),
  priority = "speed",
  ncpus = opt$ncpus,
  vc.method = "BRENT",
  maxLine = 10000,
  method = methods,
  file.output = c("pmap", "pmap.signal", "plot"),
  memo = opt$out
)

message("rMVP finished for methods: ", paste(methods, collapse = ","))
