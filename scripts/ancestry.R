#!/usr/bin/env Rscript

library(optparse)
library(EthSEQ)

option_list <- list(
  make_option(c("-v", "--vcf"), type = "character", default = NULL, 
              help = "Path to the target VCF file", metavar = "character"),
  make_option(c("-m", "--model"), type = "character", default = NULL, 
              help = "Path to the reference GDS model", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./", 
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

ethseq.Analysis(
  target.vcf = opt$vcf,
  model.gds = opt$model,
  out.dir = opt$out,
  cores = 1,
  verbose = TRUE,
  composite.model.call.rate = 0.99,
  space = "3D"
)
