#!/usr/bin/env Rscript

# STEP 1: CREATE MATRIX OF 30K INDIVIDUALS AND 20K SNPS

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','QC_DIR'))])
library("crayon",lib.loc=R_LIB)
library("cli",lib.loc=R_LIB)
library("dplyr",lib.loc=R_LIB)

set.seed(1)
n=300000  # 30k sample size

# binonmial sampling of genotype for all indivudals based on allele freqency
get_genotypes <- function(snp_freq) {
  allele1 <- rbinom(n=n, size=1, prob=snp_freq)
  allele2 <- rbinom(n=n, size=1, prob=snp_freq)
  genotypes <- allele1 + allele2
  return(genotypes)
}

# obtain allele frequencies for 20k SNPs
setwd(QC_DIR)
gwas_snps <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; gwas_snps <- gwas_snps$x

setwd(GWAS_DIR)
genotype_matrix_k <- NULL
snp_list <- c(5000,10000,15000,20000)
# create matrix of 30k individuals and 20k genotypes
# split into 4 matrices of 5k genotypes
for (j in snp_list) {
for (i in (j-4999):(j)) {
  snp <- get_genotypes(gwas_snps[i])    # get allele count for all individuals
  genotype_matrix_k <- suppressMessages(bind_cols(genotype_matrix_k, snp))
}
save(genotype_matrix_k, file=paste0("simulation_matrix_k_",j/1000,"k.RData"))
}
