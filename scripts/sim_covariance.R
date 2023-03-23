#!/usr/bin/env Rscript

# load libraries
source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','QC_DIR','SIM_COV'))])
library(optparse, lib.loc=R_LIB)
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)
library(MASS, lib.loc=R_LIB)

# parameters
option_list = list(
  make_option(c("-i","--snps"), type="integer", default=1000,help="i: num of causal SNPs",metavar="character"),
  make_option(c("-e","--environment"), type="numeric", default=1.2,help="e: environment",metavar="character"),
  make_option(c("-a","--amplification"), type="numeric", default=14,help="a: amplification weight",metavar="character"),  
  make_option(c("-s","--seed"), type="numeric", default=1,help="s: seed",metavar="character"),
  make_option(c("-m","--matrix"), type="character", default="1000",help="m: matrix",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (nchar(opt$matrix) != 4) {
  print_help(opt_parser)
  stop("-m flag must be 4 digits", call.=FALSE)
}
snp_num <- opt$snps; print(snp_num)
E_ratio <- opt$environment; print(E_ratio)
m <- opt$matrix; print(m)
a_weight <- (opt$amplification)/100; print(a_weight)
s <- opt$seed; print(s)

# create sample of individuals -- random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

# set seed
set.seed(s)

# genotype matrix for 300k
setwd(QC_DIR)
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric')
snp_freqs <- snp_freqs$x

setwd(GWAS_DIR)
load("simulation_matrix_k_5k.RData")

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]

# Beta (effect size)
print("create matrix of Betas sampled from prespecified covariance matrices")
### PRE-SPECIFIED COV MATRICES ###
sigma_1 <- matrix(c(1,1,1,1), 2, 2)
sigma_2 <- matrix(c(substr(m,1,1), substr(m,2,2), substr(m,3,3), substr(m,4,4)), 2, 2)
class(sigma_2) <- "numeric"
###
mu <- c(0,0)
Beta <- mvrnorm((snp_num*(1-a_weight)), mu=mu, Sigma=sigma_1)
Beta <- rbind(Beta, mvrnorm((snp_num*a_weight), mu=mu, Sigma=sigma_2))

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G_m <- (Beta[,2]^2) * 2 * snp_freqs * (1-snp_freqs)   # genetic variance
  var_E_m <- (var_G_m * (1-h2)) / h2          # environmental variance (males)
  var_E_f <- var_E_m * E_ratio              # environmental variance for female based on proportion
  E <- NULL
  for (i in 1:length(sex)) {
    if (sex[i] == 0) {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_f))
    } else {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_m))
    }
  }
  return(E)
}

# environment and heritability combination
print("get vector of environmental effect")
E <- get_environment(0.5, E_ratio)
rm(snp_freqs)
print("get phenotype vector")
genotype_matrix_i <- t(genotype_matrix_i)
genotype_matrix_i[,m_id] <- genotype_matrix_i[,m_id]*Beta[,2]
genotype_matrix_i[,f_id] <- genotype_matrix_i[,f_id]*Beta[,1]
pheno <- colSums(genotype_matrix_i) + E
print(length(pheno))
rm(genotype_matrix_i)

# GWAS
GWAS <- function(genotype, pheno, E) {
  f_model <- lm(pheno[f_id] ~ genotype[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ genotype[m_id] + E[m_id])
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

# k snps, LD blocks
print("perform sex-specific GWAS")
mash_BETA <- NULL
mash_SE <- NULL
for (j in c(1,10,15,20)) {
  if (j != 1){
    load(paste0("simulation_matrix_k_",j,"k.RData"))
  }
  for (i in 1:length(colnames(genotype_matrix_k))) {
    gwas <- GWAS(genotype_matrix_k[,i], pheno, E)
    mash_BETA <- rbind(mash_BETA, data.frame(female= gwas[1,1], male= gwas[2,1]))
    mash_SE <- rbind(mash_SE, data.frame(female= gwas[1,2], male= gwas[2,2]))
  }
  print(j)
  rm(genotype_matrix_k)
}
print("write to file")

## CHANGE FILE NAME IF NEEDED
if (a_weight == 0) {
    file_name <- paste0("mash_",a_weight,"_",s,".RData")
} else {
    file_name <- paste0("mash_",m,"_",s,".RData")
}

setwd(SIM_COV)
save(mash_BETA, mash_SE, file=file_name)
print(paste0("File, ", file_name, " , in directory: ", SIM_COV))
