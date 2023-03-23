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
  make_option(c("-p","--pheno"), type="character", default="height",help="phenotype",metavar="character"),
  make_option(c("-i","--snps"), type="integer", default=1000,help="i: num of causal SNPs",metavar="character"),
  make_option(c("-e","--environment"), type="numeric", default=1.2,help="e: environment",metavar="character"),
  make_option(c("-s","--seed"), type="numeric", default=1,help="s: seed",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno_name <- opt$pheno; print(pheno_name)
snp_num <- opt$snps; print(snp_num)
E_ratio <- opt$environment; print(E_ratio)
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
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; snp_freqs <- snp_freqs$x

setwd(GWAS_DIR)
load("simulation_matrix_k_5k.RData")

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]

# Beta (effect size)
print("create matrix of Betas sampled from prespecified covariance matrices")

### PRE-SPECIFIED COV MATRICES ###
if (pheno == "height") {
    a <- c(0.184, 0.124, 0.07, 0.042, 0.027, 0.017, 0.015, 0.007)
    null <- 0.439
    sigma_1 <- matrix(c(1, 1, 1, 1), 2, 2)
    sigma_2 <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
    sigma_3 <- matrix(c(1, 1.5, 1.5, 2.25), 2, 2)
    sigma_4 <- matrix(c(1, 1.125, 1.125, 2.25), 2, 2)
    sigma_5 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    sigma_6 <- matrix(c(2.25, 1.5, 1.5, 1), 2, 2)
    sigma_7 <- matrix(c(1, 0.75, 0.75, 2.25), 2, 2)
    sigma_8 <- matrix(c(4,2,2,1), 2, 2)
} else if (pheno == "bmi") {
    a <- c(0.129, 0.111, 0.108, 0.076, 0.027, 0.019, 0.008, 0.007)
    null <- 0.446
    sigma_1 <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
    sigma_2 <- matrix(c(1, 1, 1, 1), 2, 2)
    sigma_3 <- matrix(c(2.25, 1.5, 1.5 ,1), 2, 2)
    sigma_4 <- matrix(c(2.25, 1.125, 1.125, 1), 2, 2)
    sigma_5 <- matrix(c(2.25, 0.75, 0.75, 1), 2, 2)
    sigma_6 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    sigma_7 <- matrix(c(1, 1.5, 1.5, 2.25), 2, 2)
    sigma_8 <- matrix(c(9, 3, 3, 1), 2, 2)
} else if (pheno == "systolicBP_auto") {
    a <- c(0.097, 0.085, 0.04, 0.025, 0.024, 0.015, 0.011, 0.008)
    null <- 0.612
    sigma_1 <- matrix(c(2.25, 1.5, 1.5, 1), 2, 2)
    sigma_2 <- matrix(c(1, 1, 1, 1), 2, 2)
    sigma_3 <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
    sigma_4 <- matrix(c(2.25, 1.125, 1.125, 1), 2, 2)
    sigma_5 <- matrix(c(4, 2, 2, 1), 2, 2)
    sigma_6 <- matrix(c(9, 3, 3, 1), 2, 2)
    sigma_7 <- matrix(c(1, 1.5, 1.5, 2.25), 2, 2)
    sigma_8 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
} else if (pheno == "creatinine") {
    a <- c(0.078, 0.035, 0.031, 0.014, 0.013, 0.01, 0.01, 0.007)
    null <- 0.721
    sigma_1 <- matrix(c(1, 1, 1, 1), 2, 2)
    sigma_2 <- matrix(c(1, 1.5, 1.5, 2.25), 2, 2)
    sigma_3 <- matrix(c(2.25, 1.5, 1.5, 1), 2, 2)
    sigma_4 <- matrix(c(4, 2, 2, 1), 2, 2)
    sigma_5 <- matrix(c(1, 1.125, 1.125, 2.25), 2, 2)
    sigma_6 <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
    sigma_7 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    sigma_8 <- matrix(c(1, 2, 2, 4), 2, 2)
} else if (pheno == "IGF1") {
    a <- c(0.083, 0.046, 0.038, 0.037, 0.027, 0.024, 0.017, 0.014)
    null <- 0.592
    sigma_1 <- matrix(c(1, 1, 1, 1), 2, 2)
    sigma_2 <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
    sigma_3 <- matrix(c(2.25, 1.5, 1.5, 1), 2, 2)
    sigma_4 <- matrix(c(1, 1.5, 1.5, 2.25), 2, 2)
    sigma_5 <- matrix(c(1, 1.125, 1.25, 2.25), 2, 2)
    sigma_6 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
    sigma_7 <- matrix(c(1, 0.75, 0.75, 2.25), 2, 2)
    sigma_8 <- matrix(c(9, 3, 3, 1), 2, 2)
}
null_sigma <- matrix(c(0,0,0,0),2,2)
null_weight <- (1-null) / sum(a)
a = round(a * null_weight,3)
mu <- c(0,0)
Beta <- mvrnorm((snp_num*(a[1])), mu=mu, Sigma=sigma_1)
Beta <- rbind( rbind( rbind( rbind( rbind( rbind( rbind( rbind( Beta, mvrnorm((snp_num*a[2]), mu=mu, Sigma=sigma_2)) ,
                    mvrnorm((snp_num*a[3]), mu=mu, Sigma=sigma_3) ),
                    mvrnorm((snp_num*a[4]), mu=mu, Sigma=sigma_4) ),
                    mvrnorm((snp_num*a[5]), mu=mu, Sigma=sigma_5) ),
                    mvrnorm((snp_num*a[6]), mu=mu, Sigma=sigma_6) ),
                    mvrnorm((snp_num*a[7]), mu=mu, Sigma=sigma_7) ),
                    mvrnorm((snp_num*a[8]), mu=mu, Sigma=sigma_8) ),
                    mvrnorm((snp_num*null), mu=mu, Sigma=null_sigma) )

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
E <- get_environment(0.05, E_ratio)
rm(snp_freqs)
print("get phenotype vector")
genotype_matrix_i <- t(genotype_matrix_i)
genotype_matrix_i[,m_id] <- genotype_matrix_i[,m_id]*Beta[,2]
genotype_matrix_i[,f_id] <- genotype_matrix_i[,f_id]*Beta[,1]
pheno <- colSums(genotype_matrix_i) + E
# standardize by units of std dev
pheno[m_id] <- pheno[m_id] / sd(pheno[m_id])
pheno[f_id] <- pheno[f_id]  / sd(pheno[f_id])
print(length(pheno))
rm(genotype_matrix_i)

# GWAS ## TODo: more snps per gwas
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
setwd(SIM_COV)
file_name <- paste0("mash_", pheno_name, "_", s, ".RData")
save(mash_BETA, mash_SE, file=file_name)
print(paste0("File, ", file_name, " , in directory: ", SIM_COV))