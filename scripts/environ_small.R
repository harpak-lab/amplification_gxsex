#!/usr/bin/env Rscript

# STEP 2: CREATE MATRIX EFFECT ESTIMATES AND SE FOR MASH INPUT (100, 1K snps)

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','QC_DIR'))])
library(optparse, lib.loc=R_LIB)
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)
setwd(QC_DIR)

# parameters
option_list = list(
  make_option(c("-i","--snps"), type="integer", default=100,help="i: num of causal SNPs",metavar="character"), 
  make_option(c("-g","--heritability"), type="numeric", default=0.05,help="g: h2",metavar="character"), 
  make_option(c("-e","--environment"), type="numeric", default=1,help="e: environment",metavar="character") 
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$snps) | is.null(opt$heritability) | is.null(opt$environment)) {
  print_help(opt_parser)
  stop("Missing argument for snp #, heritability, or environmental variance ratio", call.=FALSE)
}
snp_num <- opt$snps; print(snp_num)
h2 <- opt$heritability; print(h2)
E_ratio <- opt$environment; print(E_ratio)

# create sample of individuals -- random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

# set seed
set.seed(1)
# genotype matrix for 300k
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; snp_freqs <- snp_freqs$x
setwd(GWAS_DIR)
load("simulation_matrix_k_5k.RData")

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]
genotype_matrix_i <- genotype_matrix_k[,(1:snp_num)]

# Beta (effect size)
Beta <- rnorm(snp_num, mean=0, sd=1)

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G <- (Beta^2) * 2 * snp_freqs * (1-snp_freqs)   # genetic variance
  var_E_m <- (var_G * (1-h2)) / h2          # environmental variance (males)
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
#h2 in c(0.5, 0.05)   E_ratio in c(1, 1.5, 5)
E <- get_environment(h2, E_ratio)
rm(snp_freqs)
pheno <- rowSums(t(t(genotype_matrix_i)*Beta)) + E
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
file_name <- paste0("mash_",snp_num,"_",h2,"_",E_ratio,".RData")
save(mash_BETA, mash_SE, file=file_name)
print(paste0("File, ", file_name, ", in directory: " GWAS_DIR))