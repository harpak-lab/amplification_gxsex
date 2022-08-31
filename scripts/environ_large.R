#!/usr/bin/env Rscript

# STEP 2: CREATE MATRIX EFFECT ESTIMATES AND SE FOR MASH INPUT (10K snps)

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','QC_DIR'))])
library(optparse, lib.loc=R_LIB)
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)
setwd(QC_DIR)

# parameters
option_list = list(
  make_option(c("-g","--heritability"), type="numeric", default=100,help="g: h2",metavar="character"), 
  make_option(c("-e","--environment"), type="numeric", default=100,help="e: environment",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$heritability) | is.null(opt$environment)) {
  print_help(opt_parser)
  stop("Missing argument for heritability or environmental variance ratio", call.=FALSE)
}
h2 <- opt$heritability; print(h2)
E_ratio <- opt$environment; print(E_ratio)
snp_num <- 10000

# create sample of individuals with random sex
# 0=female  1=male
n=300000
sex <- rbinom(n,1,0.5)   
f_id <- which(sex==0) ; m_id <- which(sex==1)

# set seed
set.seed(1)
# genotype matrix for 300k individuals
snp_freqs <- read.csv("maf_sample_20k.txt", colClasses = 'numeric') ; snp_freqs <- snp_freqs$x    # load 20k total snp frequencies
setwd(GWAS_DIR)
load("simulation_matrix_k_5k.RData")  # load matrix of 300k by 5k snp genotypes 

# i causal SNPs
snp_freqs <- snp_freqs[1:snp_num]   # subset 10k causal snp frequencies

######################################################################################

# Beta (effect size)
Beta <- rnorm(snp_num, mean=0, sd=1)    # get vector of 10k Beta effect sizes

# get environmental effect
get_environment <- function(h2, E_ratio) {
  var_G <- (Beta^2) * 2 * snp_freqs * (1-snp_freqs)   # genetic variance
  var_E_m <- (var_G * (1-h2)) / h2          # environmental variance (males)
  var_E_f <- var_E_m * E_ratio              # environmental variance for female based on proportion
  E <- NULL
  for (i in 1:length(sex)) {    # assign environmental effect based on gender
    if (sex[i] == 0) {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_f))
    } else {
      E[i] <- rnorm(1, mean=0, sd=sqrt(var_E_m))
    }
  }
  return(E)
}

get_pheno_split <- function(addon) {
    pheno_sum <- rep(0,250)
    for (i in 1:(5000/250)) {
        start <- ((250*i) - 249)
        end <- (250 * i)
        pheno_sum <- pheno_sum + rowSums(t(t(genotype_matrix_k[start:end])*Beta[start:end+addon]))
    }
    return(pheno_sum)
}

# environment and heritability combination; get phenotype
E <- get_environment(h2, E_ratio)   # call function to get vector of 300k environmental effects
rm(snp_freqs)
pheno <- get_pheno_split(0) # vector of 5k phenotypes
rm(genotype_matrix_k)
load("simulation_matrix_k_10k.RData")
pheno <- pheno + get_pheno_split(5000) + E
rm(genotype_matrix_k)
load("simulation_matrix_k_5k.RData")
print(length(pheno))

# GWAS function for each SNP
GWAS <- function(genotype, pheno, E) {
  # linear regression of phenotype on genotype + environmental effect
  f_model <- lm(pheno[f_id] ~ genotype[f_id] + E[f_id])
  m_model <- lm(pheno[m_id] ~ genotype[m_id] + E[m_id])
  # grab Beta and SE coefficients for the SNP
  f_gwas <- summary(f_model)$coefficients[2,1:2]
  m_gwas <- summary(m_model)$coefficients[2,1:2]
  return(rbind(f_gwas,m_gwas))
}

# k snps, LD blocks
mash_BETA <- NULL
mash_SE <- NULL
for (j in c(1,10,15,20)) {    # for each of the premade genotype matrices
  if (j != 1){
    load(paste0("simulation_matrix_k_",j,"k.RData"))    # go through each SNP
  }
  for (i in 1:length(colnames(genotype_matrix_k))) {
    gwas <- GWAS(genotype_matrix_k[,i], pheno, E)     # call GWAS function
    mash_BETA <- rbind(mash_BETA, data.frame(female= gwas[1,1], male= gwas[2,1]))
    mash_SE <- rbind(mash_SE, data.frame(female= gwas[1,2], male= gwas[2,2]))
  }
  print(j)
  rm(genotype_matrix_k)
}
file_name <- paste0("mash_",snp_num,"_",h2,"_",E_ratio,".RData")
save(mash_BETA, mash_SE, file=file_name)
print(paste0("File, ", file_name, ", in directory: " GWAS_DIR))