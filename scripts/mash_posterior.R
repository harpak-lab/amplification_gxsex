#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library(optparse, lib.loc=R_LIB)
library(crayon,lib.loc=R_LIB)
library(dplyr,lib.loc=R_LIB)
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-s","--set"), type="character", default="additive",help="input set # if using PGS pipeline",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
mode <- opt$set; print(mode)

# add PGS suffix if using PGS method 
if (mode == "additive") {
    suffix <- "" ; wd <- paste0(GWAS_DIR,"/",pheno,"/mash")

} else { 
    suffix <- "_pgs" ; wd <- paste0(GWAS_DIR,"/",pheno,"/PGS_",mode)
} 
setwd(wd)
load(file= paste0(pheno,"_mash",suffix,".RData"))      # setup file
load(file= paste0(pheno,"_mash_100g",suffix,".RData"))      # fitted file

# adjust table for mixture weights
weight_col <- function(df) {
    colnames(df) <- gsub("^(.*)[.].*", "\\1",colnames(df))
    df <- t(rowsum(t(df), group = colnames(df)))
    return(df)
}

# posterior summaries for all
header <- c("female", "male")
pm_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(pm_all) <- header
psd_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(psd_all) <- header
lfsr_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(lfsr_all) <- header

# split calculation into chunks of 30k
interval <- 30000
num <- floor(nrow(data$Bhat) / interval)
for (i in 0:num) {
    print(paste0("progress: ", i,"/",num))
    start <- (i*interval) + 1
    end <- (i+1)*interval
    if (end > nrow(data$Bhat)) {
        end <- nrow(data$Bhat)
    } 
    datasub=mash_set_data(data$Bhat[start:end,], data$Shat[start:end,])
    msub = mash(datasub, g=g_ave, fixg=TRUE)    
    if (i == 0) {
        cov_names <- colnames(msub$posterior_weights)
        weights_all <- data.frame(matrix(ncol = length(cov_names), nrow = 0)) ; colnames(weights_all) <- cov_names
    }
    pm = get_pm(msub)
    psd = get_psd(msub)
    lfsr = get_lfsr(msub)
    pm_all <- rbind(pm_all, pm)
    psd_all <- rbind(psd_all, psd)
    lfsr_all <- rbind(lfsr_all, lfsr)
    weights <- weight_col(msub$posterior_weights)
    weights <- round(weights, 10)
    weights_all <- rbind(weights_all, weights)
}

# write posterior estimates to table
write.table(pm_all, file=paste0(pheno,"_mash_pm",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(psd_all, file=paste0(pheno,"_mash_psd",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(lfsr_all, file=paste0(pheno,"_mash_lfsr",suffix,".txt"), sep="\t", row.names=FALSE)
write.table(weights_all, file=paste0(pheno,"_mash_weights",suffix,".txt"), sep="\t", row.names=FALSE)

print(paste0(pheno," - done"))