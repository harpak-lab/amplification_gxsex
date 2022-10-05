#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library(optparse, lib.loc=R_LIB)

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-s", "--set"), type = "character", default = "1", help = "input set # of PGS", metavar = "character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
mode <- opt$set; print(mode)

# load mash lfsr and pm results
wd <- paste0(GWAS_DIR, "/", pheno, "/PGS_", mode)
setwd(wd)
lfsr_df <- read.table(paste0(pheno, "_mash_lfsr_pgs.txt"), sep="\t", head=TRUE, colClasses=c(rep("numeric", 2)))
pm_df <- read.table(paste0(pheno, "_mash_pm_pgs.txt"), sep="\t", head=TRUE, colClasses=c(rep("numeric", 2)))

get_pseudo_p <- function(sex) {
    # load gwas p-values, 2 columns: ID and P
    file_name <- paste0(sex, "_train.", pheno, ".glm.linear")
    gwas_df <- read.table(file_name, sep = "\t", head = FALSE,
    col.names = c("CHROM", "POS", "ID", "REF", "ALT", "A1", "AX", "TEST", "OBS_CT", "BETA", "SE", "TSTAT", "P"),
    colClasses = c(rep("NULL", 2), "character", rep("NULL", 2), "character", rep("NULL",6), "numeric"))

    # sort p-values
    sorted_p <- gwas_df[order(gwas_df$P), 3]
    # concatenate ids, A1, lfsr and pm
    if (sex == "female") {
        mash_ids <- data.frame(gwas_df$ID, gwas_df$A1, pm_df$female, lfsr_df$female)
    } else {
        mash_ids <- data.frame(gwas_df$ID, gwas_df$A1, pm_df$male, lfsr_df$male)
    }
    colnames(mash_ids) <- c("ID", "A1", "pm", "lfsr")

    # sort by lfsr
    mash_ids <- mash_ids[order(mash_ids$lfsr), ]

    # concatenate sorted ids with pseudo pvalues
    mash_ids$P <- sorted_p

    # write to text file
    write.table(mash_ids, file=paste0(sex,"_pseudoP_pgs.",pheno,".txt"), sep="\t", row.names=FALSE, quote=FALSE)

    return(mash_ids)
}
