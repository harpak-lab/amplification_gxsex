#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','PHENO_DIR'))])
library(optparse, lib.loc=R_LIB)

# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)

setwd(PHENO_DIR)
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", head=TRUE, col.names=c("FID","IID",pheno), colClasses=c("integer","integer","numeric"))
df_sex <- read.table("sex_ids.txt", sep="\t", head=TRUE, col.names=c("IID","sex"))        

# merge sex, and pheno
df <- merge(df_pheno, df_sex, by="IID")
female_df <- df[df$sex == 0,]
male_df <- df[df$sex == 1,]

# standardized within sex
female_df[,3] <- (female_df[,3] - mean(female_df[,3])) / sd(female_df[,3])
male_df[,3] <- (male_df[,3] - mean(male_df[,3])) / sd(male_df[,3])

# combine df and write as text file
df <- rbind(female_df[,1:3], male_df[,1:3])
colnames(df) <- c("#FID","IID",pheno)
setwd("/scratch/08005/cz5959/Phenotypes/standardized")
write.table(df, file=paste0("pheno_",pheno,"_std.txt"), sep="\t", row.names=FALSE, quote=FALSE)
