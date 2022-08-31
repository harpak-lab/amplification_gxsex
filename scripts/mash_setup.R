#!/usr/bin/env Rscript


source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','LD_SCORE'))])
library(optparse, lib.loc=R_LIB)
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
pheno
pheno <- opt$pheno

# load LD scores and set working directory
load(file=paste(LD_SCORE,"/LD_groups.RData"))
wd <- paste0(GWAS_DIR,"/",pheno)
setwd(wd)

setup_sex_df <- function(sex) {
    # load results - summstats
    if (opt$set == "additive") {
        file_name <- paste0(sex,"_all.",pheno,".glm.linear")
    } else {
        file_name <- paste0(wd,"/PGS_",opt$set,"/",sex,"_train.",pheno,".glm.linear")
    }
    gwas_df <- read.table(file_name,sep="\t",head=FALSE, 
    col.names=c("CHROM","POS","ID","REF","ALT","A1","AX","TEST","OBS_CT","BETA","SE","TSTAT","P"), 
    colClasses = c(rep("integer",2), rep("character", 3), rep("NULL", 4), rep("numeric",2), "NULL", "numeric")) 
    # add new ID column CHROM:POS:REF:ALT:ID, and index column
    gwas_df$VAR <- paste(gwas_df$CHROM, gwas_df$POS, gwas_df$REF, gwas_df$ALT, gwas_df$ID, sep=":")
    gwas_df$index <- seq.int(nrow(gwas_df))

    return(gwas_df)
}

# format sex-specific summary statistics
female_df <- setup_sex_df("female")
male_df <- setup_sex_df("male")

# take union of SNPs if file lengths do not match
if (nrow(female_df) != nrow(male_df)) {
    print('ERROR: different gwas lengths!')
    if (nrow(female_df) > nrow(male_df)) {
        print( female_df[!female_df$ID %in% male_df$ID,] )
        female_df <- female_df[female_df$ID %in% male_df$ID,]
    } else {
        print( male_df[!male_df$ID %in% female_df$ID,] )
        male_df <- male_df[male_df$ID %in% female_df$ID,]
    }
}

# merge LD groups with p-values; pos_groups from loaded LD scores
LD_groups <- merge(pos_groups, female_df[c("index","P")], by="index")
LD_groups <- merge(LD_groups, male_df[c("index","P")], by="index")


# create matrix of BETA and SE
conditions <- c("female", "male")
r <- nrow(female_df)
BETA <- matrix(c(female_df$BETA, male_df$BETA), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))
SE <- matrix(c(female_df$SE, male_df$SE), nrow=r, ncol=2, dimnames=list(c(female_df$VAR),conditions))

# create mash data object
data = mash_set_data(BETA, SE)

# save variables into Rdata
if (opt$set != "additive") {
    save(data, LD_groups, file= paste0(wd,"/PGS_",opt$set,"/",pheno,"_mash_pgs.RData"))
} else {
    dir.create(file.path(wd, "mash"))
    summstat_pos <- female_df$POS
    save(summstat_pos, LD_groups, data, file= paste0(wd,"/mash/",pheno,"_mash.RData"))
}