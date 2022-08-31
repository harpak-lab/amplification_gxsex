#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR',"PHENO_DIR"))])
library(optparse, lib.loc=R_LIB)
# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-s","--set"), type="character", default="1",help="set number",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
set <- opt$set; print(set)

p.threshold <- c('1', '0.01', '1e-5', '1e-8')

# load phenotype and covariate file
setwd(PHENO_DIR)
pheno_file <- paste0("pheno_",pheno,".txt")
phenotype <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("integer","integer","numeric"))
covariates <- read.table("covariates.txt", header=FALSE, colClasses=c(rep("integer",2), rep("numeric", 12), rep("NULL", 3)))
colnames(covariates)=c("FID","IID",paste0("PC",1:10),"sex","birthyear")
# merge phenotype and covariate file
pheno_covar <- merge(phenotype, covariates, by=c("FID", "IID"))

# test set
wd <- file.path(GWAS_DIR,pheno,paste0("PGS_",set))
setwd(wd)
female_ids <- read.table(paste0(pheno,"_female_testIIDs.txt"), sep="\t", col.names=c("FID","IID"), colClasses=c("NULL","integer"))
male_ids <- read.table(paste0(pheno,"_male_testIIDs.txt"), sep="\t", col.names=c("FID","IID"), colClasses=c("NULL","integer"))
both_sex_ids <- rbind(female_ids,male_ids)
test_ids <- list(female = female_ids, male = male_ids, both_sex = both_sex_ids)

### FUNCTIONS
pgs_null <- function(sex) {
    # keep test set ids
    pheno_covar_null <- pheno_covar[pheno_covar$IID %in% test_ids[[sex]]$IID,]
    # null model
    null.model <- lm(paste0(pheno," ~ ."), data=pheno_covar_null[,!colnames(pheno_covar_null) %in% c("FID","IID")])
    null.r2 <- summary(null.model)$r.squared
    print(paste0("Null R2 for ", sex, ": ", null.r2))
    return(null.r2)
}

pgs_prediction <- function(sex, type = "additive", null.r2, mode = 'specific') {  
    pgs.result <- NULL
    for(i in p.threshold){
        print(i)
        # go through each p-value threshold
        if (mode == 'combined') {
            pgs_f <- read.table(paste0("female_",type,"_",pheno,".",i,".profile"), header=T)
            pgs_m <- read.table(paste0("male_",type,"_",pheno,".",i,".profile"), header=T)
            pgs_f <- pgs_f[pgs_f$IID %in% test_ids[['female']]$IID,]
            pgs_m <- pgs_m[pgs_m$IID %in% test_ids[['male']]$IID,]
            pgs <- rbind(pgs_f,pgs_m)
        } else if (mode == 'both_sex_additive') {
            pgs <- read.table(paste0("both_sex_additive_",pheno,".",i,".profile"), header=T)
        } else {
            pgs <- read.table(paste0(sex,"_",type,"_",pheno,".",i,".profile"), header=T)
        }
        pgs <- pgs[pgs$IID %in% test_ids[[sex]]$IID,]
        
        pgs.snps <- mean(pgs$CNT)/2
        # Merge pgs with phenotype matrix - only take FID, IID and pgs 
        pheno.pgs <- merge(pheno_covar, pgs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
        # logistic regression
        model <- lm(paste0(pheno," ~ ."), data=pheno.pgs[,!colnames(pheno.pgs)%in%c("FID","IID")])
        model.r2 <- summary(model)$r.squared        # model R2
        pgs.r2 <- model.r2-null.r2                  # incremental R2
        # other summary results                                           
        pgs.coef <- summary(model)$coefficients["SCORE",]
        pgs.beta <- as.numeric(pgs.coef[1])
        pgs.se <- as.numeric(pgs.coef[2])
        pgs.p <- as.numeric(pgs.coef[4])
        # store results
        pgs.result <- rbind(pgs.result, data.frame(Threshold=i, R2=model.r2, incR2=pgs.r2, P=pgs.p, BETA=pgs.beta, SE=pgs.se, SNP=pgs.snps))        
    }
    # Best threshold based on R2 is:
    print(paste("Best result for", sex, type, mode, sep=", "))
    print(pgs.result[which.max(pgs.result$R2),])
    # return results
    #pgs_summary <- list("results" = pgs.result, "residuals" = pgs.residuals, "fitted" = pgs.fitted)
    #return(pgs_summary)
}

### MAIN

# null
female_null <- pgs_null("female"); male_null <- pgs_null("male"); both_sex_null <- pgs_null("both_sex"); 

print("Order for each method: \n 1) prediction on both-sex \n 2) prediction on females 3) prediction on males")

# additive
print("COVARIANCE-NAIVE, SEX-SPECIFIC")
pgs_prediction("both_sex", "additive", both_sex_null, "combined")
pgs_prediction("female", "additive", female_null)
pgs_prediction("male", "additive", male_null)

# mash
print("COVARIANCE-AWARE, SEX-SPECIFIC")
pgs_prediction("both_sex", "mash", both_sex_null, "combined")
pgs_prediction("female", "mash", female_null)
pgs_prediction("male", "mash", male_null)

# both-sex additive
print("ADDITIVE, BOTH-SEX")
pgs_prediction("both_sex","additive", both_sex_null, "both_sex_additive")
pgs_prediction("female","additive", female_null, "both_sex_additive")
pgs_prediction("male","additive", male_null, "both_sex_additive")


