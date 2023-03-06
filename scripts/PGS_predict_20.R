#!/usr/bin/env Rscript

# load libraries
source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR',"PHENO_DIR"))])
library(optparse, lib.loc=R_LIB)
# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)

p.threshold <- c('1', '0.01','1e-5','1e-8')

# load phenotype and covariate file
setwd(PHENO_DIR)
pheno_file <- paste0("pheno_",pheno,".txt")
phenotype <- read.table(pheno_file, sep="\t", head=FALSE, col.names=c("FID","IID",pheno), colClasses=c("integer","integer","numeric"))
covariates <- read.table("covariates.txt", header=FALSE, colClasses=c(rep("integer",2), rep("numeric", 12), rep("NULL", 3)))
colnames(covariates)=c("FID","IID",paste0("PC",1:10),"sex","birthyear")
# merge phenotype and covariate file
pheno_covar <- merge(phenotype, covariates, by=c("FID", "IID"))

### FUNCTIONS

pgs_null <- function(sex) {
    # keep test set ids
    pheno_covar_null <- pheno_covar[pheno_covar$IID %in% test_ids[[sex]]$IID,]
    # null model
    null.model <- lm(paste0(pheno," ~ ."), data=pheno_covar_null[,!colnames(pheno_covar_null) %in% c("FID","IID")])
    null.r2 <- summary(null.model)$r.squared
    return(null.r2)
}

# sex: male or female
# type: using additive or covariance-aware (mash)
# mode: sex-specific or standardized or nothing
pgs_prediction <- function(sex, type = "additive", mode = "base") {
    pgs.result <- NULL
    # go through each p-value threshold
    for (i in p.threshold){

        # load .profile based on model 
        if (mode == "sex-specific") {
            pgs_f <- read.table(paste0("female_",type,"_",pheno,".",i,".profile"), header=T)
            pgs_m <- read.table(paste0("male_",type,"_",pheno,".",i,".profile"), header=T)
            pgs_f <- pgs_f[pgs_f$IID %in% test_ids[['female']]$IID,]
            pgs_m <- pgs_m[pgs_m$IID %in% test_ids[['male']]$IID,]
            pgs <- rbind(pgs_f,pgs_m)
        } else if (mode == "standardized") {
            pgs <- read.table(paste0("both_sex_std_additive_",pheno,".",i,".profile"), header=T)
        } else {
            pgs <- read.table(paste0("both_sex_additive",pheno,".",i,".profile"), header=T)
        }
        pgs <- pgs[pgs$IID %in% test_ids[[sex]]$IID,]

        pgs.snps <- mean(pgs$CNT)/2
        # Merge pgs with phenotype matrix - only take FID, IID and pgs 
        pheno.pgs <- merge(pheno_covar, pgs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
        # regression
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
    print(paste("Best threshold: ", pgs.result[which.max(pgs.result$R2),1]))
    return(pgs.result[which.max(pgs.result$R2),c(2,3)])
}

### MAIN
add_result_list <- NULL
as_result_list <- NULL
ss_result_list <- NULL
m_result_list <- NULL
for (i in c(1:20)) {
    print(i)
    wd <- paste0(GWAS_DIR, "/", pheno, "/PGS_", set)
    setwd(wd)

    # get test ids
    female_ids <- read.table(paste0(pheno,"_female_testIIDs.txt"), sep="\t", col.names=c("FID","IID"), colClasses=c("NULL","integer"))
    male_ids <- read.table(paste0(pheno,"_male_testIIDs.txt"), sep="\t", col.names=c("FID","IID"), colClasses=c("NULL","integer"))
    both_sex_ids <- rbind(female_ids,male_ids)
    test_ids <- list(female = female_ids, male = male_ids, both_sex = both_sex_ids)

    # obtain null model r2
    null.r2 <- pgs_null("both_sex")

    # both-sex additive
    add_result <- pgs_prediction("both_sex","additive", "base")
    add_result_list <- rbind(add_result_list, add_result)

     # both-sex additive, standardized by sex
    as_result <- pgs_prediction("both_sex","additive", "standardized")
    as_result_list <- rbind(as_result_list, as_result)

    # sex-specific additive
    ss_result <- pgs_prediction("both_sex", "additive", "sex-specific")
    ss_result_list <- rbind(ss_result_list, ss_result)

    # sex-specific covariance-aware
    m_result <- pgs_prediction("both_sex", "mash", "sex-specific")
    m_result_list <- rbind(m_result_list, m_result)
}

#TODO# combine the above lists and then write as table
add_result_list$type <- rep("add", nrow(add_result_list))
as_result_list$type <- rep("as", nrow(as_result_list))
ss_result_list$type <- rep("ss", nrow(add_result_list))
m_result_list$type <- rep("mash", nrow(add_result_list))
all_results <- rbind(rbind(rbind(add_result_list, ss_result_list), m_result_list), as_result_list)

wd <- paste0(GWAS_DIR, "/", pheno, "/PGS_1")
setwd(wd)
file_name <- paste0(pheno, "_PGS20.txt")
write.table(all_results, file=file_name, sep="\t", row.names=FALSE, quote=FALSE)

print(paste0("File, ", file_name, ", in directory: ", wd))




