#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)
library(optparse, lib.loc=R_LIB)

# arguments and set working directory
# argument parser
option_list = list(
    make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
    make_option(c("-m","--mode"), type="character", default="",help="input '_same' if using same size subset",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
pheno <- opt$pheno; print(pheno)
method <- opt$mode

wd <- paste0(GWAS_DIR,"/",pheno,"/mash")
setwd(wd)
load(file= paste(pheno,"_mash.RData",sep=""))
set.seed(1)

#### MASH ####
# random subset
random_subset <- function(method="") {
    random <- numeric(0)
    if (method=="") {
        # METHOD3: LD blocks
        for (i in unique(LD_groups_subset$group)) {
        sample_subset <- LD_groups_subset[LD_groups_subset$group == i, 'index']
        random[[(length(random) + 1)]] <- sample(sample_subset,1)
        }
    } else {
        # METHOD4: LD blocks, same sample size
        unique_groups <- sample(unique(LD_groups_subset$group))
        while (length(random) <= 1703) {
            for (i in unique_groups) {
                if (length(random) = 1703) {
                    break
                }
                sample_subset <- LD_groups_subset[LD_groups_subset$group == i, 'index']
                random[[(length(random) + 1)]] <- sample(sample_subset,1)
            }
        }  
    }
    print(paste0("random subset ", length(random)))
    return(random)
}

fit_mash <- function(random) {
    #correlation structure
    data.temp = mash_set_data(data$Bhat[random,],data$Shat[random,])
    Vhat = estimate_null_correlation_simple(data.temp)
    rm(data.temp)
    data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)

    # set up canoncial covar matrices and add hypothesis
    U.c = cov_canonical(data.random)
    # add hypothesis
    corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
    effect = c(1.5,2,3)
    for (c in corr) {
        for (e in effect) {
            U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
            U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
        }
    }
    U.c[['equal_-0.25_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
    U.c[['equal_-0.5_1']] <- matrix(c(1,-0.5,-0.5,1),2,2)
    U.c[['equal_-0.75_1']] <- matrix(c(1,-0.75,-0.75,1),2,2)
    U.c[['equal_-1_1']] <- matrix(c(1,-1,-1,1),2,2)
    names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")


    # fit mash model 
    m = mash(data.random, Ulist= U.c, outputlevel = 1)
    # mixture proportions
    mixture_prop <- get_estimated_pi(m)
    return(mixture_prop)
}

# p-value thresholds for LD_groups
thresholds <- [1, 5e-2, 1e-5, 5e-8]
for (p in thresholds) {
    LD_groups_subset <- LD_groups[LD_groups$P.x < p | LD_groups$P.y < p,]
    for (method %in% c("", "_same")){
        random <- random_subset(method=method)
        m <- fit_mash(random)
        write.table(m, file=paste0(pheno,"_",as.character(p),method,".txt"), sep="\t", row.names=FALSE, quotes=FALSE)
    }
}

print(paste0("Results in: ", wd))