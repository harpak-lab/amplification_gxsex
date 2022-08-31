#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
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
pheno <- opt$pheno; print(pheno)
mode <- opt$set; print(mode)

# load mash data from setup
if (mode == "additive") {
    suffix <- "" ; wd <- paste0(GWAS_DIR,"/",pheno,"/mash")

} else { 
    suffix <- "_pgs" ; wd <- paste0(GWAS_DIR,"/",pheno,"/PGS_",mode)
} 
setwd(wd)
load(file= paste0(pheno,"_mash",suffix,".RData")) 

#### MASH ####
# random subset
random_subset <- function(seed=1) {
    set.seed(seed); print(paste0("Seed #: ", seed))
    # METHOD3: LD blocks
    random <- numeric(0)
    for (i in unique(LD_groups$group)) {
        sample_subset <- LD_groups[LD_groups$group == i, 'index']
        random[[(length(random) + 1)]] <- sample(sample_subset,1)
    }      
    print(paste0("random subset length: ", length(random)))
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
    # mixture model
    mixture_prop <- get_estimated_pi(m)
    g <- get_fitted_g(m)
    return(list(mixture_prop, g))
}

# repeat 100 times to get average of fitted model
rep=100
g_list <- list()
for (i in 1:rep) {
    random <- random_subset(i)
    # sample w/o replacement
    LD_groups <- LD_groups[ ! LD_groups$index %in% random,]
    # mash
    results <- fit_mash(random)
    mix <- results[[1]]
    g <- results[[2]]
    if (i==1) {
        mixture <- matrix(names(mix), ncol=1)
        g_ave <- g
    } else {
        # combine all g pi and grid
        g_ave$pi <- g_ave$pi + g$pi
        g_ave$grid <- g_ave$grid + g$grid
    }
    mixture <- cbind(mixture, mix)
    g_list[[i]] <- g
}
# get average for g_all by dividing by number of repetitions
g_ave$pi <- g_ave$pi / rep
g_ave$grid <- g_ave$grid / rep

# save fitted results
colnames(mixture) <- cbind(paste0("mix_",0:rep))
if (mode == "additive") {
    write.table(mixture, file=paste0(pheno,"mixprop_100_all.txt"), sep="\t", row.names=FALSE)
    save(g_list, g_ave, file= paste0(pheno,"_mash_100g.RData"))
} else {
    save(g_ave, file= paste0(pheno,"_mash_100g_pgs.RData"))
}
