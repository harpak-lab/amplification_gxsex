#!/usr/bin/env Rscript

# STEP 3: TEST EVIDENCE OF ENVIRONMENTAL VARIANCE IN MASH
source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','QC_DIR'))])
library(optparse, lib.loc=R_LIB)
library(ashr, lib.loc=R_LIB)
library(mashr, lib.loc=R_LIB)

setwd(GWAS_DIR)

snp_nums <- c(100, 1000, 10000)
E_ratios <- c(1,1.5,5)
h2s <- c(0.05,0.1,0.5)
for (snp_num in snp_nums) {
    for (E_ratio in E_ratios) {
        for (h2 in h2s) {
            load(paste0("mash_",snp_num,"_",h2,"_",E_ratio,".RData"))
            ## mash
            mash_BETA <- as.matrix(mash_BETA) ; mash_SE <- as.matrix(mash_SE)
            data <- mash_set_data(mash_BETA, mash_SE, zero_Shat_reset = .Machine$double.eps)
            # set up covariance matrices
            U.c = cov_canonical(data)
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
            m.c <- mash(data,U.c)
            mixture <- get_estimated_pi(m.c)
            write.table(mixture, file=paste0(snp_num,"_",h2,"_",E_ratio,".txt"),row.names = FALSE, sep="\t")
        }
    }
}

print(paste0("Files (.txt) in directory: " GWAS_DIR))


