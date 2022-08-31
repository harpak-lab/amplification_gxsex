#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("tidyverse", lib.loc=R_LIB)
library("reshape2", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("ggsci", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)

# user input phenotype
option_list = list(
  make_option(c("-m","--mode"), type="character", default="",help="input '_same' if using same size subset",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method <- opt$mode

wd <- paste0(GWAS_DIR)
setwd(wd)

all_null <- read.csv(paste0("noeffect_weight",method,".txt", sep="\t")

# NO EFFECT PLOT
pdf(file=paste0(pheno,"_null",method,".pdf"), width=4, height=3)
colors <- c("#d1b724", "#563f61", "#b0464f", "#2b62d9")
ggplot(all_null, aes(x=pval, y=Weight, group=Phenotype, color=Phenotype)) +
  geom_point() +
  geom_line() +
  labs(x="P-value Threshold") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(axis.title=element_text(size=10), axis.text=element_text(size=9))

print(paste0("Results in: ", wd))