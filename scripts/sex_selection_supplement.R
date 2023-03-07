#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','SEL_DIR'))])
library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("reshape2", lib.loc=R_LIB)
library("matrixStats", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("ggsci", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)

setwd(SEL_DIR)

# P-VALUE PLOTS
pval_plot <- function(pvalue) {

# load RData files
load(paste0("asj_zscore_plot.",pvalue,".all.Rdata"))
asj <- results; asj_i <- resultsi
load(paste0("nfe_zscore_plot.",pvalue,".all.Rdata"))
nfe <- results; nfe_i <- resultsi
load(paste0("fin_zscore_plot.",pvalue,".all.Rdata"))
fin <- results; fin_i <- resultsi
load(paste0("afr_zscore_plot.",pvalue,".all.Rdata"))
afr <- results; afr_i <- resultsi
load(paste0("amr_zscore_plot.",pvalue,".all.Rdata"))
amr <- results; amr_i <- resultsi
load(paste0("eas_zscore_plot.",pvalue,".all.Rdata"))
eas <- results; eas_i <- resultsi

# merge ancestry results (new)
results <- rbind(asj, nfe, fin, afr, amr, eas)
resultsi <- rbind(asj_i, nfe_i, fin_i, afr_i, amr_i, eas_i)
resultsi$ANC <- c("Ashkenazi Jewish", "Non-Finnish European", "Finnish", "African", "Latino/American Admixed", "East Asian")

### ZSCORE PLOT - 6 ancestries
pdf(file=paste0("selection_",pvalue,".pdf"), width=6, height=10)
print(
ggplot(results, aes(x=ZMEAN, y=TRAIT)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=CILOWER, xmax=CIUPPER), height=0.3, size=0.3) +
  geom_vline(xintercept = 0, linetype="longdash", color="#563f61", alpha=0.5, size=0.3) +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN), linetype="dashed", color="#b02e38", alpha=0.5, size=0.3) +
  #coord_cartesian(xlim=c(-2.9,4.4)) +
  scale_x_continuous(breaks = c(-2,0,2,4), limits=c(-2.9,4.4)) +
  facet_wrap(~ANC,ncol=3) +
  theme_classic() +
  xlab("Z-score for Sexually-Antagonistic Selection") +
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=12), axis.text=element_text(size=9),
        panel.grid.major.y = element_line(color="gray95", size=0.5))
)
dev.off()

}

# Call function for 3 pvalue thresholds
pval_plot("1e-8")
pval_plot("1e-5")
pval_plot("1e-3")


### UKB ZSCORE PLOT 
load("UKB_zscore_plot.1e-05.all.Rdata")
pdf(file="selection_UKB.pdf", width=4, height=5)
print(
ggplot(results, aes(x=ZMEAN, y=TRAIT)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=CILOWER, xmax=CIUPPER), height=0.3, size=0.3) +
  geom_vline(xintercept = 0, linetype="longdash", color="#563f61", alpha=0.5, size=0.3) +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN), linetype="dashed", color="#b02e38", alpha=0.5, size=0.3) +
  #coord_cartesian(xlim=c(-2.9,4.4)) +
  scale_x_continuous(breaks = c(-2,0,2,4), limits=c(-3.4,4.4)) +
  facet_wrap(~ANC,ncol=3) +
  theme_classic() +
  xlab("Z-score for Sexually-Antagonistic Selection") +
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=12), axis.text=element_text(size=9),
        panel.grid.major.y = element_line(color="gray95", size=0.5))
)
dev.off()

print(paste0("Files, selection_*.pdf, in directory: ", SEL_DIR))