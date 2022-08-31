source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','SEL_FILE'))])

library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("reshape2", lib.loc=R_LIB)
library("matrixStats", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("ggsci", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)

# load RData files (from Matt)
setwd(SEL_FILE)
load("fst_plot_testosterone.1e-05.RData"); testosterone <- pointsf
load("fst_plot_protein_total.1e-05.RData"); protein <- pointsf
load("zscore_plot.1e-05.RData")

# split results by ancestry
azj <- results[results$ANC == "Ashkenazi Jewish",]
fin <- results[results$ANC == "Finnish",]
nfe <- results[results$ANC == "Non-Finnish European",]

# lm line
t_model <- lm(FST~V, testosterone, weight=w)
t_B <- summary(t_model)$coefficients[2]; t_yi <- summary(t_model)$coefficients[1]
p_model <- lm(FST~V, protein, weight=w)
p_B <- summary(p_model)$coefficients[2]; p_yi <- summary(p_model)$coefficients[1]

### FST PLOT
pdf(file="FIG7C.selection.pdf", width=3, height=4)
f1 <- ggplot(testosterone, aes(x=V, y=FST, weight=w, size=w)) +
  geom_point(color = "black", alpha= 0.2) +
  geom_abline(slope=t_B, intercept=t_yi, size=0.5, color="#3a53a4") +
  theme_classic() +
  scale_y_continuous(breaks=c(0,5e-5,1e-4), labels=c("0","5e-5","1e-4"), limits = c(0,1e-4)) +
  scale_x_continuous(breaks=c(0,5e-4,1e-3, 1.5e-3), labels=c("0","5e-4","1e-3", "1.5e-3"), limits=c(0,1.5e-3)) +
  xlab("VGxSex (nmol/L)") +
  theme(axis.title.x=element_text(size=10), axis.text=element_text(size=8), axis.title.y=element_blank(),
        legend.position="none") +
  scale_size(range = c(0.5,3)) 

f2 <- ggplot(protein, aes(x=V, y=FST, weight=w,size=w)) +
  geom_point(color = "black", alpha= 0.2) +
  geom_abline(slope=p_B, intercept=p_yi, size=0.5, color="=#6f3c96") +
  theme_classic() +
  scale_y_continuous(breaks=c(0,4e-5, 8e-5), labels=c("0","4e-5","8e-5"), limits = c(0,9e-5)) +
  scale_x_continuous(breaks=c(0,5e-4,1e-3), labels=c("0","5e-4","1e-3"), limits=c(0,1.1e-3)) +
  xlab("VGxSex (g/L)") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), axis.text=element_text(size=8),
        legend.position="none") +
  scale_size(range = c(0.5,3)) 

fa <- ggplotGrob(f1) ; fb <- ggplotGrob(f2)
maxWidth = grid::unit.pmax(fa$widths[2:5], fb$widths[2:5])
fa$widths[2:5] <- as.list(maxWidth) ; fb$widths[2:5] <- as.list(maxWidth)
fplot <- grid.arrange(fb, fa, ncol=1, nrow=2)
annotate_figure(fplot, 
                left = text_grob("FST Between Males and Females", size=10, rot=90))
dev.off()

### ZSCORE PLOT
pdf(file="FIG7D.selection.pdf", width=6, height=8)
ggplot(results, aes(x=ZMEAN, y=TRAIT)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=CILOWER, xmax=CIUPPER), height=0.3, size=0.3) +
  geom_vline(xintercept = 0, linetype="longdash", color="#563f61", alpha=0.5, size=0.3) +
  geom_vline(data = resultsi, aes(xintercept = ANCMEAN), linetype="dashed", color="#b02e38", alpha=0.5, size=0.3) +
  coord_cartesian(xlim=c(-2.6,5.35)) +
  facet_wrap(~ANC,ncol=3) +
  theme_classic() +
  xlab("Z-score for Sexually-Antagonistic Selection") +
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=12), axis.text=element_text(size=9),
        panel.grid.major.y = element_line(color="gray95", size=0.5))
dev.off()


