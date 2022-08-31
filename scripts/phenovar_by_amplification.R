#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','PHENO_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(dplyr, lib.loc=R_LIB)
library(gridExtra, lib.loc=R_LIB)
library(grid, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)

# get mash weights
setwd(GWAS_DIR)
mash_weights <- read.csv("mash_weights.txt", sep="\t")

# get phenotype names and sex
setwd(PHENO_DIR)
pheno_names <- read.csv("pheno_names.txt", sep="\t", colClasses=c(rep("character",2), "NULL"))
pheno_vars <- read.csv("pheno_meanvar.txt", sep="\t")

# estimate phenotypic variance and mean ratio
pheno_vars <- pheno_vars %>%
  merge(pheno_names, by.x="pheno", by.y="Code") %>%
  mutate(var_ratio = m_var / f_var) %>%                 # pheno var ratio
  mutate(mean_ratio = m_mean / f_mean) %>%              # pheno mean ratio
  select(c(1,8,9,10))

# create column of difference of male and female mash weights
mash_weights <- mash_weights %>%
  mutate(diff_mf = sum_weight_m - sum_weight_f) %>%     # amplification difference
  select(c(1,6))

# combine phenotypic variance and mash weights df
df <- merge(pheno_vars, mash_weights, by.x = "pheno", by.y = "phenotype")

#### PLOT ####
# split df by outliers (arm_fatfree_mass and testosterone)
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()
anno_size=2.8   # annotation text size
axis_text = 9   # axis text size

pdf(file="FIG4b_phenovar_amplification.pdf", width=5.5, height=3.5)
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=axis_text),
        axis.ticks.x=element_blank(), axis.line.x.bottom = element_blank())

g2 <- ggplot(data=df3, aes(x= diff_mf, y= var_ratio)) +
  geom_point(size=2, color='#2b62d9') + 
  geom_hline(yintercept = 27.6) + geom_vline(xintercept = 100) +
  scale_x_continuous(breaks=c(85,90,95), limits=c(80,100), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

g3 <- ggplot(data=df1, aes(x= diff_mf, y= var_ratio)) +
  geom_smooth(method="lm", data=df_corr, aes(x=diff_mf, y=var_ratio), inherit.aes = F,
              color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=seq(-60,60,20), limits=c(-75,65), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3,seq(0.5,3.5,0.25)), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text)) +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)

g4 <- ggplot(df2, aes(x= diff_mf, y= var_ratio)) + 
  geom_smooth(method="lm", data=df_corr, aes(x=diff_mf, y=var_ratio), inherit.aes = F,
              color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  geom_vline(xintercept = 100) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(breaks=c(85,90,95), limits=c(80,100), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(size=axis_text),
        axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size) 

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g3)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, gB, g4, ncol=2, nrow=2, widths=c(5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Difference between fractions of variants with M>F and M<F effect", size=11),    # x-axis label
                left = text_grob("Ratio of Male to Female Phenotypic Variances", size=12, rot=90),                  # y-axis label
                top = text_grob("Phenotypic Variance Strongly Correlated with Amplification", size = 14))           # title
dev.off()

print("Pearson correlation between amplification difference and phenotypic variance ratio: ")
df_corr <- df[df$pheno != 'testosterone',]    # remove outlier - testosterone
model <- cor.test(df_corr$diff_mf, df_corr$var_ratio)
print(model)

print("Pearson correlation between amplification difference and phenotypic mean ratio: ")
df_corr <- df[!df$pheno %in% c('testosterone','wth_bmi_adj'),]     # remove outlier - testosterone and bmi adjusted waist:hip
model <- cor.test(df_corr$mean_ratio, df_corr$diff_mf)
print(model)
