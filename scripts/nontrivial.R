#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','LDSC_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggbreak, lib.loc=R_LIB)
library(dplyr, lib.loc=R_LIB)
library(gridExtra, lib.loc=R_LIB)
library(grid, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)

# get correlation
setwd(LDSC_FILE)
ldsc <- read.csv("relative_h2.txt", sep="\t")
ldsc_both <- ldsc[ldsc$Sex == "both_sex", c(1,2,6)]
# get ratio of relative h2
ldsc_ratio <- ldsc[ldsc$Sex != "both_sex", c(1:4,8)]
ldsc_ratio <- ldsc_ratio %>%
  dcast(Code + Phenotype ~ Sex, value.var = "relative_h2") %>%
  mutate(relative_h2_ratio = ifelse(male > female, male / female, female / male)) %>%
  select(c(1,5))

# get mash weights
setwd(GWAS_DIR)
mash <- read.csv("mash_weights.txt", sep="\t")

df <- merge(merge(ldsc_both, mash, by.x = "Code", by.y = "phenotype"), ldsc_ratio, by = "Code")
# estimate difference in weights, weight on non-perfect correlation and unequal magnitude
df <- df %>%
  mutate(diff_weight = sum_weight_m - sum_weight_f) %>%
  mutate(nonperf_corr = 1 - (perf_corr/100)) %>%
  mutate(unequal_mag = 1 - (sum_weight_e/100)) %>%
  mutate(facet = ifelse(Phenotype == "Testosterone", 1, 2)) %>%
  select(c(2,8,10,11,12))

### SCATTER PLOT ###
# x axis is (1 - equal magnitude)
# y axis is (1 - perfect correlation) 
df_poly <- data.frame(x = c(0.37, 0.37, 0.96, 1.03, 1.03), 
                      y = c(0.29, 0.37, 0.96, 0.96, 0.29))
pdf(file="nontrivial_weight.pdf", width=6, height=4)
ggplot(data= df, aes(x= unequal_mag, y= nonperf_corr)) +
  geom_point(aes(size = relative_h2_ratio), color = "#2b62d9") +
  geom_abline(slope = 1, intercept = 0, color = "gray") +
  geom_polygon(data = df_poly, aes(x,y), fill="gray", alpha = 0.2) +
  scale_x_continuous(breaks = seq(0.4,1,0.1), limit = c(0.37,1.03), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0.3,0.9,0.1), limit = c(0.29,0.96), expand = c(0,0)) +
  theme_classic() +
  theme(strip.text.x = element_blank(), plot.title = element_text(size=14), axis.title = element_text(size=11),
        axis.text = element_text(size = 9), legend.position = "none", 
        plot.margin = margin(10,10,10,10)) +
  labs(title = "Weight on Nontrivial Correlation vs Nontrivial Magnitude Matrices", 
       x = "Proportion of Weight on Unequal Magnitude Matrices", 
       y = "Proportion of Weight on \nImperfect Correlation Matrices",
       caption = "* heritability ratio between sexes = \nlarger heritability / smaller heritability") +
  geom_text_repel(aes(label = Phenotype), box.padding = 0.7, max.overlaps = Inf, size=3) +
  annotate("text", x = 0.4, y = 0.9, hjust = 0, 
           label = "point size proportional to heritability ratio") 
dev.off()



