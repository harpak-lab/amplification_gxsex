#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','PHENO_DIR','LDSC_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggbreak, lib.loc=R_LIB)
library(dplyr, lib.loc=R_LIB)
library(gridExtra, lib.loc=R_LIB)
library(grid, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)

# get phenotypic variance
setwd(PHENO_DIR)
pheno_var <- read.csv("pheno_meanvar.txt", sep="\t")

# get heritability
setwd(LDSC_DIR)
ldsc_df <- read.csv("ldsc_results.txt", sep="\t")
ldsc_df <- ldsc_df[c(1,2,3,4,5)] ; ldsc_df <- ldsc_df[ldsc_df$Sex != 'both_sex',]
ldsc_se <- dcast(ldsc_df, Code + Phenotype ~ Sex, value.var = c("h2.std.error"))
ldsc_df <- dcast(ldsc_df, Code ~ Sex, value.var = c("Heritability"))
colnames(ldsc_df) <- c("pheno","f_h2","m_h2") ; colnames(ldsc_se) <- c("pheno","pheno_name","f_h2_se","m_h2_se")
ldsc_df <- ldsc_df %>%
  mutate(f_h2e = 1-f_h2) %>%
  mutate(m_h2e = 1-m_h2)

# merge variance and heritability dataframe
df <- merge(merge(pheno_var, ldsc_df, by='pheno'), ldsc_se, by='pheno')
# estimate genetic and environmental variance, and their se
df <- df %>%
  mutate(geno_var_m = m_h2*m_var) %>% mutate(geno_var_f = f_h2*f_var) %>%
  mutate(env_var_m = m_h2e*m_var) %>% mutate(env_var_f = f_h2e*f_var) %>%
  mutate(geno_se_m = sqrt((m_varse^2*m_h2_se^2)+(m_h2^2*m_varse^2)+(m_var^2*m_h2_se^2))) %>% 
  mutate(geno_se_f = sqrt((f_varse^2*f_h2_se^2)+(f_h2^2*f_varse^2)+(f_var^2*f_h2_se^2))) %>%
  mutate(env_se_m = sqrt((m_varse^2*m_h2_se^2)+(m_h2e^2*m_varse^2)+(m_var^2*m_h2_se^2))) %>% 
  mutate(env_se_f = sqrt((f_varse^2*f_h2_se^2)+(f_h2e^2*f_varse^2)+(f_var^2*f_h2_se^2))) %>%
  # calculate ratio
  mutate(geno_var_ratio = geno_var_m/geno_var_f) %>% 
  mutate(env_var_ratio = env_var_m/env_var_f) %>%
  # estimate standard error of ratio
  mutate(geno_var_ratio_se = 
           (geno_var_ratio)*sqrt((geno_se_m^2/geno_var_m^2)+(geno_se_f^2/geno_var_f^2)) ) %>%
  mutate(env_var_ratio_se = 
           (env_var_ratio)*sqrt((env_se_m^2/env_var_m^2)+(env_se_f^2/env_var_f^2))) %>%
  select(c(1,12,23,24,25,26))

# 90% confidence interval
x <- df$geno_var_ratio; x.se <- df$geno_var_ratio_se
y <- df$env_var_ratio; y.se <- df$env_var_ratio_se
se_ratio <- (x/y)*sqrt((x.se^2/x^2)+(y.se^2/y^2))
df$ratio_se <- se_ratio
CI <- df$ratio_se * 1.645    
CI_side <- CI/sqrt(2)
df <- df %>%
  mutate(ci_x1 = geno_var_ratio - CI_side) %>%
  mutate(ci_x2 = geno_var_ratio + CI_side) %>%
  mutate(ci_y1 = env_var_ratio - CI_side) %>%
  mutate(ci_y2 = env_var_ratio + CI_side) 

############################################
# SCATTER PLOT
# 90% CI
pcolor = '#2b62d9'
# determine what phenotypes deviate from 1:1 relationship 
df <- mutate(df, on = ifelse(abs(env_var_ratio-geno_var_ratio) <= (ci_y2-ci_y1) | 
                               abs(geno_var_ratio-env_var_ratio) <= =(ci_x2-ci_x2), 2,1))
# split dataframe by outliers
df1 <- df[! df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R', 'testosterone'),]
df1_1 <- df1[(df1$env_var_ratio/df1$geno_var_ratio) >=1,]
df1_2 <- df1[df1$env_var_ratio/df1$geno_var_ratio <1,]
df2 <- df[df$pheno %in% c('arm_fatfree_mass_L', 'arm_fatfree_mass_R'),]
df3 <- df[df$pheno == 'testosterone',]
empty_df <- data.frame()
anno_size = 2.8

pdf(file="genetic_environment_amplification.pdf", width=7, height=5)
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 25) +
  scale_y_continuous(breaks=seq(23.5,24.5,0.5), limits=c(23,25), expand=c(0,0)) + 
  scale_x_continuous(breaks=seq(0, 2, 0.5), limits=c(-0.1,2.2), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y = element_text(size=10), axis.title=element_blank(), 
        axis.ticks.x=element_blank(), axis.line.x = element_blank())

g2 <- ggplot(empty_df) + geom_point() +
  geom_hline(yintercept = 25) +
  scale_x_continuous(breaks=seq(3.1,3.7,0.2), limits=c(3,3.8), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(23.5,24.5,0.5), limits=c(23,25), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
        axis.line = element_blank())

g3 <- ggplot(data=df3, aes(x= geno_var_ratio, y= env_var_ratio)) +
  geom_hline(yintercept = 25) + geom_vline(xintercept = 79) +
  geom_point(size=2, alpha=0.5) + 
  geom_segment(aes(x=ci_x2, y=ci_y1, xend=ci_x1, yend=ci_y2)) +
  scale_y_continuous(breaks=seq(23.5,24.5,0.5), limits=c(23,25), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(77,78,1), limits=c(76.5,79), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text(aes(label=pheno_name), size=anno_size,x=75,y=24) +
  coord_cartesian(clip = "off")

g4 <- ggplot(data=df1, aes(x= geno_var_ratio, y= env_var_ratio, color=on)) +
  geom_point(size=2, alpha=0.5) + 
  geom_abline(slope = 1, intercept = 0, color=pcolor, size=0.8) +
  geom_segment(aes(x=ci_x2, y=ci_y1, xend=ci_x1, yend=ci_y2)) +
  scale_y_continuous(breaks=seq(0,3,1), limits=c(-0.24,3.4), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 2, 0.5), limits=c(-0.1,2.2), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=10), axis.title=element_blank(), legend.position = "none") +
  scale_color_continuous(low = "black", high = pcolor) +
  geom_text_repel(data=(df1_1), aes(label=pheno_name), size=anno_size, max.overlaps=Inf, segment.alpha=0.5, 
                  ylim=c(1.5,NA), box.padding=1.5, segment.linetype=2) +
  geom_text_repel(data=(df1_2), aes(label=pheno_name), size=anno_size, max.overlaps=Inf, segment.alpha=0.5, 
                  ylim=c(NA,0.9), box.padding=1.5, segment.linetype=2)

g5 <- ggplot(data=df2, aes(x= geno_var_ratio, y= env_var_ratio)) +
  geom_point(size=2, alpha=0.5, color="black") + 
  geom_abline(slope = 1, intercept = 0, color=pcolor, size= 0.8) +
  geom_segment(aes(x=ci_x2, y=ci_y1, xend=ci_x1, yend=ci_y2)) +
  scale_x_continuous(breaks=seq(3.1,3.7,0.2), limits=c(3,3.8), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,3,1), limits=c(-0.24,3.4), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.y=element_blank(), axis.text.x = element_text(size=10), axis.title=element_blank(), 
        axis.ticks.y=element_blank(), axis.line.y = element_blank(), legend.position = "none") +
  geom_text_repel(aes(label=pheno_name), size=anno_size, max.overlaps=Inf, 
                  segment.alpha=0, segment.linetype=2, color="black")

g6 <- ggplot(empty_df) + geom_point() +
  geom_vline(xintercept = 79) +
  scale_x_continuous(breaks=seq(77,78,1), limits=c(76.5,79), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,3,1), limits=c(-0.24,3.4), expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.y=element_blank(), axis.text.x = element_text(size=10), axis.title=element_blank(),
        axis.ticks.y=element_blank(), axis.line.y = element_blank())

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g4)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, g3, gB, g5, g6, ncol=3, nrow=2, widths=c(6,2,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Genetic Variance Male:Female Ratio", size=12),
                left = text_grob("Environmental Variance Male:Female Ratio", size=12, rot=90))
dev.off()

