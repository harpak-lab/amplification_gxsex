#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','PHENO_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(dplyr, lib.loc=R_LIB)
library(gridExtra, lib.loc=R_LIB)
library(grid, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)

# get phenotype names
setwd(PHENO_DIR)
pheno_names <- read.csv("pheno_names.txt", sep="\t", colClasses=c(rep("character",2), "NULL"))

#### PHENOTYPE MEAN AND VARIANCE ####
set.seed(1)
pheno_stats <- NULL
pheno_list <- pheno_names$Code    # obtain list of phenotypes
df_sex <- read.csv("sex_ids.txt", sep="\t")   # obtain list of sex
for (pheno in pheno_list) {
   df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
   df_pheno <- merge(df_pheno, df_sex, by='IID')    # merge with sex column
   m <- df_pheno[df_pheno$sex == 1, 2] ; f <- df_pheno[df_pheno$sex == 0, 2]    # split by sex
   # estimate variance and mean by sex
   male_var <- var(m)
   female_var <- var(f)
   male_mean <- mean(m)
   female_mean <- mean(f)

   # bootstrap to obtain variance standard errors
   male_varse <- mean(replicate(100, sd(sample(m,replace=T))/sqrt(length(m))))
   female_varse <- mean(replicate(100, sd(sample(f,replace=T))/sqrt(length(f))))

   pheno_stats <- rbind(pheno_stats, data.frame(pheno=pheno, m_mean=male_mean, f_mean=female_mean, m_var=male_var, f_var=female_var,
                                            m_varse=male_varse, f_varse=female_varse))
}
# write as table
write.table(pheno_stats, file="pheno_meanvar.txt", sep="\t", row.names=FALSE, quote=FALSE)

# merge phenotype names with phenotype variance and mean estimates
df <- merge(pheno_names, pheno_var, by.x="Code", by.y="pheno")

# estimate  M:F phenotype mean ratio and variance ratio
df <- df %>%
  mutate(mean_ratio = m_mean / f_mean) %>%
  mutate(var_ratio = m_var / f_var) %>%
  mutate(mean_diff = m_mean - f_mean)

#### PLOT ####
# remove outliers
df1 <- df[!df$Code %in% c("testosterone","wth_bmi_adj"),]
df2 <- df[df$Code == "testosterone",]
df3 <- df[df$Code == "wth_bmi_adj",]
empty_df <- data.frame()
# text sizes
anno_size=2.8
axis_text=9

# scatter plot
# row1
pdf(file="FIG4a_phenovar_phenomean.pdf", width=5.5, height=3.5)
g1 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  scale_x_continuous(breaks=c(-1.2,-1.1), limits=c(-1.25,-1.1), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=axis_text),
        axis.ticks.x=element_blank(), axis.line.x.bottom = element_blank())

g2 <- ggplot(empty_df) + geom_point() + 
  geom_hline(yintercept = 27.6) +
  geom_vline(xintercept = 1, linetype="dashed", alpha = 0.5) +
  scale_x_continuous(trans="log10", breaks=seq(0.6,1.5,0.1), limits=c(0.55,1.7), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(), axis.line= element_blank())

g3 <- ggplot(data=df2, aes(x= mean_ratio, y= var_ratio)) +
  geom_point(size=2, color='#2b62d9') + 
  geom_hline(yintercept = 27.6) + geom_vline(xintercept = 11) +
  scale_x_continuous(breaks=c(10.4,10.8), limits=c(10,11), expand=c(0,0)) +
  scale_y_continuous(breaks=c(27.2,27.4), limits=c(27,27.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        axis.line=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

# row2
g4 <- ggplot(data=df3, aes(x= mean_ratio, y= var_ratio)) + 
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_point(size=2, color='#2b62d9') + 
  scale_x_continuous(breaks=c(-1.2,-1.1), limits=c(-1.25,-1.0), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3,seq(0.5,3.5,0.25)), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text)) +
  geom_text_repel(aes(label=Phenotype), size=anno_size)

g5 <- ggplot(data=df1, aes(x= mean_ratio, y= var_ratio)) +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_smooth(method="lm", color="gray", size=0.8, se= F, fullrange=T) +
  geom_point(size=2, color='#2b62d9') +
  scale_x_continuous(trans="log10", breaks=seq(0.6,1.5,0.1), limits=c(0.55,1.7), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic2() +
  theme(axis.title=element_blank(), legend.position = "none", axis.text.x = element_text(size=axis_text),
        axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_text_repel(aes(label=Phenotype), size=anno_size, max.overlaps = Inf, box.padding=0.7)

g6 <- ggplot(empty_df) + geom_point() +
  geom_hline(yintercept = 1, linetype="dashed", alpha = 0.5) +
  geom_vline(xintercept = 11) +
  geom_smooth(method="lm", data=df1, aes(mean_ratio, var_ratio), inherit.aes = F, color="gray", size=0.8, se= F, fullrange=T) +
  scale_x_continuous(breaks=c(10.4,10.8), limits=c(10,11), expand=c(0,0)) +
  scale_y_continuous(trans="log10", breaks = c(0.3, 1, 2, 3), limits=c(0.28,3.6), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title=element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=axis_text),
        axis.ticks.y=element_blank(), axis.line.y = element_blank())

gA <- ggplotGrob(g1) ; gB <- ggplotGrob(g4)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth) ; gB$widths[2:5] <- as.list(maxWidth)

plot <- grid.arrange(gA, g2, g3, gB, g5, g6, ncol=3, nrow=2, widths=c(1.5,5,1), heights=c(1,4))
annotate_figure(plot, 
                bottom = text_grob("Ratio of Male to Female Phenotypic Mean", size=11),             # x axis title
                left = text_grob("Ratio of Male to Female Phenotypic Variance", size=11, rot=90))   # y axis title
dev.off()
