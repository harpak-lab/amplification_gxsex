#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','PHENO_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)

# phenotype names for formatting
setwd(PHENO_DIR)
n_df <- read.csv("pheno_names.txt", sep="\t")

# read in pgs file 
setwd(GWAS_DIR)
df <- read.csv("PGS_20.txt", sep="\t")
df <- merge(n_df, df, by.x="Code", by.y="phenotype")
colnames(df)

# keep columns related to incremental r2
df <- df[c(-1,-3,-4,-5,-8,-9,-12,-13,-16,-17)]

melt_out <- function(df) {
  df1 <- melt(df[c(1,2,4,6,8)], id.vars = 'Phenotype',
              variable.name = "r2_var", value.name = "r2_val")
  df2 <- melt(df[c(1,3,5,7,9)], id.vars = 'Phenotype',
              variable.name = "r2_var", value.name = "r2_se_val")
  df2$r2_var <- df1$r2_var
  df <- merge(df1,df2, by=c('Phenotype','r2_var'))
  df <- df[order(df$Phenotype, decreasing = FALSE),]
  return(df)
}

df_combined <- melt_out(df)

# order models
df_combined$r2_var <- factor(df_combined$r2_var, levels=c("ss_incr2","as_incr2","add_incr2","mash_incr2"))
df_combined <- df_combined[order(df_combined$r2_var),]


# plot, 6 SE error bars
pdf(file="PGS_20_plot.pdf", width=5, height=6)
print(
ggplot(data=df_combined, aes(x=Phenotype, y=r2_val, fill=(r2_var))) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-3*r2_se_val, ymax=r2_val+3*r2_se_val), alpha= 0.8, 
                show.legend = FALSE, position='dodge', stat='identity', size=0.2) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title='Performace of GxSex Aware Polygenic Scores', y="Incremental R2") +
  coord_flip() + 
  theme_classic() + 
  scale_fill_manual(values = c("#2B62D9", "#D1B724","#574061","#B0464F")) +
  theme(axis.text = element_text(size=9), axis.title.x = element_text(size=11),
        plot.title=element_text(size=14), axis.title.y = element_blank(),
        legend.position = "none") +
  annotate("text", x = 25, y=0.12, label = "Sex-specific covariance aware", hjust = 1, color="#B0464F", size=3.4 ) +
  annotate("text", x = 24, y=0.12, label = "Additive", hjust = 1, color="#574061", size=3.4 ) +
  annotate("text", x = 23, y=0.12, label = "Additive, sex-standardized", hjust = 1, color="#D1B724", size=3.4 ) +
  annotate("text", x = 22, y=0.12, label = "Sex-specific additive", hjust = 1, color="#2B62D9", size=3.4)
)
dev.off()

print(paste0("File, PGS_20_plot.pdf, in directory: ", GWAS_DIR))