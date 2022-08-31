#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)

setwd(GWAS_DIR)

# read in pgs file 
df <- read.csv("pgs_linear_results_five.txt", sep="\t")
# split df into male and female
df <- df[c(-3,-4,-7,-8,-11,-12)]
df$Sex[df$Sex == 'male'] <- 'm'
df$Sex[df$Sex == 'female'] <- 'f'

df_combined <- df[df$Sex == 'combined',]
df_sep <- df[df$Sex != 'combined',]
df_sep$Phenotype <- paste0(df_sep$Phenotype,"_",df_sep$Sex)

melt_out <- function(df) {
  df1 <- melt(df[c(1,2,3,5,7)], id.vars = c('Phenotype','Sex'),
             variable.name = "r2_var", value.name = "r2_val")
  df2 <- melt(df[c(1,2,4,6,8)], id.vars = c('Phenotype','Sex'),
             variable.name = "r2_var", value.name = "r2_se_val")
  df2$r2_var <- df1$r2_var
  df <- merge(df1,df2, by=c('Phenotype','Sex','r2_var'))
  df <- df[order(df$Phenotype, decreasing = FALSE),]
  return(df)
}

df_sep <- melt_out(df_sep)
df_combined <- melt_out(df_combined)

## combined plot, 3SE
pdf(file="pgs_comparison_five_combined.pdf",width=4,height=7)
ggplot(data=df_combined, aes(x=Phenotype, y=r2_val, fill=factor(r2_var))) +
  geom_bar(position='dodge',stat='identity') +
  geom_errorbar(aes(ymin=r2_val-3*r2_se_val, ymax=r2_val+3*r2_se_val), alpha= 0.8, 
                show.legend = FALSE, position='dodge', stat='identity', size=0.2) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title='PGS Comparison', subtitle = "Combined Models", y="Incremental R2", fill="Model") +
  coord_flip() + 
  theme_classic() + 
  scale_fill_manual(values = c("#b0464f", "#d1b724", "#2b62d9")) +
  theme(axis.text = element_text(size=9), axis.title = element_text(size=11), plot.title=element_text(size=14),
      legend.position = "none") +
  annotate("text", x = 25, y=0.105, label = "Sex-specific covariance aware", hjust = 1, color="#2b62d9", size=3.4 ) +
  annotate("text", x = 24, y=0.105, label = "Sex-specific additive", hjust = 1, color="#d1b724", size=3.4) +
  annotate("text", x = 23, y=0.105, label = "Additive", hjust = 1, color="#b0464f", size=3.4 )
dev.off()


## format df_combined for analysis
df_combined = dcast(df_combined, Phenotype ~ r2_var, value.var="r2_val")
write.table(df_combined, file = "pgs_combined_r2.txt", sep="\t", row.names=FALSE, quote = FALSE)
