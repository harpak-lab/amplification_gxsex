#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','LDSC_FILE'))])
library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("reshape2", lib.loc=R_LIB)
library("matrixStats", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("ggsci", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)

# get covariance matrice names
setwd(GWAS_DIR)
df_names <- read.csv("matrice_names.txt", sep="\t")
mix <- df_names$Name; remove(df_names)

### SIMULATION 
# get simulation mixture weight dataframe
get_df <- function(snps,h,e_ratio) {
  df <- read.csv(paste0(snps,"_",h,"_",e_ratio,".txt"), sep="\t")
  df <- data.frame(cbind(mix, df$x))   # combine names with dataframe

  # split matrice names
  df <- df %>%
    rename(Name=mix, Mean=V2) %>%
    separate(Name, c("sex","correlation","effect"), sep="[_]", fill="right") %>%
    mutate(effect = paste0(sex, effect))
  return(df)
}

# order magnitude into factors, order by correlation and magnitude
prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}

##################################
envlist <- c("1","1.5","5")

for (snp in c(100,1000,10000)) {

for (h in c("0.05","0.1","0.5")) {

for (l in envlist) {
  setwd(GWAS_DIR)
  df_values <- get_df(snp, h, l)
  ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4)])
  null <- prepare_df(df_values[1,c(2,3,4)])
  null$Mean <- as.numeric(null$Mean) ; ave$Mean <- as.numeric(ave$Mean)
  if (l == list[1]) {
    df_ave <- ave ; df_null <- null
  } else {
    df_ave <- cbind(df_ave, ave$Mean) 
    df_null <- rbind(df_null, null) 
  }
}

# rename dataframe columns
colnames(df_ave) <- c('correlation','effect',sprintf("Environmental Variance Ratio: %s",list))

# filter out weights that are 0 throughout snp numbers
df_ave <- df_ave %>% 
  melt(id.vars=c('correlation','effect')) %>%
  mutate(label = ifelse(value<0.0005,"0%",sprintf("%.1f%%", (value*100)) )) %>%
  mutate(text_color = ifelse(value<0.0005, "a", "b"))

df_null$parameter <- list

effect_labels <- c('female-\nspecific','female x3', 'female x2', 'female x1.5','equal','male x1.5','male x2','male x3','male-\nspecific')

## PLOT
# plot name 
setwd(GWAS_DIR)
pname <- paste0("env_s",snp,"_h",h)

# big 8 x 7
pdf(file=paste0(pname,".pdf"), width=7, height=8)
big <- ggplot(df_ave, aes(x= effect, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=label, color=text_color), size= 3.3) +
  scale_y_continuous(breaks=seq(-1,1,0.25)) +
  scale_x_discrete(labels= effect_labels) +
  theme_classic() +
  theme(axis.text=element_text(size=10), legend.position = "none", 
        plot.title = element_text(size=14), axis.title = element_text(size=12),
        axis.line = element_blank(), 
        strip.text = element_text(hjust=0.1, size=11), strip.background = element_rect(linetype = 0)) +
  labs(title="Simulated Weights on Hypothesis Matrices", x="Magnitude", y = "Correlation") +
  scale_color_manual(values = c("gray70","black")) +
  scale_fill_gradient(low="gray98",high="#829ed9") + 
  facet_wrap(~variable,ncol=1) 
print(big)
dev.off()

# small 3 x 3
pdf(file=paste0(pname,"_noeffect.pdf"), width=3, height=3)
small <- ggplot(df_null, aes(x= 0, y= 0, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1, fill="#829ed9") +
  geom_text(aes(label=sprintf("%.0f%%", (Mean*100))), color= "black", size= 3.4) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title="Weight of No Effect Matrice") +
  theme_void() +
  theme(axis.text=element_blank(), axis.title=element_blank(), legend.position = "none",
        axis.ticks = element_blank(), plot.title = element_text(size=10), plot.margin = margin(10,20,10,20)) +
  scale_fill_continuous(limits=c(0,100)) +
  facet_wrap(~parameter, ncol=1) +
  labs(caption= paste0("# causal SNPs = ",snp,"\n",
                       "heritability = ", h, "\n"
                       ))
print(small)
dev.off()

}
}

print(paste0("Files (.pdf) in directory: " GWAS_DIR))

