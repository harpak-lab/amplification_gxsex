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
  make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
  make_option(c("-n","--name"), type="character", default=NULL,help="formatted phenotype name",metavar="character"),
  make_option(c("-m","--mode"), type="character", default="",help="input '_same' if using same size subset",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno) | is.null(opt$name)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code or name", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
pheno_name <- opt$name
method <- opt$mode

# load mixture names
wd <- paste0(GWAS_DIR,"/",pheno,"/mash")
setwd(wd)
mix <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t")
mix <- mix$mix_0
p_values <- factor(c("5e-08", "1e-05", "0.05", "1"))

############################################################################################

# initialize list for null weight
null <- NULL
for ( p_val in p_values ) {

# load mixture proportions by p-value
df <- read.csv(paste0(pheno,"_",p_val,method,".txt"), sep="\t")
df <- cbind(mix, df)

# ORGANIZE
# split matrice names
df <- df %>%
  separate(mix, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
  mutate(magnitude = paste0(sex, magnitude))

# organize by correlation and magnitude
prepare_df <- function(df) {
  df$magnitude <- factor(df$magnitude, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, magnitude)
  return(df)
}

# split between null and values
df_ave <- prepare_df(df[2:nrow(df),c(2,3,4)])
df_null <- prepare_df(df[1,c(2,3,4)])
null <- append(null, df_null$x)

# CONDENSE
nan_weight <- 1 / (1 - df_null$x[1])
df_ave$x = df_ave$x * nan_weight
df_ave$magnitude <- as.character(df_ave$magnitude) ; df_ave$correlation <- as.numeric(df_ave$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_ave %>% filter(substr(magnitude,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(x)) %>%
    as.data.frame()
  return(df_sex)
}
for (s in c('f','m','e')) {
  assign(s, group_sex(s))
}

# group by correlation
df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2])); colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
  c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
  c('partial,\npositive', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
  c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
  c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))

# edit col names and data type
colnames(df_small) <- c('corr', 'female > male', 'female = male', 'female < male')
df_small <- melt(df_small, id.vars=c('corr'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)
colnames(df_small) <- c("correlation", "magnitude", "value")

# sum across correlation and magnitude groups
sum_mag <- df_small %>%
  group_by(magnitude) %>%
  summarise(sum= sum(value))
sum_corr <- df_small %>%
  group_by(correlation) %>%
  summarise(sum = sum(value))

# assign df to variable
assign(paste0(p_val,"_mag"), as.data.frame(sum_mag))
assign(paste0(p_val,"_corr"), as.data.frame(sum_corr))

}

############################################################################################
# ORGANIZE FOR PLOT
# merge p-value mixture weight tables
corr_list <- list(`0.05_corr`, `1_corr`,`1e-05_corr`,`5e-08_corr`)
corr_all <- corr_list %>% 
  reduce(full_join, by='correlation') %>%
  rename("5e-8"=sum.y.y, "1e-5"=sum.x.x, "0.05"=sum.x, "1"=sum.y) %>%
  gather(threshold, weight, 2:5) %>%
  mutate(threshold = factor(threshold, levels=c("5e-8", "1e-5", "0.05", "1")))
  
corr_all$correlation <- factor(corr_all$correlation, 
                               levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))

mag_list <- list(`0.05_mag`, `1_mag`,`1e-05_mag`,`5e-08_mag`)
mag_all <- mag_list %>% 
  reduce(full_join, by='magnitude') %>%
  rename("5e-8"=sum.y.y, "1e-5"=sum.x.x, "0.05"=sum.x, "1"=sum.y) %>%
  gather(threshold, weight, 2:5) %>%
  mutate(threshold = factor(threshold, levels=c("5e-8", "1e-5", "0.05", "1")))


# PLOT
colors <- c("#d1b724", "#563f61", "#b0464f", "#2b62d9")
pdf(file=paste0(pheno,"_corr",method,".pdf"), width=4, height=3)
corr_plot <- 
  ggplot(corr_all, aes(x=threshold, y=weight*100, fill=correlation)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(title=pheno_name, x="P-value Threshold", y="Weight", fill="Correlation") +
  scale_fill_manual(values = colors) +
  theme(axis.title=element_text(size=11), axis.text=element_text(size=9))
print(corr_plot)
dev.off()

colors <- c("#d67629", "#829ed9", "#207335")
pdf(file=paste0(pheno,"_mag",method,".pdf"), width=4, height=3)
mag_plot <- 
  ggplot(mag_all, aes(x=threshold, y=weight*100, fill=magnitude)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(title=pheno_name, x="P-value Threshold", y="Weight", fill="Magnitude") +
  scale_fill_manual(values = colors) +
  theme(axis.title=element_text(size=11), axis.text=element_text(size=9))
print(mag_plot)
dev.off()

print(paste0("Results in: ", wd))




