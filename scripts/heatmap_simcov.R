#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','SIM_COV','GWAS_DIR'))])
library("optparse", lib.loc=R_LIB)
require("dplyr", lib.loc=R_LIB)
require("tidyr", lib.loc=R_LIB)
require("reshape2", lib.loc=R_LIB)
require("matrixStats", lib.loc=R_LIB)
require("ggplot2", lib.loc=R_LIB)
require("ggsci", lib.loc=R_LIB)
require("ggpubr", lib.loc=R_LIB)
require("gridExtra", lib.loc=R_LIB)

option_list = list(
  make_option(c("-n","--name"), type="character", default="",help="file name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
name <- opt$name; print(name)

# get list of matrices
setwd(GWAS_DIR)
m_names <- read.csv("matrice_names.txt", sep="\t")

setwd(SIM_COV)
df <- read.csv(paste0(name, ".txt"), sep="\t")
mean <- rowMeans(df)
se <- rowSds(as.matrix(df)) / sqrt(length(colnames(df)) - 1)

df <- data.frame(cbind(m_names[1:4], mean, se))
df$effect <- paste0(df$sex, df$magnitude)
prepare_df <- function(df) {
  df$effect <- factor(df$effect, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, effect)
  return(df)
}
df_nonnull <- prepare_df(df[2:nrow(df),c(3,5,6,7)])
df_null <- prepare_df(df[1,c(3,4,5,6,7)])

# format labels
df_nonnull <- df_nonnull %>% 
  mutate(mean = as.numeric(mean)) %>%
  mutate(mean_lab = ifelse(mean < 0.0005, "0%", sprintf("%.1f%%", round(mean*100,1)) )) %>%
  mutate(se = as.numeric(se)) %>%
  mutate(se_lab = ifelse(se < 0.0005, "0%", sprintf("%.1f%%", round(se*100,1)) )) 
df_null <- df_null %>% 
  mutate(mean = as.numeric(mean)) %>%
  mutate(mean_lab = ifelse(mean < 0.0005, "0%", sprintf("%.1f%%", round(mean*100,1)) )) %>%
  mutate(se = as.numeric(se)) %>%
  mutate(se_lab = ifelse(se < 0.0005, "0%", sprintf("%.1f%%", round(se*100,1)) )) 

# CONDENSE PLOT
nan_weight <- 1 / (1 - df_null$mean[1])
df_nonnull$mean = df_nonnull$mean * nan_weight
df_nonnull$effect <- as.character(df_nonnull$effect) ; df_nonnull$correlation <- as.numeric(df_nonnull$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_nonnull %>% filter(substr(effect,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(mean)) %>%
    as.data.frame()
  return(df_sex)
}
for (s in c('f','m','e')) {
  assign(s, group_sex(s))
}

df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2]))
colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
  c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
  c('partial,\npositive', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
  c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
  c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))
colnames(df_small) <- c('corr', 'female > male', 'female = male', 'female < male')
df_small <- melt(df_small, id.vars=c('corr'))

df_small$corr <- factor(df_small$corr, levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))
df_small$variable <- factor(df_small$variable, levels = c('female > male', 'female = male', 'female < male'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)

colnames(df_small) <- c("correlation", "magnitude", "value")

sum_mag <- df_small %>%
  group_by(magnitude) %>%
  summarise(sum= sum(value))
sum_corr <- df_small %>%
  group_by(correlation) %>%
  summarise(sum = sum(value))

# PLOT 
#png(file=paste0(name,"_small.png"), width=3.5, height=2.3, units="in", res=200)
pdf(file=paste0(name,"_small.pdf"), width=3.5, height=2.3)

small <- ggplot(df_small, aes(x=magnitude, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=sprintf("%.0f%%", (value*100))), color= "black", size= 2.5) +
  labs(title="Covariance of Genetic Effects: Compact Representation") +
  xlab("Magnitude") + ylab("Correlation") +
  theme_pubclean() +
  # theme(axis.text=element_text(size=9), axis.title = element_text(size=10, hjust=0.5),
  #        plot.title = element_text(size=11, hjust=1), legend.position = "none") +
  theme(axis.text=element_text(size=8), axis.title = element_text(size=9, hjust=0.5),
        plot.title = element_blank(), legend.position = "none") +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low="gray98",high="#829ed9") + 
  annotate("text", x=1:3, y = 0.4, size = 2.8, label = sprintf("%.0f%%", (sum_mag$sum*100))) +
  annotate("text", x = 3.65, y=1:4, size = 2.8, label = sprintf("%.0f%%", (sum_corr$sum*100))) +
  coord_cartesian(xlim=c(1,3.15), ylim=c(0.9,3.9), clip="off")
print(small)
dev.off()

print(paste0("Files (.pdf) in directory: ", SIM_COV))
