#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("reshape2", lib.loc=R_LIB)
library("matrixStats", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("ggsci", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)

# user input phenotype
option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
  make_option(c("-n","--name"), type="character", default=NULL,help="formatted phenotype name",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno) | is.null(opt$name)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code or name", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
pheno_name <- opt$name

# load mixture proportion file
setwd(paste0(GWAS_DIR,pheno,"/mash"))
df <- read.csv(paste0(pheno,"mixprop_100_all.txt"), sep="\t")

# get mean and sem (NOT SIM)
df_values <- data.frame( Name = df$mix_0, Mean = rowMeans(df[2:101]), 
                         SE = rowSds(as.matrix(df[2:101])) / sqrt(length(colnames(df)) - 1) )
# split matrice names
df_values <- df_values %>%
  separate(Name, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
  mutate(magnitude = paste0(sex, magnitude))

prepare_df <- function(df) {
  df$magnitude <- factor(df$magnitude, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, magnitude)
  return(df)
}

# split between null and values
df_ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4,5)])
df_null <- prepare_df(df_values[1,c(2,3,4,5)])
df_null <- df_null %>% 
  mutate(Mean = as.numeric(Mean)) %>%
  mutate(mean_lab = ifelse(Mean < 0.0005, "0%", sprintf("%.1f%%", round(Mean*100,1)) )) %>%
  mutate(se_lab = ifelse(SE < 0.0005, "0%", sprintf("%.1f%%", round(SE*100,1)) ))

df_ave <- df_ave %>% 
  mutate(Mean = as.numeric(Mean)) %>%
  mutate(mean_lab = ifelse(Mean < 0.0005, "0%", sprintf("%.1f%%", round(Mean*100,1)) )) %>%
  mutate(se_lab = ifelse(SE < 0.0005, "0%", sprintf("%.1f%%", round(SE*100,1)) ))

effect_labels <-  c('female-\nspecific','female x3', 'female x2', 'female x1.5','equal','male x1.5','male x2','male x3','male-\nspecific')

##########################################################################

### LARGE HEATMAP ###
pdf(file=paste0(pheno,"_mash_large.pdf"), width=6.5, height=4.8)
big <- ggplot(df_ave, aes(x= magnitude, y= correlation, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=mean_lab), color= "black", size= 2.7, vjust=-0.1) +
  geom_text(aes(label=paste("\u00B1",se_lab)), color= "black", size= 2.2, vjust=1.5) +
  scale_y_continuous(breaks=seq(-1,1,0.25), expand=c(0,0)) +
  scale_x_discrete(labels= effect_labels) +
  labs(title="Weights on Hypothesis Covariance Matrices - ",pheno_name) +
  xlab("Magnitude") + ylab("Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=9), axis.title = element_text(size=11), plot.title = element_text(size=12), 
        legend.position = "none") +
  scale_fill_gradient(low="gray98",high="#829ed9")

small <- ggplot(df_null, aes(x= 0, y= 0, fill= Mean)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=mean_lab), color= "white", size= 2.7, vjust=-0.1) +
  geom_text(aes(label=paste("\u00B1",se_lab)), color= "white", size= 2.2, vjust=1.5) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Weight of \nNo Effect Matrix") +
  theme_pubclean() +
  theme(axis.text=element_blank(), axis.title=element_blank(), 
        axis.title.y = element_text(size=10, angle=360, vjust=0.5),
        legend.position = "none",axis.ticks = element_blank(), plot.title = element_blank()) +
  scale_fill_gradient(low="#2b62d9",high="#2b62d9")

lay <- rbind( c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(1,1,1,1,1), c(2,2,3,3,3))
p <- gridExtra::grid.arrange(big, small, ncol=1, layout_matrix=lay)
print(p)
dev.off()

### SMALL HEATMAP ###
nan_weight <- 1 / (1 - df_null$Mean[1])
df_ave$Mean = df_ave$Mean * nan_weight
df_ave$magnitude <- as.character(df_ave$magnitude) ; df_ave$correlation <- as.numeric(df_ave$correlation)

# group by sex
group_sex <- function(sex){
  df_sex <- df_ave %>% filter(substr(magnitude,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(Mean)) %>%
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
colnames(df_small) <- c('corr', 'female > male', 'female = male', 'female < male')
df_small <- melt(df_small, id.vars=c('corr'))

# format table names
df_small$corr <- factor(df_small$corr, levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))
df_small$variable <- factor(df_small$variable, levels = c('female > male', 'female = male', 'female < male'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
  arrange(corr, variable)
colnames(df_small) <- c("correlation", "magnitude", "value")

# calculate sums 
sum_mag <- df_small %>%
  group_by(magnitude) %>%
  summarise(sum= sum(value))
sum_corr <- df_small %>%
  group_by(correlation) %>%
  summarise(sum = sum(value))

pdf(file=paste0(pheno,"_mash_small.pdf"), width=3.5, height=2)
print(
ggplot(df_small, aes(x=magnitude, y= correlation, fill= value)) +
  geom_tile(color= "white", lwd= 1.5, linetype= 1) +
  geom_text(aes(label=sprintf("%.0f%%", (value*100))), color= "black", size= 2.5) +
  labs(title=pheno_name) +
  xlab("Magnitude") + ylab("Correlation") +
  theme_pubclean() +
  theme(axis.text=element_text(size=8), axis.title = element_text(size=9, hjust=0.5),
        plot.title = element_text(size=10,hjust=1), legend.position = "none") +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low="gray98",high="#829ed9") + 
  annotate("text", x=1:3, y = 0.4, size = 2.8, label = sprintf("%.0f%%", (sum_mag$sum*100))) +
  annotate("text", x = 3.65, y=1:4, size = 2.8, label = sprintf("%.0f%%", (sum_corr$sum*100))) +
  coord_cartesian(xlim=c(1,3.15), ylim=c(0.9,3.9), clip="off")
)
dev.off()