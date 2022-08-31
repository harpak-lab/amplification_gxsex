#!/usr/bin/env Rscript

### REGRESSION OF GENETIC EFFECT ON TESTOSTERONE LEVEL - INDIVIDUAL PHENOTYPE ###
source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR', 'PHENO_DIR'))])
library("optparse", lib.loc=R_LIB)
library(ggplot2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)

# user input phenotype and mode
option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL,help="phenotype name",metavar="character"),
  make_option(c("-n","--name"), type="character", default=NULL,help="formatted phenotype name",metavar="character"),
  make_option(c("-m","--mode"), type="character", default="both-sex",help="both-sex or sex-specific",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$pheno) | is.null(opt$name)) {
  print_help(opt_parser)
  stop("Missing argument for phenotype code or name", call.=FALSE)
}
pheno <- opt$pheno; print(pheno)
title <- opt$name
mode <- opt$mode; print(mode)

# get PGS scores
setwd(file.path(GWAS_DIR,pheno))
if (mode == 'both-sex') {
  # for both-sex
  file_name <- list.files(pattern="both_sex_additive_")
  df_both <- read.csv(file_name,sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
} else {
  # for sex-specific
  file_name <- list.files(pattern="male_additive_")
  df_male <- read.csv(file_name[2], sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  df_female <- read.csv(file_name[1], sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
  df_both <- merge(df_male,df_female, by="IID")
}

# get testosterone levels, phenotype values, and sex
setwd(PHENO_DIR)
df_testosterone <- read.csv("pheno_testosterone.txt", sep="\t", colClasses = c("NULL","integer","numeric"))
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
df_sex <- read.csv("sex_ids.txt", sep="\t")

# merge dataframes - testosterone, sex, pheno, pgs scores
df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
df <- merge(df,df_both,by='IID')
if (mode != 'both_sex') {
  ## for sex-specific
  df$SCORE <- ifelse(df$sex == 1, df$SCORE.x, df$SCORE.y)    
  df <- df[-c(5,6)]
}
colnames(df) <- c("IID", "testosterone", "sex", "pheno", "SCORE")

# order by testosterone
df <- df[order(df$testosterone),]
# label then split by sex
df$sex[df$sex == 1] <- 'male'
df$sex[df$sex == 0] <- 'female'
df_m <- df[df$sex == 'male',]
df_f <- df[df$sex == 'female',]

# find intercept for overlaps
f_inter = mean(df_f$testosterone) + 2*sd(df_f$testosterone)
m_inter = mean(df_m$testosterone) - 2*sd(df_m$testosterone)

# create dataframe with the overlaps
overlap <- df_m[df_m$testosterone <= m_inter,]
overlap <- rbind(overlap, df_f[df_f$testosterone >= f_inter,])
# remove overlapping from the non-overlapping df
df_m <- df_m[! df_m$IID %in% overlap$IID,]
df_f <- df_f[! df_f$IID %in% overlap$IID,]

# get lm results from each testosterone bin
bin_fun <- function(data, n, sex, pgs_sd) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  for (i in 1:n) {
    # linear regression
    bin <- data[cuts == levels(cuts)[i],]
    model <- lm(paste0("pheno ~ SCORE"), data = bin)
    beta <- model$coefficients[2]
    stderror <- summary(model)$coefficients[2,2]
    beta <- beta * pgs_sd; stderror <- stderror * pgs_sd
    T_mean <- mean(bin$testosterone)
    # correlation
    corr <- cor.test(bin$pheno, unlist(bin["SCORE"]))
    corr_est <- corr$estimate
    corr_error <- corr$conf.int[2] - corr_est
    results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Error=stderror, Corr=corr_est, Corr_err=corr_error, Sex=sex))
  }
  return(results)
}

# call function for overlap and nonoverlap
pgs_sd <- sd(df$SCORE)
m_results <- bin_fun(df_m,10,'male', pgs_sd)
f_results <- bin_fun(df_f,10,'female', pgs_sd)
overlap_results_m <- bin_fun(overlap[overlap$sex=='male',],1,'male', pgs_sd)
overlap_results_f <- bin_fun(overlap[overlap$sex=='female',],1,'female', pgs_sd)
overlap_results <- rbind(overlap_results_m, overlap_results_f)
results <- rbind(m_results, f_results)

# trendlines x range
m_over <- overlap_results$Testosterone[1]
f_over <- overlap_results$Testosterone[2]
x_max <- max(results$Testosterone)

# trendline boundary points
trend_y <- function(m_model, f_model) {
  trendline <- NULL
  m_y1 <- m_model$coefficients[1] + (m_model$coefficients[2] * m_over)
  m_y2 <- m_model$coefficients[1] + (m_model$coefficients[2] * x_max)
  f_y1 <- f_model$coefficients[1] + (f_model$coefficients[2] * 0)
  f_y2 <- f_model$coefficients[1] + (f_model$coefficients[2] * f_over)
  mp <- summary(m_model)$coefficients[2,4]
  fp <- summary(f_model)$coefficients[2,4]
  trendline <- rbind(trendline,data.frame(x1=m_over, x2=x_max, y1=m_y1, y2=m_y2, p=mp, Sex='male'))
  trendline <- rbind(trendline,data.frame(x1=0, x2=f_over, y1=f_y1, y2=f_y2, p=fp, Sex='female'))
  return(trendline)
}

# linear regression for BETA
trend <- NULL
result_sub <- results[results$Sex == 'male',]
male_model <- lm(Beta ~ Testosterone, data=result_sub)
result_sub <- results[results$Sex == 'female',]
female_model <- lm(Beta ~ Testosterone, data=result_sub)
trend <- rbind(trend, trend_y(male_model, female_model))
text_size = 2.56
sex_label_y = 5

setwd(file.path(GWAS_DIR,pheno))
pdf(file=paste0(pheno,"_pgssd_testosterone.pdf"), width=3.5, height=2.4)
g <- ggplot(results, aes(x=Testosterone, y=Beta, color=Sex)) +
  geom_point(size=1.5) +
  geom_point(data=overlap_results, aes(x=Testosterone, y=Beta), shape=1, stroke=1, size=1.5) +
  geom_errorbar(aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  geom_errorbar(data= overlap_results, aes(ymin=Beta-Error, ymax=Beta+Error), alpha= 0.4) +
  labs(title=title, x="Testosterone Level (nmol/L)", y="Effect of 1 PGS SD on Phenotype") +
  geom_segment(data=trend, aes(x=x1,xend=x2,y=y1,yend=y2)) +
  theme_classic() + 
  theme(axis.text = element_text(size=8.5), axis.title = element_text(size=9), 
        plot.title=element_text(size=11), legend.position= "none", plot.margin = margin(5,5,5,5)) +
  scale_color_manual(values=c("#d67629","#207335")) +
  stat_cor(method='pearson', p.accuracy=0.001, label.x.npc=0.55, label.y.npc=0.95, 
           show.legend=FALSE, size=text_size) +
  annotate("text", x=min(results$Testosterone), y=max(results$Beta), vjust=sex_label_y,
           label="female", color="#d67629", size=text_size) +
  annotate("text", x=min(results$Testosterone), y=max(results$Beta), vjust=sex_label_y, hjust=-3,
           label="male", color="#207335", size=text_size)
print(g)
dev.off()

