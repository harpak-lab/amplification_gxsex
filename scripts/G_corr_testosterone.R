#!/usr/bin/env Rscript

### REGRESSION OF GENETIC EFFECT ON TESTOSTERONE LEVEL - OVERALL ###
source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR', 'PHENO_DIR'))])
library("optparse", lib.loc=R_LIB)
library(ggplot2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)
library(ggrepel, lib.loc=R_LIB)

# phenotype names
setwd(PHENO_DIR)
pheno_list <- read.csv("pheno_names.txt", sep="\t")

# get lm results from each testosterone level bin
bin_fun <- function(data, n, sex) {
  intervals = seq(0,nrow(data),nrow(data)/n)
  cuts <- cut(1:nrow(data), breaks = intervals)
  results <- NULL
  for (i in 1:n) {
    # linear regression
    bin <- data[cuts == levels(cuts)[i],]
    model <- lm(paste0("pheno ~ SCORE"), data = bin)
    beta <- model$coefficients[2]
    stderror <- summary(model)$coefficients[2,2]
    T_mean <- mean(bin$testosterone)
    results <- rbind(results, data.frame(Testosterone=T_mean, Beta=beta, Error=stderror, Sex=sex))
  }
  return(results)
}

get_corr <- function(data) {
  # correlation between Beta and testosterone
  corrs_row <- NULL
  for (sex in c("male", "female")) {
    data_sub <- data[data$Sex == sex,]
    corr <- cor.test(data_sub$Testosterone, data_sub$Beta)
    corr_est <- corr$estimate
    corr_err <- corr$conf.int[2] - corr_est
    corrs_row <- rbind(corrs_row, data.frame(Pheno=pheno, Est=corr_est, Err=corr_err, Sex=sex))
  }
  corrs_row$est_diff <- abs(corrs_row[1,2] - corrs_row[2,2])
  corrs_row$err_sum <- abs(corrs_row[1,3] + corrs_row[2,3])
  return(corrs_row)
}

corrs_result <- NULL
for (pheno in pheno_list$Code) {
  # phenotype value, testosterone level, and sex
  setwd(PHENO_DIR)
  df_testosterone <- read.csv("pheno_testosterone.txt", sep="\t", colClasses = c("NULL","integer","numeric"))
  df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
  df_sex <- read.csv("sex_ids.txt", sep="\t")
  colnames(df_pheno) <- c('IID','pheno')

  # both-sex PGS scores
  setwd(file.path(GWAS_DIR,pheno))
  file_name <- list.files(pattern="both_sex_additive_")
  df_both <- read.csv(file_name,sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))

  # merge dataframes
  df <- merge(merge(df_testosterone, df_sex, by='IID'), df_pheno, by='IID')
  df <- merge(df,df_both,by='IID')

  # order by testosterone
  df <- df[order(df$testosterone),]
  # label then split by sex
  df$sex[df$sex == 1] <- 'male'
  df$sex[df$sex == 0] <- 'female'
  df_m <- df[df$sex == 'male',]
  df_f <- df[df$sex == 'female',]

  # call function to get betas for each bin
  m_results <- bin_fun(df_m,10,'male')
  f_results <- bin_fun(df_f,10,'female')
  results <- rbind(m_results, f_results)

  # correlation between Beta and testosterone
  corrs_result <- rbind(corrs_result, get_corr(results))
}

# save correlation result file
setwd(GWAS_DIR)
write.table(corrs_result, file="pgs_testosterone_corr.txt", sep="\t", row.names=FALSE)

# merge correlation values with phenotype names
corrs_result <- corrs_result[order(corrs_result$Pheno),]
corrs_result <- merge(corrs_result, pheno_list, by.x='Pheno', by.y='Code')

#####################################################
# SCATTER PLOT
# 90% confidence interval (mean +/- 1.645*SE)
pdf(file="G_testosterone_corr.pdf", width=4.5, height=6)
rects <- data.frame(xstart = seq(0.5,26.5,1), xend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
rects <- rects[1:26,]
plot <- ggplot(corrs_result, aes(x= reorder(Phenotype, est_diff), y=Est, color=Sex)) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size=2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=Est-(1.645*Err), ymax=Est+(1.645*Err)), alpha= 0.6, width=0.5, position=position_dodge(width=0.7)) +
  geom_rect(data=rects, aes(xmin=xstart,xmax=xend,ymin=-1.7,ymax=1.15), 
            inherit.aes = FALSE, alpha=0.2, fill = c(rep(c("grey","white"),13))) +
  labs(y="R - Correlation between Testosterone Level and \nEffect of Polygenic Score on Phenotype") +
  scale_y_continuous(expand=c(0,0), breaks=seq(-1.5,1,0.5)) +
  theme_classic() + 
  theme(axis.text = element_text(size=10), axis.title.x = element_blank(), axis.title.y = element_text(size=11),
            plot.title=element_blank(), legend.position = "none") +
  scale_color_manual(values=c("#d67629","#207335")) +
  coord_flip()

annotate_figure(plot, 
                bottom = text_grob("R - Correlation between Testosterone Level and \nEffect of Polygenic Score on Phenotype", size=12))
dev.off()

