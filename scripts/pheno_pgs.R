#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR','PHENO_DIR'))])
library(ggplot2, lib.loc=R_LIB)
library(ggpubr, lib.loc=R_LIB)
library(reshape2, lib.loc=R_LIB)
library(ggsci, lib.loc=R_LIB)
library(gridExtra, lib.loc=R_LIB)
library(grid, lib.loc=R_LIB)
library(dplyr, lib.loc=R_LIB)
library(boot, lib.loc=R_LIB)

# 1 male; 0 female
set.seed(1)

# phenotype names
setwd(PHENO_DIR)
pheno_names <- read.csv("pheno_names.txt", sep="\t")

# bootstrap correlation between phenotype and pgs
bootstrap_r2 <- function(df, v) {
  df_sample <- df[v,]
  model <- lm(pheno ~ SCORE, df_sample)
  bootstrap_R2 <- summary(model)$r.squared
  return(bootstrap_R2)
}

# function for retrieving regression coefficients and values
get_lm_coefficients <- function(df, sex) {
  results <- NULL
  for (i in 0:1) {
    sex_df <- df[df$sex == i,]
    model <- lm(pheno ~ SCORE, sex_df); model <- summary(model)   # lm regression phenotype on pgs
    # obtain slope, error, intercept, r2, and r2 se
    b <- model$coefficients[2]
    err <- model$coefficients[2,2]
    y_inter <- model$coefficients[1]
    boot_r2 <- boot(sex_df, bootstrap_r2, R=100)    # bootstrap (100)
    r2.se <- summary(boot_r2)$bootSE
    r2 <- model$r.squared
    
    results <- rbind(results, data.frame(Phenotype = pheno, Beta=b, Err=err, Y0=y_inter, R2=r2, R2.Err=r2.se, Sex=i))
  }
  return( results )
}

# get phenotype means from each bin
bin_fun <- function(data, n) {
  results <- NULL
  for (j in 0:1) {
    sex_df <- data[data$sex == j,]
    # split by bins
    intervals = seq(0,nrow(sex_df),nrow(sex_df)/n)   
    cuts <- cut(1:nrow(sex_df), breaks = intervals)
    for (i in 1:n) {
      # bin
      bin <- sex_df[cuts == levels(cuts)[i],]
      # get mean
      pheno_mean <- mean(bin$pheno)
      pheno_se <- sd(bin$pheno) / sqrt(nrow(bin))
      pgs_mean <- mean(bin$SCORE)
      results <- rbind(results, data.frame(Pheno=pheno_mean, Pheno_SE=pheno_se, Score=pgs_mean, Sex=j))
    }
  }
  return(results)
}

### PLOT FOR EACH PHENOTYPE ###
slope_list <- NULL
for (i in 1:nrow(pheno_names)) {

# phenotype name and unit
pheno <- pheno_names$Code[i]
pheno_title <- pheno_names$Phenotype[i]
unit <- pheno_names$Unit[i]
print(pheno)

# load phenotype files
setwd(PHENO_DIR)
df_pheno <- read.csv(paste0("pheno_",pheno,".txt"), sep="\t", colClasses = c("NULL","integer","numeric"))
df_sex <- read.csv("sex_ids.txt", sep="\t")
df_pheno <- merge(df_pheno, df_sex, by="IID")
colnames(df_pheno) <- c("IID", "pheno","sex")

# load pgs file
setwd(file.path(GWAS_DIR,pheno))
file_name <- list.files(pattern="male_additive_")[2]
df_male <- read.csv(file_name, sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))
file_name <- list.files(pattern="female_additive_")
df_female <- read.csv(file_name, sep="", colClasses= c("NULL","integer",rep("NULL",3),"numeric"))

# merge phenotype and pgs df
df_m <- merge(df_pheno, df_male, by="IID")
df_f <- merge(df_pheno, df_female, by="IID")

# order by pgs score
df_m <- df_m[order(df_m$SCORE),]
df_f <- df_f[order(df_f$SCORE),]

# convert to score sd
df_m$SCORE <- df_m$SCORE / sd(df_m$SCORE)
df_f$SCORE <- df_f$SCORE / sd(df_f$SCORE)

# get phenotype mean for each of 10 bins
m_results <- bin_fun(df_m, 10)
f_results <- bin_fun(df_f, 10)

# lm regression of phenotype over pgs
m_lm <-  get_lm_coefficients(df_m, "male")
f_lm <-  get_lm_coefficients(df_f, "female")

# add to slope list
slope_list <- rbind(rbind(slope_list, m_lm[-c(4)]), f_lm[-c(4)])

### PLOT ###
vdown = 1
label_m <- m_lm %>%
  mutate(Beta = sprintf("%.2f",Beta)) %>%
  mutate(Err = sprintf("%.2f",Err)) %>%
  mutate(R2 = sprintf("%.1f",R2*100)) %>%
  mutate(R2.Err = sprintf("%.1f",R2.Err*100))
label_f <- f_lm %>%
  mutate(Beta = sprintf("%.2f",Beta)) %>%
  mutate(Err = sprintf("%.2f",Err)) %>%
  mutate(R2 = sprintf("%.1f",R2*100)) %>%
  mutate(R2.Err = sprintf("%.1f",R2.Err*100))

pdf(file=paste0(pheno,"_pheno_pgs.pdf"), width=2, height=3)
p1 <- ggplot(data = m_results, aes(x=Score, y=Pheno, color=as.character(Sex))) +
  geom_point() +
  geom_abline(data=m_lm, aes(slope=Beta, intercept = Y0, color=as.character(Sex))) +
  theme_classic() +
  xlab("Male PGS SD") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), 
        axis.text = element_text(size=8.5), legend.position="none") + 
  scale_color_manual(values=c("#d67629", "#207335")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=0.05, vjust=vdown, size=3,color="#207335",
           label= bquote(slope == .(label_m$Beta[2])*"\u00b1"*.(label_m$Err[2])*", ")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=0.05, vjust=vdown+1, size=3,color="#207335",
           label=bquote( R^2 == .(label_m$R2[2])*"%\u00b1"*.(label_m$R2.Err[2])*"%" )) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=-1.1, vjust=vdown+9.8, size=3,color="#d67629",
           label= bquote(slope == .(label_m$Beta[1])*"\u00b1"*.(label_m$Err[1])*", ")) +
  annotate("text", x=min(m_results$Score), y=max(m_results$Pheno), hjust=-1.3, vjust=vdown+10.3, size=3,color="#d67629",
           label=bquote( R^2 == .(label_m$R2[1])*"%\u00b1"*.(label_m$R2.Err[1])*"%" )) 

p2 <- ggplot(data = f_results, aes(x=Score, y=Pheno, color=as.character(Sex))) +
  geom_point() +
  geom_abline(data=f_lm, aes(slope=Beta, intercept = Y0, color=as.character(Sex))) +
  theme_classic() +
  xlab("Female PGS SD") +
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=10), 
        axis.text = element_text(size=8.5), legend.position="none") + 
  scale_color_manual(values=c("#d67629", "#207335")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=0.05, vjust=vdown, size=3,color="#207335",
         label= bquote(slope == .(label_f$Beta[2])*"\u00b1"*.(label_f$Err[2])*", ")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=0.05, vjust=vdown+1, size=3,color="#207335",
           label=bquote( R^2 == .(label_f$R2[2])*"%\u00b1"*.(label_f$R2.Err[2])*"%" )) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=-1.1, vjust=vdown+9.8, size=3,color="#d67629",
           label= bquote(slope == .(label_f$Beta[1])*"\u00b1"*.(label_f$Err[1])*", ")) +
  annotate("text", x=min(f_results$Score), y=max(f_results$Pheno), hjust=-1.3, vjust=vdown+10.3, size=3,color="#d67629",
           label=bquote( R^2 == .(label_f$R2[1])*"%\u00b1"*.(label_f$R2.Err[1])*"%" )) 

p <- grid.arrange(p1, p2, ncol=1, nrow=2,
                  left=textGrob(paste0(pheno_title, " (",unit,")"), gp=gpar(fontsize=10), rot=90))
print(p)
dev.off()

}

# write table with lm results for every phenotype
pheno <- rep(pheno_names$Code,each=2)
pheno_title <- rep(pheno_names$Phenotype,each=2)
slope_list <- cbind(pheno, pheno_title, slope_list)
setwd(GWAS_DIR)
write.table(slope_list, file = "sexspecific_pheno_pgs_lm.txt", sep="\t", row.names = F, quote = F)
  


