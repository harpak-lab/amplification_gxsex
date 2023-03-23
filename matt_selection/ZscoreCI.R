#Z score confidence interval

p <- commandArgs(trailingOnly = TRUE)[1]
anclist <- commandArgs(trailingOnly = TRUE)[2]
cat("Pval Threshold:",p)

df <- read.csv(paste0("results.",p,".z.NEW.csv"),header = TRUE)
df$BOOTSE <- rep(NA,81)
df$CILOWER <- rep(NA,81)
df$CIUPPER <- rep(NA,81)

# anclist <- c("asj","fin","nfe")
traitlist <- c("calcium","albumin","arm_fatfree_mass_L","arm_fatfree_mass_R",
               "bmi","creatinine","diastolicBP_auto","eosinophil_perc",
               "FVC_best","HbA1c","height","hip_circ","IGF1","lymphocyte_perc",
               "protein_total","pulse_rate","RBC_count","SHBG","systolicBP_auto",
               "testosterone","urate","urea","waist_circ","waist_to_hip","weight",
               "whole_body_fat_mass","wth_bmi_adj")

set.seed(123)
# anclist <- c("asj")
# traitlist <- c("calcium")
for(t in traitlist){
  cat("\n\nTrait:",t)
  for(a in anclist){
    cat("\nAnc:",a)
    fi <- read.csv(paste0("./",t,"/",a,"/",a,".",t,".slopes.",p,".z.NEW.csv"),
                        header = TRUE)
    zscores <- fi$Zscore
    # print(head(zscores))
    bootz <- rep(NA,10000)
    for(i in 1:10000){
      # print(mean(sample(zscores,replace = TRUE)))
      bootz[i] <- mean(sample(zscores,replace = TRUE))
      # print(bootz[i])
    }
    # print(head(bootz))
    # hist(bootz,breaks = 20)
    se <- sd(bootz)
    print(se)
    cis <- quantile(bootz,c(0.05,0.95))
    df$BOOTSE[which(df$TRAIT == t & df$ANC == a)] <- se
    df$CILOWER[which(df$TRAIT == t & df$ANC == a)] <- cis[1]
    df$CIUPPER[which(df$TRAIT == t & df$ANC == a)] <- cis[2]
  }
}

cat("\n")

write.table(df,file = paste0("results.",p,".z.CI.csv"),
            append = FALSE,row.names = FALSE,sep = ",",quote = FALSE,
            col.names = TRUE)