trait <- commandArgs(trailingOnly = TRUE)[1]
n <- commandArgs(trailingOnly = TRUE)[2]

# AncList <- c("afr","amr","asj","eas","fin","nfe","sas","oth")

f <- paste0("./results.",n,".csv")
df <- read.csv(f,header = TRUE)

a <- "UKB"

Amax <- df$MAXVSLOPE[which(df$TRAIT == trait & df$ANC == a)]

fi <- paste0("./",trait,"/",a,"/",a,".",trait,".slopes.",n,".csv")
outputz <- read.csv(fi,header = TRUE)

Avals <- outputz$slopes
SEvals <- outputz$SEs

outputz$Zscore <- Avals/SEvals
outfile <- paste0("./",trait,"/",a,"/",a,".",trait,".slopes.",n,".z.csv")
write.table(outputz,file = fi,sep = ",",quote = FALSE,
            row.names = FALSE,col.names = TRUE,append = FALSE)

rowz <- data.frame(trait,a,
                   mean(Avals,na.rm = TRUE),
                   mean(SEvals,na.rm = TRUE),
                   sd(Avals,na.rm = TRUE),
                   Amax,mean(outputz$Zscore),sd(outputz$Zscore))
resultsfile <- paste0("./results.",n,".z.csv")
write.table(rowz,file = resultsfile,append = TRUE,sep = ",",quote = FALSE,
            col.names = FALSE,row.names = FALSE)

cat("\n")


