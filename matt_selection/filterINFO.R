#Filter the AllChrINFO.txt data

#First, get AllChrINFO.txt as dataframe
# chr <- commandArgs(trailingOnly = TRUE)[1]
Anc <- commandArgs(trailingOnly = TRUE)[1]
splitid <- commandArgs(trailingOnly = TRUE)[2]
# chr <- paste0("chr",chr)
# print(chr)
print("Reading Table")
df <- read.table(paste0("/scratch/08312/mjm8356/data/",Anc,"/",Anc,"INFO.txt.",splitid))
# df <- read.table(paste0("/scratch/08312/mjm8356/data/",chr,"/AllAnc.37.",chr,".INFO"),header = TRUE)


# names(df) <- c("CHROM","POS","REF","ALT","AN_asj_XX","AN_fin_XX","AN_amr_XX",
#                "AN_eas_XX","AN_afr_XX","AN_nfe_XX","AN_oth_XX","AN_asj_XY",
#                "AN_fin_XY","AN_amr_XY","AN_eas_XY","AN_afr_XY","AN_nfe_XY",
#                "AN_oth_XY","AF_asj_XX","AF_fin_XX","AF_amr_XX","AF_eas_XX",
#                "AF_afr_XX","AF_nfe_XX","AF_oth_XX","AF_asj_XY","AF_fin_XY",
#                "AF_amr_XY","AF_eas_XY","AF_afr_XY","AF_nfe_XY","AF_oth_XY",
#                "AC_asj_XX","AC_fin_XX","AC_amr_XX","AC_eas_XX","AC_afr_XX",
#                "AC_nfe_XX","AC_oth_XX","AC_asj_XY","AC_fin_XY","AC_amr_XY",
#                "AC_eas_XY","AC_afr_XY","AC_nfe_XY","AC_oth_XY")

names(df) <- c("CHROM","POS","REF","ALT","AN_XX","AN_XY","AF_XX","AF_XY","AC_XX",
               "AC_XY","AF_TOT","FST","AF_DIFF")

#Filter by removing all rows where there are NAs. Write this table to a new
#file called AllChrINFO_filt.txt
# print("Filtering NA")
# df <- df[,-c(grep("sas",colnames(df)),grep("ami",colnames(df)),grep("mid",colnames(df)))]
# df2 <- df[-which(apply(df,1,function(x) any(x %in% c("?")))),]
# df2[,grepl("AF",names(df2))] <- sapply(df2[,grepl("AF",names(df2))],as.numeric)
# # write.table(df2,file = paste0("./",chr,"/AllAnc_Filt.",chr,".txt"),sep = "\t",row.names = FALSE)

boot_allele_se <- function(p,n,m = 1000){
  c <- p*2*n
  pop <- c(rep(1,c),rep(0,((2*n)-c)))
  # print(pop)
  sampDist <- rep(NA,m)
  for(i in 1:m){
    repl <- sample(pop,replace = TRUE)
    sampDist[i] <- (sum(repl))/length(repl)
  }
  # print(sampDist)
  return(sd(sampDist))
}
vecboot <- Vectorize(boot_allele_se)

print("Getting Male SE")

mse <- vecboot(df$AF_XY,(df$AN_XY/2))

print("Getting Female SE")

fse <- vecboot(df$AF_XX,(df$AN_XX/2))

print("Getting New FST Estimate")

a <- ((df$AF_XY-df$AF_XX)^2) - mse^2 - fse^2
b <- (4 * df$AF_TOT * (1-df$AF_TOT)) - mse^2 - fse^2

df$FSTNEW <- a/b

print("Saving File")

filename <- paste0("/scratch/08312/mjm8356/Revision/New_Estimator/",Anc,"/",Anc,"INFO.new.",splitid,".txt")
write.table(df,file = filename,sep = "\t",row.names = FALSE,col.names = FALSE)

#Split up every Ancestry Group into its own dataframe, then filter these
#Ancestry Groups by only loci with >1000 alleles sampled.
#Calculate FST and Allele Frequency Difference (FST/4pq)
#Remove any sites with an Na FST (which would indicate No XX or XY individuals
#had the given allele)
#Save filtered Ancestry Group dataframe as new file called [ancgroup]INFO.txt
# AncList <- c("afr","amr","asj","eas","fin","nfe","oth")
# AncList <- c("asj","fin","nfe")
# for(i in AncList){
#   print(paste("Anc",i))
#   dfi <- cbind(df2[,1:4],df2[,grepl(i,names(df2))])
#   colnames(dfi) <- sub("_.*_","_",colnames(dfi))
#   dfi$AF_TOT <- (dfi$AC_XX + dfi$AC_XY)/(dfi$AN_XX + dfi$AN_XY)
#   # dfi$FST <- ((dfi$AF_female - dfi$AF_male)^2)/(4 * dfi$AF_TOT * (1 - dfi$AF_TOT))
#   
#   mse <- vecboot(dfi$AF_XY,(dfi$AN_XY/2))
#   fse <- vecboot(dfi$AF_XX,(dfi$AN_XX/2))
#   
#   a <- ((dfi$AF_XY-dfi$AF_XX)^2) - mse^2 - fse^2
#   b <- (4 * dfi$AF_TOT * (1-dfi$AF_TOT)) - mse^2 - fse^2
#   
#   dfi$FST <- a/b
#   
#   # dfi$AF_DIFF <- ((dfi$AF_female - dfi$AF_male)/(4 * dfi$AF_TOT * (1 - dfi$AF_TOT)))^2
#   dfi <- dfi[-which(is.na(dfi$FST)),]
#   filename <- paste0("./",chr,"/",i,"INFO.",chr,".txt")
#   write.table(dfi,file = filename,sep = "\t",row.names = FALSE,col.names = FALSE)
#   
#   # if(length(which(dfi[,paste("AN","XX",sep = "_")] >= 1000)) == 0 |
#   #    length(which(dfi[,paste("AN","XY",sep = "_")] >= 1000)) == 0){next}
#   # dfi2 <- df2[-union(which(dfi[,paste("AN","XX",sep = "_")] < 1000),
#   #                   which(dfi[,paste("AN","XY",sep = "_")] < 1000)),]
#   # colnames(dfi2) <- sub("_.*_","_",colnames(dfi2))
#   # dfi2$AF_TOT <- (dfi2$AC_female + dfi2$AC_male)/(dfi2$AN_female + dfi2$AN_male)
#   # dfi2$FST <- ((dfi2$AF_female - dfi2$AF_male)^2)/(4 * dfi2$AF_TOT * (1 - dfi2$AF_TOT))
#   # dfi2$AF_DIFF <- ((dfi2$AF_female - dfi2$AF_male)/(4 * dfi2$AF_TOT * (1 - dfi2$AF_TOT)))^2
#   # dfi2 <- dfi2[-which(is.na(dfi2$FST)),]
#   # filename <- paste0("./",chr,"/",i,"INFO.",chr,".1000.txt")
#   # write.table(dfi2,file = filename,sep = "\t",row.names = FALSE,col.names = FALSE)
# }
