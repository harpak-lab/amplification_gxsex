#Analyze GWAS data

#First, get the current trait (from the bash loop)
trait <- commandArgs(trailingOnly = TRUE)[1]
print(paste("Current Trait:",trait))

#Move to the correct directory for the trait of interest
# getwd()

#Get the file containing the Male and Female Beta values and generate a data-
#frame with this data
filename <- paste0(trait,"/combined_tab.",trait,".txt")
GWAS_data <- read.table(filename,header = TRUE)

#Convert the two BETA columns to numerics (rather than character strings)
GWAS_data$MBETA <- as.numeric(GWAS_data$MBETA)
GWAS_data$FBETA <- as.numeric(GWAS_data$FBETA)
Bz <- list(GWAS_data$MBETA,GWAS_data$FBETA)
sz <- list(GWAS_data$MSE,GWAS_data$FSE)
GWAS_data$VAR <- 2*(sz[[1]]^4 + sz[[2]]^4) + 4*sz[[1]]^2*sz[[2]]^2 + 
  4*(Bz[[1]]^2*sz[[1]]^2 + Bz[[2]]^2*sz[[2]]^2) + 
  4*(sz[[1]]^2*Bz[[2]]^2 + sz[[2]]^2*Bz[[1]]^2) - 8*Bz[[1]]*Bz[[2]]*(sz[[1]]^2 + sz[[2]]^2)

GWAS_data <- GWAS_data[order(nchar(GWAS_data[,"CHROM"]),
                                   GWAS_data[,"CHROM"],
                                   GWAS_data[,"POS"]),]

AncList <- c("afr","amr","asj","eas","fin","nfe","sas","oth")
for(i in AncList){
  # print(paste0("Current Ancestry: ",i))
  pathname <- paste0("/scratch/08312/mjm8356/data/",i,"INFO.txt")
  INFO_data <- read.table(pathname,header = TRUE)
  INFOloci <- paste(INFO_data$CHROM,INFO_data$POS,
                    sep = ":")

  # print(head(s))
  # print(s$CHROM[1:5])
  # print(s$POS[1:5])
  GWASloci <- paste(GWAS_data$CHROM,GWAS_data$POS,sep = ":")
  # print(head(GWASloci))
  GWASkeep <- which(GWASloci %in% INFOloci)
  # print(length(GWASkeep))
  INFOloci[which(duplicated(INFOloci))] <- "duplicate"
  INFOkeep <- which(INFOloci %in% GWASloci[GWASkeep])
  
  p <- INFO_data$AF_TOT[GWASkeep]
  q <- 1-p
  b <- (GWAS_data$MBETA[GWASkeep] - GWAS_data$FBETA[GWASkeep])^2
  v <- ((2*p*q)^2)*b
  dfi <- data.frame(CHROM = GWAS_data$CHROM[GWASkeep],
                   POS = GWAS_data$POS[GWASkeep],
                   V = v,
                   FST = INFO_data$FST[INFOkeep],
                   varV = GWAS_data$VAR[GWASkeep])
  # print(head(df))
  tabname <- paste0("./",trait,"/",i,"/",i,".table.",trait,".txt")
  # print(tabname)
  write.table(dfi,file = tabname,sep = "\t",row.names = FALSE,append = FALSE)
}
