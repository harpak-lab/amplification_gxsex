# Analyze all traits for all ancestry groups

# Get the slope of a regression line for FST vs V (function of beta values and p)
# Using
# 1) A random sample of points, one point per LD block
# 2) The maximum V value for each LD block

# Run #1 100 times, getting a distribution of slopes, as well as a variance

# Add the max-V slope, mean slope, and variance of slopes to a csv
# This csv should contain info for all traits and all ancestries collected in one
# place

# Also store each vector of slopes and variances in a dataframe per trait per anc

#Analyze GWAS data

#First, get the current trait (from the bash loop)
library(tictoc)

set.seed(123)
AncList <- commandArgs(trailingOnly = TRUE)[1]
trait <- commandArgs(trailingOnly = TRUE)[2]
pstr <- commandArgs(trailingOnly = TRUE)[3]
cat("\n\nCurrent P-val Threshold:",pstr)
cat("\nCurrent Trait:",trait)

#Move to the correct directory for the trait of interest
# getwd()

#Get the file containing the Male and Female Beta values and generate a data-
#frame with this data

# print("Read in GWAS data file")
filename <- paste0("/scratch/08312/mjm8356/GWAS_results/",trait,"/pstats.",trait,".txt")
pstats <- read.table(filename,header = TRUE)

# pstats <- pstats[order(nchar(pstats[,"CHROM"]),
#                              pstats[,"CHROM"],
#                              pstats[,"POS"]),]


p <- as.numeric(pstr)
pFilt <- which(pstats$MP <= p | pstats$FP <= p)

blocks <- read.table("/scratch/08312/mjm8356/GWAS_results/Pickrell_breakpoints_EUR.bed",header = TRUE)

# print("Loop through ancestries")

for(a in AncList){
  cat("\nCurrent ancestry:",a,"\n")
  # Get the dataframe for the trait of interest
  pathname <- paste0("./",trait,"/NewFST.",a,".table.",trait,".txt")
  dfi <- read.table(pathname,header = TRUE)
  
  # Get loci which pass the p-val threshold
  ploci <- paste(pstats$CHROM[pFilt],pstats$POS[pFilt],sep = ":")
  dfiloci <- paste(dfi$CHROM,dfi$POS,sep = ":")
  # Check for overlap of p-val threshold-passing loci
  dKeep <- which(dfiloci %in% ploci)
  
  dfi <- dfi[dKeep,]
  
  # Initialize dataframes for slopes (Avals), Standard Errors (SEvals), and
  # Z-scores (BootZ)
  Avals <- rep(NA,1000)
  SEvals <- rep(NA,1000)
  BootZ <- rep(NA,1000)
  
  # tic()
  
  # Iterate 1000 times for the bootstrapping
  for(n in 1:1000){
    # print(n)
    
    # Initialize a dataframe with the Variance (x-axis), FST (y-axis) and weights
    # (varV) for each block
    points <- data.frame(V = rep(NA,1703),
                         FST = rep(NA,1703),
                         varV = rep(NA,1703),
                         BLOCK = rep(NA,1703))
    
    
    for(i in 1:length(blocks$chr)){
      # cat("\rRep ",i,"/",1703,sep = "")
      pos <- which(dfi$CHROM == blocks$chr[i] & 
                     dfi$POS >= blocks$start[i] &
                     dfi$POS <= blocks$stop[i])
      if(length(pos) == 0){next}
      # print("Getting random sample points")
      samp <- sample(pos,size = 1)
      points[i,1:3] <- dfi[samp,3:5]
      points$BLOCK[i] <- i

    }
    # print("Getting slope of random points")
    na.omit(points)
    
    # Get the sites that fall within each block
    bpoints <- data.frame(V = rep(NA,length(points[[1]])),
                          FST = rep(NA,length(points[[1]])),
                          varV = rep(NA,length(points[[1]])))
    
    # Choose one point at random
    bootp <- sample(points$BLOCK,replace = TRUE)
    for(b in 1:length(bootp)){
      blockn <- bootp[b]
      pos <- which(dfi$CHROM == blocks$chr[blockn] & 
                     dfi$POS >= blocks$start[blockn] &
                     dfi$POS <= blocks$stop[blockn])
      if(length(pos) == 0){next}
      # print("Getting random sample points")
      samp <- sample(pos,size = 1)
      
      # Get the values of V, FST, and varV
      bpoints[b,1:3] <- dfi[samp,3:5]
    }
    
    # Calculate regression weighting by 1/varV
    points$w <- 1/points$varV
    
    # An if statement to skip any blocks with no sites in them
    if(class(try(lm(formula = FST ~ V,data = points,weights = w),
                 silent = TRUE)) == "try-error"){next}
    
    # Calculate the slope of a weighted regression A
    A <- lm(formula = FST ~ V,data = points,weights = w)
    S <- summary(A)
    Avals[n] <- A$coefficients[2]
    
    # Similarly, get the SE of the regression
    if(class(try(S$coefficients[2,2],silent = TRUE)) == "try-error"){next
    }else{
      SEvals[n] <- S$coefficients[2,2]
    }
    
    # Save the regression slope and SE for calculating a Z-score
    bpoints$w <- 1/bpoints$varV
    if(class(try(lm(formula = FST ~ V,data = points,weights = w),
                 silent = TRUE)) == "try-error"){next}
    bA <- lm(formula = FST ~ V,data = bpoints,weights = w)
    bS <- summary(bA)
    Aboot <- bA$coefficients[2]
    if(class(try(bS$coefficients[2,2],silent = TRUE)) == "try-error"){next
    }else{
      Sboot <- bS$coefficients[2,2]
    }
    Zboot <- Aboot/Sboot
    
    # Add the z-score to the initialized zscore vector
    
    BootZ[n] <- Zboot
    # cat("\n")
  }
  print(paste("ABoot",Aboot))
  print(paste("SBoot",Sboot))
  print(paste("Zboot",Zboot))
  print("Head BootZ")
  print(head(BootZ))
  t = toc()
  
  # Get the values of slopes and SEs from all 1000 iterations
  outputdf <- data.frame(slopes = Avals,SEs = SEvals)
  outputdf <- na.omit(outputdf)
  if(length(outputdf$slopes) == 0){next}
  # Save in a table
  fi <- paste0("./",trait,"/",a,".",trait,".slopes.",pstr,".NEW.all.csv")
  write.table(outputdf,file = fi,sep = ",",quote = FALSE,
              row.names = FALSE,col.names = TRUE,append = FALSE)
  
  # cat("\nCreating table of results")
  
  # Append the mean slope and SE in a table that will contain info from all traits
  rowi <- data.frame(trait,a,
                     mean(Avals,na.rm = TRUE),
                     mean(SEvals,na.rm = TRUE),
                     sd(Avals,na.rm = TRUE))
  resultsfile <- paste0("./results.",pstr,".NEW.all.csv")
  write.table(rowi,file = resultsfile,append = TRUE,sep = ",",quote = FALSE,
              col.names = FALSE,row.names = FALSE)
  
  outputz <- outputdf
  
  # Get a table that includes z-score info
  outputz$Zscore <- outputz$slopes/outputz$SEs
  fi <- paste0("./",trait,"/",a,".",trait,".slopes.",pstr,".z.NEW.all.csv")
  write.table(outputz,file = fi,sep = ",",quote = FALSE,
              row.names = FALSE,col.names = TRUE,append = FALSE)
  
  rowz <- data.frame(trait,a,
                     mean(Avals,na.rm = TRUE),
                     mean(SEvals,na.rm = TRUE),
                     sd(Avals,na.rm = TRUE),
                     mean(outputz$Zscore),sd(outputz$Zscore))
  resultsfile <- paste0("./results.",pstr,".z.NEW.all.csv")
  write.table(rowz,file = resultsfile,append = TRUE,sep = ",",quote = FALSE,
              col.names = FALSE,row.names = FALSE)
  
  # Get a table that includes Confidence Interval info
  rowCI <- data.frame(trait,a,
                      mean(Avals,na.rm = TRUE),
                      mean(SEvals,na.rm = TRUE),
                      sd(Avals,na.rm = TRUE),
                      mean(outputz$Zscore),
                      sd(outputz$Zscore),sd(BootZ),
                      quantile(BootZ,0.05),quantile(BootZ,0.95))
  resultsfile <- paste0("./results.",pstr,".z.CI.all.csv")
  write.table(rowCI,file = resultsfile,append = TRUE,sep = ",",quote = FALSE,
              col.names = FALSE,row.names = FALSE)
  
  cat("\n")
}
