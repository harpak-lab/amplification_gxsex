# Getting final version of plots for Main Text of manuscript

# Initialize Environment ########################################

# We will use the ggplot2 package to generate figures
suppressPackageStartupMessages(library(ggplot2))
# We will use the ggpubr package to combine figures
suppressPackageStartupMessages(library(ggpubr))
# We will use the dplyr package to wrangle data.frames
suppressPackageStartupMessages(library(dplyr))

# We will draw the p-value threshold from a command-line variable passed into
# our Rscript, which will allow us to analyze all p-value thresholds (1e-03,
# 1e-05, and 1e-08) but we will focus on 1e-05 for our main text
num <- commandArgs(trailingOnly = TRUE)[2]
AncList <- commandArgs(trailingOnly = TRUE)[1]
p <- as.numeric(num)

cat("P-value Threshold:",num)

# Get list of ancestry mean and median sample sizes for XX and XY individuals
AncSizes <- read.csv("/scratch/08312/mjm8356/data/AvgSampleSizes.csv",header = TRUE)
# Use a list of GWAS traits (as abbreviations), and store corresponding "Full"
# trait names
TraitList <- data.frame(abr = c("albumin","arm_fatfree_mass_L","arm_fatfree_mass_R","bmi",
               "calcium","creatinine","diastolicBP_auto","eosinophil_perc",
               "FVC_best","HbA1c","height","hip_circ","IGF1","lymphocyte_perc",
               "protein_total","pulse_rate","RBC_count","SHBG","systolicBP_auto",
               "testosterone","urate","urea","waist_circ","waist_to_hip","weight",
               "whole_body_fat_mass","wth_bmi_adj"),
               full = c("Albumin","Arm fat free mass (L)","Arm fat free mass (R)",
                        "BMI","Calcium","Creatinine","Diastolic BP",
                        "Eosinophil percentage","Forced vital capacity","HbA1c",
                        "Height","Hip circumference","IGF-1","Lymphocyte percentage",
                        "Total protein","Pulse rate","RBC count","SHBG",
                        "Systolic BP","Testosterone","Urate","Urea",
                        "Waist circumference","Waist to hip ratio","Weight",
                        "Whole body fat mass","Waist:Hip (BMI adjusted)"))
# Concatenate all Ancestry and Trait abbreviated and full names in one data.frame
FullList <- data.frame(abr = c(TraitList$abr,"nfe","asj","fin","afr","amr","sas","eas","oth"),
                       full = c(TraitList$full,"Non-Finnish European",
                                "Ashkenazi Jewish","Finnish","African","Latino/American Admixed",
                                "South Asian","East Asian","Other"))

# Setting up universal variables used in generating all FST vs. V plots
blocks <- read.table("/scratch/08312/mjm8356/GWAS_results/Pickrell_breakpoints_EUR.bed",header = TRUE)

# Create a helper function to get the "Full" name from a given abbreviation for
# Ancestry or Trait
abrToFull <- function(str){
  # Argument: str = string containing abbreviation for Ancestry or Trait
  # Output: corresponding full name of Ancestry or Trait
  
  out <- FullList %>% filter(abr %in% str) %>% 
    arrange(factor(abr,levels = str)) %>% pull(full)
  return(out)
}

# AncList <- "nfe"

# FST vs. V plot + Regression #################################

# Create a helper function to get the full regression line equation for points
# from a given sampling method
LMeq <- function(samp,df){
  # Function generates regression lines for FST vs V weighted by 1/var[V], then 
  # saves the equation as "FST = [slope]V + [intercept] (SE slope = [SEslope],
  # SE Int = [SEint]"
  
  # Arguments: samp = Sampling method ("R" for "Random" or "M" for "Max V Value")
  # output: "Y = mx + b" form of equation as a character string
  
  pe <- df[which(df$Samp == samp),]
  r <- lm(formula = FST ~ V,data = pe,weights = w)
  se <- summary(r)$coefficients
  a <- formatC(se[2,1],digits = 4)
  b <- formatC(se[1,1],digits = 3)
  SEa <- formatC(se[2,2],format = "e",digits = 2)
  SEb <- formatC(se[1,2],format = "e",digits = 2)
  eq <- paste0("FST = ",a,"V + ",b," (SE slope = ",SEa,", SE Int = ",SEb,")")
  return(eq)
}

# Create a function to generate representative FST vs. V_{GxSex} plots
FST_Vplot <- function(a,t,u){
  # Function generates a plot of FST vs. V (where V = 2p(1-p)(beta_m-beta_f)^2)
  # and gets a regression line for FST vs. V weighted by 1/variance[V]
  # with slope, y-intercept, and SEs of each
  
  # 1) Filter loci by a GWAS p-value cutoff
  # 2) Move through a list of LD blocks (defined by Berisa and
  # Pickrell, 2015) and randomly select one of the loci that falls within
  # that LD block range, as well as the locus with the highest V value per
  # LD block
  # 3) Find a weighted regression of FST vs. V with weights 
  # determined by 1/var[V]
  # 4) Repeat the above steps 100 times
  
  # For the iteration with the median slope
  # 4) Plot the FST and V values for all loci selected this way,
  # and plot them on a graph with V on x-axis and FST on y-axis; opacity of
  # each point is determined by 1/var[V]
  # 5) Print the slope, y-intercept, and Standard Error (SE) of each,
  # then plot this regression line
  
  # Arguments:
  # a = Ancestry group
  # t = Trait of interest
  # u = Units (for axis label)
  
  # outputs:
  # Figure of FST vs. V
  # Slope, Slope SE, Y-Intercept, Y-Intercept SE
  # cat("This is my wd:",wd)
  # setwd(wd)
  # cat("\n",getwd(),sep = "")
  
  if(t == "protein_total"){set.seed(123)}
  else{set.seed(321)}# Allow for reproducible results when re-running script
  tFull <- abrToFull(t)
  # aFull <- abrToFull(a)
  
  cat("\n\nTrait:",tFull)
  # cat("\nAncestry:",aFull)
  
  # Read in p-value from GWAS data
  cat("\n\nReading p-val Table")
  # pstats <- read.table(paste0("./pstats.",t,".txt"),header = TRUE)
  pstats <- read.table(paste0("/scratch/08312/mjm8356/GWAS_results/",t,"/pstats.",t,".txt"),header = TRUE)
  
  cat("\n")
  print(head(pstats))
  
  cat("\nFiltering p-vals by p <",num)
  pFilt <- which(pstats$MP <= p | pstats$FP <= p)
  pl <- paste(pstats$CHROM[pFilt],pstats$POS[pFilt],sep = ":")
  
  # Read in data.frame of FSTs, V values, and var[V]
  cat("\n\nReading Data Table")
  # di <- read.table(paste0("./",a,".table.",t,".txt"),header = TRUE)
  di <- read.table(paste0("./",t,"/NewFST.",a,".table.",t,".txt"),header = TRUE)
  
  cat("\n")
  print(head(di))
  
  # Define chromosomal loci in data as Chr:Pos (e.g., chr1:1234)
  cat("\nFiltering data")
  dl <- paste(di$CHROM,di$POS,sep = ":")
  
  # mismap <- read.table(paste0("/scratch/08312/mjm8356/data/GRCH37/nfe_allMisMatch.tsv"),
  #                      fill = TRUE)
  # mmloci <- paste0(mismap[[2]],":",mismap[[1]])
  
  # Find which loci in data are also found in set of all GWAS loci filtered by 
  # p-value, and filter data to only include these loci
  # dk <- which(dl %in% pl & !(dl %in% mmloci))
  dk <- which(dl %in% pl)
  d <- di[dk,]
  
  cat("\nAfter filtering, number of sites decreases from",length(di$CHROM),
      "to",length(d$CHROM))
  
  # Initialize data.frame which will store V, FST, and var[V] for each of 
  # 1) one randomly sampled locus per LD block (Samp = "R")
  # 2) the locus with the highest V value per LD block (Samp = "M")
  points <- data.frame(V = rep(NA,1703),
                       FST = rep(NA,1703),
                       varV = rep(NA,1703),
                       Samp = rep(NA,1703))
  
  # Loop over LD blocks
  cat("\n")
  for(i in 1:length(blocks$chr)){
    cat("\rLooping over LD block ",i,"/",1703,sep = "")
    # Which data loci fall within LD block
    pos <- which(d$CHROM == blocks$chr[i] & 
                   d$POS >= blocks$start[i] &
                   d$POS <= blocks$stop[i])
    if(length(pos) == 0){next} # If there are no data loci within LD block, move
                               # to next block
    samp <- sample(pos,size = 1) # Sample one locus
    R <- d[samp,3:5]
    R$Samp <- "R"
    points[i,] <- R  # Add that locus' V, FST, and var[V] to points
    sampm <- pos[which(d$V[pos] == max(d$V[pos])[1])] # Get single locus w/ highest
                                                  # V value in LD block
    # M <- d[sampm,3:5]
    # M$Samp <- "M"
    # print(M)
    # points[i*2,] <- M # Add that locus' V, FST, and var[V] to mpoints
  }
  
  cat("\nGetting Equations")
  # Filter out any rows in points and mpoints (i.e., any LD blocks) which didn't
  # have any data loci in them
  points <- na.omit(points)
  
  # Define weights as 1/var[V]
  points$w <- 1/points$varV
  
  eqR <- LMeq("R",points)

  # eqM <- LMeq("M",points)
  
  cat("\nCreating Graph")
  # Plot points and weighted linear regression on FST vs. V graph, with FST on 
  # y-axis and V on x-axis, with point opacity determined by 1/var[V]
  pointsf <- points[which(points$Samp == "R"),]
  
  # Get color of regression line and CI shading based on ancestry/trait combo
  invisible(ifelse(a == a1 & t == t1,lcol <- "magenta2",lcol <- "#00BFC4"))
  # invisible(ifelse(a == a1 & t == t1,fcol <- "orchid1",fcol <- "skyblue"))
  
  mxf <- sort(pointsf$V,decreasing = TRUE)[3]
  myf <- sort(pointsf$FST,decreasing = TRUE)[2]
  
  eqdf <- data.frame(
    X = 0.5*mxf,
    Y = myf - 0.05*myf,
    lab = eqR
  )
  
  hascarrot <- which(strsplit(u,"")[[1]] == "^")
  
  if(length(hascarrot) > 0){
    varstr <- substr(u,1,hascarrot-1)
    xax <- bquote(V[GxS]~ "\u2014 Phenotypic Variance due to GxSex (" 
                  *.(varstr)^2* ")")
  }else if(length(hascarrot) == 0){
    xax <- bquote(V[GxS]~ "\u2014 Phenotypic Variance due to GxSex (" *.(u)* ")")
  }

  save(pointsf, file=paste0("nfe_fst_plot_",t,".",p,".Rdata")) ### cz
  
  p1 <- ggplot(pointsf) + 
    geom_smooth(mapping = aes(V,FST,weight = w),
                method = "lm",se = FALSE,color = lcol) +
    geom_point(aes(V,FST,size = w),color = "black",alpha = 0.1) +
    # geom_text(data = eqdf,aes(x = X,y = Y,label = lab),
    #           color = "black",size = 3,family = "Ariel") +
    coord_cartesian(xlim = c(0,mxf),ylim = c(0,myf)) +
    ggtitle(paste("LST vs V for",tFull,"\nin the nfe Population")) +
    xlab(xax) +
    ylab(bquote("Male\u2013Female " ~L[ST])) + 
    scale_y_continuous(breaks = seq(0,(0.75*myf),length.out = 4),
                       labels = function(FST) ifelse(FST == 0,"0",formatC(FST,format = "e",digits = 2))) +
    scale_x_continuous(labels = function(V) ifelse(V == 0,"0",formatC(V,format = "e",digits = 1))) +
    theme(plot.title = element_text(size = 16,hjust = 0.5),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x.bottom = element_line(color = "black"),
          axis.line.y.left = element_line(color = "black"),
          axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 14))
  of1 <- paste0(t,".",a,".",num,".FST_V.GAM.jpeg")
  ggsave(of1,plot = p1,dpi = 300,path = "./Figs/",device = "jpeg")
  return(p1)
}

# Z-score Plot ###################################################

# Read in appropriate table of Z-scores based on p-value filtering level
z <- read.csv(paste0("./results.",num,".z.CI.all.csv"),header = TRUE)

# Create a function to generate a plot comparing Z-scores across traits by
# ancestry
ZscorePlot <- function(Ancs){
  # Function will create a wrapped set of plots
  results <- z[which(z$ANC %in% Ancs),]
  resultsi <- results %>% group_by(ANC) %>% summarise(ANCMEAN = mean(ZMEAN))
  resultsi$LAB <- "Mean Z-score for Ancestry"
  
  texts <- AncSizes[which(AncSizes$anc %in% Ancs),]
  for(n in names(texts)[2:5]){
    texts[[n]] <- paste(gsub("(n)(*)", "\\1 \\2", n),"AC:",
                        round(texts[[n]],digits = 1))
  }
  
  # Set up a "COL" (colors) column to allow for coloring two focal bars/points
  results$COLS <- ifelse(results$TRAIT == t1,"col1",
                         ifelse(results$TRAIT == t2,"col2",
                                "b"))
  
  # Organize the TRAIT and ANC columns so they appear alphabetically in the fig
  ord <- results %>% group_by(ANC) %>% arrange(desc(ZMEAN)) %>% 
    filter(ANC == Ancs[1]) %>% pull(TRAIT)
  
  cat(rev(abrToFull(ord)),"\n\n")
  
  results$TRAIT <- factor(results$TRAIT,
                          levels = rev(ord),
                          labels = rev(abrToFull(ord)))
  results$ANC <- factor(results$ANC,levels = Ancs,
                        labels = abrToFull(Ancs))

  # Format table containing mean information so that it also is in alphabetical
  # order by ancestry, and then add in a column to print the median size of
  # each ancestry
  resultsi$ANC <- factor(resultsi$ANC,levels = Ancs,
                         labels = "nfe")
  texts$anc <- factor(texts$anc,levels = Ancs,
                      labels = "nfe")
  resultsi$LAB2 <- paste(texts$medianmale,texts$medianfemale,sep = "\n")

  save(resultsi, results, data, 
       file=paste0(Ancs,"_zscore_plot.",num,".all.Rdata")) ### cz
  
  print(resultsi$ANCMEAN)
  
  results$ZSENORM <- results$ZSE/sqrt(length(results$ZSE)-1)
  
  p2 <- ggplot() + 
    geom_vline(data = resultsi,mapping = aes(xintercept = ANCMEAN),
               linetype = "dashed",size = 0.5,show.legend = F,color = "red") +
    geom_vline(data = resultsi,mapping = aes(xintercept = 0),
               linetype = "dashed",size = 0.5,show.legend = F,color = "blue") +
    # geom_errorbar(data = results,
    #               aes(ZMEAN,TRAIT,xmin = ZMEAN - ZSENORM,xmax = ZMEAN + ZSENORM),
    #               color = "black",linetype = "dashed",width = 0,size = 0.5) +
    geom_errorbar(data = results,
                  aes(ZMEAN,TRAIT,xmin = CILOWER,xmax = CIUPPER,
                      color = COLS),
                  width = 0,size = 1) +
    geom_point(data = results,aes(ZMEAN,TRAIT,color = COLS),size = 2) +
    scale_color_manual(name = NULL,
                       values = c("magenta2","#00BFC4","black"),
                       breaks = c("col1","col2","b"),
                       guide = NULL) +
    # facet_wrap(~ANC,ncol = 4) +
    ggtitle(paste("Testing for Sexually-Antagonistic Selection")) +
    xlab("Z-score for Sexually-Antagonistic Selection") +
    geom_text(data = resultsi,mapping = aes(x = Inf,y = 39,label = LAB2),
              size = 6,hjust = 1,vjust = 0.5,family = "Ariel") +
    geom_text(data = resultsi,mapping = aes(x = (ANCMEAN + 0.2),y = 34,label = LAB),
              angle = 270,size = 6,color = "red",family = "Ariel") +
    geom_text(data = resultsi,mapping = aes(x = -0.2,y = 34,label = "zero"),
              angle = 90,size = 6,color = "blue",family = "Ariel") +
    coord_cartesian(xlim = c(-2,5.5),ylim = c(0,40)) +
    theme(plot.title = element_text(size = 16,hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 25),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "grey92"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x.bottom = element_line(color = "black"),
          axis.line.y.left = element_line(color = "black"),
          axis.text = element_text(size = 20),
          strip.text = element_text(size = 20))
  of2 <- paste0("zscores.",num,".ALL.",Ancs,".jpeg")
  ggsave(of2,plot = p2,dpi = 300,path = "./Figs/",device = "jpeg",
         height = 20,width = 15)
  return(p2)
}

# Run Code ###################################################

# Whole Body Fat Mass in the Finnish population
a1 <- AncList
t1 <- "protein_total"
u1 <- "g/L"
# plot1 <- FST_Vplot(a1,t1,u1) +
#   theme(plot.title = element_blank(),
#         axis.title.y.left = element_blank())

cat("\n")

# Waist Circumference in the Ashkenazi Jewish population
a2 <- AncList
t2 <- "testosterone"
u2 <- "nmol/L"
# plot2 <- FST_Vplot(a2,t2,u2) +
#   theme(plot.title = element_blank(),
#         axis.text.x.top = element_blank(),
#         axis.ticks.x.top = element_blank(),
#         axis.text.y.right = element_blank(),
#         axis.ticks.y.right = element_blank(),
#         axis.title.y.left = element_blank(),
#         axis.title.y.right = element_blank())

# Z-score plot for Non-Finnish European, Finnish, and Ashkenazi Jewish
# populations
a3 <- c(AncList)
# a3 <- c("asj","fin","nfe")
plot3 <- ZscorePlot(a3) + theme(plot.title = element_blank())

# Put all figures together #########################################

# Collect all plots together with zoomed-in plots (scatter plots) on top of
# z-score (error bar) plot

# text1 <- "Point size\nproportional to\nregression weight"
# 
# figtopr <- ggarrange(plot1,plot2,ncol = 1,nrow = 2,
#                     heights = c(1,1))
# figtopr <- annotate_figure(figtopr,
#                            left = text_grob(bquote(F[ST]~ "Between Males and Females"),
#                                              rot = 90,size = 20,x = 0),
#                            right = text_grob(bquote(atop("Points are SNPs\nwith size\nproportional to\nregression weight\n\nLines Fit to\nTheoretical Model\n",
#                                                          F[ST] %~~% "A" ~V[GxS])),
#                                              rot = 0,size = 18,
#                                              hjust = 0,x = 0,vjust = 1))
# bl <- ggplot() + theme_void()
# figtop <- ggarrange(bl,figtopr,widths = c(1.2,1.2),labels = c("A","B"),
#                     font.label = list(color = "black",size = 18))
# 
# 
# 
# figfull <- ggarrange(figtop,plot3,ncol = 1,nrow = 2,labels = c("","C"),
#                      heights = c(0.7,1.2),
#                      font.label = list(color = "black",size = 18))
# 
# 
# figfinal <- annotate_figure(figfull,
#                              top = text_grob("Testing for Sexually-Antagonistic Selection",
#                                              hjust = 0.5,vjust = 0,size = 35,
#                                              y = 0.25,
#                                              face = "bold"))
# 
# # ojpeg <- paste0("FullFig.",num,".jpeg")
# opdf <- paste0("FullFig.",num,".ALL.pdf")
# 
# # ggsave(ojpeg,plot = figfinal,path = "./Figs/Main/",device = "jpeg",
# #        width = 20,height = 18)
# ggsave(opdf,plot = figfinal,path = "./Figs/Main/",device = cairo_pdf,
#        width = 20,height = 18,dpi = 300)

# df <- read.csv("results.1e-05.z.CI.csv",header = TRUE)
# 
# for(i in c(which(df$TRAIT %in% c("whole_body_fat_mass",
#                                  "waist_circ",
#                                  "hip_circ") & df$ANC == "asj"),
#            which(df$TRAIT %in% c("whole_body_fat_mass",
#                                  "waist_circ",
#                                  "bmi") & df$ANC == "fin"))){
#   cat("\n",df$ANC[i],df$TRAIT[i])
#   cat("\nZ =",sprintf(signif(df$ZMEAN[i],digits = 4),fmt = "%#.4g"),
#       "\u00B1",
#       sprintf(signif(abs(df$CILOWER[i] - df$ZMEAN[i]),digits = 4),
#               fmt = "%#.4g"),"\n")
# }

