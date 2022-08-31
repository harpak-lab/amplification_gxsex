#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','GWAS_DIR'))])
library("optparse", lib.loc=R_LIB)
library("crayon",lib.loc=R_LIB)
library("dplyr",lib.loc=R_LIB)
library("tidyr",lib.loc=R_LIB)
library("withr",lib.loc=R_LIB)
library("ggplot2",lib.loc=R_LIB)
library("grid",lib.loc=R_LIB)
library("gridExtra",lib.loc=R_LIB)
library("labeling",lib.loc=R_LIB)
library("farver",lib.loc=R_LIB)
library("digest",lib.loc=R_LIB)
library("backports",lib.loc=R_LIB)
library("ggpubr",lib.loc=R_LIB)
library("ggsci",lib.loc=R_LIB)

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

setwd( paste0(GWAS_DIR,"/",pheno) )
male_name <- paste0("male_all.",pheno,".glm.linear")
female_name <- paste0("female_all.",pheno,".glm.linear")
male_df <- read.csv(male_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ))
female_df <- read.csv(female_name, sep="\t", colClasses=c(rep("integer",2),"character",rep("NULL",9),"numeric" ))

edit_df <- function(df) {
  for (p in c(5e-30, 5e-20, 5e-10, 5e-5, 0.001, 0.05, 0.75, 0.1, 0.5, 0.7, 0.8)) {
    df <- bind_rows(filter(df, P<=p), sample_frac(filter(df, P>p), 0.25))
  }
  plot_df <- df %>%
  drop_na() %>%
  group_by(X.CHROM) %>%
  summarize(CHR_LEN=max(POS)) %>%
  mutate(TOT=cumsum(as.numeric(CHR_LEN))-CHR_LEN) %>%
  select(-CHR_LEN) %>%
  left_join(df, ., by=c("X.CHROM"="X.CHROM")) %>%
  arrange(X.CHROM, POS) %>%
  mutate(POS_CUM=POS+TOT) %>%
  mutate(COLOR=ifelse(X.CHROM %% 2, 1,2))

  return(plot_df)
}

x_axis <- function(plot_df) {
  axis_df <- plot_df %>% 
    group_by(X.CHROM) %>%
    summarize(CENTER=(max(POS_CUM)+min(POS_CUM))/2 )
  return(axis_df)
}

female_df <- edit_df(female_df)
male_df <- edit_df(male_df)
axis_df <- x_axis(female_df)
axis_df$X.CHROM <- c(1:18,"",20,"",22)

pdf(file=paste0(pheno,"_miami.pdf"), width=7, height=3)

female_plot <- ggplot(female_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=0.1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_continuous(expand=c(0,0)) +
  
  theme_pubclean() +
  theme(legend.position = "none", axis.title = element_blank(), plot.margin = margin(20,5.5,0,5.5),
        axis.text = element_text(size=9), panel.grid.major.y=element_blank(), plot.title=element_text(size=16, hjust=0.5)) +
  scale_color_manual(values=c("#d67629","#1d47a1")) +
  labs(title=pheno_name) +
  annotation_custom(grobTree(textGrob("female", x=0.1, y=0.9, gp = gpar(col="#d67629", fontsize=11) )))

male_plot <- ggplot(male_df, aes(x=POS_CUM,y=-log10(P))) +
  geom_point( aes(color=as.factor(COLOR)), alpha=0.7, size=0.1) +
  
  scale_x_continuous(label=axis_df$X.CHROM, breaks=axis_df$CENTER) +
  scale_y_reverse(expand=c(0,0)) +
  
  theme_pubclean() + 
  theme(legend.position = "none",  axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
        plot.margin =  margin(0,5.5,40,5.5), panel.grid.major.y=element_blank(), 
        axis.text.y = element_text(size=9), axis.title.y = element_blank(), axis.title.x = element_text(size=11)) +
  scale_color_manual(values=c("#d67629","#1d47a1")) +
  labs(x="Chromosome") +
  annotation_custom(grobTree(textGrob("male", x=0.1, y=0.1, gp = gpar(col="#207335", fontsize=11) )))

p <- grid.arrange(female_plot, male_plot, nrow = 2, left=textGrob(bquote(-log[10] ~P), rot=90, vjust=1, gp=gpar(fontsize=11)))
p

dev.off()

