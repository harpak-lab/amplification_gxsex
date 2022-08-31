#!/usr/bin/env Rscript

source("config.R")
rm(list= ls()[!(ls() %in% c('R_LIB','LDSC_FILE'))])
library("dplyr", lib.loc=R_LIB)
library("tidyr", lib.loc=R_LIB)
library("ggplot2", lib.loc=R_LIB)
library("grid", lib.loc=R_LIB)
library("gridExtra", lib.loc=R_LIB)
library("ggpubr", lib.loc=R_LIB)

setwd(LDSC_FILE)
df <- read.csv("ldsc_results.txt", sep="\t")

both_h2 <- rep(df[df$Sex == "both_sex", "Heritability"], each=3)
asterick_pheno <- c("arm_fatfree_mass_R", "weight", "arm_fatfree_mass_L", "bmi",
                    "whole_body_fat_mass", "waist_circ", "hip_circ", "waist_to_hip_ratio")

df <- df %>%
  mutate(relative_h2 = Heritability / both_h2) %>%
  mutate(relative_h2_se = ((Heritability + h2.std.error) / both_h2) - relative_h2 ) %>%
  mutate(Sex = factor(Sex, levels=c("female","both_sex","male"))) %>%
  mutate(star = ifelse((Code %in% asterick_pheno) & (Sex == 'both_sex'), "*", NA)) %>%
  arrange(Correlation, Phenotype, Sex) %>%
  mutate(Phenotype = factor(Phenotype, levels=unique(Phenotype)))

# write to table
write.table(df, file = "relative_h2.txt", quote=FALSE, sep="\t", row.names=FALSE) 

pdf(file="FIG1.r2_by_h2.pdf", width=7, height=8)

rects <- data.frame(ystart = seq(0.5,26.5,1), yend = seq(1.5,27.5,1), col = c(1,rep(c(2,1),13)))
p1 <- ggplot(df, aes(x=relative_h2, y=Phenotype, col=(Sex))) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.5) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(y=Phenotype, xmin=relative_h2-relative_h2_se, xmax=relative_h2+relative_h2_se), 
                 height=0, position=position_dodge(width=0.5)) +
  geom_rect(data=rects, aes(ymin=ystart,ymax=yend,xmin=0.6,xmax=2.3), 
            inherit.aes = FALSE, alpha=0.2, fill = c("white",rep(c("grey","white"),13))) +
  geom_point(aes(x=relative_h2-0.4, shape=as.factor(star)), color="#b01924", size=4) +
  scale_shape_identity() +
  scale_x_continuous(breaks=c(0.5,1,1.5,2), limits=c(0.6,2.3), expand=c(0,0)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=10), axis.text.y = element_text(hjust=0), 
        axis.title.y = element_blank(), axis.title.x = element_text(size=12),
        plot.margin = margin(5.5,5.5,5.5,0), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(labels = c("female", "both sex", "male"), values=c("#d67629","#563f61","#207335")) +
  labs(x="Heritability Relative to Heritability of Both-sex Sample", col="Sex") + 
  # legend
  geom_rect(inherit.aes = FALSE, xmin=1.75,xmax=2.1, ymin=20.5,ymax=23.5, fill="white") +
  annotate("text", x=1.8, y=23, label="male", hjust=0, color="#207335") +
  annotate("text", x=1.8, y=22, label="both sex", hjust=0, color="#563f61") +
  annotate("text", x=1.8, y=21, label="female", hjust=0, color="#d67629")

corr_df <- df %>%
  select(c(2,6)) %>%
  distinct() %>%
  mutate(pheno_point = seq(-0.03,1.03,by=1.06/(length(Phenotype)-1)))

p2 <- ggplot(corr_df, aes(x=0,y=Correlation)) +
  geom_segment(aes(x = 0, y = -0.05, xend = 0, yend = 1.05),
               arrow = arrow(length = unit(0.3, "cm"), end="both"), size=1, color="#2b62d9") +
  geom_segment(aes(x=0, y=Correlation, xend=1, yend=pheno_point), alpha=0.3) +
  geom_point(shape=1, size=3) +
  geom_point(aes(x=1,y=pheno_point),size=0.1) +
  scale_y_continuous(limits=c(-0.05,1.05), expand=c(0,0)) +
  scale_x_continuous(limits=c(-0.1,1.01)) +
  theme(plot.margin = margin(10,0,19,5), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_blank(), axis.title.y = element_text(size=12), axis.text.y = element_text(size=10)) +
  labs(x="", y="Genetic Correlation between Males and Females")

lay <- rbind( c(1,2,2,2,2,2))
p <- grid.arrange(p2, p1, ncol = 2, layout_matrix=lay, 
                  top=textGrob("", gp=gpar(fontsize=16)))

p
dev.off()

## STAT
# relative heritability diff by correlation
stat_df <- df[df$Sex != 'both_sex', c(1,3,6,8,9)]
f <- stat_df[stat_df$Sex == 'female',]
m <- stat_df[stat_df$Sex == 'male',]
stat_df <- data.frame(Code = f$Code, Correlation = f$Correlation, h2_diff = abs(f$relative_h2 - m$relative_h2))
# heritability difference
model <- cor.test(stat_df$Correlation, stat_df$h2_diff)
print(model)
