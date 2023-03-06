library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggrepel)

# names
setwd("~/Documents/Harpak/Phenotypes/")
n_df <- read.csv("names.txt", sep="\t")
tail(n_df)

# mash weights
setwd("~/Documents/Harpak/GxSex/mash")
mash_df <- read.csv("mash_weights.txt", sep="\t")
head(mash_df)
mash_df <- mash_df %>%
  merge(n_df, by.x="phenotype", by.y="Code") %>%
  mutate(amp = abs(sum_weight_f-sum_weight_m)) %>%
  select(c(1,6,8))


# polygenic score inc r2
setwd("~/Documents/Harpak/GxSex/PGS")
pgs_df <- read.csv("pgs_20_all.txt", sep="\t")
pgs_add <- pgs_df[pgs_df$model == "add",c(2:28)]
pgs_mash <- pgs_df[pgs_df$model == "mash",c(2:28)]
# summarize pgs dataframe
pgs_all <- pgs_mash/pgs_add
# pgs average
pgs_sum <- colSums(pgs_all)/20
pgs_sum <- cbind(names(pgs_sum), data.frame(pgs_sum, row.names=NULL))
colnames(pgs_sum)[1] <- "phenotype"
# pgs 20
pgs_all <- t(pgs_all); pgs_all <- cbind(rownames(pgs_all), data.frame(pgs_all, row.names=NULL))
colnames(pgs_all)[1] <- "phenotype"
pgs_all <- merge(merge(pgs_all, pgs_sum, by="phenotype"), mash_df, by = "phenotype")
pgs_all$Phenotype <- factor(pgs_all$Phenotype, levels = pgs_all$Phenotype[order(pgs_all$pgs_sum)])
tail(pgs_all)

psize=0.3
ggplot(data=pgs_all, aes(x=Phenotype, color=amp)) +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.5) +
  geom_vline(aes(xintercept=Phenotype), alpha=0.1) +
  geom_point(aes(y=X21), size=psize) + geom_point(aes(y=X22), size=psize) + geom_point(aes(y=X23), size=psize) + geom_point(aes(y=X24), size=psize) + geom_point(aes(y=X25), size=psize) +
  geom_point(aes(y=X26), size=psize) + geom_point(aes(y=X27), size=psize) + geom_point(aes(y=X28), size=psize) + geom_point(aes(y=X29), size=psize) + geom_point(aes(y=X30), size=psize) +
  geom_point(aes(y=X31), size=psize) + geom_point(aes(y=X32), size=psize) + geom_point(aes(y=X33), size=psize) + geom_point(aes(y=X34), size=psize) + geom_point(aes(y=X35), size=psize) +
  geom_point(aes(y=X36), size=psize) + geom_point(aes(y=X37), size=psize) + geom_point(aes(y=X38), size=psize) + geom_point(aes(y=X39), size=psize) + geom_point(aes(y=X40), size=psize) +
  geom_point(aes(y=pgs_sum, color=amp), size=3.5) +
  xlab("") + ylab("Relative Prediction Accuracy (GxSex PGS : Additive PGS)") +
  scale_y_continuous(trans="log2", breaks = seq(0.25,2,0.25)) +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_text(size=11), axis.title.y = element_blank(), axis.text = element_text(size=9)) +
  scale_color_continuous(low="grey", high="#2b62d9")
  
# correlation
cor.test(pgs_all$amp, pgs_all$pgs_sum, method="pearson")