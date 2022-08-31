#####################
# configuration file
#####################

# process paramerters
mem=64000
threads=16

#####################
# FILE PATHS

### R library files ###
R_LIB = "~/R/x86_64-pc-linux-gnu-library/4.0/"

### General ###
QC_DIR="~/QC"                   # quality controlled pgen sample and variant files
PHENO_DIR="~/Phenotypes"        # phenotype files
GWAS_DIR="~/GWAS_Results"       # summary statistics, results, intermediate tables
# mash and PGS files directories are subdirectories of $GWAS_DIR/{phenotype}/ 

### ENSEMBL annotation files ###
ENSEMBL="~/ensembl-vep"         # ensembl vep github (https://github.com/Ensembl/ensembl-vep.git)
ENSEMBL_CACHE="~/cache"         # homo_sapiens_CRCh38.vcf cache

### LDSC files ###
LDSC_DIR="~/ldsc"               # ldsc github (https://github.com/bulik/ldsc.git)
LDSC_FILE="~/LD_files"          # results from LDSC regression
LD_SCORE="~/LD_scores"          # LD scores
LD_1000G="~/1000G/all_phase3"   # 1000 genomes phase 3 bfile (all_phase3 is a prefix for the bfiles, not folder)

### SELECTION ###
SEL_DIR="~/selection"           # sexually-antagonistic selection results
