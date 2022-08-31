#!/bin/sh

# set up environment
cd $LDSC_DIR
source activate ldsc

source config.R
# cmd line argument - phenotype
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO

# file paths
GWAS_DIR=${GWAS_DIR}/$PHENO
LD_FILE=${LDSC_FILE}/$PHENO

declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    # reformat summary statistics
    munge_sumstats.py --sumstats $GWAS_DIR/${sex}_all.${PHENO}.glm.linear \
    --snp ID --a2 AX --N-col OBS_CT --merge-alleles $LD_SCORE/w_hm3.snplist \
    --chunksize 500000 --out $LD_FILE/${PHENO}_${sex}

    # calculate heritability
    ldsc.py --h2 $LD_FILE/${PHENO}_${sex}.sumstats.gz \
    --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ \
    --out $LD_FILE/${PHENO}_${sex}_h2
done

# calculate correlation
ldsc.py --rg $LD_FILE/${PHENO}_female.sumstats.gz,$LD_FILE/${PHENO}_male.sumstats.gz \
--w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ \
--out $LD_FILE/${PHENO}_male_female

#### grep get heritability and genetic correlation
grep "h2: " $LD_FILE/*h2.log
grep "Genetic Correlation: " $LD_FILE/*male_female.log

