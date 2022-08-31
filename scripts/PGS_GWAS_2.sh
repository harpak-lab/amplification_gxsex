#!/bin/sh

source config.R
while getopts p:m: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $MODE
PGS_DIR=$GWAS_DIR/$PHENO/PGS_$SET
mkdir -p $PGS_DIR/{both_sex,female,male}

# GWAS for both sex, female, and male specfic
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p
cd $PHENO_DIR
declare -a arr=("female" "male")
for i in {1..22}
do
    # both sex
    plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_DIR/ukb_imp_chr${i}_v3_11 \
    --remove $PGS_DIR/${PHENO}_female_testIIDs.txt $PGS_DIR/${PHENO}_male_testIIDs.txt \
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PGS_DIR/both_sex/both_sex_${i} \
    --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize
    # sex-specific
    for sex in "${arr[@]}"
    do 
        plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
        --pfile $QC_DIR/ukb_imp_chr${i}_v3_11 --keep-${sex}s \
        --remove $PGS_DIR/${PHENO}_female_testIIDs.txt $PGS_DIR/${PHENO}_male_testIIDs.txt \
        --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PGS_DIR/${sex}/${sex}_${i} \
        --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize 
    done
done

# create test bfile
plink2 --memory 64000 --threads 16 --pfile $QC_DIR/ukb_imp_all_v3_11 --keep $PGS_DIR/${PHENO}_female_testIIDs.txt $PGS_DIR/${PHENO}_male_testIIDs.txt \
--make-bed --out $PGS_DIR/${PHENO}_test

# combine GWAS results
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    cd $PGS_dir/$sex
    cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $PGS_dir/${sex}_train.${PHENO}.glm.linear
done