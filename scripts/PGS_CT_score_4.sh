#!/bin/sh

while getopts p:s: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $SET

PGS_DIR=$GWAS_DIR/$PHENO/PGS_$SET

# thresholds to use
echo "1 0 1" > range_list 
echo "0.01 0 0.01" >> range_list
echo "1e-5 0 1e-5" >> range_list
echo "1e-8 0 1e-8" >> range_list

### ADDITIVE ###
declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    # clumping of base data
    plink \
        --bfile $LD_1000G \
        --clump-p1 1 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump $PGS_DIR/${sex}_train.${PHENO}.glm.${TYPE} \
        --clump-snp-field ID \
        --clump-field P \
        --out $PGS_DIR/${sex}_additive.${PHENO}
    
    # extract SNP id from clumping results
    awk 'NR!=1{print $3}' $PGS_DIR/${sex}_additive.${PHENO}.clumped >  $PGS_DIR/${sex}_additive.${PHENO}.clumped.snpid
    # get snp and p-values from training set
    awk '{print $3,$13}' $PGS_DIR/${sex}_train.${PHENO}.glm.${TYPE} > $PGS_DIR/${sex}_additive.${PHENO}.snp_pvalue
    
    # calculate PRS
    ## 3:SNP    6:effective allele (A1)  12: effect size
    # score is a sum across SNPs of the # of ref alleles (0,1,or 2) at that SNP multiplied by the score for that SNP
    plink \
        --bfile $PGS_DIR/${PHENO}_test \
        --score $PGS_DIR/${sex}_train.${PHENO}.glm.${TYPE} 3 6 10 header \
        --q-score-range range_list $PGS_DIR/${sex}_additive.${PHENO}.snp_pvalue \
        --extract $PGS_DIR/${sex}_additive.${PHENO}.clumped.snpid \
        --out $PGS_DIR/${sex}_additive_${PHENO}
done

### Covariance aware (mash) ###
declare -a arr=("female" "male")
for sex in "${arr[@]}"
do
    plink \
        --bfile $LD_1000G \
        --clump-p1 1 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump $PGS_DIR/${sex}_pseudoP_pgs.${PHENO}.txt \
        --clump-snp-field ID \
        --clump-field P \
        --out $PGS_DIR/${sex}_mash.${PHENO}
    
    awk 'NR!=1{print $3}' $PGS_DIR/${sex}_mash.${PHENO}.clumped >  $PGS_DIR/${sex}_mash.${PHENO}.clumped.snpid
    awk '{print $1,$5}' $PGS_DIR/${sex}_pseudoP_pgs.${PHENO}.txt > $PGS_DIR/${sex}_mash.${PHENO}.snp_pvalue

    plink \
        --bfile $PGS_DIR/${PHENO}_test \
        --score $PGS_DIR/${sex}_pseudoP_pgs.${PHENO}.txt 1 2 3 header \
        --q-score-range range_list $PGS_DIR/${sex}_mash.${PHENO}.snp_pvalue \
        --extract $PGS_DIR/${sex}_mash.${PHENO}.clumped.snpid \
        --out $PGS_DIR/${sex}_mash_${PHENO}
done

