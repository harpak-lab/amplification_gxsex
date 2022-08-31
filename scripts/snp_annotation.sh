#!/bin/sh

# get phenotype
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

source config.R
GWAS_DIR=${GWAS_DIR}/$PHENO
mkdir -p $GWAS_DIR/annotation
cd $SCRATCH/ensembl-vep

# CLUMP
declare -a arr=("female" "male")
for sex in "${arr[@]}"
do
plink \
    --bfile $LD_1000G \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $GWAS_DIR/${sex}_all.${PHENO}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $GWAS_DIR/annotation/${sex}_$PHENO
# format clumped file
sed 's/ \+/\t/g' $GWAS_DIR/annotation/${sex}_${PHENO}.clumped | cut -f4 > $GWAS_DIR/annotation/${sex}_${PHENO}.clumped.ids

# ENSEMBL annotation
./vep -i $GWAS_DIR/annotation/${sex}_${PHENO}.clumped.ids \
--dir_cache $ENSEMBL_CACHE --cache -o $GWAS_DIR/annotation/${sex}_${PHENO}_annotation.txt \
--force_overwrite --tab --symbol --nearest symbol --pick \
--fields "Uploaded_variation,Location,Allele,SYMBOL,Gene,NEAREST,Feature,Feature_type,Consequence" --stats_text
done

