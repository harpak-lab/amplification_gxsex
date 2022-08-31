# Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits
Carrie Zhu, Matthew J. Ming, Jared M. Cole, Mark Kirkpatrick, Arbel Harpak

Provided below are instructions and details for scripts used to generate the results and figures in ["Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits"](https://www.biorxiv.org/content/10.1101/2022.05.06.490973v1).

## Outline:
1. Download GWAS summary statistics files from [UTBox](https://utexas.box.com/s/ef25198jq6owpcq5j2wq6najovlq75b8)
2. Download Final_scripts directory which contains scripts and files used to replicate the results and figures
2. Install software: [plink 1.9](https://www.cog-genomics.org/plink/), [plink 2.0](https://www.cog-genomics.org/plink/2.0/), [ldsc](https://github.com/bulik/ldsc), [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/script/index.html), [mashr](https://github.com/stephenslab/mashr)
3. Follow file path outline shown below or in directory_outline, (directory names in config.R)
4. Update configuration file with own file paths: config.R
5. For each section in Documentation, it is best to follow the code in order

## Software
*plink v1.9 beta* (Purcell, S. & Chang, C 2021)  
*plink v2.0 alpha* (Purcell, S. & Chang, C 2020)  
- Download both to the directory containing all the scripts
*LDSC v1.0.1*  (Bulik-Sullivan et al. 2015)  
Ensembl command line *variant effect predictor (VEP) v106* (McLaren et al. 2016)  
- We used the command line VEP tool to annotate SNPs, following documentation listed on the website 
- We downloaded cache files for human genome assembly GRCh37 using INSTALL.pl to $ENSEMBL_CACHE
- perl module Set::IntervalTree also needs to be installed to use the --nearest flag for VEP  
*mashr* package in R (Urbut, et al. 2019)  

## Documentation
### General Flags
- ```-p``` or ```--pheno``` flag indicates the phenotype code
- ```-n``` or ```--name``` flag indicates a formated phenotype name, often used for a title of a plot

### Phenotype files
Phenotype files are obtained from UK Biobank and renamed pheno_(phenotype code).txt. A list of phenotype codes and formatted names for labels are provided in pheno_names.txt

### Single snp analysis
##### Miami plots from GWAS summary statistics estimated in males and females only
Download sex-specific summary statistics from [UTBox](https://utexas.box.com/s/ef25198jq6owpcq5j2wq6najovlq75b8) in the corresponding $GWAS_DIR/[phenotype code] folder
Code Example: ```./manhattan.R -p arm_fatfree_mass_L -n "Arm fat-free mass (L)"```
<br> <br/>
##### SNP annotation for list of SNPs after clumping and thresholding, removing SNPs with p-value>5e-8, pairwise LD threshold r<sup>2</sup>>0.1, or within 250kb  
In the $LD_1000G directory, download 1000G phase 3 genotype data (all_phase3 files), which were created using the following code:  
```plink2 --pfile all_phase3 --chr 1-22 --max-alleles 2 --keep eur_ids.txt --rm-dup exclude-all --king-cutoff 0.0442 --make-bed --out all_phase3```  
eur_ids.txt contained subpopulation codes CEU and GBR on separate lines to be kept in the sample  

Code Example: ```./snp_annotation.sh -p height```  

### LD Score regression
##### Estimate heritability and genetic correlation from sex-specific GWAS summary statistics
Download the pre-computed LD scores from the [ldsc tutorial](https://github.com/bulik/ldsc) (Bulik-Sullivan et al. 2015) to $LD_SCORE.  
Code Example: ```./ldsc_basic.sh -p height``` 

#### Create plot for Figure 1
Download ldsc_results.txt to $LDSC_FILE, which contains sex-specific heritability estimates and male-female genetic correlations, estimated in the previous step.  
relative_h2.txt, which is used in nontrivial.R, is created in this step and placed in $LD_FILE
Code: ```./r2_by_h2.R```

### Multivariate adaptive shrinkage (mashr)
The ```-s``` or ```--set``` flag may be used in the three scripts below. The flag indicates whether or not to produce results specifically for the polygenic score pipeline (**Text S1**), which uses a smaller sample size for the cross-validation procedure. Use the flag and input the set number [1-5] only if performing the PGS pipeline, otherwise do not use the flag. 

Read in the data to create a matrix of effect estimates and standard errors by sex from sex-specific GWAS summary statistics.  
Code Example: ```./mash_setup.R -p height```

Incorporate pre-specified hypothesis covariance matrices and estimate mixture proportions.  
Code Example: ```./mash_100.R -p height```

Compute posterior estimates: posterior mean, standard deviation, weight, and lfsr.  
Code Example: ```./mash_posterior.R -p height```

#### Plots using mixture weights
Create a overall and compact heatmaps of mixture weights. Plots for **Fig. S4**.   
Code Example: ```./mash_heatmap.R -p height -n Height```

Plot for **Fig. 4A**.  Download pheno_names.txt to $PHENO_DIR. Create a text file names sex_ids.txt with two columns [IID, sex], which contains the IIDs and corresponding sex in (0,1) format.
Code: ```./phenovar_by_phenomean.R```

Plot for **Fig. 4B**. Download mash_weights.txt to $GWAS_DIR, which summarizes weights from all traits.  
Code: ```./phenovar_by_amplification.R```

Plot for **Fig. S7**  
Code: ```nontrivial.R```

#### Test different p-value threshold (Methods and Fig. S3)
The ```-m``` or ```--mode``` flag may be used in the three scripts below. To keep the same random sample size for input to *mash*, use the flag, with the parameter '_same'. Otherwise, do not use the flag. 

```./mash_setup.R``` needs to be run first for the same trait.  
Code Example: ```mash_p_threshold.R -p height``` or ```mash_p_threshold.R -p height -m _same```  

Plots for **Fig. S3B**.  
Code Example: ```mash_pvalue_plot.R -p height -n Height``` or ```mash_pvalue_plot.R -p height -n Height -m _same```  

Plots for **Fig. S3A**. In $GWAS_DIR, Download noeffect_weight.txt and noeffect_weight_same.txt, which has the weight on the no effect matrix for each phenotype and p-value threshold used in this particular analysis.  
Code Example: ```mash_pvalue_null_plot.R``` or ```mash_pvalue_null_plot.R -m _same```  

#### Environmental variance simulation for *mash*
##### Create a matrix of 30K individuals and 20K genotypes. 
In $QC_DIR, Download maf_sample_20k.txt, which contained a random sample of 20K mean allele frequencies from UK Biobank [Resource 1967](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1967).  
Code: ```environ_matrix.R```

Estimate effect sizes and standard errors for input into mash. Flags ```-i```, ```-g```. and ```-e```, indicate parameters for number of causal SNPs, heritability, and female to male environmental variance ratio. We used the following parameters:  
```-i``` or ```--snps```: [100, 1000]  
```-g``` or ```--heritability```: [0.05, 0.5]  
```-e``` or ```--environment```: [1, 1.5, 5]  
Code Example: ```environ_small.R -i 100 -g 0.5 -e 1.5```

Same as the one above, but only for a causal SNP sample size of 10k. Flags ```-g``` and ```-e``` are the same, and there is no ```-i```.   
Code Example: ```environ_large.R -g 0.05 -e 5```

*mash* fitting procedure to test if differences in environmental variance is captured. Results for all parameters from the previous two scripts are required before running this script.  
Code: ```environ_mash.R```

Produce heatmap plots as depicted in **Fig. S6**. Download matrice_names.txt, which is a list of all hypothesis covariance matrices used, to $GWAS_DIR.  
Code Example: ```environ_heatmap.R```

### Polygenic Score  
Create five test sets to be used for cross validation. This step will create folders PGS_1-5 and place a test set in each folder.  Download QC_ids.txt to $PHENO_DIR.  
Code Example: ```PGS_testset_1.R -p height```

For the following four scripts, provide the set number [1-5] with the ```-s``` flag.  

Generate sex-specific GWAS summary statistics estimated in both sexes and sexes separately.  
Code Example: ```PGS_GWAS_2.sh -p height -s 1```

Continue procedure using the following scripts adding the set number to the ```-s``` flag:  
Code Example: ```./mash_setup.R -p height -s 2```  
Code Example: ```./mash_100.R -p height -s 2```  
Code Example: ```./mash_posterior.R -p height -s 2```  

Clumping and thresholding procedure.  
Code Example: ```PGS_CT_SCORE_4.sh -p height -s 2```

Print out predictions using PGS based on three methods described in **Text S1**.   
Code Example: ```PGS_predict_5_linear.R -p height -s 2```

#### Plots
Plot for **Fig. S14**. Download pgs_linear_results_five.txt, which provides a summary of correlations from the script before. pgs_combined_r2.txt is outputted.  
Code: ```PGS_plot_final.R```

Plot for **Fig. S16**. Move the best .profile file (female_additive..., male_additive..., male_mash..., female_mash..., both_sex_additive...) from set 1 to $GWAS_dir/[pheno_code] as printed by PGS_predict_5_linear.R. There should be a total of 5 .profile files for each phenotype in pheno_names.txt.    
Code: ```pheno_pgs.R```

Plot for **Fig. 2I,J**. This script uses sexspecific_pheno_pgs_lm.txt which is produced in the script before. The Spearman correlation between the male and females panels is printed out.  
Code: ```pheno_pgs_overall.R```

### Testosterone as an modulator of amplification
Plots for **Fig. S9,10**. If using PGS estimated from sex-specific summary statistics (**Fig. S10**), input 'sex-specific' for the ```-m``` or ```--mode``` flag. Otherwise, do not use that flag.  
Code Example: ```G_testosterone.R -p height -n Height```

Plots for **Fig. S11A**.  
Code Example: ```G_testosterone_pgs.R -p height -n Height```

Plot for **Fig. 5B**.  
Code: ```G_corr_testosterone.R```

Plot for **Fig. S11B**.  
Code: ```G_corr_testosterone_pgs.R```

Plot for **Fig. S12**.  
Code: ```G_corr_testosterone_age.R```

### Model of shared amplification
Plot for **Fig. 6**.  This script uses pheno_meanvar.txt, which was produce in phenovar_by_phenomean.R. Download ldsc_results.txt, if not already.  
Code: ```gen_env_bootstrap.R```

### Sexually-Antagonistic Selection
Plot for **Fig. 7C,D**.  Download the following RData files:  
fst_plot_testosterone.1e-05.RData  
fst_plot_protein_total.1e-05.RData  
zscore_plot.1e-05.RData  
Code: ```sex_selection_plot.R```

Analysis for sexually-antagonistic selection was done by Matthew J. Ming. Code can be found [here](https://github.com/MattJMing)

