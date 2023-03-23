# Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits
Carrie Zhu, Matthew J. Ming, Jared M. Cole, Michael D. Edge, Mark Kirkpatrick, Arbel Harpak

Provided below are instructions and details for scripts used to generate the results and figures in ["Amplification is the Primary Mode of Gene-by-Sex Interaction in Complex Human Traits"](https://www.biorxiv.org/content/10.1101/2022.05.06.490973v1).
<br> <br/>
## Outline:
1. Download GWAS summary statistics files from Zenodo - [additive (10.5281/zenodo.7508246)](https://doi.org/10.5281/zenodo.7508246) and [sex-specific (10.5281/zenodo.7222725)](https://doi.org/10.5281/zenodo.7222725)
2. Download [scripts](/scripts) directory which contains scripts and files used to replicate the results and figures
2. Install software: [plink 1.9](https://www.cog-genomics.org/plink/), [plink 2.0](https://www.cog-genomics.org/plink/2.0/), [ldsc](https://github.com/bulik/ldsc), [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/script/index.html), [mashr](https://github.com/stephenslab/mashr)
3. Follow file path outline shown in [directory_outline](/directory_outline/), (directory names in [config.R](/scripts/config.R))
4. Update configuration file with own file paths: [config.R](/scripts/config.R)
5. Scripts often require items generated from previous scripts, so it is best to follow the documentation in order
<br> <br/>
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

Other *R* packages used are listed in [R_packages](/R_packages.txt)
<br> <br/>
## Documentation
### General Flags
- ```-p``` or ```--pheno``` flag indicates the phenotype code; ex. hip_circ  
- ```-n``` or ```--name``` flag indicates a formated phenotype name, often used for a title of a plot; ex. Hip circumference  
<br> <br/>
### Phenotype files
Phenotype files are obtained from UK Biobank and renamed pheno_(phenotype code).txt. A list of phenotype codes and formatted names for labels are provided in pheno_names.txt. A list of phenotype codes and UKBB data fields are provided in pheno_ids.txt
<br> <br/>
### Single snp analysis
##### Miami plots from GWAS summary statistics estimated in males and females only
Download sex-specific summary statistics from [Zenodo (10.5281/zenodo.7222725)](https://doi.org/10.5281/zenodo.7222725) to the corresponding $GWAS_DIR/[phenotype code] folder. 
Create miami plot.  
Code Example: ```./manhattan.R -p arm_fatfree_mass_L -n "Arm fat-free mass (L)"```

##### SNP annotation for list of SNPs after clumping and thresholding, removing SNPs with p-value>5e-8, pairwise LD threshold r<sup>2</sup>>0.1, or within 250kb  
###### 1000G Data  
In the $LD_1000G directory, download 1000G phase 3 GrCh 37 genotype data ([all_phase3 bfiles](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg)) to $LD_1000G  
Download [eur_ids.txt](/intermediate_files/eur_ids.txt), which contains subpopulation codes CEU and GBR to be kept in the sample  
Filter the bfiles using the following
```plink2 --pfile all_phase3 --chr 1-22 --max-alleles 2 --keep eur_ids.txt --rm-dup exclude-all --king-cutoff 0.0442 --make-bed --out all_phase3```  
###### Annotation  
Generate list of annotations after clumping  
Code Example: ```./snp_annotation.sh -p height```  
<br> <br/>
### LD Score regression
##### Estimate heritability and genetic correlation from sex-specific GWAS summary statistics
Download the pre-computed LD scores from the [ldsc tutorial](https://github.com/bulik/ldsc) (Bulik-Sullivan et al. 2015) to $LD_SCORE.  
Estimate heritability and genetic correlation.  
Code Example: ```./ldsc_basic.sh -p height``` 

#### Create plot for Figure 1
Download [ldsc_results.txt](/intermediate_files/ldsc_results.txt) to $LDSC_FILE, which contains all the sex-specific heritability estimates and male-female genetic correlations, estimated in the previous step.  
relative_h2.txt, which is used in nontrivial.R, is created in this step and placed in $LD_FILE.  
Code: ```./r2_by_h2.R```
<br> <br/>
### Multivariate adaptive shrinkage (mashr)
The ```-s``` or ```--set``` flag may be used in the three scripts below. The flag indicates whether or not to produce results specifically for the polygenic score pipeline (**Text S8**), which uses a smaller sample size for the cross-validation procedure. Use the flag and input the set number [1-20] only if using the PGS pipeline, otherwise do not use the flag.  
The examples below generate mixture weights and posterior estimates.  

Read in the data to create a matrix of effect estimates and standard errors by sex from sex-specific GWAS summary statistics.  
Code Example: ```./mash_setup.R -p height```

Incorporate pre-specified hypothesis covariance matrices and estimate mixture proportions.  
Code Example: ```./mash_100.R -p height```

Compute posterior estimates: posterior mean, standard deviation, weight, and lfsr.  
Code Example: ```./mash_posterior.R -p height```

#### Plots using mixture weights
Create overall and compact heatmaps of mixture weights. Plots for **Data S1-27**.   
Code Example: ```./mash_heatmap.R -p height -n Height```

Plot for **Fig. S8**.  
Download [pheno_names.txt](/intermediate_files/pheno_names.txt) to $PHENO_DIR. Create a text file named sex_ids.txt with two columns, ["IID", "sex"], which contains all the sample IIDs and corresponding sex in (0,1) format.  
Code: ```./phenovar_by_phenomean.R```

Plot for **Fig. 4A**.  
Download [mash_weights.txt](/intermediate_files/mash_weights.txt) to $GWAS_DIR, which summarizes weights from all traits.  
Code: ```./phenovar_by_amplification.R```

Plot for **Fig. S7**  
This script uses mash_weights.txt and pheno_names.txt.  
Code: ```nontrivial.R```

#### Test different p-value threshold (Methods and Fig. S3)
The ```-m``` or ```--mode``` flag may be used in the three scripts below. To keep the same random sample size for input to *mash*, use the flag, with the parameter '_same'. Otherwise, do not use the flag. 

```./mash_setup.R``` needs to be run first for the trait.  
Code Example: ```mash_p_threshold.R -p height``` or ```mash_p_threshold.R -p height -m _same```  

Plots for **Fig. S4**.  
Code Example: ```mash_pvalue_plot.R -p height -n Height``` or ```mash_pvalue_plot.R -p height -n Height -m _same```  

Plots for **Fig. S4**.  
In $GWAS_DIR, Download noeffect_weight.txt and noeffect_weight_same.txt, which has the weight on the no effect matrix for each phenotype and p-value threshold used in this particular analysis.  
Code Example: ```mash_pvalue_null_plot.R``` or ```mash_pvalue_null_plot.R -m _same```  

#### Environmental variance simulation for *mash*
##### Create a matrix of 30K individuals and 20K genotypes. 
In $QC_DIR, Download [maf_sample_20k.txt](/intermediate_files/maf_sample_20k.txt), which contains a random sample of 20K mean allele frequencies from UK Biobank [Resource 1967](https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1967).  
This script creates simulation_matrix_k_5k.RData, simulation_matrix_k_10k.RData, simulation_matrix_k_15k.RData, and simulation_matrix_k_20k.RData.  
Code: ```environ_matrix.R```

Estimate effect sizes and standard errors for input into mash. Flags ```-i```, ```-g```. and ```-e```, indicate parameters for number of causal SNPs, heritability, and female to male environmental variance ratio. We used the following parameters:  
```-i``` or ```--snps```: [100, 1000]  
```-g``` or ```--heritability```: [0.05, 0.5]  
```-e``` or ```--environment```: [1, 1.5, 5]  
Code Example: ```environ_small.R -i 100 -g 0.5 -e 1.5```

Same as the one above, but only for a causal SNP sample size of 10k. Flags ```-g``` and ```-e``` are the same, and there is no ```-i```.   
Code Example: ```environ_large.R -g 0.05 -e 5```

*mash* fitting procedure to test if differences in environmental variance are captured. Results for all parameters from the previous two scripts are required before running this script. The results we obtained are posted in [intermediate_files](/intermediate_files/environ_var_sim/) and are used for plotting.  
Code: ```environ_mash.R```

Produce heatmap plots as depicted in **Fig. S5**.  
Download matrice_names.txt, which is a list of all hypothesis covariance matrices used, to $GWAS_DIR.  
Code Example: ```environ_heatmap.R```
<br> <br/>
### Polygenic Score  
Create 20 test sets to be used for cross validation. This step will create folders PGS_1-20 and place a test set in each folder.  
Create a file, QC_ids.txt, containing all the sample IDs in one column, with column title 'IID', and move to $PHENO_DIR. The sample IDs are obtained after performing sample quality checks on UK Biobank data. The quality checks are further detailed in [QC.sh](/scripts/QC.sh).  
Code Example: ```PGS_testset_1.R -p height```

For the additive, standardized by sex model, first perform within sex standardization for each phenotype. These files will have a _std suffix.  
Code Example: ```standardize_pheno.R -p height```

For the following scripts, provide the set (cross-fold) number [1-20] with the ```-s``` flag.  

Generate sex-specific GWAS summary statistics estimated in both sexes and sexes separately.
You will need quality checked genotype files (from [QC.sh](/scripts/QC.sh)) labeled ukb_imp_chr(1-22)_v3_11. The covariate file contains the FID, IID, first 10 PCs, sex, and birth year. This requires UK Biobank access.  
Code Example: ```PGS_GWAS_2.sh -p height -s 1```

Continue procedure using the following scripts adding the set number to the ```-s``` flag:  
Code Example: ```./mash_setup.R -p height -s 2```  
Code Example: ```./mash_100.R -p height -s 2```  
Code Example: ```./mash_posterior.R -p height -s 2```  
Code Example: ```./lfsr_to_pvalue.R -p height -s 2```  

Clumping and thresholding procedure. Download range_list.txt to the directory your scripts are located. These contain the p-value ranges for the CT procedure.  
Code Example: ```PGS_CT_score_4.sh -p height -s 2```

Print out predictions using PGS based on four models described in **Text S8**. The result will output in the $GWAS_DIR/[phenotype code]/PGS_1 folder. It will include columns for the R2, incremental R2, and model (add=additive, both-sex; as=additive, standardized by sex; ss=sex-specific, covariance-naive; and mash=sex-specific, covariance-aware). Along with the .profile files, you will need the phenotype and covariate files.  
Code Example: ```PGS_predict_20.R -p height -s 2```

#### Plots
Plot for **Fig. S12**.  
Download [pgs_20.txt](/intermediate_files/PGS_20.txt) to $GWAS_DIR, which provides a summary of means and standard errors of the PGS results for all phenotypes from the script above.   
Code: ```PGS_plot.R```

Plot for **Fig. 2A-H and S11**.  
Move the best .profile file (female_additive..., male_additive..., male_mash..., female_mash..., both_sex_additive...) from set 1 to $GWAS_DIR/[phenotype code]. There should be a total of 5 .profile files for each phenotype in pheno_names.txt.    
Code: ```pheno_pgs.R```  
Note: We used a 5 fold cross validation for this script and the scripts in "Testosterone as an modulator of amplification" section, with a 50K total test set. 

Plot for **Fig. 2I,J**.  
This script uses sexspecific_pheno_pgs_lm.txt which is produced in the script before. The Spearman correlation between the male and females panels is printed out.  
Code: ```pheno_pgs_overall.R```

Plot for **Fig. 4B**.  
Download [PGS_20_all.txt](/intermediate_files/pgs_20_all.txt), which contains all R2 from the 20 folds for the 27 phenotypes.  
Code: ```pgs_20_utility.R```
<br> <br/>
### Testosterone as an modulator of amplification
Plots for **Fig. 5A and S13,14**. If using PGS estimated from sex-specific summary statistics (**Fig. S14**), input 'sex-specific' for the ```-m``` or ```--mode``` flag. Otherwise, do not use that flag.  
Code Example: ```G_testosterone.R -p height -n Height```

Plots for **Fig. S15A**.  
Code Example: ```G_testosterone_pgs.R -p height -n Height```

Plot for **Fig. 5B**.  
Code: ```G_corr_testosterone.R```

Plot for **Fig. S15B**.  
Code: ```G_corr_testosterone_pgs.R```

Plot for **Fig. S16**.  
Code: ```G_corr_testosterone_age.R```
<br> <br/>
### Model of shared amplification
Plot for **Fig. 6**.  This script uses pheno_meanvar.txt, which was produced in phenovar_by_phenomean.R. You will also need [ldsc_results.txt](/intermediate_files/ldsc_results.txt)  
Code: ```gen_env_bootstrap.R```
<br> <br/>
### Sexually-Antagonistic Selection
Plot for **Fig. 7C,D**.   
Download [sex_selection](/intermediate_files/sex_selection/) files in $SEL_DIR.  
Code: ```sex_selection_plot.R```

Plot for **Fig. S19,20**.  
Code Example: ```sex_selection_supplement.R```

Analysis for sexually-antagonistic selection was done by [Matthew J. Ming](https://github.com/MattJMing). Code can be found in the [matt_selection](/matt_selection/) directory. Follow the associated README.  
<br> <br/>
### More Supplementary Scripts
#### Simulation of covariance structure
The script is similar to that in [Environmental variance simulation for *mash*](#environmental-variance-simulation-for-mash). We used the defaults for flags: -i and -e. The flag for the pre-specified matrix, -m, must be 4 digits (only single digit matrices, but you can directly customize the matrices in the script). The -a flag refers to the weight of the nonnull matrix. We repeated this step 100 times, using random seed, -s 1-100.  
 
For Fig. S6A, use ```-a 0```. For Fig. S6B, use ```-a 14 -m 4221```. For Fig. S6C, use ```-a 14 -m 1000```.  
Code Example: ```sim_covariance.R -a 14 -s 1 -m 4221```  

Run mash. Use RData file names generated from previous script as the input.  
Code Example:  ```mash_simcov.R -n mash_4221_1```  

Plot for **Fig. S6**.  
Concatenate all 100 .txt files generated from the script above by column and rename. Use the new name as the input. Our results are located in [intermediate_files](/intermediate_files/simulation%20covariance/).  
Code Example: ```heatmap_simcov.R -n mash_4221_all```

#### Characterizing GxSex based on independent analysis of individual sex-heterogenous SNPs 
Use one the following phenotypes as the flag: [height, bmi, creatinine, IGF1, systolicBP_auto]. We used the defaults for flags: -i and -e. We repeated this step 10 times, using random seed, -s 1-10.       
Code Example: ```sim_sexhet.R -p height -s 1```  

Results for **Table S6**.  
Input the phenotype code and name of RData file generated from the previous step.  
Code Example: ```sexhet_p.R -p height -n mash_height_1```

#### Competing models for sex differences in trait variance
Plot for **Fig. S18B**. 
Similar to script in [Model of shared amplification](#model-of-shared-amplification).  
Code: ```gen_env_models.R```




