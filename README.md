# BVG_7003_GWAS

This page is created by Jiaxu Wu

There are some R scripts for Genome Wide Association Study (GWAS) in in BVG_7003 @ Université Laval.

All the scripts used in this project were collected in `script` folder. The recoreds and results were saved in the `GAPIT_result` and `rMVP_result` folder.
We also attached test files for **GWAS** in the `data_input` folder, which includes `geno.hmp.txt` and `pheno.txt`.

**GWAS** is a research approach used to identify genomic varients that are statistically associated with a risk for the particular trait. The methods involves surverying genomes of large numbers of individuals. 


## Requirement

### Environment

**R** _version 4.1.1_ --see https://cran.r-project.org/

### IDE 

**RStudio Desktop** _version 2022.07.2_ --see https://posit.co/products/open-source/rstudio/

### Package

**GAPIT** _version 3.1.0_ --see https://github.com/jiabowang/GAPIT3

**rMVP** _version 1.0.0 --see https://github.com/xiaolei-lab/rMVP

## GAPIT

GAPIT (Genomic Association and Prediction Integrated Tool) 
implemented a series of methods for GWAS. The GWAS models include General Linear Model (GLM), Mixed Linear Model (MLM or Q+K), Compressed MLM (CMLM), Enriched CMLM, SUPPER, Multiple Loci Mixed Model (MLMM), FarmCPU and BLINK. 
GAPIT accepts multiple input data formats, including both numeric, hapmap, and PLINK genotype formats. GAPIT will produce comprehensive reports to interpret data and results in publication ready formats.

## rMVP

