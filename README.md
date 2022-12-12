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

**rMVP** _version 1.0.0_ --see https://github.com/xiaolei-lab/rMVP

## Data

In this project, the Haplotype Map file `geno.hmp.txt` consists of 11 columns plus the information associated with a single SNP of 768 genotypes. 

The phenotype file `pheno.txt` includes the protein contains of 768 genotypes.

## GAPIT

**GAPIT (Genomic Association and Prediction Integrated Tool)** implemented a series of methods for GWAS. The GWAS models include General Linear Model (GLM), Mixed Linear Model (MLM or Q+K), Compressed MLM (CMLM), Enriched CMLM, SUPPER, Multiple Loci Mixed Model (MLMM), FarmCPU and BLINK. 

GAPIT accepts multiple input data formats, including both numeric, hapmap, and PLINK genotype formats. GAPIT will produce comprehensive reports to interpret data and results in publication ready formats.

### Usage

*Create a working directory `GAPIT_result` and copy `geno.hmp.txt` and `pheno.txt` to the folder before run this script.*

#### Installation and load the package
    
    install.packages("devtools")
    devtools::install_github("jiabowang/GAPIT3",force=TRUE)
    library(GAPIT3)
    
#### Change the working directory

    setwd("/path/to/GAPIT_result/")
    
#### Import the data

    hapmap_geno<-read.table(file.choose(), head=F) # make sure header=F
    pheno<-read.table(file.choose(), head=TRUE)
    
#### The first analysis (no compression)

    dir.create("No-compression")
    setwd("/path/to/GAPIT_result/No-compression/")
    
    analysis1<-GAPIT(
      Y=pheno,
      G=hapmap_geno,
      SNP.impute="Major",
      PCA.total=3,
      Major.allele.zero=T,
      group.from=768,
      group.to=768,
      group.by=1)
      
*The result output files were saved in the `No-compression` folder.*

#### THe second analysis (with compression)

    setwd("path/to/GAPIT_result/")
    dir.create("with-compression")
    setwd("path/to/GAPIT_result/with-compression/")
    
    analysis2<-GAPIT(
      Y=pheno,
      G=hapmap_geno,
      SNP.impute="Major",
      PCA.total=3,
      Major.allele.zero=T)

*The result output files were saved in the `with-compression` folder.*

#### The third analysis (use different kinship clustering methods to group individuals based on their kinship)

    setwd("path/to/GAPIT_result/")
    dir.create("parameters")
    setwd("path/to/GAPIT_result/parameters/")

    analysis3<-GAPIT(
      Y=pheno,
      G=hapmap_geno,
      SNP.impute="Major",
      kinship.cluster=c("complete","ward"),
      kinship.group=c("Mean","Max","Median"),
      PCA.total=3,
      Major.allele.zero=T,
      Model.selection=T)
      
*The result output files were saved in the `parameters` folder.*

#### The fourth analysis (Multiple Loci Mixied linear Model)

    setwd("path/to/GAPIT_result/")
    dir.create("MLMM")
    setwd("path/to/GAPIT_result/MLMM/")
    
    analysis4 <- GAPIT(
      Y=pheno[,c(1,2)],
      G=hapmap_geno,
      model="MLMM",
      PCA.total=3,
      file.output=T)
      
*The result output files were saved in the `MLMM` folder.*

#### The fifth analysis (Fixed and random model Circulating Probability Unification)

    setwd("path/to/GAPIT_result/")
    dir.create("FarmCPU")
    setwd("path/to/GAPIT_result/FarmCPU/")
    
    analysis5 <- GAPIT(
      Y=pheno[,c(1,2)],
      G=hapmap_geno,
      model="FarmCPU",
      PCA.total=3,
      file.output=T)

*The result output files were saved in the `FarmCPU` folder.*

## rMVP

**rMVP** is a Memory-efficient, Visualization-enhanced, and Parallel-accelerated Tool for GWAS.

### Usage

*Create a working directory `rMVP_result` and copy `geno.hmp.txt` and `pheno.txt` to the folder before run this script.*

#### Installation and load the package

    install.package ("rMVP")
    library (rMVP)

#### Change the working directory

    setwd("/path/to/rMVP_result")
    
#### Import data

    MVP.Data(fileHMP="geno.hmp.txt",
         filePhe="pheno.txt",
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.hmp",
         #priority="speed",
         #maxLine=10000
    )


*By using the `MVP.Data` function, we can generate the a log file, genotype (mvp.hmp.geno.desc), phenotype (mvp.hmp.phe) and map mvp.hmp.geno.map data, which are saved in the `rMVP_result` folder.*

#### Input data

    genotype <- attach.big.matrix("mvp.hmp.geno.desc")
    phenotype <- read.table("mvp.hmp.phe",head=TRUE)
    map <- read.table("mvp.hmp.geno.map" , head = TRUE)
    
#### Run the GWAS 

*Three models are included in MVP package: General Linear Model (GLM), Mixed Linear Model (MLM), and FarmCPU.**

    for(i in 2:ncol(phenotype)){
    imMVP <- MVP(
      phe=phenotype[, c(1, i)],
      geno=genotype,
      map=map,
      #K=Kinship,
      #CV.GLM=Covariates,  ##if you have additional covariates, please keep there open.
      #CV.MLM=Covariates,
      #CV.FarmCPU=Covariates,
      #nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
      #nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
      nPC.FarmCPU=3,
      priority="speed",   ##for Kinship construction
      ncpus=5,
      vc.method="BRENT",  ##only works for MLM
      maxLoop=10,
      method.bin="FaST-LMM", #"static" (#only works for FarmCPU)
      #permutation.threshold=TRUE,
      #permutation.rep=100,
      threshold=10,
      method=c("FarmCPU")
      )
      gc()
    }
    
### Output

After finishing the analysis, we can generate *phenotype distribution plot*, *SNP-density plot*, *PCA plot (2D and 3D)*, *manhattan plot in Circular fashion*, *Q-Q plot*, and *manhattan plot in Rectangular fashion for single trait or method*, we also generated a log file, you can find all of them in the `rMVP_result` folder.
