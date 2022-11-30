install.packages("rMVP")
library(rMVP)

### Change working directory

getwd()

setwd("/Users/jiaxuwu/BVG_7003/GWAS/")

dir.create("rMVP")
setwd("/Users/jiaxuwu/BVG_7003/GWAS/rMVP")


### Import genotypic and phenotypic data
MVP.Data(fileHMP="geno.hmp.txt",
         filePhe="pheno.txt",
         sep.hmp="\t",
         sep.phe="\t", # separator of phenotype file
         SNP.effect="Add",
         fileKin=TRUE, # TRUE or FALSE, if TRUE, kinship matrix represents relationship among individuals will be calculated
         filePC=TRUE, # TRUE or FALSE, if TRUE, principal component analysis will be performed
         out="mvp.hmp",
         #priority="speed",
         #maxLine=10000
)


###Kinship
MVP.Data.Kin("mvp.hmp.geno.desc", TRUE, mvp_prefix='mvp', out='mvp')
Kinship <- attach.big.matrix("mvp.kin.desc")

###PCA
MVP.Data.PC("mvp.hmp.geno.desc", TRUE, out='mvp', pcs.keep=5)


### Data input
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
map <- read.table("mvp.hmp.geno.map" , head = TRUE)


####Run GWAS
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



