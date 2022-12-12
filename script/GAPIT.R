# to run GAPIT we need to access some tools that are located in other packages. The following lines of 
# code will install the required packages in to your R library and then will import them so that GAPIT
# may access them. The easiest way to do this is to highlight all of the following code and select run. 
# Once you have done this, you will only need to use the library functions to access them at future dates.


install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT3)

### Change working directory
getwd()
setwd("/path/to/GAPIT_result/")


#these lines of code install the actual GAPIT functions from the maize genetics website
#source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

# Now we need to import the data files that we will be using for doing the analysis. Best to have
# data saved in a text format. If you are using R studio you can import the files using the import 
# wizard located under the tools tab at the top.This may be the easiest as it walks you through the 
# process. You may perfer the command line method as it gives you more controll over the import of the 
# data. Using the following command, file.choose(), opens a browser allowing you to select the name of 
# the file in the approriate directory. Select the file and the name of the file will appear in the
# consol window. Copy and paste this name into the read.table command.


# If your data is in the HapMap format, you will only need the two files-the phenotypic and the hapmap
# file. We will do this analysis first for the sake of time. 

hapmap_geno<-read.table(file.choose(), head=F) # make sure header=F
pheno<-read.table(file.choose(), head=TRUE)

#myCV <- read.table("Q3_ADMIX.txt", head=TRUE) ## K = 3 provided a good assessment of population structure.

# good idea to check our pehontype data to make sure the file strucutre is correct and how the data
# is distributed, checking for outliers

str(pheno) # gives us information on the object, in this case a data frame and other information
hist(pheno$protein) # creates a histogram plot of our data, things look pretty good, we have a
 #lines that are bit extreme but with 775 lines we would expect about 27 to be at least 3 sd's away from
 # the mean. 

#some basic statistics to look at 
mean(pheno$protein) # 12.23
range(pheno$protein) # 9.5 to 15.8
sd(pheno$protein) # .88
which(is.na(pheno$protein)) # look for lines with missing data, there should be none which is confimred
  # by the result in the console "integer(0) meaning the number of NA values was 0

# the following is a very basic analysis. Make sure to change directory to where the results
# will be saved. In R studio go to session at the top and select "set working directory"
# and select "choose directory"- this will allow you to browse to the correct folder in which you
# would like to save the results. Once you locate the correct folder, highlight it and hit select.

# or you can do it manually using the something like the following
dir.create("No-compression")
setwd("/path/to/GAPIT/No-compression/")

# first analysis where compression is not used in the model, this is done by setting the group.from 
# and group.to equal to the size of your population and then setting the group.by to 1 so that each 
# entry is consider its own "group." Again, refer to manual for complete description.
analysis1<-GAPIT(
  Y=pheno,
  G=hapmap_geno,
  SNP.impute="Major",
  PCA.total=3,
  Major.allele.zero=T,
  group.from=768,
  group.to=768,
  group.by=1)


# now we use compression to see how that effects the outcome of the analysis. The default settings 
# for compression are to group by 10. You can change this to whatever value you feel like. It is
# important to note that not all values will work depending on your population. You will get the error
# that "matrix is singular, select another level." You can keep selecting levels until you find one 
# that works. However, compression may not make that big of a difference for your analyses. To read
# more about it see "Mixed linear model approach for genome-wide association studies" Zhang et al. 
# Nature Genetics 2010

# if you want to use the default settings for compression, you do not need to specify anything.

# change the directory where the output will be stored
setwd("path/to/GAPIT/")
dir.create("with-compression")
setwd("path/to/GAPIT/with-compression/")

analysis2<-GAPIT(
  Y=pheno,
  G=hapmap_geno,
  SNP.impute="Major",
  PCA.total=3,
  Major.allele.zero=T)

# the differences were pretty small
# no compression used 
# SNP    FDR_Adjusted_P-values
# 12_10811  	1.74E-06
# 12_10199		1.74E-06
# 12_10575		1.74E-06
# 12_20685		1.75E-06
# 12_11437		1.17E-05
# 12_30301		0.000439948

# with compression
# SNP  	FDR_Adjusted_P-values
# 12_10811	1.60E-06
# 12_10199	1.60E-06
# 12_10575  1.60E-06
# 12_20685	2.36E-06
# 12_11437	1.40E-05
# 12_30301	0.000599332

# change to a location that works for you
setwd("path/to/GAPIT/")
dir.create("parameters")
setwd("path/to/GAPIT/parameters/")

analysis3<-GAPIT(
  Y=pheno,
  G=hapmap_geno,
  SNP.impute="Major",
  kinship.cluster=c("complete","ward"),
  kinship.group=c("Mean","Max","Median"),
  PCA.total=3,
  Major.allele.zero=T,
  Model.selection=T)

# this analysis used different kinship clustering methods to group individuals based on their kinship.
# For kinship.group in used three different methods to derive kinship among the groups. The final line
# of code implements the model selection feature of GAPIT. This feature uses the Bayesian information
# criteria (BIC) to select the optimal number of principal components to use in the model. This is done
# because the degree of population structure can vary form trait to trait, therefore you don't always need
# to use the same number of principal components for each trait. 

# Below is the output from this analysis, note the difference in FDR p-values. The optimal model
# did not use any PCs, the best method for deriving kinship between groups was "Max", and the 
# kinship clustering method that was used was "Complete". The number of groups did not change, neither
# did the estimate of heritability.

# GWAS output
# SNP  	FDR_Adjusted_P-values
# 12_10811	7.67E-07
# 12_10199	7.67E-07
# 12_10575	7.67E-07
# 12_20685	1.54E-06
# 12_11437	9.17E-06
# 12_30301	0.000428139

# the FDR p-values that we got were a bit smaller than our analysis using compression with the 
# default settings. This example illustrates the value of trying a few analyses.

#### Multi-locus analysis
# change to a location that works for you

setwd("path/to/GAPIT/")
dir.create("MLMM")
setwd("path/to/GAPIT/MLMM/")

analysis4 <- GAPIT(
  Y=pheno[,c(1,2)],
  G=hapmap_geno,
  model="MLMM",
  PCA.total=3,
  file.output=T)



# change to a location that works for you
setwd("path/to/GAPIT/")
dir.create("FarmCPU")
setwd("path/to/GAPIT/FarmCPU/")

analysis5 <- GAPIT(
  Y=pheno[,c(1,2)],
  G=hapmap_geno,
  model="FarmCPU",
  PCA.total=3,
  file.output=T)











