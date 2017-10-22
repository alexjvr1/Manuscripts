# Gradient Forest CHall

### Datasets:

2 Datasets for each regions: Neutral, Potential Adaptive Loci (identified by at least 2 programs)

1. CHN

2. CZ

3. CHS

4. CHS.VS

5. CHS.TI


# Gradient Forest Analysis

Here I aim to: 

1)determine the environmental and spatial determinants of variation in allele frequency across the Switzerland. 

2)compare global vs local effects

For thic chapter I just want to see how the adaptive loci are changing in allele frequency across each environmental variable, and determine which of the variables are more important. So for that i am using the same environmental variables that I used in the EAA tests. 

For Chapter 4 I will map the adaptive loci back to BioClim data. 


I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013


## Data

I'm using 2 different datasets: 

1. Potential adaptive loci (LFMM, BayENV2, XtX, PCAdapt, Bayescan)

2. "Neutral" loci (1000 loci randomly sampled from all - Potential Adaptive loci)

It's important to subsample the data, since the datasets run for a long time if they have more than 1000 loci. 


## Input files

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017

### 1. MEM variables

#### CHN
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHN <- read.csv("AllEnv.Data_CHN_20161010.csv", header=T)
head(env.data.CHN)
#extract x and y
CHN.xy <- env.data.CHN[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CHN.dxy <- distm(CHN.xy)
CHN.dxy <- as.dist(CHN.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CHN.th <- give.thresh(CHN.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CHN.nb1 <- mst.nb(CHN.dxy)
CHN.wh1 <- which(as.matrix(CHN.dxy)==CHN.th,arr.ind=TRUE)
plot(CHN.nb1,CHN.xy,pch=20,cex=2,lty=3)
lines(CHN.xy[CHN.wh1[1,],1],CHN.xy[CHN.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CHN.th 
#[1] 50457.46
CHN.nb1
Neighbour list object:
Number of regions: 19 
Number of nonzero links: 36 
Percentage nonzero weights: 9.972299 
Average number of links: 1.894737 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CHN.listw=nb2listw(CHN.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CHN.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
CHN.scores=scores.listw(CHN.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(CHN.scores,CHN.listw,nsim=999)
       stat  pval
1 1.0475936 0.001
2 0.9817148 0.001
3 0.8595811 0.001
4 0.7500000 0.001
5 0.6971080 0.001
6 0.4930407 0.013
7 0.3301269 0.037
8 0.1124090 0.242

#7 significant MEM eigenfunctions with positive correlations. OBS. use first half (4 in our case)

write.table (CHN.scores$vectors[,1], "CHN.scores_MEM1.txt") 
write.table (CHN.scores$vectors[,2], "CHN.scores_MEM2.txt") 
write.table (CHN.scores$vectors[,3], "CHN.scores_MEM3.txt") 
write.table (CHN.scores$vectors[,4], "CHN.scores_MEM4.txt") 
```

## Great circle distance between all sampling points

```
##Geographic distance between all CHN populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017")

#rdist.earth (in fields package) wants only long & lat
CHN_lon.lat <- cbind(env.data.CHN$long, env.data.CHN$lat)
CHN_lon.lat

#calculate great circle distances
distance.matrix.CHN <- rdist.earth(CHN_lon.lat)
summary(distance.matrix.CHN)
dim(distance.matrix.CHN)

#and use only the lower half of the matrix
upper.tri(distance.matrix.CHN)
distance.matrix.CHN[lower.tri(distance.matrix.CHN)]<-NA
distance.matrix.CHN

#change from matrix to dataframe
CHN.bli <- as.data.frame(distance.matrix.CHN)
head(CHN.bli)
colnames(CHN.bli) <- env.data.CHN$site
rownames(CHN.bli) <- env.data.CHN$site

CHN.bli[lower.tri(CHN.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
CHN.bli=as.data.frame(as.table(as.matrix(CHN.bli)))  #Turn into a 3-column table
CHN.bli
CHN.bli=na.omit(CHN.bli)  #Get rid of the junk we flagged above
CHN.bli
colnames(CHN.bli)<-c("site1", "site2", "dist(km)")
head(CHN.bli)

CHN.bli2 <- CHN.bli[sort(CHN.bli$site2),]

head(CHN.bli2)


##write to csv
write.csv(CHN.bli, file="distance.CHN.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHN.AdaptiveLoci <- read.table("CHN.duplicated.outliers.20171020")
head(CHN.AdaptiveLoci)

CHN.AdaptiveLoci <- gsub("\\.", ":", CHN.AdaptiveLoci$V1)
CHN.AdaptiveLoci <- gsub("X:", "", CHN.AdaptiveLoci)
head(CHN.AdaptiveLoci)
write.table(CHN.AdaptiveLoci, "CHN.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.AdaptiveLoci.names --recode --recode-INFO-all --out CHN.229.158

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.158.recode.vcf --plink --out CHN.229.158.plink

plink --file CHN.229.158.plink --noweb --recode --recodeA --out CHN.229.158.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster

```
##in R

CHN.nosex <- read.table("CHN.229.158.plink.nosex", header=F)
head(CHN.nosex)
CHN.pop <- gsub("_\\d+", "", CHN.nosex$V1)
CHN.pop
CHN.nosex$V3 <- CHN.pop
head(CHN.nosex)
write.table(CHN.nosex, "CHN.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHN.229.158.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.158

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.Adaptive.MAF <- read.table("CHN.229.158.frq.strat", header=T)
head(CHN.Adaptive.MAF)

CHN.Adaptive.MAF <- CHN.Adaptive.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.Adaptive.MAF2 <- melt(CHN.Adaptive.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.Adaptive.MAF2)
head(CHN.Adaptive.MAF2)


CHN.Adaptive.MAF3 <- dcast(CHN.Adaptive.MAF2, formula= CLST ~ SNP)
head(CHN.Adaptive.MAF3)
colnames(CHN.Adaptive.MAF3) <- paste("X", colnames(CHN.Adaptive.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.Adaptive.MAF3, file="CHN.229.158.Adaptive.MAF.csv")
```


##### 2. Neutral loci

I'm keeping the same neutral file as before, but note that I have changed the number of loci in the adaptive file. So technically this should be rerun with the new file names. 


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --exclude CHN.AdaptiveLoci.names --recode --recode-INFO-all --out CHN.229.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.Neutral.recode.vcf --plink --out CHN.229.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHN.Neutral.loci.names.map <- read.table("CHN.229.Neutral.plink.map", header=F)
CHN.Neutral.loci.names <-  CHN.Neutral.loci.names.map$V2
CHN.Neutral.loci.names <- as.data.frame(CHN.Neutral.loci.names)
CHN.1000.loci.names <- CHN.Neutral.loci.names[sample(nrow(CHN.Neutral.loci.names),1000),]
CHN.1000.loci.names <- as.data.frame(CHN.1000.loci.names)
summary(CHN.1000.loci.names)
write.table(CHN.1000.loci.names, "CHN.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHN.229.Neutral.recode.vcf --snps CHN.1000.loci.names --recode --recode-INFO-all --out CHN.1000NeutralLoci

vcftools --vcf CHN.1000NeutralLoci.recode.vcf --plink --out CHN.1000NeutralLoci.plink
plink --file CHN.1000NeutralLoci.plink --noweb --recode --recodeA --out CHN.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHN.1000NeutralLoci.plink --within CHN.PlinkCluster --noweb --freq --out CHN.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.Neutral.MAF <- read.table("CHN.1000NeutralLoci.frq.strat", header=T)
head(CHN.Neutral.MAF)

CHN.Neutral.MAF <- CHN.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.Neutral.MAF2 <- melt(CHN.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.Neutral.MAF2)
head(CHN.Neutral.MAF2)


CHN.Neutral.MAF3 <- dcast(CHN.Neutral.MAF2, formula= CLST ~ SNP)
head(CHN.Neutral.MAF3)
colnames(CHN.Neutral.MAF3) <- paste("X", colnames(CHN.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.Neutral.MAF3, file="CHN.229.Neutral.MAF.csv")
```

#### CZ
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CZ <- read.csv("AllEnv.Data_CZ_20161010.csv", header=T)
head(env.data.CZ)
#extract x and y
CZ.xy <- env.data.CZ[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CZ.dxy <- distm(CZ.xy)
CZ.dxy <- as.dist(CZ.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CZ.th <- give.thresh(CZ.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CZ.nb1 <- mst.nb(CZ.dxy)
CZ.wh1 <- which(as.matrix(CZ.dxy)==CZ.th,arr.ind=TRUE)
plot(CZ.nb1,CZ.xy,pch=20,cex=2,lty=3)
lines(CZ.xy[CZ.wh1[1,],1],CZ.xy[CZ.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CZ.th 
#[1] 122285.5
CZ.nb1
Neighbour list object:
Number of regions: 38 
Number of nonzero links: 74 
Percentage nonzero weights: 5.124654 
Average number of links: 1.947368 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CZ.listw=nb2listw(CZ.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CZ.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
CZ.scores=scores.listw(CZ.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(CZ.scores,CZ.listw,nsim=999)
          stat  pval
1  1.037362823 0.001
2  1.021372356 0.001
3  0.996965797 0.001
4  0.963918447 0.001
5  0.923493554 0.001
6  0.889534383 0.001
7  0.799794898 0.001
8  0.780295687 0.001
9  0.750000000 0.001
10 0.701384093 0.001
11 0.621955851 0.001
12 0.589551689 0.001
13 0.462419831 0.001
14 0.389652397 0.007
15 0.248021336 0.044
16 0.235486531 0.045
17 0.003281208 0.446

#16 significant MEM eigenfunctions with positive correlations. OBS. use first half (8 in our case)

write.table (CZ.scores$vectors[,1], "CZ.scores_MEM1.txt") 
write.table (CZ.scores$vectors[,2], "CZ.scores_MEM2.txt") 
write.table (CZ.scores$vectors[,3], "CZ.scores_MEM3.txt") 
write.table (CZ.scores$vectors[,4], "CZ.scores_MEM4.txt")
write.table (CZ.scores$vectors[,5], "CZ.scores_MEM5.txt")
write.table (CZ.scores$vectors[,6], "CZ.scores_MEM6.txt")
write.table (CZ.scores$vectors[,7], "CZ.scores_MEM7.txt")
write.table (CZ.scores$vectors[,8], "CZ.scores_MEM8.txt") 
```

## Great circle distance between all sampling points

```
##Geographic distance between all CZ populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017")

#rdist.earth (in fields package) wants only long & lat
CZ_lon.lat <- cbind(env.data.CZ$long, env.data.CZ$lat)
CZ_lon.lat

#calculate great circle distances
distance.matrix.CZ <- rdist.earth(CZ_lon.lat)
summary(distance.matrix.CZ)
dim(distance.matrix.CZ)

#and use only the lower half of the matrix
upper.tri(distance.matrix.CZ)
distance.matrix.CZ[lower.tri(distance.matrix.CZ)]<-NA
distance.matrix.CZ

#change from matrix to dataframe
CZ.bli <- as.data.frame(distance.matrix.CZ)
head(CZ.bli)
colnames(CZ.bli) <- env.data.CZ$site
rownames(CZ.bli) <- env.data.CZ$site

CZ.bli[lower.tri(CZ.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
CZ.bli=as.data.frame(as.table(as.matrix(CZ.bli)))  #Turn into a 3-column table
CZ.bli
CZ.bli=na.omit(CZ.bli)  #Get rid of the junk we flagged above
CZ.bli
colnames(CZ.bli)<-c("site1", "site2", "dist(km)")
head(CZ.bli)

CZ.bli2 <- CZ.bli[sort(CZ.bli$site2),]

head(CZ.bli2)


##write to csv
write.csv(CZ.bli, file="distance.CZ.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

###Change Duplicated Adaptive Loci names to vcf formate

##R

CZ.AdaptiveLoci <- read.table("CZ.duplicated.outliers.20171020")
head(CZ.AdaptiveLoci)
CZ.AdaptiveLoci <- gsub("\\.", ":", CZ.AdaptiveLoci$V1)
CZ.AdaptiveLoci <- gsub("X:", "", CZ.AdaptiveLoci)

head(CZ.AdaptiveLoci)
write.table(CZ.AdaptiveLoci, "CZ.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --snps CZ.AdaptiveLoci.names --recode --recode-INFO-all --out CZ.404.198

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.198.recode.vcf --plink --out CZ.404.198.plink

plink --file CZ.404.198.plink --noweb --recode --recodeA --out CZ.404.198.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster

```
##in R

CZ.nosex <- read.table("CZ.404.198.plink.nosex", header=F)
head(CZ.nosex)
CZ.pop <- gsub("_\\d+", "", CZ.nosex$V1)
CZ.pop
CZ.nosex$V3 <- CZ.pop
head(CZ.nosex)
write.table(CZ.nosex, "CZ.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CZ.404.198.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.198

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.Adaptive.MAF <- read.table("CZ.404.198.frq.strat", header=T)
head(CZ.Adaptive.MAF)

CZ.Adaptive.MAF <- CZ.Adaptive.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.Adaptive.MAF2 <- melt(CZ.Adaptive.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.Adaptive.MAF2)
head(CZ.Adaptive.MAF2)


CZ.Adaptive.MAF3 <- dcast(CZ.Adaptive.MAF2, formula= CLST ~ SNP)
head(CZ.Adaptive.MAF3)
colnames(CZ.Adaptive.MAF3) <- paste("X", colnames(CZ.Adaptive.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.Adaptive.MAF3, file="CZ.404.198.Adaptive.MAF.csv")
```


##### 2. Neutral loci



```

#Filter the Neutral loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --exclude CZ.AdaptiveLoci.names --recode --recode-INFO-all --out CZ.404.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.Neutral.recode.vcf --plink --out CZ.404.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CZ.Neutral.loci.names.map <- read.table("CZ.404.Neutral.plink.map", header=F)
CZ.Neutral.loci.names <-  CZ.Neutral.loci.names.map$V2
CZ.Neutral.loci.names <- as.data.frame(CZ.Neutral.loci.names)
CZ.1000.loci.names <- CZ.Neutral.loci.names[sample(nrow(CZ.Neutral.loci.names),1000),]
CZ.1000.loci.names <- as.data.frame(CZ.1000.loci.names)
summary(CZ.1000.loci.names)
write.table(CZ.1000.loci.names, "CZ.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CZ.404.Neutral.recode.vcf --snps CZ.1000.loci.names --recode --recode-INFO-all --out CZ.1000NeutralLoci

vcftools --vcf CZ.1000NeutralLoci.recode.vcf --plink --out CZ.1000NeutralLoci.plink
plink --file CZ.1000NeutralLoci.plink --noweb --recode --recodeA --out CZ.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CZ.1000NeutralLoci.plink --within CZ.PlinkCluster --noweb --freq --out CZ.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.Neutral.MAF <- read.table("CZ.1000NeutralLoci.frq.strat", header=T)
head(CZ.Neutral.MAF)

CZ.Neutral.MAF <- CZ.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.Neutral.MAF2 <- melt(CZ.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.Neutral.MAF2)
head(CZ.Neutral.MAF2)


CZ.Neutral.MAF3 <- dcast(CZ.Neutral.MAF2, formula= CLST ~ SNP)
head(CZ.Neutral.MAF3)
colnames(CZ.Neutral.MAF3) <- paste("X", colnames(CZ.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.Neutral.MAF3, file="CZ.404.Neutral.MAF.csv")
```

#### CHS
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHS <- read.csv("AllEnv.Data_CHS_20161010.csv", header=T)
env.data.CHS <- env.data.CHS[-23,]  ##remove stba
head(env.data.CHS)
#extract x and y
CHS.xy <- env.data.CHS[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CHS.dxy <- distm(CHS.xy)
CHS.dxy <- as.dist(CHS.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CHS.th <- give.thresh(CHS.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CHS.nb1 <- mst.nb(CHS.dxy)
CHS.wh1 <- which(as.matrix(CHS.dxy)==CHS.th,arr.ind=TRUE)
plot(CHS.nb1,CHS.xy,pch=20,cex=2,lty=3)
lines(CHS.xy[CHS.wh1[1,],1],CHS.xy[CHS.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CHS.th 
#[1] 47654.6
CHS.nb1
Neighbour list object:
Number of regions: 24 
Number of nonzero links: 46 
Percentage nonzero weights: 7.986111 
Average number of links: 1.916667 


#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CHS.listw=nb2listw(CHS.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CHS.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
CHS.scores=scores.listw(CHS.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(CHS.scores,CHS.listw,nsim=999)
          stat  pval
1  1.055515864 0.001
2  1.000816527 0.001
3  0.935993728 0.001
4  0.854928765 0.001
5  0.750000000 0.001
6  0.741297080 0.001
7  0.683044717 0.001
8  0.493745692 0.003
9  0.304139751 0.051
10 0.009410945 0.401

#8 significant MEM eigenfunctions with positive correlations. OBS. use first half (4 in our case)

write.table (CHS.scores$vectors[,1], "CHS.scores_MEM1.txt") 
write.table (CHS.scores$vectors[,2], "CHS.scores_MEM2.txt") 
write.table (CHS.scores$vectors[,3], "CHS.scores_MEM3.txt") 
write.table (CHS.scores$vectors[,4], "CHS.scores_MEM4.txt")
```

## Great circle distance between all sampling points

```
##Geographic distance between all CHS populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017")

#rdist.earth (in fields package) wants only long & lat
CHS_lon.lat <- cbind(env.data.CHS$long, env.data.CHS$lat)
CHS_lon.lat

#calculate great circle distances
distance.matrix.CHS <- rdist.earth(CHS_lon.lat)
summary(distance.matrix.CHS)
dim(distance.matrix.CHS)

#and use only the lower half of the matrix
upper.tri(distance.matrix.CHS)
distance.matrix.CHS[lower.tri(distance.matrix.CHS)]<-NA
distance.matrix.CHS

#change from matrix to dataframe
CHS.bli <- as.data.frame(distance.matrix.CHS)
head(CHS.bli)
colnames(CHS.bli) <- env.data.CHS$site
rownames(CHS.bli) <- env.data.CHS$site

CHS.bli[lower.tri(CHS.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
CHS.bli=as.data.frame(as.table(as.matrix(CHS.bli)))  #Turn into a 3-column table
CHS.bli
CHS.bli=na.omit(CHS.bli)  #Get rid of the junk we flagged above
CHS.bli
colnames(CHS.bli)<-c("site1", "site2", "dist(km)")
head(CHS.bli)

CHS.bli2 <- CHS.bli[sort(CHS.bli$site2),]

head(CHS.bli2)


##write to csv
write.csv(CHS.bli, file="distance.CHS.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

###Change Duplicated Adaptive Loci names to vcf formate

##R

CHS.AdaptiveLoci <- read.table("CHS.duplicated.outliers.20171020")
head(CHS.AdaptiveLoci)
CHS.AdaptiveLoci <- gsub("\\.", ":", CHS.AdaptiveLoci$V1)
CHS.AdaptiveLoci <- gsub("X:", "", CHS.AdaptiveLoci)

head(CHS.AdaptiveLoci)
write.table(CHS.AdaptiveLoci, "CHS.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --snps CHS.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.275.355

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.355.recode.vcf --plink --out CHS.275.355.plink

plink --file CHS.275.355.plink --noweb --recode --recodeA --out CHS.275.355.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.PlinkCluster

```
##in R

CHS.nosex <- read.table("CHS.275.355.plink.nosex", header=F)
head(CHS.nosex)
CHS.pop <- gsub("_\\d+", "", CHS.nosex$V1)
CHS.pop
CHS.nosex$V3 <- CHS.pop
head(CHS.nosex)
write.table(CHS.nosex, "CHS.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHS.275.355.plink --within CHS.PlinkCluster --freq --noweb --out CHS.275.355

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.Adaptive.MAF <- read.table("CHS.275.355.frq.strat", header=T)
head(CHS.Adaptive.MAF)

CHS.Adaptive.MAF <- CHS.Adaptive.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.Adaptive.MAF2 <- melt(CHS.Adaptive.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.Adaptive.MAF2)
head(CHS.Adaptive.MAF2)


CHS.Adaptive.MAF3 <- dcast(CHS.Adaptive.MAF2, formula= CLST ~ SNP)
head(CHS.Adaptive.MAF3)
colnames(CHS.Adaptive.MAF3) <- paste("X", colnames(CHS.Adaptive.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.Adaptive.MAF3, file="CHS.275.355.Adaptive.MAF.csv")
```


##### 2. Neutral loci



```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --exclude CHS.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.275.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.Neutral.recode.vcf --plink --out CHS.275.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.Neutral.loci.names.map <- read.table("CHS.275.Neutral.plink.map", header=F)
CHS.Neutral.loci.names <-  CHS.Neutral.loci.names.map$V2
CHS.Neutral.loci.names <- as.data.frame(CHS.Neutral.loci.names)
CHS.1000.loci.names <- CHS.Neutral.loci.names[sample(nrow(CHS.Neutral.loci.names),1000),]
CHS.1000.loci.names <- as.data.frame(CHS.1000.loci.names)
summary(CHS.1000.loci.names)
write.table(CHS.1000.loci.names, "CHS.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.275.Neutral.recode.vcf --snps CHS.1000.loci.names --recode --recode-INFO-all --out CHS.1000NeutralLoci

vcftools --vcf CHS.1000NeutralLoci.recode.vcf --plink --out CHS.1000NeutralLoci.plink
plink --file CHS.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.1000NeutralLoci.plink --within CHS.PlinkCluster --noweb --freq --out CHS.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.Neutral.MAF <- read.table("CHS.1000NeutralLoci.frq.strat", header=T)
head(CHS.Neutral.MAF)

CHS.Neutral.MAF <- CHS.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.Neutral.MAF2 <- melt(CHS.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.Neutral.MAF2)
head(CHS.Neutral.MAF2)


CHS.Neutral.MAF3 <- dcast(CHS.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.Neutral.MAF3)
colnames(CHS.Neutral.MAF3) <- paste("X", colnames(CHS.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.Neutral.MAF3, file="CHS.275.Neutral.MAF.csv")
```



#### CHS.TI
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHS.TI <- read.csv("AllEnv.Data_CHS.TI_20161010.csv", header=T)
head(env.data.CHS.TI)
#extract x and y
CHS.TI.xy <- env.data.CHS.TI[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CHS.TI.dxy <- distm(CHS.TI.xy)
CHS.TI.dxy <- as.dist(CHS.TI.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CHS.TI.th <- give.thresh(CHS.TI.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CHS.TI.nb1 <- mst.nb(CHS.TI.dxy)
CHS.TI.wh1 <- which(as.matrix(CHS.TI.dxy)==CHS.TI.th,arr.ind=TRUE)
plot(CHS.TI.nb1,CHS.TI.xy,pch=20,cex=2,lty=3)
lines(CHS.TI.xy[CHS.TI.wh1[1,],1],CHS.TI.xy[CHS.TI.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CHS.TI.th 
#[1] 19167.21
CHS.TI.nb1
Neighbour list object:
Number of regions: 15 
Number of nonzero links: 28 
Percentage nonzero weights: 12.44444 
Average number of links: 1.866667 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CHS.TI.listw=nb2listw(CHS.TI.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CHS.TI.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
CHS.TI.scores=scores.listw(CHS.TI.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(CHS.TI.scores,CHS.TI.listw,nsim=999)
         stat  pval
1 1.035461987 0.001
2 0.925047840 0.001
3 0.750000000 0.001
4 0.735948033 0.002
5 0.426558214 0.029
6 0.008914643 0.370

#5 significant MEM eigenfunctions with positive correlations. OBS. use first half (3 in our case)

write.table (CHS.TI.scores$vectors[,1], "CHS.TI.scores_MEM1.txt") 
write.table (CHS.TI.scores$vectors[,2], "CHS.TI.scores_MEM2.txt") 
write.table (CHS.TI.scores$vectors[,3], "CHS.TI.scores_MEM3.txt") 
```

## Great circle distance between all sampling points

```
##Geographic distance between all CHS.TI populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017")

#rdist.earth (in fields package) wants only long & lat
CHS.TI_lon.lat <- cbind(env.data.CHS.TI$long, env.data.CHS.TI$lat)
CHS.TI_lon.lat

#calculate great circle distances
distance.matrix.CHS.TI <- rdist.earth(CHS.TI_lon.lat)
summary(distance.matrix.CHS.TI)
dim(distance.matrix.CHS.TI)

#and use only the lower half of the matrix
upper.tri(distance.matrix.CHS.TI)
distance.matrix.CHS.TI[lower.tri(distance.matrix.CHS.TI)]<-NA
distance.matrix.CHS.TI

#change from matrix to dataframe
CHS.TI.bli <- as.data.frame(distance.matrix.CHS.TI)
head(CHS.TI.bli)
colnames(CHS.TI.bli) <- env.data.CHS.TI$site
rownames(CHS.TI.bli) <- env.data.CHS.TI$site

CHS.TI.bli[lower.tri(CHS.TI.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
CHS.TI.bli=as.data.frame(as.table(as.matrix(CHS.TI.bli)))  #Turn into a 3-column table
CHS.TI.bli
CHS.TI.bli=na.omit(CHS.TI.bli)  #Get rid of the junk we flagged above
CHS.TI.bli
colnames(CHS.TI.bli)<-c("site1", "site2", "dist(km)")
head(CHS.TI.bli)

CHS.TI.bli2 <- CHS.TI.bli[sort(CHS.TI.bli$site2),]

head(CHS.TI.bli2)


##write to csv
write.csv(CHS.TI.bli, file="distance.CHS.TI.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

###Change Duplicated Adaptive Loci names to vcf formate

##R

CHS.TI.AdaptiveLoci <- read.table("CHS.TI.duplicated.outliers.20171020")
head(CHS.TI.AdaptiveLoci)
CHS.TI.AdaptiveLoci <- gsub("\\.", ":", CHS.TI.AdaptiveLoci$V1)
CHS.TI.AdaptiveLoci <- gsub("X.", "", CHS.TI.AdaptiveLoci)
head(CHS.TI.AdaptiveLoci)
write.table(CHS.TI.AdaptiveLoci, "CHS.TI.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.TI.140.260

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.260.recode.vcf --plink --out CHS.TI.140.260.plink

plink --file CHS.TI.140.260.plink --noweb --recode --recodeA --out CHS.TI.140.260.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster

```
##in R

CHS.TI.nosex <- read.table("CHS.TI.140.260.plink.nosex", header=F)
head(CHS.TI.nosex)
CHS.TI.pop <- gsub("_\\d+", "", CHS.TI.nosex$V1)
CHS.TI.pop
CHS.TI.nosex$V3 <- CHS.TI.pop
head(CHS.TI.nosex)
write.table(CHS.TI.nosex, "CHS.TI.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHS.TI.140.260.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.260

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.Adaptive.MAF <- read.table("CHS.TI.140.260.frq.strat", header=T)
head(CHS.TI.Adaptive.MAF)

CHS.TI.Adaptive.MAF <- CHS.TI.Adaptive.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.Adaptive.MAF2 <- melt(CHS.TI.Adaptive.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.Adaptive.MAF2)
head(CHS.TI.Adaptive.MAF2)


CHS.TI.Adaptive.MAF3 <- dcast(CHS.TI.Adaptive.MAF2, formula= CLST ~ SNP)
head(CHS.TI.Adaptive.MAF3)
colnames(CHS.TI.Adaptive.MAF3) <- paste("X", colnames(CHS.TI.Adaptive.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.Adaptive.MAF3, file="CHS.TI.140.260.Adaptive.MAF.csv")
```


##### 2. Neutral loci



```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --exclude CHS.TI.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.TI.140.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.Neutral.recode.vcf --plink --out CHS.TI.140.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.TI.Neutral.loci.names.map <- read.table("CHS.TI.140.Neutral.plink.map", header=F)
CHS.TI.Neutral.loci.names <-  CHS.TI.Neutral.loci.names.map$V2
CHS.TI.Neutral.loci.names <- as.data.frame(CHS.TI.Neutral.loci.names)
CHS.TI.1000.loci.names <- CHS.TI.Neutral.loci.names[sample(nrow(CHS.TI.Neutral.loci.names),1000),]
CHS.TI.1000.loci.names <- as.data.frame(CHS.TI.1000.loci.names)
summary(CHS.TI.1000.loci.names)
write.table(CHS.TI.1000.loci.names, "CHS.TI.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.TI.140.Neutral.recode.vcf --snps CHS.TI.1000.loci.names --recode --recode-INFO-all --out CHS.TI.1000NeutralLoci

vcftools --vcf CHS.TI.1000NeutralLoci.recode.vcf --plink --out CHS.TI.1000NeutralLoci.plink
plink --file CHS.TI.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.TI.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.TI.1000NeutralLoci.plink --within CHS.TI.PlinkCluster --noweb --freq --out CHS.TI.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.Neutral.MAF <- read.table("CHS.TI.1000NeutralLoci.frq.strat", header=T)
head(CHS.TI.Neutral.MAF)

CHS.TI.Neutral.MAF <- CHS.TI.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.Neutral.MAF2 <- melt(CHS.TI.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.Neutral.MAF2)
head(CHS.TI.Neutral.MAF2)


CHS.TI.Neutral.MAF3 <- dcast(CHS.TI.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.TI.Neutral.MAF3)
colnames(CHS.TI.Neutral.MAF3) <- paste("X", colnames(CHS.TI.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.Neutral.MAF3, file="CHS.TI.140.Neutral.MAF.csv")
```


#### CHS.VS
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHS.VS <- read.csv("AllEnv.Data_CHS.VS_20161010.csv", header=T)
head(env.data.CHS.VS)
#extract x and y
CHS.VS.xy <- env.data.CHS.VS[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CHS.VS.dxy <- distm(CHS.VS.xy)
CHS.VS.dxy <- as.dist(CHS.VS.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CHS.VS.th <- give.thresh(CHS.VS.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CHS.VS.nb1 <- mst.nb(CHS.VS.dxy)
CHS.VS.wh1 <- which(as.matrix(CHS.VS.dxy)==CHS.VS.th,arr.ind=TRUE)
plot(CHS.VS.nb1,CHS.VS.xy,pch=20,cex=2,lty=3)
lines(CHS.VS.xy[CHS.VS.wh1[1,],1],CHS.VS.xy[CHS.VS.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CHS.VS.th 
#[1] 41004.58
CHS.VS.nb1
Neighbour list object:
Number of regions: 10 
Number of nonzero links: 18 
Percentage nonzero weights: 18 
Average number of links: 1.8  

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CHS.VS.listw=nb2listw(CHS.VS.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CHS.VS.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
CHS.VS.scores=scores.listw(CHS.VS.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(CHS.VS.scores,CHS.VS.listw,nsim=999)
         stat  pval
1 0.998929647 0.001
2 0.750000000 0.001
3 0.523833554 0.043
4 0.008291494 0.350

#3 significant MEM eigenfunctions with positive correlations. OBS. use first half (2 in our case)

write.table (CHS.VS.scores$vectors[,1], "CHS.VS.scores_MEM1.txt") 
write.table (CHS.VS.scores$vectors[,2], "CHS.VS.scores_MEM2.txt") 
```

## Great circle distance between all sampling points

```
##Geographic distance between all CHS.VS populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017")

#rdist.earth (in fields package) wants only long & lat
CHS.VS_lon.lat <- cbind(env.data.CHS.VS$long, env.data.CHS.VS$lat)
CHS.VS_lon.lat

#calculate great circle distances
distance.matrix.CHS.VS <- rdist.earth(CHS.VS_lon.lat)
summary(distance.matrix.CHS.VS)
dim(distance.matrix.CHS.VS)

#and use only the lower half of the matrix
upper.tri(distance.matrix.CHS.VS)
distance.matrix.CHS.VS[lower.tri(distance.matrix.CHS.VS)]<-NA
distance.matrix.CHS.VS

#change from matrix to dataframe
CHS.VS.bli <- as.data.frame(distance.matrix.CHS.VS)
head(CHS.VS.bli)
colnames(CHS.VS.bli) <- env.data.CHS.VS$site
rownames(CHS.VS.bli) <- env.data.CHS.VS$site

CHS.VS.bli[lower.tri(CHS.VS.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
CHS.VS.bli=as.data.frame(as.table(as.matrix(CHS.VS.bli)))  #Turn into a 3-column table
CHS.VS.bli
CHS.VS.bli=na.omit(CHS.VS.bli)  #Get rid of the junk we flagged above
CHS.VS.bli
colnames(CHS.VS.bli)<-c("site1", "site2", "dist(km)")
head(CHS.VS.bli)

CHS.VS.bli2 <- CHS.VS.bli[sort(CHS.VS.bli$site2),]

head(CHS.VS.bli2)


##write to csv
write.csv(CHS.VS.bli, file="distance.CHS.VS.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

###Change Duplicated Adaptive Loci names to vcf formate

##R

CHS.VS.AdaptiveLoci <- read.table("CHS.VS.duplicated.outliers.20171020")
head(CHS.VS.AdaptiveLoci)
CHS.VS.AdaptiveLoci <- gsub("\\.", ":", CHS.VS.AdaptiveLoci$V1)
CHS.VS.AdaptiveLoci <- gsub("X:", "", CHS.VS.AdaptiveLoci)
head(CHS.VS.AdaptiveLoci)
write.table(CHS.VS.AdaptiveLoci, "CHS.VS.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.VS.135.408

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.408.recode.vcf --plink --out CHS.VS.135.408.plink

plink --file CHS.VS.135.408.plink --noweb --recode --recodeA --out CHS.VS.135.408.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster

```
##in R

CHS.VS.nosex <- read.table("CHS.VS.135.408.plink.nosex", header=F)
head(CHS.VS.nosex)
CHS.VS.pop <- gsub("_\\d+", "", CHS.VS.nosex$V1)
CHS.VS.pop
CHS.VS.nosex$V3 <- CHS.VS.pop
head(CHS.VS.nosex)
write.table(CHS.VS.nosex, "CHS.VS.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHS.VS.135.408.plink --within CHS.VS.PlinkCluster --freq --noweb --out CHS.VS.135.408

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.Adaptive.MAF <- read.table("CHS.VS.135.408.frq.strat", header=T)
head(CHS.VS.Adaptive.MAF)

CHS.VS.Adaptive.MAF <- CHS.VS.Adaptive.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.Adaptive.MAF2 <- melt(CHS.VS.Adaptive.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.Adaptive.MAF2)
head(CHS.VS.Adaptive.MAF2)


CHS.VS.Adaptive.MAF3 <- dcast(CHS.VS.Adaptive.MAF2, formula= CLST ~ SNP)
head(CHS.VS.Adaptive.MAF3)
colnames(CHS.VS.Adaptive.MAF3) <- paste("X", colnames(CHS.VS.Adaptive.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.Adaptive.MAF3, file="CHS.VS.135.408.Adaptive.MAF.csv")
```


##### 2. Neutral loci



```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --exclude CHS.VS.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.Neutral.recode.vcf --plink --out CHS.VS.135.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.VS.Neutral.loci.names.map <- read.table("CHS.VS.135.Neutral.plink.map", header=F)
CHS.VS.Neutral.loci.names <-  CHS.VS.Neutral.loci.names.map$V2
CHS.VS.Neutral.loci.names <- as.data.frame(CHS.VS.Neutral.loci.names)
CHS.VS.1000.loci.names <- CHS.VS.Neutral.loci.names[sample(nrow(CHS.VS.Neutral.loci.names),1000),]
CHS.VS.1000.loci.names <- as.data.frame(CHS.VS.1000.loci.names)
summary(CHS.VS.1000.loci.names)
write.table(CHS.VS.1000.loci.names, "CHS.VS.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.VS.135.Neutral.recode.vcf --snps CHS.VS.1000.loci.names --recode --recode-INFO-all --out CHS.VS.1000NeutralLoci

vcftools --vcf CHS.VS.1000NeutralLoci.recode.vcf --plink --out CHS.VS.1000NeutralLoci.plink
plink --file CHS.VS.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.VS.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.VS.1000NeutralLoci.plink --within CHS.VS.PlinkCluster --noweb --freq --out CHS.VS.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.Neutral.MAF <- read.table("CHS.VS.1000NeutralLoci.frq.strat", header=T)
head(CHS.VS.Neutral.MAF)

CHS.VS.Neutral.MAF <- CHS.VS.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.Neutral.MAF2 <- melt(CHS.VS.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.Neutral.MAF2)
head(CHS.VS.Neutral.MAF2)


CHS.VS.Neutral.MAF3 <- dcast(CHS.VS.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.VS.Neutral.MAF3)
colnames(CHS.VS.Neutral.MAF3) <- paste("X", colnames(CHS.VS.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.Neutral.MAF3, file="CHS.VS.135.Neutral.MAF.csv")
```




#### ENV data

Variables used: 

Same as for the EAA: 

1. sol.rad.60d

2. temp.laying.date

3. pcpt.60d

4. shadow.days

5. day10cm


Create a .csv input file with

1. pop names

2. distance

Use the last population (trpa) as distance 0. So all distances are measured from here.

##### Input Files

CHN.GF.AdaptiveLoci.Input.csv

CHN.GF.NeutralLoci.Input.csv

CZ.GF.AdaptiveLoci.Input.csv

CZ.GF.NeutralLoci.Input.csv

CHS.TI.GF.AdaptiveLoci.Input.csv

CHS.TI.GF.NeutralLoci.Input.csv

CHS.VS.GF.AdaptiveLoci.Input.csv

CHS.VS.GF.NeutralLoci.Input.csv


### Run GF

1. Adaptive Loci
```
library(gradientForest)

###CHN

gf.CHN.Adaptive <- read.csv("CHN.GF.AdaptiveLoci.Input.csv", header=T)
envGF.CHN.Adaptive <- gf.CHN.Adaptive[,2:11]
colnames(envGF.CHN.Adaptive)

CHN.AdaptiveSNPS <- gf.CHN.Adaptive[,grep("X.", colnames(gf.CHN.Adaptive))]
CHN.AdaptiveSNPS <- CHN.AdaptiveSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHN.Adaptive)/2)
maxLevel

gf.CHN.Adaptive.SNPs <- gradientForest(cbind(envGF.CHN.Adaptive, CHN.AdaptiveSNPS), predictor.vars=colnames(envGF.CHN.Adaptive), response.vars=colnames(CHN.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##4 warnings (loci variable in less than 5 populations)


###CZ

gf.CZ.Adaptive <- read.csv("CZ.GF.AdaptiveLoci.Input.csv", header=T)
colnames(gf.CZ.Adaptive)
envGF.CZ.Adaptive <- gf.CZ.Adaptive[,2:15]
colnames(envGF.CZ.Adaptive)

CZ.AdaptiveSNPS <- gf.CZ.Adaptive[,grep("X.", colnames(gf.CZ.Adaptive))]
CZ.AdaptiveSNPS <- CZ.AdaptiveSNPS[,-1]
colnames(CZ.AdaptiveSNPS)
maxLevel <- log2(0.368*nrow(envGF.CZ.Adaptive)/2)
maxLevel

gf.CZ.Adaptive.SNPs <- gradientForest(cbind(envGF.CZ.Adaptive, CZ.AdaptiveSNPS), predictor.vars=colnames(envGF.CZ.Adaptive), response.vars=colnames(CZ.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##0 warnings (loci variable in less than 5 populations)

###CHS.TI

gf.CHS.TI.Adaptive <- read.csv("CHS.TI.GF.AdaptiveLoci.Input.csv", header=T)
colnames(gf.CHS.TI.Adaptive)
envGF.CHS.TI.Adaptive <- gf.CHS.TI.Adaptive[,2:10]
colnames(envGF.CHS.TI.Adaptive)

CHS.TI.AdaptiveSNPS <- gf.CHS.TI.Adaptive[,grep("X.", colnames(gf.CHS.TI.Adaptive))]
CHS.TI.AdaptiveSNPS <- CHS.TI.AdaptiveSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.Adaptive)/2)
maxLevel

gf.CHS.TI.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.TI.Adaptive, CHS.TI.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.TI.Adaptive), response.vars=colnames(CHS.TI.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##29 warnings - loci variable in less than 5 pops

###CHS.VS

gf.CHS.VS.Adaptive <- read.csv("CHS.VS.GF.AdaptiveLoci.Input.csv", header=T)
colnames(gf.CHS.VS.Adaptive)
envGF.CHS.VS.Adaptive <- gf.CHS.VS.Adaptive[,2:9]
colnames(envGF.CHS.VS.Adaptive)

CHS.VS.AdaptiveSNPS <- gf.CHS.VS.Adaptive[,grep("X.", colnames(gf.CHS.VS.Adaptive))]
CHS.VS.AdaptiveSNPS <- CHS.VS.AdaptiveSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Adaptive)/2)
maxLevel

gf.CHS.VS.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.VS.Adaptive, CHS.VS.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.VS.Adaptive), response.vars=colnames(CHS.VS.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

###more than 50 warnings about non-polymorphic loci


###CHS

gf.CHS.Adaptive <- read.csv("CHS.GF.AdaptiveLoci.Input.csv", header=T)
colnames(gf.CHS.Adaptive)
envGF.CHS.Adaptive <- gf.CHS.Adaptive[,2:11]
colnames(envGF.CHS.Adaptive)

CHS.AdaptiveSNPS <- gf.CHS.Adaptive[,grep("X.", colnames(gf.CHS.Adaptive))]
CHS.AdaptiveSNPS <- CHS.AdaptiveSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.Adaptive)/2)
maxLevel

gf.CHS.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.Adaptive, CHS.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.Adaptive), response.vars=colnames(CHS.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##9 warnings

```

2. Neutral Loci

```
###CHN

gf.CHN.Neutral <- read.csv("CHN.GF.NeutralLoci.Input.csv", header=T)
envGF.CHN.Neutral <- gf.CHN.Neutral[,2:11]
colnames(envGF.CHN.Neutral)

CHN.Neutral.SNPs <- gf.CHN.Neutral[,12:ncol(gf.CHN.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CHN.Neutral)/2)
maxLevel

gf.CHN.Neutral.SNPs <- gradientForest(cbind(envGF.CHN.Neutral, CHN.Neutral.SNPs), predictor.vars=colnames(envGF.CHN.Neutral), response.vars=colnames(CHN.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##3 warnings

###CZ

gf.CZ.Neutral <- read.csv("CZ.GF.NeutralLoci.Input.csv", header=T)
envGF.CZ.Neutral <- gf.CZ.Neutral[,2:15]
colnames(envGF.CZ.Neutral)

CZ.Neutral.SNPs <- gf.CZ.Neutral[,16:ncol(gf.CZ.Neutral)]
colnames(CZ.Neutral.SNPs)
maxLevel <- log2(0.368*nrow(envGF.CZ.Neutral)/2)
maxLevel

gf.CZ.Neutral.SNPs <- gradientForest(cbind(envGF.CZ.Neutral, CZ.Neutral.SNPs), predictor.vars=colnames(envGF.CZ.Neutral), response.vars=colnames(CZ.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

gf.CHS.TI.Neutral <- read.csv("CHS.TI.GF.NeutralLoci.Input.csv", header=T)
envGF.CHS.TI.Neutral <- gf.CHS.TI.Neutral[,2:10]
colnames(envGF.CHS.TI.Neutral)

CHS.TI.Neutral.SNPs <- gf.CHS.TI.Neutral[,11:ncol(gf.CHS.TI.Neutral)]
colnames(CHS.TI.Neutral.SNPs)
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.Neutral)/2)
maxLevel

gf.CHS.TI.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.TI.Neutral, CHS.TI.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.TI.Neutral), response.vars=colnames(CHS.TI.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##25 warnings


###CHS.VS

gf.CHS.VS.Neutral <- read.csv("CHS.VS.GF.NeutralLoci.Input.csv", header=T)
envGF.CHS.VS.Neutral <- gf.CHS.VS.Neutral[,2:9]
colnames(envGF.CHS.VS.Neutral)

CHS.VS.Neutral.SNPs <- gf.CHS.VS.Neutral[,10:ncol(gf.CHS.VS.Neutral)]
colnames(CHS.VS.Neutral.SNPs)
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Neutral)/2)
maxLevel

gf.CHS.VS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.VS.Neutral, CHS.VS.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.VS.Neutral), response.vars=colnames(CHS.VS.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##50 warnings

###CHS

gf.CHS.Neutral <- read.csv("CHS.GF.NeutralLoci.Input.csv", header=T)
envGF.CHS.Neutral <- gf.CHS.Neutral[,2:11]
colnames(envGF.CHS.Neutral)

CHS.Neutral.SNPs <- gf.CHS.Neutral[,12:ncol(gf.CHS.Neutral)]
colnames(CHS.Neutral.SNPs)
maxLevel <- log2(0.368*nrow(envGF.CHS.Neutral)/2)
maxLevel

gf.CHS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.Neutral, CHS.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.Neutral), response.vars=colnames(CHS.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##3 warnings

```



## Presenting Results

As in the Fitzpatric paper.

1. Table1 = summary of the datasets

2. Fig1 = R2 heatmap of reletive importance of different predictor variables

3. Graph of the average allele turnover for each of the three datasets (on one Figure)

4. Graph of the SNPs for BIO5 for ENV & Fst datasets

5. PCA of env variables with sample sites on (PCA centered but not scaled)

6. Geographic map of PC space: colour according to PC space

7. Map difference between datasets (Neutral vs Fst, Neutral vs Env)

8. Pop-level vulnerability to climate change: Map mean genetic offset between current and future climates. 


I had a lot of problems creating the geographic maps. I got scripts from Karina (from Victoria Sork's lab) which created "greengrid" raster files with data extracted from 2-10k points. However, the resolution for this was always really low. It turns out I can create the maps using the scripts from Fitzpatrick & Keller 2015 by just using the raster layers directly (i.e. all the points rather than a subset of random points). 


### Table 1

The R2 can be obtained from the importance of each "Species" in the gf model:  see getAnywhere(performance.plot)
```
###CHS.VS

perf.NEUTRAL.CHS.VS <- importance(gf.CHS.VS.Neutral.SNPs, type="Species")
perf.Fst.CHS.VS <- importance(gf.CHS.VS.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL.CHS.VS)
summary(perf.Fst.CHS.VS)

perf.NEUT.CHS.VS.df <- as.data.frame(perf.NEUTRAL.CHS.VS)
perf.Fst.CHS.VS.df <- as.data.frame(perf.Fst.CHS.VS)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.CHS.VS.df[which(perf.NEUT.CHS.VS.df$perf.NEUTRAL>0.5),])
length(perf.Fst.CHS.VS.df[which(perf.Fst.CHS.VS.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHS.VS.Neutral.SNPs
gf.CHS.VS.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.CHS.VS.df$perf.NEUT.CHS.VS)
hist(perf.Fst.CHS.VS.df$perf.Fst.CHS.VS)



###CHS.TI

perf.NEUTRAL.CHS.TI <- importance(gf.CHS.TI.Neutral.SNPs, type="Species")
perf.Fst.CHS.TI <- importance(gf.CHS.TI.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL.CHS.TI)
summary(perf.Fst.CHS.TI)

perf.NEUT.CHS.TI.df <- as.data.frame(perf.NEUTRAL.CHS.TI)
perf.Fst.CHS.TI.df <- as.data.frame(perf.Fst.CHS.TI)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.CHS.TI.df[which(perf.NEUT.CHS.TI.df$perf.NEUTRAL>0.5),])
length(perf.Fst.CHS.TI.df[which(perf.Fst.CHS.TI.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHS.TI.Neutral.SNPs
gf.CHS.TI.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)

###CHS

perf.NEUTRAL.CHS <- importance(gf.CHS.Neutral.SNPs, type="Species")
perf.Fst.CHS <- importance(gf.CHS.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL.CHS)
summary(perf.Fst.CHS)

perf.NEUT.CHS.df <- as.data.frame(perf.NEUTRAL.CHS)
perf.Fst.CHS.df <- as.data.frame(perf.Fst.CHS)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.CHS.df[which(perf.NEUT.CHS.df$perf.NEUTRAL>0.5),])
length(perf.Fst.CHS.df[which(perf.Fst.CHS.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHS.Neutral.SNPs
gf.CHS.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)


###CZ

perf.NEUTRAL.CZ <- importance(gf.CZ.Neutral.SNPs, type="Species")
perf.Fst.CZ <- importance(gf.CZ.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL.CZ)
summary(perf.Fst.CZ)

perf.NEUT.CZ.df <- as.data.frame(perf.NEUTRAL.CZ)
perf.Fst.CZ.df <- as.data.frame(perf.Fst.CZ)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.CZ.df[which(perf.NEUT.CZ.df$perf.NEUTRAL>0.5),])
length(perf.Fst.CZ.df[which(perf.Fst.CZ.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CZ.Neutral.SNPs
gf.CZ.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)


###CHN

perf.NEUTRAL.CHN <- importance(gf.CHN.Neutral.SNPs, type="Species")
perf.Fst.CHN <- importance(gf.CHN.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL.CHN)
summary(perf.Fst.CHN)

perf.NEUT.CHN.df <- as.data.frame(perf.NEUTRAL.CHN)
perf.Fst.CHN.df <- as.data.frame(perf.Fst.CHN)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.CHN.df[which(perf.NEUT.CHN.df$perf.NEUTRAL>0.5),])
length(perf.Fst.CHN.df[which(perf.Fst.CHN.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHN.Neutral.SNPs
gf.CHN.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)


```




### Fig1: R2 heatmap
```
#R

###CHN

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHN.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHN.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHN.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHN") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()

###CHS

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHS.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHS.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()


###CZ
#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CZ.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CZ.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CZ.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CZ") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###CHS.TI

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHS.TI.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHS.TI.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.TI.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.TI") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()


###CHS.VS

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHS.VS.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHS.VS.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.VS.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.VS") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()
```

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/31220345-07976010-a9b8-11e7-83de-bf481fa4cc4e.png



### Fig2: Average allelic turnover across all three datasets. 

It took me a while to find the code for the gf plots so that I could figure out what was being plotted, and so combine the plots. 
```
leg.txt <- c("Neutral", "Adaptive")  ##define the text that will be added to the legend
```






##### CHS.TI Cumulative plot
```
##sol.rad.60d

CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "sol.rad.60d") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "sol.rad.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

ymax=0.0015 #set the ymax so that it's the same for the the plots to be overlaid

pdf("sol.rad.60d.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="sol.rad.60d: mean solar radiation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##temp.laying.date
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "temp.laying.date") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "temp.laying.date")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("temp.laying.date.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="temp.laying.date", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##pcpt.60d
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "pcpt.60d") 
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "pcpt.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("pcpt.60d.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="pcpt.60d: Precipitation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##shadow.days
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "shadow.days") ##find the cumulative importance for each gf.model output. 
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "shadow.days")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("shadow.days.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="shadow.days", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##day10cm
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "day10cm") ##find the cumulative importance for each gf.model output. Correcting label (see above)
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "day10cm")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("day10cm.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="day10cm: calendar day with 10+cm snow", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##MEM1
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "MEM1") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "MEM1")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("MEM1.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM1", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM2
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "MEM2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "MEM2")


isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("MEM2.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="MEM2", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "dist")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("dist.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="Distance (km)", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##Create the combined plot in AdobeIllustrator
```



##### CHS.VS Cumulative plot
```
##sol.rad.60d

CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "sol.rad.60d") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "sol.rad.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

ymax=0.0015 #set the ymax so that it's the same for the the plots to be overlaid

pdf("sol.rad.60d.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="sol.rad.60d: mean solar radiation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##temp.laying.date
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "temp.laying.date") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "temp.laying.date")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("temp.laying.date.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="temp.laying.date", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##pcpt.60d
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "day10cm") ##I mislabeled these in the original excel sheet (now corrected). Only for CHS.VS. The rest of the transects are correct
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "day10cm")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("pcpt.60d.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="pcpt.60d: Precipitation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##shadow.days
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "pcpt.60d") ##find the cumulative importance for each gf.model output. Correcting the excel mislabelling here: see above
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "pcpt.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("shadow.days.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="shadow.days", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##day10cm
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "shadow.days") ##find the cumulative importance for each gf.model output. Correcting label (see above)
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "shadow.days")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("day10cm.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="day10cm: calendar day with 10+cm snow", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##MEM1
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "MEM1") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "MEM1")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("MEM1.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM1", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM2
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "MEM2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "MEM2")


isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("MEM2.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="MEM2", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "dist")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("dist.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="Distance (km)", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##Create the combined plot in AdobeIllustrator
```

![alt_txt][Fig2]

[Fig2]:https://cloud.githubusercontent.com/assets/12142475/22886423/13fc2d1c-f1fe-11e6-9f5b-2fbad3e0ca54.png





# 4 Individual SNPs

###First I need to redefine the species.cumulative.plot function that comes with gradient forest:
I've changed the y-limit so everythings on the same axis. And I've definede the imp.var.names. 
I've also change show.overall to FALSE so that I can see the individual SNPs

getAnywhere(species.cumulative.plot)

```
species.cumulative.plot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", 
    leg.nspecies = 10, legend = TRUE, mfrow = rev(n2mfrow(length(imp.vars) * 
        (show.species + show.overall))), show.species = TRUE, 
    show.overall = FALSE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 
        0.75, 0.1, 0.1), common.scale = F, line.ylab = 1, cex.legend = 0.75, 
    ...) 
{
    if (is.null(imp.vars)) 
        imp.vars <- imp.var.names <- c("sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm")
    par(mfrow = mfrow)
    cols <- rainbow(length(levels(obj$res.u$spec)))
    names(cols) <- levels(obj$res.u$spec)
    xaxt <- if (show.overall) 
        "n"
    else "s"
    if (show.species) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX, "Species")
            xlim <- range(sapply(CU, "[[", "x"))
            ylim <- c(0,0.02)   ###I'm editing this so that I can plot on one plot
            plot(xlim, ylim, type = "n", xlab = if (show.overall) 
                ""
            else imp.vars.names[imp.vars == varX], ylab = "", 
                xaxt = xaxt, ...)
            for (species in names(CU)) {
                isub <- seq(1, length(CU[[species]]$x), len = pmin(500, 
                  length(CU[[species]]$x)))
                lines(CU[[species]]$x[isub], CU[[species]]$y[isub], 
                  type = "s", col = cols[species])
            }
            no.species <- length(names(cols))
            imp.sp <- sapply(CU, function(cu) max(cu$y))
            best <- order(-imp.sp)[1:min(leg.nspecies, length(imp.sp))]
            if (legend) 
                legend(x = leg.posn, legend = names(cols)[best], 
                  pch = rep(1, no.species)[best], col = cols[best], 
                  bty = "n", cex = cex.legend, pt.lwd = 2)
        }
    }
    if (show.overall) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX)
            ymax <- 1.2 ##max(CU$y)
            if (varX == imp.vars[1]) 
                ymax1 <- 1.2  ##change this from ymax 
            isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
            plot(CU$x[isub], CU$y[isub], type = "s", ylab = "", 
                xlab = imp.vars.names[imp.vars == varX], ylim = c(0, 
                  if (common.scale) ymax1 else ymax), ...)
        }
    }
    mtext("Cumulative importance", side = 2, line = line.ylab, 
        outer = TRUE)
}


```



Plot of individual SNPs for a specific dataset and environmental variable
```
### CHN
pdf("CHN.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.Adaptive.SNPs)
dev.off()

pdf("CHN.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.Neutral.SNPs)
dev.off()

### CZ
pdf("CZ.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.Adaptive.SNPs)
dev.off()

pdf("CZ.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.Neutral.SNPs)
dev.off()


### CHS
pdf("CHS.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.Adaptive.SNPs)
dev.off()

pdf("CHS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.Neutral.SNPs)
dev.off()


### CHS.VS
pdf("CHS.VS.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.Adaptive.SNPs)
dev.off()

pdf("CHS.VS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.Neutral.SNPs)
dev.off()

### CHS.TI
pdf("CHS.TI.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.Adaptive.SNPs)
dev.off()

pdf("CHS.TI.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.Neutral.SNPs)
dev.off()

```
![alt_txt][Fig3]

[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png



### Best fit SNPs

Find the SNPs from the Adaptive Loci dataset that contribute the most to the Cumulative importance for each Environmental Variable. 

I'll then plot the CumImp of these loci vs the Neutral dataset for each of the variables. 




```
CU.Adapt.CHN.sol.rad.60d <- cumimp(gf.CHN.Adaptive.SNPs, "sol.rad.60d", "Species")

head(CU.Adapt.CHN.sol.rad.60d)
imp.sp.Adapt.CHN.sol.rad.60d <- sapply(CU.Adapt.CHN.sol.rad.60d, function(cu) max(cu$y))
imp.sp.Adapt.CHN.sol.rad.60d
imp.sp.Adapt.CHN.sol.rad.60d.df <- as.data.frame(imp.sp.Adapt.CHN.sol.rad.60d)
imp.sp.Adapt.CHN.sol.rad.60d.df$loci <- names(imp.sp.Adapt.CHN.sol.rad.60d)
imp.sp.Adapt.CHN.sol.rad.60d.df$isub <- 1:nrow(imp.sp.Adapt.CHN.sol.rad.60d.df)  ##index the loci
imp.sp.Adapt.CHN.sol.rad.60d.df <- imp.sp.Adapt.CHN.sol.rad.60d.df[order(imp.sp.Adapt.CHN.sol.rad.60d.df$imp.sp.Adapt.CHN.sol.rad.60d, decreasing=T),]  ##order

library(data.table)
imp.sp.Adapt.CHN.sol.rad.60d.df$cumsum <- cumsum(imp.sp.Adapt.CHN.sol.rad.60d.df$imp.sp.Adapt.CHN.sol.rad.60d)  ##add a column that is the cumulative sum of all the loci

hist(imp.sp.Adapt.CHN.sol.rad.60d.df$cumsum, breaks=69)  ##figure of the cumulative importance over all the loci

##Select only the loci that contribute to the cumulative importance. 

best.Adapt.CHN.sol.rad.60d <- na.omit(subset(imp.sp.Adapt.CHN.sol.rad.60d.df, imp.sp.Adapt.CHN.sol.rad.60d.df>0))
length(best.Adapt.CHN.sol.rad.60d$loci)  ##get the number of loci for Table

species.Adapt.CHN.sol.rad.60d.names <- best.Adapt.CHN.sol.rad.60d$loci  ## get the names
CU.Adapt.CHN.sol.rad.60d.best <- CU.Adapt.CHN.sol.rad.60d[species.Adapt.CHN.sol.rad.60d.names]

###Now I have to rerun GF using only those loci for this env variable

colnames(envGF.CHN.Adaptive)

CHN.AdaptiveSNPS.sol.rad.best <- CHN.AdaptiveSNPS[isub.CHN.sol.rad.60d]  ##select only those snps

maxLevel <- log2(0.368*nrow(envGF.CHN.Adaptive)/2)
maxLevel

gf.CHN.Adaptive.SNPs.sol.rad.60d <- gradientForest(cbind(envGF.CHN.Adaptive, CHN.AdaptiveSNPS.sol.rad.best), predictor.vars=colnames(envGF.CHN.Adaptive), response.vars=colnames(CHN.AdaptiveSNPS.sol.rad.best), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

CU.CHN.sol.rad.60d <- cumimp(gf.CHN.Adaptive.SNPs.sol.rad.60d, "sol.rad.60d")
isub.CHN.sol.rad.60d <- seq(1, length(CU.CHN.sol.rad.60d$x), len=pmin(500, length(CU.CHN.sol.rad.60d$x)))
plot(CU.CHN.sol.rad.60d$x[isub.CHN.sol.rad.60d], CU.CHN.sol.rad.60d$y[isub.CHN.sol.rad.60d], type="s", ylab="Cumulative Plot", xlab="sol.rad.60d", ylim=c(0,0.12), lty=1, lwd=1.5)
```





###define the new CU without species info, so that we can draw a cumulative plot
CU.Adapt.CHN.sol.rad.60d <- cumimp(gf.CHN.Adaptive.SNPs, "sol.rad.60d")

cols <- rainbow(length(levels(gf.CHN.Adaptive.SNPs$res.u$spec)))
names(cols) <- levels(gf.CHN.Adaptive.SNPs$res.u$spec)

isub.CHN.sol.rad.60d <- best.Adapt.CHN.sol.rad.60d$isub  ##this is the index of all the loci to be used in the plot

ymax=1.2 #set the ymax so that it's the same for the the plots to be overlaid
plot(CU.Adapt.CHN.sol.rad.60d$x[isub.CHN.sol.rad.60d], CU.Adapt.CHN.sol.rad.60d$y[isub.CHN.sol.rad.60d], type = "s", show.overall=T, ylab = "Cumulative Importance", xlab="sol.rad.60d", ylim=c(0,0.02), xlim=c(50000, 100000), lty=3, lwd=1.5)
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", xlim=c(5000, 13500), ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border


plot(CU.season.best$X.1080976.29, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1227797.33, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1310662.61, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1669736.10, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.790027.13, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1531198.75, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1723931.38, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.469873.12, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.249703.61, xlim=c(0,260), ylim=c(0,0.5), type="l")
par(new=T)
plot(CU.season.best$X.1419042.71, xlim=c(0,260), ylim=c(0,0.5), type="l")




```


