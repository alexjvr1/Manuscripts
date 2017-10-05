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

###Change Duplicated Adaptive Loci names to vcf formate

##R

CHN.AdaptiveLoci <- read.table("CHN.duplicated.outliers")
head(CHN.AdaptiveLoci)
CHN.AdaptiveLoci <- gsub("\\.", ":", CHN.AdaptiveLoci$V1)
head(CHN.AdaptiveLoci)
write.table(CHN.AdaptiveLoci, "CHN.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.AdaptiveLoci.names --recode --recode-INFO-all --out CHN.229.247

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.247.recode.vcf --plink --out CHN.229.247.plink

plink --file CHN.229.247.plink --noweb --recode --recodeA --out CHN.229.247.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster

```
##in R

CHN.nosex <- read.table("CHN.229.247.plink.nosex", header=F)
head(CHN.nosex)
CHN.pop <- gsub("_\\d+", "", CHN.nosex$V1)
CHN.pop
CHN.nosex$V3 <- CHN.pop
head(CHN.nosex)
write.table(CHN.nosex, "CHN.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHN.229.247.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.247

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.Adaptive.MAF <- read.table("CHN.229.247.frq.strat", header=T)
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
write.csv(CHN.Adaptive.MAF3, file="CHN.229.247.Adaptive.MAF.csv")
```


##### 2. Neutral loci



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

CZ.AdaptiveLoci <- read.table("CZ.duplicated.outliers")
head(CZ.AdaptiveLoci)
CZ.AdaptiveLoci <- gsub("\\.", ":", CZ.AdaptiveLoci$V1)
head(CZ.AdaptiveLoci)
write.table(CZ.AdaptiveLoci, "CZ.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --snps CZ.AdaptiveLoci.names --recode --recode-INFO-all --out CZ.404.212

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.212.recode.vcf --plink --out CZ.404.212.plink

plink --file CZ.404.212.plink --noweb --recode --recodeA --out CZ.404.212.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster

```
##in R

CZ.nosex <- read.table("CZ.404.212.plink.nosex", header=F)
head(CZ.nosex)
CZ.pop <- gsub("_\\d+", "", CZ.nosex$V1)
CZ.pop
CZ.nosex$V3 <- CZ.pop
head(CZ.nosex)
write.table(CZ.nosex, "CZ.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CZ.404.212.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.212

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.Adaptive.MAF <- read.table("CZ.404.212.frq.strat", header=T)
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
write.csv(CZ.Adaptive.MAF3, file="CZ.404.212.Adaptive.MAF.csv")
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

CHS.TI.AdaptiveLoci <- read.table("CHS.TI.duplicated.outliers")
head(CHS.TI.AdaptiveLoci)
CHS.TI.AdaptiveLoci <- gsub("\\.", ":", CHS.TI.AdaptiveLoci$V1)
head(CHS.TI.AdaptiveLoci)
write.table(CHS.TI.AdaptiveLoci, "CHS.TI.AdaptiveLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Adaptive loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.TI.140.233

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.233.recode.vcf --plink --out CHS.TI.140.233.plink

plink --file CHS.TI.140.233.plink --noweb --recode --recodeA --out CHS.TI.140.233.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster

```
##in R

CHS.TI.nosex <- read.table("CHS.TI.140.233.plink.nosex", header=F)
head(CHS.TI.nosex)
CHS.TI.pop <- gsub("_\\d+", "", CHS.TI.nosex$V1)
CHS.TI.pop
CHS.TI.nosex$V3 <- CHS.TI.pop
head(CHS.TI.nosex)
write.table(CHS.TI.nosex, "CHS.TI.PlinkCluster", quote=F, row.names=F, col.names=F)

```



And calculate MAF with Plink

```

plink --file CHS.TI.140.233.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.233

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.Adaptive.MAF <- read.table("CHS.TI.140.233.frq.strat", header=T)
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
write.csv(CHS.TI.Adaptive.MAF3, file="CHS.TI.140.233.Adaptive.MAF.csv")
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


### Run GF

1. Adaptive Loci
```
gf.CHS.VS.Adaptive <- read.csv("CHS.VS.135.358.Adaptive.MAF.csv", header=T)
envGF.CHS.VS.Adaptive <- gf.CHS.VS.Adaptive[,-1]
colnames(envGF.CHS.VS.Adaptive)

CHS.VS.AdaptiveSNPS <- CHS.VS.Adaptive.MAF3[,grep("X.", colnames(CHS.VS.Adaptive.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Adaptive)/2)
maxLevel

gf.CHS.VS.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.VS.Adaptive, CHS.VS.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.VS.Adaptive), response.vars=colnames(CHS.VS.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
```

2. Neutral Loci

```
gf.CHS.VS.Neutral <- read.csv("CHS.VS.135.Neutral.MAF.csv", header=T)
envGF.CHS.VS.Neutral <- gf.CHS.VS.Neutral[,-1]
colnames(envGF.CHS.VS.Neutral)

CHS.VS.Neutral.SNPs <- CHS.VS.Neutral.MAF3[,grep("X.", colnames(CHS.VS.Neutral.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Neutral)/2)
maxLevel

gf.CHS.VS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.VS.Neutral, CHS.VS.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.VS.Neutral), response.vars=colnames(CHS.VS.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
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
perf.NEUTRAL <- importance(gf.CHS.VS.Neutral.SNPs, type="Species")
perf.Fst <- importance(gf.CHS.VS.Adaptive.SNPs, type="Species")


##mean and range
summary(perf.NEUTRAL)
summary(perf.Fst)

perf.NEUT.df <- as.data.frame(perf.NEUTRAL)
perf.Fst.df <- as.data.frame(perf.Fst)


##count the number of loci above R2 of x (here 0.25)
length(perf.NEUT.df[which(perf.NEUT.df$perf.NEUTRAL>0.25),])
length(perf.Fst.df[which(perf.Fst.df>0.25),])



#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHS.VS.Neutral.SNPs
gf.CHS.VS.Adaptive.SNPs

#How many loci were originally included?


##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)
```




### Fig1: R2 heatmap
```
#R

#Retrieve results from the individual outputs for the three GF models: 

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
leg.txt <- c("Adaptive", "Neutral")  ##define the text that will be added to the legend

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
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "pcpt.60d") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "pcpt.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("pcpt.60d.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="pcpt.60d: Precipitation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##shadow.days
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "shadow.days") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "shadow.days")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("shadow.days.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="shadow.days", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##day10cm
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "day10cm") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "day10cm")

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
plot(gf.CHN.Adaptive.SNPs)
dev.off()

pdf("CHN.individual.NeutralSNPs.GF.pdf")
plot(gf.CHN.Neutral.SNPs)
dev.off()

### CZ
pdf("CZ.individual.AdaptiveSNPs.GF.pdf")
plot(gf.CZ.Adaptive.SNPs)
dev.off()

pdf("CZ.individual.NeutralSNPs.GF.pdf")
plot(gf.CZ.Neutral.SNPs)
dev.off()


### CHS
pdf("CHS.individual.AdaptiveSNPs.GF.pdf")
plot(gf.CHS.Adaptive.SNPs)
dev.off()

pdf("CHS.individual.NeutralSNPs.GF.pdf")
plot(gf.CHS.Neutral.SNPs)
dev.off()


### CHS.VS
pdf("CHS.VS.individual.AdaptiveSNPs.GF.pdf")
plot(gf.CHS.VS.Adaptive.SNPs)
dev.off()

pdf("CHS.VS.individual.NeutralSNPs.GF.pdf")
plot(gf.CHS.VS.Neutral.SNPs)
dev.off()

### CHS.TI
pdf("CHS.TI.individual.AdaptiveSNPs.GF.pdf")
plot(gf.CHS.TI.Adaptive.SNPs)
dev.off()

pdf("CHS.TI.individual.NeutralSNPs.GF.pdf")
plot(gf.CHS.TI.Neutral.SNPs)
dev.off()

```
![alt_txt][Fig3]

[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png



