# Gradient Forest Analysis

Aim: determine whether mean temp during development or growth season are the most important environmental drivers of adaptive divergence. 

I've identified candidate loci associated with mean temp 60 days and number of days above 6deg C. I have two sets of candidate loci for 
each of the transects: 

1. temp

any locus that was found in at least two of the following tests: 

bayenv2 for temp, lfmm for temp, pcadapt, XtX 100 top loci, and bayescan

2. season

any locus that was found in at least two of the following tests: 

bayenv2 for season, lfmm for season, pcadapt, XtX 100 top loci, and bayescan

3. Neutral loci

A set of 1000 loci from ALL - Candidate loci. 
It's important to subsample the data, since the datasets run for a long time if they have more than 1000 loci. 


Transects: 

CHN, CHS, CHS.VS, CHS.TI, CZ, and SE


# Gradient Forest CHall

### Datasets:

2 Datasets for each regions: Neutral, Potential Adaptive Loci (identified by at least 2 programs)

1. CHN

2. CZ

3. CHS

4. CHS.VS

5. CHS.TI

6. SE




I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013



## Input files

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/GradientForest

MEM and distance have been calculated for Swiss pops in Chapter 2, so can be copied directly from the .csv input files here: 

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017/*Input.csv

### CHN

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHN.TempLoci <- read.table("CHN.temp.duplicated.outliers.20171009")
head(CHN.TempLoci)
CHN.TempLoci <- gsub("\\.", ":", CHN.TempLoci$V1)
head(CHN.TempLoci)
write.table(CHN.TempLoci, "CHN.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.TempLoci.names --recode --recode-INFO-all --out CHN.229.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.TempLoci.recode.vcf --plink --out CHN.229.TempLoci.plink

plink --file CHN.229.TempLoci.plink --noweb --recode --recodeA --out CHN.229.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster

```
##in R

CHN.nosex <- read.table("CHN.229.TempLoci.plink.nosex", header=F)
head(CHN.nosex)
CHN.pop <- gsub("_\\d+", "", CHN.nosex$V1)
CHN.pop
CHN.nosex$V3 <- CHN.pop
head(CHN.nosex)
write.table(CHN.nosex, "CHN.PlinkCluster", quote=F, row.names=F, col.names=F)

```

And calculate MAF with Plink

```

plink --file CHN.229.TempLoci.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.TempLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.TempLoci.MAF <- read.table("CHN.229.TempLoci.frq.strat", header=T)
head(CHN.TempLoci.MAF)

CHN.TempLoci.MAF <- CHN.TempLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.TempLoci.MAF2 <- melt(CHN.TempLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.TempLoci.MAF2)
head(CHN.TempLoci.MAF2)


CHN.TempLoci.MAF3 <- dcast(CHN.TempLoci.MAF2, formula= CLST ~ SNP)
head(CHN.TempLoci.MAF3)
colnames(CHN.TempLoci.MAF3) <- paste("X", colnames(CHN.TempLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.TempLoci.MAF3, file="CHN.229.TempLoci.MAF.csv")
```



##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHN.SeasonLoci <- read.table("CHN.Season.duplicated.outliers.20171009")
head(CHN.SeasonLoci)
CHN.SeasonLoci <- gsub("\\.", ":", CHN.SeasonLoci$V1)
head(CHN.SeasonLoci)
write.table(CHN.SeasonLoci, "CHN.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.SeasonLoci.names --recode --recode-INFO-all --out CHN.229.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.SeasonLoci.recode.vcf --plink --out CHN.229.SeasonLoci.plink

plink --file CHN.229.SeasonLoci.plink --noweb --recode --recodeA --out CHN.229.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHN.229.SeasonLoci.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.SeasonLoci.MAF <- read.table("CHN.229.SeasonLoci.frq.strat", header=T)
head(CHN.SeasonLoci.MAF)

CHN.SeasonLoci.MAF <- CHN.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.SeasonLoci.MAF2 <- melt(CHN.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.SeasonLoci.MAF2)
head(CHN.SeasonLoci.MAF2)


CHN.SeasonLoci.MAF3 <- dcast(CHN.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHN.SeasonLoci.MAF3)
colnames(CHN.SeasonLoci.MAF3) <- paste("X", colnames(CHN.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.SeasonLoci.MAF3, file="CHN.229.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --exclude CHN.SeasonLoci.names --recode --recode-INFO-all --out CHN.229.Neutral.temporary
vcftools --vcf CHN.229.Neutral.temporary.recode.vcf --exclude CHN.TempLoci.names --recode --recode-INFO-all --out CHN.229.Neutral

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




### CZ

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CZ.TempLoci <- read.table("CZ.temp.duplicated.outliers.20171009")
head(CZ.TempLoci)
CZ.TempLoci <- gsub("\\.", ":", CZ.TempLoci$V1)
head(CZ.TempLoci)
write.table(CZ.TempLoci, "CZ.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CZ.404.5265.recode.vcf --snps CZ.TempLoci.names --recode --recode-INFO-all --out CZ.404.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.TempLoci.recode.vcf --plink --out CZ.404.TempLoci.plink

plink --file CZ.404.TempLoci.plink --noweb --recode --recodeA --out CZ.404.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster

```
##in R

CZ.nosex <- read.table("CZ.404.TempLoci.plink.nosex", header=F)
head(CZ.nosex)
CZ.pop <- gsub("_\\d+", "", CZ.nosex$V1)
CZ.pop
CZ.nosex$V3 <- CZ.pop
head(CZ.nosex)
write.table(CZ.nosex, "CZ.PlinkCluster", quote=F, row.names=F, col.names=F)

```

And calculate MAF with Plink

```

plink --file CZ.404.TempLoci.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.TempLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.TempLoci.MAF <- read.table("CZ.404.TempLoci.frq.strat", header=T)
head(CZ.TempLoci.MAF)

CZ.TempLoci.MAF <- CZ.TempLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.TempLoci.MAF2 <- melt(CZ.TempLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.TempLoci.MAF2)
head(CZ.TempLoci.MAF2)


CZ.TempLoci.MAF3 <- dcast(CZ.TempLoci.MAF2, formula= CLST ~ SNP)
head(CZ.TempLoci.MAF3)
colnames(CZ.TempLoci.MAF3) <- paste("X", colnames(CZ.TempLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.TempLoci.MAF3, file="CZ.404.TempLoci.MAF.csv")
```


##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CZ.SeasonLoci <- read.table("CZ.Season.duplicated.outliers.20171009")
head(CZ.SeasonLoci)
CZ.SeasonLoci <- gsub("\\.", ":", CZ.SeasonLoci$V1)
head(CZ.SeasonLoci)
write.table(CZ.SeasonLoci, "CZ.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --snps CZ.SeasonLoci.names --recode --recode-INFO-all --out CZ.404.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.SeasonLoci.recode.vcf --plink --out CZ.404.SeasonLoci.plink

plink --file CZ.404.SeasonLoci.plink --noweb --recode --recodeA --out CZ.404.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CZ.404.SeasonLoci.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.SeasonLoci.MAF <- read.table("CZ.404.SeasonLoci.frq.strat", header=T)
head(CZ.SeasonLoci.MAF)

CZ.SeasonLoci.MAF <- CZ.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.SeasonLoci.MAF2 <- melt(CZ.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.SeasonLoci.MAF2)
head(CZ.SeasonLoci.MAF2)


CZ.SeasonLoci.MAF3 <- dcast(CZ.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CZ.SeasonLoci.MAF3)
colnames(CZ.SeasonLoci.MAF3) <- paste("X", colnames(CZ.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.SeasonLoci.MAF3, file="CZ.404.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --exclude CZ.SeasonLoci.names --recode --recode-INFO-all --out CZ.404.Neutral.temporary
vcftools --vcf CZ.404.Neutral.temporary.recode.vcf --exclude CZ.TempLoci.names --recode --recode-INFO-all --out CZ.404.Neutral

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


### CHS

##### Dist and MEM


```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHS <- read.csv("AllEnv.Data_CHS_20161010.csv", header=T)
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
Number of regions: 25 
Number of nonzero links: 48 
Percentage nonzero weights: 7.68 
Average number of links: 1.92 

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
1  1.055611750 0.001
2  1.001597791 0.001
3  0.953101581 0.001
4  0.886104818 0.001
5  0.750000000 0.001
6  0.741690676 0.001
7  0.735557055 0.001
8  0.494402427 0.003
9  0.426721469 0.008
10 0.009058062 0.383

#9 significant MEM eigenfunctions with positive correlations. OBS. use first half (5 in our case)

write.table (CHS.scores$vectors[,1], "CHS.scores_MEM1.txt") 
write.table (CHS.scores$vectors[,2], "CHS.scores_MEM2.txt") 
write.table (CHS.scores$vectors[,3], "CHS.scores_MEM3.txt") 
write.table (CHS.scores$vectors[,4], "CHS.scores_MEM4.txt")
write.table (CHS.scores$vectors[,5], "CHS.scores_MEM5.txt")

```

######## Distance

```
##Geographic distance between all CHS populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)



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





##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TempLoci <- read.table("CHS.temp.duplicated.outliers.20171009")
head(CHS.TempLoci)
CHS.TempLoci <- gsub("\\.", ":", CHS.TempLoci$V1)
head(CHS.TempLoci)
write.table(CHS.TempLoci, "CHS.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --snps CHS.TempLoci.names --recode --recode-INFO-all --out CHS.275.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.TempLoci.recode.vcf --plink --out CHS.275.TempLoci.plink

plink --file CHS.275.TempLoci.plink --noweb --recode --recodeA --out CHS.275.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.PlinkCluster

```
##in R

CHS.nosex <- read.table("CHS.275.TempLoci.plink.nosex", header=F)
head(CHS.nosex)
CHS.pop <- gsub("_\\d+", "", CHS.nosex$V1)
CHS.pop
CHS.nosex$V3 <- CHS.pop
head(CHS.nosex)
write.table(CHS.nosex, "CHS.PlinkCluster", quote=F, row.names=F, col.names=F)

```

And calculate MAF with Plink

```

plink --file CHS.275.TempLoci.plink --within CHS.PlinkCluster --freq --noweb --out CHS.275.TempLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TempLoci.MAF <- read.table("CHS.275.TempLoci.frq.strat", header=T)
head(CHS.TempLoci.MAF)

CHS.TempLoci.MAF <- CHS.TempLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TempLoci.MAF2 <- melt(CHS.TempLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TempLoci.MAF2)
head(CHS.TempLoci.MAF2)


CHS.TempLoci.MAF3 <- dcast(CHS.TempLoci.MAF2, formula= CLST ~ SNP)
head(CHS.TempLoci.MAF3)
colnames(CHS.TempLoci.MAF3) <- paste("X", colnames(CHS.TempLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TempLoci.MAF3, file="CHS.275.TempLoci.MAF.csv")
```


##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.SeasonLoci <- read.table("CHS.Season.duplicated.outliers.20171009")
head(CHS.SeasonLoci)
CHS.SeasonLoci <- gsub("\\.", ":", CHS.SeasonLoci$V1)
head(CHS.SeasonLoci)
write.table(CHS.SeasonLoci, "CHS.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --snps CHS.SeasonLoci.names --recode --recode-INFO-all --out CHS.275.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.SeasonLoci.recode.vcf --plink --out CHS.275.SeasonLoci.plink

plink --file CHS.275.SeasonLoci.plink --noweb --recode --recodeA --out CHS.275.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.275.SeasonLoci.plink --within CHS.PlinkCluster --freq --noweb --out CHS.275.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.SeasonLoci.MAF <- read.table("CHS.275.SeasonLoci.frq.strat", header=T)
head(CHS.SeasonLoci.MAF)

CHS.SeasonLoci.MAF <- CHS.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.SeasonLoci.MAF2 <- melt(CHS.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.SeasonLoci.MAF2)
head(CHS.SeasonLoci.MAF2)


CHS.SeasonLoci.MAF3 <- dcast(CHS.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.SeasonLoci.MAF3)
colnames(CHS.SeasonLoci.MAF3) <- paste("X", colnames(CHS.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.SeasonLoci.MAF3, file="CHS.275.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --exclude CHS.SeasonLoci.names --recode --recode-INFO-all --out CHS.275.Neutral.temporary
vcftools --vcf CHS.275.Neutral.temporary.recode.vcf --exclude CHS.TempLoci.names --recode --recode-INFO-all --out CHS.275.Neutral

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




### CHS.TI

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TI.TempLoci <- read.table("CHS.TI.temp.duplicated.outliers.20171009")
head(CHS.TI.TempLoci)
CHS.TI.TempLoci <- gsub("\\.", ":", CHS.TI.TempLoci$V1)
head(CHS.TI.TempLoci)
write.table(CHS.TI.TempLoci, "CHS.TI.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.TempLoci.names --recode --recode-INFO-all --out CHS.TI.140.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.TempLoci.recode.vcf --plink --out CHS.TI.140.TempLoci.plink

plink --file CHS.TI.140.TempLoci.plink --noweb --recode --recodeA --out CHS.TI.140.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster

```
##in R

CHS.TI.nosex <- read.table("CHS.TI.140.TempLoci.plink.nosex", header=F)
head(CHS.TI.nosex)
CHS.TI.pop <- gsub("_\\d+", "", CHS.TI.nosex$V1)
CHS.TI.pop
CHS.TI.nosex$V3 <- CHS.TI.pop
head(CHS.TI.nosex)
write.table(CHS.TI.nosex, "CHS.TI.PlinkCluster", quote=F, row.names=F, col.names=F)

```

And calculate MAF with Plink

```

plink --file CHS.TI.140.TempLoci.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.TempLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.TempLoci.MAF <- read.table("CHS.TI.140.TempLoci.frq.strat", header=T)
head(CHS.TI.TempLoci.MAF)

CHS.TI.TempLoci.MAF <- CHS.TI.TempLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.TempLoci.MAF2 <- melt(CHS.TI.TempLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.TempLoci.MAF2)
head(CHS.TI.TempLoci.MAF2)


CHS.TI.TempLoci.MAF3 <- dcast(CHS.TI.TempLoci.MAF2, formula= CLST ~ SNP)
head(CHS.TI.TempLoci.MAF3)
colnames(CHS.TI.TempLoci.MAF3) <- paste("X", colnames(CHS.TI.TempLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.TempLoci.MAF3, file="CHS.TI.140.TempLoci.MAF.csv")
```


##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TI.SeasonLoci <- read.table("CHS.TI.Season.duplicated.outliers.20171009")
head(CHS.TI.SeasonLoci)
CHS.TI.SeasonLoci <- gsub("\\.", ":", CHS.TI.SeasonLoci$V1)
head(CHS.TI.SeasonLoci)
write.table(CHS.TI.SeasonLoci, "CHS.TI.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.SeasonLoci.names --recode --recode-INFO-all --out CHS.TI.140.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.SeasonLoci.recode.vcf --plink --out CHS.TI.140.SeasonLoci.plink

plink --file CHS.TI.140.SeasonLoci.plink --noweb --recode --recodeA --out CHS.TI.140.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.TI.140.SeasonLoci.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.SeasonLoci.MAF <- read.table("CHS.TI.140.SeasonLoci.frq.strat", header=T)
head(CHS.TI.SeasonLoci.MAF)

CHS.TI.SeasonLoci.MAF <- CHS.TI.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.SeasonLoci.MAF2 <- melt(CHS.TI.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.SeasonLoci.MAF2)
head(CHS.TI.SeasonLoci.MAF2)


CHS.TI.SeasonLoci.MAF3 <- dcast(CHS.TI.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.TI.SeasonLoci.MAF3)
colnames(CHS.TI.SeasonLoci.MAF3) <- paste("X", colnames(CHS.TI.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.SeasonLoci.MAF3, file="CHS.TI.140.SeasonLoci.MAF.csv")
```

##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --exclude CHS.TI.SeasonLoci.names --recode --recode-INFO-all --out CHS.TI.140.Neutral.temporary
vcftools --vcf CHS.TI.140.Neutral.temporary.recode.vcf --exclude CHS.TI.TempLoci.names --recode --recode-INFO-all --out CHS.TI.140.Neutral

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


### SE


##### Dist and MEM


```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.SE <- read.csv("SE.Derived.EnvData_20171009.csv", header=T)
head(env.data.SE)
#extract x and y
SE.xy <- env.data.SE[, c("Long","Lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
SE.dxy <- distm(SE.xy)
SE.dxy <- as.dist(SE.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
SE.th <- give.thresh(SE.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
SE.nb1 <- mst.nb(SE.dxy)
SE.wh1 <- which(as.matrix(SE.dxy)==SE.th,arr.ind=TRUE)
plot(SE.nb1,SE.xy,pch=20,cex=2,lty=3)
lines(SE.xy[SE.wh1[1,],1],SE.xy[SE.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
SE.th 
#[1] 493707.1
SE.nb1
Neighbour list object:
Number of regions: 15 
Number of nonzero links: 28 
Percentage nonzero weights: 12.44444 
Average number of links: 1.866667  

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
SE.listw=nb2listw(SE.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(SE.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
SE.scores=scores.listw(SE.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(SE.scores,SE.listw,nsim=999)
       stat  pval
1 0.9934033 0.001
2 0.9164982 0.001
3 0.7945233 0.001
4 0.6336452 0.003
5 0.4408883 0.022
6 0.2268045 0.138

#5 significant MEM eigenfunctions with positive correlations. OBS. use first half (3 in our case)

write.table (SE.scores$vectors[,1], "SE.scores_MEM1.txt") 
write.table (SE.scores$vectors[,2], "SE.scores_MEM2.txt") 
write.table (SE.scores$vectors[,3], "SE.scores_MEM3.txt") 
```

######## Distance

```
##Geographic distance between all SE populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)



#rdist.earth (in fields package) wants only long & lat
SE_lon.lat <- cbind(env.data.SE$Long, env.data.SE$Lat)
SE_lon.lat

#calculate great circle distances
distance.matrix.SE <- rdist.earth(SE_lon.lat)
summary(distance.matrix.SE)
dim(distance.matrix.SE)

#and use only the lower half of the matrix
upper.tri(distance.matrix.SE)
distance.matrix.SE[lower.tri(distance.matrix.SE)]<-NA
distance.matrix.SE

#change from matrix to dataframe
SE.bli <- as.data.frame(distance.matrix.SE)
head(SE.bli)
colnames(SE.bli) <- env.data.SE$pop
rownames(SE.bli) <- env.data.SE$pop

SE.bli[lower.tri(SE.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
SE.bli=as.data.frame(as.table(as.matrix(SE.bli)))  #Turn into a 3-column table
SE.bli
SE.bli=na.omit(SE.bli)  #Get rid of the junk we flagged above
SE.bli
colnames(SE.bli)<-c("site1", "site2", "dist(km)")
head(SE.bli)

SE.bli2 <- SE.bli[sort(SE.bli$site2),]

head(SE.bli2)


##write to csv
write.csv(SE.bli, file="distance.SE.csv",row.names=F)

```


##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.TempLoci <- read.table("SE.temp.duplicated.outliers.20171009")
head(SE.TempLoci)
SE.TempLoci <- gsub("\\.", ":", SE.TempLoci$V1)
head(SE.TempLoci)
write.table(SE.TempLoci, "SE.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.TempLoci.names --recode --recode-INFO-all --out SE.132.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.TempLoci.recode.vcf --plink --out SE.132.TempLoci.plink

plink --file SE.132.TempLoci.plink --noweb --recode --recodeA --out SE.132.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster

For SE I created this in excel because the names were formatted strangely. 


And calculate MAF with Plink

```

plink --file SE.132.TempLoci.plink --within SE.PlinkCluster --freq --noweb --out SE.132.TempLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.TempLoci.MAF <- read.table("SE.132.TempLoci.frq.strat", header=T)
head(SE.TempLoci.MAF)

SE.TempLoci.MAF <- SE.TempLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.TempLoci.MAF2 <- melt(SE.TempLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.TempLoci.MAF2)
head(SE.TempLoci.MAF2)


SE.TempLoci.MAF3 <- dcast(SE.TempLoci.MAF2, formula= CLST ~ SNP)
head(SE.TempLoci.MAF3)
colnames(SE.TempLoci.MAF3) <- paste("X", colnames(SE.TempLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.TempLoci.MAF3, file="SE.132.TempLoci.MAF.csv")
```


##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.SeasonLoci <- read.table("SE.Season.duplicated.outliers.20171009")
head(SE.SeasonLoci)
SE.SeasonLoci <- gsub("\\.", ":", SE.SeasonLoci$V1)
head(SE.SeasonLoci)
write.table(SE.SeasonLoci, "SE.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.SeasonLoci.recode.vcf --plink --out SE.132.SeasonLoci.plink

plink --file SE.132.SeasonLoci.plink --noweb --recode --recodeA --out SE.132.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file SE.132.SeasonLoci.plink --within SE.PlinkCluster --freq --noweb --out SE.132.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.SeasonLoci.MAF <- read.table("SE.132.SeasonLoci.frq.strat", header=T)
head(SE.SeasonLoci.MAF)

SE.SeasonLoci.MAF <- SE.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.SeasonLoci.MAF2 <- melt(SE.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.SeasonLoci.MAF2)
head(SE.SeasonLoci.MAF2)


SE.SeasonLoci.MAF3 <- dcast(SE.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(SE.SeasonLoci.MAF3)
colnames(SE.SeasonLoci.MAF3) <- paste("X", colnames(SE.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.SeasonLoci.MAF3, file="SE.132.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --exclude SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.Neutral.temporary
vcftools --vcf SE.132.Neutral.temporary.recode.vcf --exclude SE.TempLoci.names --recode --recode-INFO-all --out SE.132.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.Neutral.recode.vcf --plink --out SE.132.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
SE.Neutral.loci.names.map <- read.table("SE.132.Neutral.plink.map", header=F)
SE.Neutral.loci.names <-  SE.Neutral.loci.names.map$V2
SE.Neutral.loci.names <- as.data.frame(SE.Neutral.loci.names)
SE.1000.loci.names <- SE.Neutral.loci.names[sample(nrow(SE.Neutral.loci.names),1000),]
SE.1000.loci.names <- as.data.frame(SE.1000.loci.names)
summary(SE.1000.loci.names)
write.table(SE.1000.loci.names, "SE.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf SE.132.Neutral.recode.vcf --snps SE.1000.loci.names --recode --recode-INFO-all --out SE.1000NeutralLoci

vcftools --vcf SE.1000NeutralLoci.recode.vcf --plink --out SE.1000NeutralLoci.plink
plink --file SE.1000NeutralLoci.plink --noweb --recode --recodeA --out SE.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file SE.1000NeutralLoci.plink --within SE.PlinkCluster --noweb --freq --out SE.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.Neutral.MAF <- read.table("SE.1000NeutralLoci.frq.strat", header=T)
head(SE.Neutral.MAF)

SE.Neutral.MAF <- SE.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.Neutral.MAF2 <- melt(SE.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.Neutral.MAF2)
head(SE.Neutral.MAF2)


SE.Neutral.MAF3 <- dcast(SE.Neutral.MAF2, formula= CLST ~ SNP)
head(SE.Neutral.MAF3)
colnames(SE.Neutral.MAF3) <- paste("X", colnames(SE.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.Neutral.MAF3, file="SE.132.Neutral.MAF.csv")
```






### CHS.VS

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.VS.TempLoci <- read.table("CHS.VS.temp.duplicated.outliers.20171009")
head(CHS.VS.TempLoci)
CHS.VS.TempLoci <- gsub("\\.", ":", CHS.VS.TempLoci$V1)
head(CHS.VS.TempLoci)
write.table(CHS.VS.TempLoci, "CHS.VS.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.TempLoci.names --recode --recode-INFO-all --out CHS.VS.135.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.TempLoci.recode.vcf --plink --out CHS.VS.135.TempLoci.plink

plink --file CHS.VS.135.TempLoci.plink --noweb --recode --recodeA --out CHS.VS.135.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster

```
##in R

CHS.VS.nosex <- read.table("CHS.VS.135.TempLoci.plink.nosex", header=F)
head(CHS.VS.nosex)
CHS.VS.pop <- gsub("_\\d+", "", CHS.VS.nosex$V1)
CHS.VS.pop
CHS.VS.nosex$V3 <- CHS.VS.pop
head(CHS.VS.nosex)
write.table(CHS.VS.nosex, "CHS.VS.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.VS.SeasonLoci <- read.table("CHS.VS.Season.duplicated.outliers.20171009")
head(CHS.VS.SeasonLoci)
CHS.VS.SeasonLoci <- gsub("\\.", ":", CHS.VS.SeasonLoci$V1)
head(CHS.VS.SeasonLoci)
write.table(CHS.VS.SeasonLoci, "CHS.VS.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.SeasonLoci.names --recode --recode-INFO-all --out CHS.VS.135.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.SeasonLoci.recode.vcf --plink --out CHS.VS.135.SeasonLoci.plink

plink --file CHS.VS.135.SeasonLoci.plink --noweb --recode --recodeA --out CHS.VS.135.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.VS.135.SeasonLoci.plink --within CHS.VS.PlinkCluster --freq --noweb --out CHS.VS.135.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.SeasonLoci.MAF <- read.table("CHS.VS.135.SeasonLoci.frq.strat", header=T)
head(CHS.VS.SeasonLoci.MAF)

CHS.VS.SeasonLoci.MAF <- CHS.VS.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.SeasonLoci.MAF2 <- melt(CHS.VS.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.SeasonLoci.MAF2)
head(CHS.VS.SeasonLoci.MAF2)


CHS.VS.SeasonLoci.MAF3 <- dcast(CHS.VS.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.VS.SeasonLoci.MAF3)
colnames(CHS.VS.SeasonLoci.MAF3) <- paste("X", colnames(CHS.VS.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.SeasonLoci.MAF3, file="CHS.VS.135.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --exclude CHS.VS.SeasonLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral.temporary
vcftools --vcf CHS.VS.135.Neutral.temporary.recode.vcf --exclude CHS.VS.TempLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral

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





### SE

```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.SE <- read.csv("SE.Derived.EnvData_20171009.csv", header=T)
head(env.data.SE)
#extract x and y
SE.xy <- env.data.SE[, c("Long","Lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
SE.dxy <- distm(SE.xy)
SE.dxy <- as.dist(SE.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
SE.th <- give.thresh(SE.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
SE.nb1 <- mst.nb(SE.dxy)
SE.wh1 <- which(as.matrix(SE.dxy)==SE.th,arr.ind=TRUE)
plot(SE.nb1,SE.xy,pch=20,cex=2,lty=3)
lines(SE.xy[SE.wh1[1,],1],SE.xy[SE.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
SE.th 
#[1] 493707.1
SE.nb1
Neighbour list object:
Number of regions: 15 
Number of nonzero links: 28 
Percentage nonzero weights: 12.44444 
Average number of links: 1.866667 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
SE.listw=nb2listw(SE.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(SE.listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
SE.scores=scores.listw(SE.listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(SE.scores,SE.listw,nsim=999)
       stat  pval
1 0.9934033 0.001
2 0.9164982 0.001
3 0.7945233 0.001
4 0.6336452 0.002
5 0.4408883 0.018
6 0.2268045 0.127

#5 significant MEM eigenfunctions with positive correlations. OBS. use first half (3 in our case)

write.table (SE.scores$vectors[,1], "SE.scores_MEM1.txt") 
write.table (SE.scores$vectors[,2], "SE.scores_MEM2.txt") 
write.table (SE.scores$vectors[,3], "SE.scores_MEM3.txt") 

```

## Great circle distance between all sampling points

```
##Geographic distance between all SE populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/GradientForest")

#rdist.earth (in fields package) wants only long & lat
SE_lon.lat <- cbind(env.data.SE$Long, env.data.SE$Lat)
SE_lon.lat

#calculate great circle distances
distance.matrix.SE <- rdist.earth(SE_lon.lat)
summary(distance.matrix.SE)
dim(distance.matrix.SE)

#and use only the lower half of the matrix
upper.tri(distance.matrix.SE)
distance.matrix.SE[lower.tri(distance.matrix.SE)]<-NA
distance.matrix.SE

#change from matrix to dataframe
SE.bli <- as.data.frame(distance.matrix.SE)
head(SE.bli)
colnames(SE.bli) <- env.data.SE$pop
rownames(SE.bli) <- env.data.SE$pop

SE.bli[lower.tri(SE.bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
SE.bli=as.data.frame(as.table(as.matrix(SE.bli)))  #Turn into a 3-column table
SE.bli
SE.bli=na.omit(SE.bli)  #Get rid of the junk we flagged above
SE.bli
colnames(SE.bli)<-c("site1", "site2", "dist(km)")
head(SE.bli)

SE.bli2 <- SE.bli[sort(SE.bli$site2),]

head(SE.bli2)


##write to csv
write.csv(SE.bli, file="distance.SE.csv",row.names=F)
```

### SE

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.TempLoci <- read.table("SE.temp.duplicated.outliers.20171009")
head(SE.TempLoci)
SE.TempLoci <- gsub("\\.", ":", SE.TempLoci$V1)
head(SE.TempLoci)
write.table(SE.TempLoci, "SE.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.TempLoci.names --recode --recode-INFO-all --out SE.132.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.TempLoci.recode.vcf --plink --out SE.132.TempLoci.plink

plink --file SE.132.TempLoci.plink --noweb --recode --recodeA --out SE.132.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster

```
##in R

SE.nosex <- read.table("SE.132.TempLoci.plink.nosex", header=F)
head(SE.nosex)
SE.pop <- gsub("_\\d+", "", SE.nosex$V1)
SE.pop
SE.nosex$V3 <- SE.pop
head(SE.nosex)
write.table(SE.nosex, "SE.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.SeasonLoci <- read.table("SE.Season.duplicated.outliers.20171009")
head(SE.SeasonLoci)
SE.SeasonLoci <- gsub("\\.", ":", SE.SeasonLoci$V1)
head(SE.SeasonLoci)
write.table(SE.SeasonLoci, "SE.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.SeasonLoci.recode.vcf --plink --out SE.132.SeasonLoci.plink

plink --file SE.132.SeasonLoci.plink --noweb --recode --recodeA --out SE.132.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file SE.132.SeasonLoci.plink --within SE.PlinkCluster --freq --noweb --out SE.132.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.SeasonLoci.MAF <- read.table("SE.132.SeasonLoci.frq.strat", header=T)
head(SE.SeasonLoci.MAF)

SE.SeasonLoci.MAF <- SE.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.SeasonLoci.MAF2 <- melt(SE.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.SeasonLoci.MAF2)
head(SE.SeasonLoci.MAF2)


SE.SeasonLoci.MAF3 <- dcast(SE.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(SE.SeasonLoci.MAF3)
colnames(SE.SeasonLoci.MAF3) <- paste("X", colnames(SE.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.SeasonLoci.MAF3, file="SE.132.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --exclude SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.Neutral.temporary
vcftools --vcf SE.132.Neutral.temporary.recode.vcf --exclude SE.TempLoci.names --recode --recode-INFO-all --out SE.132.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.Neutral.recode.vcf --plink --out SE.132.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
SE.Neutral.loci.names.map <- read.table("SE.132.Neutral.plink.map", header=F)
SE.Neutral.loci.names <-  SE.Neutral.loci.names.map$V2
SE.Neutral.loci.names <- as.data.frame(SE.Neutral.loci.names)
SE.1000.loci.names <- SE.Neutral.loci.names[sample(nrow(SE.Neutral.loci.names),1000),]
SE.1000.loci.names <- as.data.frame(SE.1000.loci.names)
summary(SE.1000.loci.names)
write.table(SE.1000.loci.names, "SE.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf SE.132.Neutral.recode.vcf --snps SE.1000.loci.names --recode --recode-INFO-all --out SE.1000NeutralLoci

vcftools --vcf SE.1000NeutralLoci.recode.vcf --plink --out SE.1000NeutralLoci.plink
plink --file SE.1000NeutralLoci.plink --noweb --recode --recodeA --out SE.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file SE.1000NeutralLoci.plink --within SE.PlinkCluster --noweb --freq --out SE.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.Neutral.MAF <- read.table("SE.1000NeutralLoci.frq.strat", header=T)
head(SE.Neutral.MAF)

SE.Neutral.MAF <- SE.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.Neutral.MAF2 <- melt(SE.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.Neutral.MAF2)
head(SE.Neutral.MAF2)


SE.Neutral.MAF3 <- dcast(SE.Neutral.MAF2, formula= CLST ~ SNP)
head(SE.Neutral.MAF3)
colnames(SE.Neutral.MAF3) <- paste("X", colnames(SE.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.Neutral.MAF3, file="SE.132.Neutral.MAF.csv")
```




#### ENV data

Variables used: 

Same as for the EAA: 

1. mean.temp

2. season.length



Create a .csv input file with

1. pop names

2. distance

Use the last population (trpa) as distance 0. So all distances are measured from here.


##### Input Files

CHN.GF.TempLoci.Input.csv

CHN.GF.SeasonLoci.Input.csv

CHN.GF.NeutralLoci.Input.csv


CZ.GF.TempLoci.Input.csv

CZ.GF.NeutralLoci.Input.csv

CZ.GF.SeasonLoci.Input.csv


CHS.TI.GF.TempLoci.Input.csv

CHS.TI.GF.NeutralLoci.Input.csv

CHS.TI.GF.SeasonLoci.Input.csv


CHS.VS.GF.TempLoci.Input.csv

CHS.VS.GF.NeutralLoci.Input.csv

CHS.VS.GF.SeasonLoci.Input.csv


SE.GF.TempLoci.Input.csv

SE.GF.NeutralLoci.Input.csv

SE.GF.SeasonLoci.Input.csv


### Run GF

#### 1.Temp Loci
```
library(gradientForest)

###CHN

gf.CHN.TempLoci <- read.csv("CHN.229.TempLoci.MAF.csv", header=T)
envGF.CHN.TempLoci <- gf.CHN.TempLoci[,2:8]
colnames(envGF.CHN.TempLoci)

CHN.TempLoci.SNPS <- gf.CHN.TempLoci[,9:ncol(gf.CHN.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHN.TempLoci)/2)
maxLevel

gf.CHN.TempLoci.SNPs <- gradientForest(cbind(envGF.CHN.TempLoci, CHN.TempLoci.SNPS), predictor.vars=colnames(envGF.CHN.TempLoci), response.vars=colnames(CHN.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##4 warnings (loci variable in less than 5 populations)



###CHS

#library(gradientForest)
gf.CHS.TempLoci <- read.csv("CHS.275.TempLoci.MAF.csv", header=T)
colnames(gf.CHS.TempLoci)
envGF.CHS.TempLoci <- gf.CHS.TempLoci[,2:9]
colnames(envGF.CHS.TempLoci)

CHS.TempLoci.SNPS <- gf.CHS.TempLoci[,10:ncol(gf.CHS.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.TempLoci)/2)
maxLevel

gf.CHS.TempLoci.SNPs <- gradientForest(cbind(envGF.CHS.TempLoci, CHS.TempLoci.SNPS), predictor.vars=colnames(envGF.CHS.TempLoci), response.vars=colnames(CHS.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##4 warnings (loci variable in less than 5 populations)



###CZ

#library(gradientForest)
gf.CZ.TempLoci <- read.csv("CZ.404.TempLoci.MAF.csv", header=T)
colnames(gf.CZ.TempLoci)
envGF.CZ.TempLoci <- gf.CZ.TempLoci[,2:12]
colnames(envGF.CZ.TempLoci)

CZ.TempLoci.SNPS <- gf.CZ.TempLoci[,13:ncol(gf.CZ.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.CZ.TempLoci)/2)
maxLevel

gf.CZ.TempLoci.SNPs <- gradientForest(cbind(envGF.CZ.TempLoci, CZ.TempLoci.SNPS), predictor.vars=colnames(envGF.CZ.TempLoci), response.vars=colnames(CZ.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

#library(gradientForest)
gf.CHS.TI.TempLoci <- read.csv("CHS.TI.140.TempLoci.MAF.csv", header=T)
colnames(gf.CHS.TI.TempLoci)
envGF.CHS.TI.TempLoci <- gf.CHS.TI.TempLoci[,2:7]
colnames(envGF.CHS.TI.TempLoci)

CHS.TI.TempLoci.SNPS <- gf.CHS.TI.TempLoci[,8:ncol(gf.CHS.TI.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.TempLoci)/2)
maxLevel

gf.CHS.TI.TempLoci.SNPs <- gradientForest(cbind(envGF.CHS.TI.TempLoci, CHS.TI.TempLoci.SNPS), predictor.vars=colnames(envGF.CHS.TI.TempLoci), response.vars=colnames(CHS.TI.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##10 warnings (loci variable in less than 5 populations)



###CHS.VS

#library(gradientForest)
gf.CHS.VS.TempLoci <- read.csv("CHS.VS.135.TempLoci.MAF.csv", header=T)
colnames(gf.CHS.VS.TempLoci)
envGF.CHS.VS.TempLoci <- gf.CHS.VS.TempLoci[,2:6]
colnames(envGF.CHS.VS.TempLoci)

CHS.VS.TempLoci.SNPS <- gf.CHS.VS.TempLoci[,7:ncol(gf.CHS.VS.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.TempLoci)/2)
maxLevel

gf.CHS.VS.TempLoci.SNPs <- gradientForest(cbind(envGF.CHS.VS.TempLoci, CHS.VS.TempLoci.SNPS), predictor.vars=colnames(envGF.CHS.VS.TempLoci), response.vars=colnames(CHS.VS.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##32 warnings (loci variable in less than 5 populations)




###SE

#library(gradientForest)
gf.SE.TempLoci <- read.csv("SE.132.TempLoci.MAF.csv", header=T)
colnames(gf.SE.TempLoci)
envGF.SE.TempLoci <- gf.SE.TempLoci[,2:7]
colnames(envGF.SE.TempLoci)

SE.TempLoci.SNPS <- gf.SE.TempLoci[,8:ncol(gf.SE.TempLoci)]
maxLevel <- log2(0.368*nrow(envGF.SE.TempLoci)/2)
maxLevel

gf.SE.TempLoci.SNPs <- gradientForest(cbind(envGF.SE.TempLoci, SE.TempLoci.SNPS), predictor.vars=colnames(envGF.SE.TempLoci), response.vars=colnames(SE.TempLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##6 warnings (loci variable in less than 5 populations)

```


#### 1.Season Loci
```
library(gradientForest)

###CHN

gf.CHN.SeasonLoci <- read.csv("CHN.229.SeasonLoci.MAF.csv", header=T)
envGF.CHN.SeasonLoci <- gf.CHN.SeasonLoci[,2:8]
colnames(envGF.CHN.SeasonLoci)

CHN.SeasonLoci.SNPS <- gf.CHN.SeasonLoci[,9:ncol(gf.CHN.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHN.SeasonLoci)/2)
maxLevel

gf.CHN.SeasonLoci.SNPs <- gradientForest(cbind(envGF.CHN.SeasonLoci, CHN.SeasonLoci.SNPS), predictor.vars=colnames(envGF.CHN.SeasonLoci), response.vars=colnames(CHN.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##8 warnings (loci variable in less than 5 populations)



###CHS

#library(gradientForest)
gf.CHS.SeasonLoci <- read.csv("CHS.275.SeasonLoci.MAF.csv", header=T)
colnames(gf.CHS.SeasonLoci)
envGF.CHS.SeasonLoci <- gf.CHS.SeasonLoci[,2:9]
colnames(envGF.CHS.SeasonLoci)

CHS.SeasonLoci.SNPS <- gf.CHS.SeasonLoci[,10:ncol(gf.CHS.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.SeasonLoci)/2)
maxLevel

gf.CHS.SeasonLoci.SNPs <- gradientForest(cbind(envGF.CHS.SeasonLoci, CHS.SeasonLoci.SNPS), predictor.vars=colnames(envGF.CHS.SeasonLoci), response.vars=colnames(CHS.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##4 warnings (loci variable in less than 5 populations)



###CZ

#library(gradientForest)
gf.CZ.SeasonLoci <- read.csv("CZ.404.SeasonLoci.MAF.csv", header=T)
colnames(gf.CZ.SeasonLoci)
envGF.CZ.SeasonLoci <- gf.CZ.SeasonLoci[,2:12]
colnames(envGF.CZ.SeasonLoci)

CZ.SeasonLoci.SNPS <- gf.CZ.SeasonLoci[,13:ncol(gf.CZ.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.CZ.SeasonLoci)/2)
maxLevel

gf.CZ.SeasonLoci.SNPs <- gradientForest(cbind(envGF.CZ.SeasonLoci, CZ.SeasonLoci.SNPS), predictor.vars=colnames(envGF.CZ.SeasonLoci), response.vars=colnames(CZ.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

#library(gradientForest)
gf.CHS.TI.SeasonLoci <- read.csv("CHS.TI.140.SeasonLoci.MAF.csv", header=T)
colnames(gf.CHS.TI.SeasonLoci)
envGF.CHS.TI.SeasonLoci <- gf.CHS.TI.SeasonLoci[,2:7]
colnames(envGF.CHS.TI.SeasonLoci)

CHS.TI.SeasonLoci.SNPS <- gf.CHS.TI.SeasonLoci[,8:ncol(gf.CHS.TI.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.SeasonLoci)/2)
maxLevel

gf.CHS.TI.SeasonLoci.SNPs <- gradientForest(cbind(envGF.CHS.TI.SeasonLoci, CHS.TI.SeasonLoci.SNPS), predictor.vars=colnames(envGF.CHS.TI.SeasonLoci), response.vars=colnames(CHS.TI.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##17 warnings (loci variable in less than 5 populations)



###CHS.VS

#library(gradientForest)
gf.CHS.VS.SeasonLoci <- read.csv("CHS.VS.135.SeasonLoci.MAF.csv", header=T)
colnames(gf.CHS.VS.SeasonLoci)
envGF.CHS.VS.SeasonLoci <- gf.CHS.VS.SeasonLoci[,2:6]
colnames(envGF.CHS.VS.SeasonLoci)

CHS.VS.SeasonLoci.SNPS <- gf.CHS.VS.SeasonLoci[,7:ncol(gf.CHS.VS.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.SeasonLoci)/2)
maxLevel

gf.CHS.VS.SeasonLoci.SNPs <- gradientForest(cbind(envGF.CHS.VS.SeasonLoci, CHS.VS.SeasonLoci.SNPS), predictor.vars=colnames(envGF.CHS.VS.SeasonLoci), response.vars=colnames(CHS.VS.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##21 warnings (loci variable in less than 5 populations)




###SE

#library(gradientForest)
gf.SE.SeasonLoci <- read.csv("SE.132.SeasonLoci.MAF.csv", header=T)
colnames(gf.SE.SeasonLoci)
envGF.SE.SeasonLoci <- gf.SE.SeasonLoci[,2:7]
colnames(envGF.SE.SeasonLoci)

SE.SeasonLoci.SNPS <- gf.SE.SeasonLoci[,8:ncol(gf.SE.SeasonLoci)]
maxLevel <- log2(0.368*nrow(envGF.SE.SeasonLoci)/2)
maxLevel

gf.SE.SeasonLoci.SNPs <- gradientForest(cbind(envGF.SE.SeasonLoci, SE.SeasonLoci.SNPS), predictor.vars=colnames(envGF.SE.SeasonLoci), response.vars=colnames(SE.SeasonLoci.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##6 warnings (loci variable in less than 5 populations)

```

#### 3. Neutral Loci
```
library(gradientForest)

###CHN

gf.CHN.Neutral <- read.csv("CHN.229.Neutral.MAF.csv", header=T)
envGF.CHN.Neutral <- gf.CHN.Neutral[,2:8]
colnames(envGF.CHN.Neutral)

CHN.Neutral.SNPS <- gf.CHN.Neutral[,9:ncol(gf.CHN.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CHN.Neutral)/2)
maxLevel

gf.CHN.Neutral.SNPs <- gradientForest(cbind(envGF.CHN.Neutral, CHN.Neutral.SNPS), predictor.vars=colnames(envGF.CHN.Neutral), response.vars=colnames(CHN.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##3 warnings (loci variable in less than 5 populations)




###CHS

#library(gradientForest)
gf.CHS.Neutral <- read.csv("CHS.275.Neutral.MAF.csv", header=T)
colnames(gf.CHS.Neutral)
envGF.CHS.Neutral <- gf.CHS.Neutral[,2:9]
colnames(envGF.CHS.Neutral)

CHS.Neutral.SNPS <- gf.CHS.Neutral[,10:ncol(gf.CHS.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CHS.Neutral)/2)
maxLevel

gf.CHS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.Neutral, CHS.Neutral.SNPS), predictor.vars=colnames(envGF.CHS.Neutral), response.vars=colnames(CHS.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##2 warnings (loci variable in less than 5 populations)


###CZ

#library(gradientForest)
gf.CZ.Neutral <- read.csv("CZ.404.Neutral.MAF.csv", header=T)
colnames(gf.CZ.Neutral)
envGF.CZ.Neutral <- gf.CZ.Neutral[,2:12]
colnames(envGF.CZ.Neutral)

CZ.Neutral.SNPS <- gf.CZ.Neutral[,13:ncol(gf.CZ.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CZ.Neutral)/2)
maxLevel

gf.CZ.Neutral.SNPs <- gradientForest(cbind(envGF.CZ.Neutral, CZ.Neutral.SNPS), predictor.vars=colnames(envGF.CZ.Neutral), response.vars=colnames(CZ.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

#library(gradientForest)
gf.CHS.TI.Neutral <- read.csv("CHS.TI.140.Neutral.MAF.csv", header=T)
colnames(gf.CHS.TI.Neutral)
envGF.CHS.TI.Neutral <- gf.CHS.TI.Neutral[,2:7]
colnames(envGF.CHS.TI.Neutral)

CHS.TI.Neutral.SNPS <- gf.CHS.TI.Neutral[,8:ncol(gf.CHS.TI.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.Neutral)/2)
maxLevel

gf.CHS.TI.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.TI.Neutral, CHS.TI.Neutral.SNPS), predictor.vars=colnames(envGF.CHS.TI.Neutral), response.vars=colnames(CHS.TI.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##33 warnings (loci variable in less than 5 populations)


###CHS.VS

#library(gradientForest)
gf.CHS.VS.Neutral <- read.csv("CHS.VS.135.Neutral.MAF.csv", header=T)
colnames(gf.CHS.VS.Neutral)
envGF.CHS.VS.Neutral <- gf.CHS.VS.Neutral[,2:6]
colnames(envGF.CHS.VS.Neutral)

CHS.VS.Neutral.SNPS <- gf.CHS.VS.Neutral[,7:ncol(gf.CHS.VS.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Neutral)/2)
maxLevel

gf.CHS.VS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.VS.Neutral, CHS.VS.Neutral.SNPS), predictor.vars=colnames(envGF.CHS.VS.Neutral), response.vars=colnames(CHS.VS.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##47 warnings (loci variable in less than 5 populations)


###SE

#library(gradientForest)
gf.SE.Neutral <- read.csv("SE.132.Neutral.MAF.csv", header=T)
colnames(gf.SE.Neutral)
envGF.SE.Neutral <- gf.SE.Neutral[,2:7]
colnames(envGF.SE.Neutral)

SE.Neutral.SNPS <- gf.SE.Neutral[,8:ncol(gf.SE.Neutral)]
maxLevel <- log2(0.368*nrow(envGF.SE.Neutral)/2)
maxLevel

gf.SE.Neutral.SNPs <- gradientForest(cbind(envGF.SE.Neutral, SE.Neutral.SNPS), predictor.vars=colnames(envGF.SE.Neutral), response.vars=colnames(SE.Neutral.SNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##6 warnings (loci variable in less than 5 populations)
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

library(reshape)
library(ggplot2)

###CHN

#Retrieve results from the individual outputs for both GF models: 

R.sq.CHN.TempLoci <- (rowMeans(gf.CHN.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHN.TempLoci <- as.data.frame(R.sq.CHN.TempLoci)
R.sq.CHN.TempLoci
colnames(R.sq.CHN.TempLoci) <- "TempLoci"

R.sq.CHN.SeasonLoci <- (rowMeans(gf.CHN.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHN.SeasonLoci <- as.data.frame(R.sq.CHN.SeasonLoci)
R.sq.CHN.SeasonLoci
colnames(R.sq.CHN.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.CHN.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 7 variables
R.sq.CHN.EnvVariables <- R.sq.CHN.TempLoci[1:7,]  #get the first 8 rows from the Adaptive dataset
R.sq.CHN.EnvVariables <- as.data.frame(R.sq.CHN.EnvVariables)
colnames(R.sq.CHN.EnvVariables) <- "TempLoci"
row.names(R.sq.CHN.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2", "MEM3", "MEM4")

R.sq.CHN.EnvVariables$SeasonLoci <- R.sq.CHN.SeasonLoci[1:7,] 
R.sq.CHN.EnvVariables$Neutral <- R.sq.Neutral[1:7,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.CHN.EnvVariables <- as.matrix(R.sq.CHN.EnvVariables)  ##turn into matrix for heatmap
R.sq.CHN.all.melt <- melt(R.sq.CHN.EnvVariables)  ##melt for ggplot heatmap
R.sq.CHN.all.melt <- R.sq.CHN.all.melt[order(R.sq.CHN.all.melt$X1),]  ##order by predictor variable

R.sq.CHN.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHN.CHP4.R2plot.pdf")
ggplot(R.sq.CHN.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHN") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###CHS

#Retrieve results from the individual outputs for both GF models: 

R.sq.CHS.TempLoci <- (rowMeans(gf.CHS.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.TempLoci <- as.data.frame(R.sq.CHS.TempLoci)
R.sq.CHS.TempLoci
colnames(R.sq.CHS.TempLoci) <- "TempLoci"

R.sq.CHS.SeasonLoci <- (rowMeans(gf.CHS.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.SeasonLoci <- as.data.frame(R.sq.CHS.SeasonLoci)
R.sq.CHS.SeasonLoci
colnames(R.sq.CHS.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.CHS.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.CHS.EnvVariables <- R.sq.CHS.TempLoci[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.CHS.EnvVariables <- as.data.frame(R.sq.CHS.EnvVariables)
colnames(R.sq.CHS.EnvVariables) <- "TempLoci"
row.names(R.sq.CHS.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2", "MEM3", "MEM4", "MEM5")

R.sq.CHS.EnvVariables$SeasonLoci <- R.sq.CHS.SeasonLoci[1:8,] 
R.sq.CHS.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.CHS.EnvVariables <- as.matrix(R.sq.CHS.EnvVariables)  ##turn into matrix for heatmap
R.sq.CHS.all.melt <- melt(R.sq.CHS.EnvVariables)  ##melt for ggplot heatmap
R.sq.CHS.all.melt <- R.sq.CHS.all.melt[order(R.sq.CHS.all.melt$X1),]  ##order by predictor variable

R.sq.CHS.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.CHP4.R2plot.pdf")
ggplot(R.sq.CHS.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###CZ

#Retrieve results from the individual outputs for both GF models: 

R.sq.CZ.TempLoci <- (rowMeans(gf.CZ.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CZ.TempLoci <- as.data.frame(R.sq.CZ.TempLoci)
R.sq.CZ.TempLoci
colnames(R.sq.CZ.TempLoci) <- "TempLoci"

R.sq.CZ.SeasonLoci <- (rowMeans(gf.CZ.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CZ.SeasonLoci <- as.data.frame(R.sq.CZ.SeasonLoci)
R.sq.CZ.SeasonLoci
colnames(R.sq.CZ.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.CZ.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 11 variables
R.sq.CZ.EnvVariables <- R.sq.CZ.TempLoci[1:11,]  #get the first 8 rows from the Adaptive dataset
R.sq.CZ.EnvVariables <- as.data.frame(R.sq.CZ.EnvVariables)
colnames(R.sq.CZ.EnvVariables) <- "TempLoci"
row.names(R.sq.CZ.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2", "MEM3", "MEM4", "MEM5", "MEM6", "MEM7", "MEM8")

R.sq.CZ.EnvVariables$SeasonLoci <- R.sq.CZ.SeasonLoci[1:11,] 
R.sq.CZ.EnvVariables$Neutral <- R.sq.Neutral[1:11,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.CZ.EnvVariables <- as.matrix(R.sq.CZ.EnvVariables)  ##turn into matrix for heatmap
R.sq.CZ.all.melt <- melt(R.sq.CZ.EnvVariables)  ##melt for ggplot heatmap
R.sq.CZ.all.melt <- R.sq.CZ.all.melt[order(R.sq.CZ.all.melt$X1),]  ##order by predictor variable

R.sq.CZ.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CZ.CHP4.R2plot.pdf")
ggplot(R.sq.CZ.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CZ") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###CHS.VS

#Retrieve results from the individual outputs for both GF models: 

R.sq.CHS.VS.TempLoci <- (rowMeans(gf.CHS.VS.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.VS.TempLoci <- as.data.frame(R.sq.CHS.VS.TempLoci)
R.sq.CHS.VS.TempLoci
colnames(R.sq.CHS.VS.TempLoci) <- "TempLoci"

R.sq.CHS.VS.SeasonLoci <- (rowMeans(gf.CHS.VS.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.VS.SeasonLoci <- as.data.frame(R.sq.CHS.VS.SeasonLoci)
R.sq.CHS.VS.SeasonLoci
colnames(R.sq.CHS.VS.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.CHS.VS.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 5 variables
R.sq.CHS.VS.EnvVariables <- R.sq.CHS.VS.TempLoci[1:5,]  #get the first 8 rows from the Adaptive dataset
R.sq.CHS.VS.EnvVariables <- as.data.frame(R.sq.CHS.VS.EnvVariables)
colnames(R.sq.CHS.VS.EnvVariables) <- "TempLoci"
row.names(R.sq.CHS.VS.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2")

R.sq.CHS.VS.EnvVariables$SeasonLoci <- R.sq.CHS.VS.SeasonLoci[1:5,] 
R.sq.CHS.VS.EnvVariables$Neutral <- R.sq.Neutral[1:5,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.CHS.VS.EnvVariables <- as.matrix(R.sq.CHS.VS.EnvVariables)  ##turn into matrix for heatmap
R.sq.CHS.VS.all.melt <- melt(R.sq.CHS.VS.EnvVariables)  ##melt for ggplot heatmap
R.sq.CHS.VS.all.melt <- R.sq.CHS.VS.all.melt[order(R.sq.CHS.VS.all.melt$X1),]  ##order by predictor variable

R.sq.CHS.VS.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.VS.CHP4.R2plot.pdf")
ggplot(R.sq.CHS.VS.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.VS") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()




###CHS.TI

#Retrieve results from the individual outputs for both GF models: 

R.sq.CHS.TI.TempLoci <- (rowMeans(gf.CHS.TI.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.TI.TempLoci <- as.data.frame(R.sq.CHS.TI.TempLoci)
R.sq.CHS.TI.TempLoci
colnames(R.sq.CHS.TI.TempLoci) <- "TempLoci"

R.sq.CHS.TI.SeasonLoci <- (rowMeans(gf.CHS.TI.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.CHS.TI.SeasonLoci <- as.data.frame(R.sq.CHS.TI.SeasonLoci)
R.sq.CHS.TI.SeasonLoci
colnames(R.sq.CHS.TI.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.CHS.TI.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 6 variables
R.sq.CHS.TI.EnvVariables <- R.sq.CHS.TI.TempLoci[1:6,]  #get the first 8 rows from the Adaptive dataset
R.sq.CHS.TI.EnvVariables <- as.data.frame(R.sq.CHS.TI.EnvVariables)
colnames(R.sq.CHS.TI.EnvVariables) <- "TempLoci"
row.names(R.sq.CHS.TI.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2", "MEM3")

R.sq.CHS.TI.EnvVariables$SeasonLoci <- R.sq.CHS.TI.SeasonLoci[1:6,] 
R.sq.CHS.TI.EnvVariables$Neutral <- R.sq.Neutral[1:6,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.CHS.TI.EnvVariables <- as.matrix(R.sq.CHS.TI.EnvVariables)  ##turn into matrix for heatmap
R.sq.CHS.TI.all.melt <- melt(R.sq.CHS.TI.EnvVariables)  ##melt for ggplot heatmap
R.sq.CHS.TI.all.melt <- R.sq.CHS.TI.all.melt[order(R.sq.CHS.TI.all.melt$X1),]  ##order by predictor variable

R.sq.CHS.TI.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.TI.CHP4.R2plot.pdf")
ggplot(R.sq.CHS.TI.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.TI") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###SE

#Retrieve results from the individual outputs for both GF models: 

R.sq.SE.TempLoci <- (rowMeans(gf.SE.TempLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.SE.TempLoci <- as.data.frame(R.sq.SE.TempLoci)
R.sq.SE.TempLoci
colnames(R.sq.SE.TempLoci) <- "TempLoci"

R.sq.SE.SeasonLoci <- (rowMeans(gf.SE.SeasonLoci.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.SE.SeasonLoci <- as.data.frame(R.sq.SE.SeasonLoci)
R.sq.SE.SeasonLoci
colnames(R.sq.SE.SeasonLoci) <- "SeasonLoci"

R.sq.Neutral <- (rowMeans(gf.SE.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 6 variables
R.sq.SE.EnvVariables <- R.sq.SE.TempLoci[1:6,]  #get the first 8 rows from the Adaptive dataset
R.sq.SE.EnvVariables <- as.data.frame(R.sq.SE.EnvVariables)
colnames(R.sq.SE.EnvVariables) <- "TempLoci"
row.names(R.sq.SE.EnvVariables) <- c("dist", "mean.temp", "season.length", "MEM1", "MEM2", "MEM3")

R.sq.SE.EnvVariables$SeasonLoci <- R.sq.SE.SeasonLoci[1:6,] 
R.sq.SE.EnvVariables$Neutral <- R.sq.Neutral[1:6,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.SE.EnvVariables <- as.matrix(R.sq.SE.EnvVariables)  ##turn into matrix for heatmap
R.sq.SE.all.melt <- melt(R.sq.SE.EnvVariables)  ##melt for ggplot heatmap
R.sq.SE.all.melt <- R.sq.SE.all.melt[order(R.sq.SE.all.melt$X1),]  ##order by predictor variable

R.sq.SE.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("SE.CHP4.R2plot.pdf")
ggplot(R.sq.SE.all.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 SE") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()
```

![alt_txt][Fig1]

[Fig1]:



### Fig2: Average allelic turnover across all three datasets. 

It took me a while to find the code for the gf plots so that I could figure out what was being plotted, and so combine the plots. 
```
leg.txt <- c("Neutral", "Mean.Temp", "Season.Length")  ##define the text that will be added to the legend
```






##### CHS.TI Cumulative plot
```
##Temp.Loci

CU.Temp.CHS.TI.Temp.Loci <- cumimp(gf.CHS.TI.TempLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Season.CHS.TI.Temp.Loci <- cumimp(gf.CHS.TI.SeasonLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.TI.Temp.Loci <- cumimp(gf.CHS.TI.Neutral.SNPs, "mean.temp.60d")

isub.Temp.CHS.TI.Temp.Loci <- seq(1, length(CU.Temp.CHS.TI.Temp.Loci$x), len = pmin(500, length(CU.Temp.CHS.TI.Temp.Loci$x)))
isub.Season.CHS.TI.Temp.Loci <- seq(1, length(CU.Season.CHS.TI.Temp.Loci$x), len = pmin(500, length(CU.Season.CHS.TI.Temp.Loci$x)))
isub.Neutral.CHS.TI.Temp.Loci <- seq(1, length(CU.Neutral.CHS.TI.Temp.Loci$x), len = pmin(500, length(CU.Neutral.CHS.TI.Temp.Loci$x)))

ymax=0.04 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.CHS.TI.pdf")
plot(CU.Neutral.CHS.TI.Temp.Loci$x[isub.Neutral.CHS.TI.Temp.Loci], CU.Neutral.CHS.TI.Temp.Loci$y[isub.Neutral.CHS.TI.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.Temp.Loci$x[isub.Temp.CHS.TI.Temp.Loci], CU.Temp.CHS.TI.Temp.Loci$y[isub.Temp.CHS.TI.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.Temp.Loci$x[isub.Season.CHS.TI.Temp.Loci], CU.Season.CHS.TI.Temp.Loci$y[isub.Season.CHS.TI.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.CHS.TI.Season.Loci <- cumimp(gf.CHS.TI.TempLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Season.CHS.TI.Season.Loci <- cumimp(gf.CHS.TI.SeasonLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.TI.Season.Loci <- cumimp(gf.CHS.TI.Neutral.SNPs, "days.above.6")

isub.Temp.CHS.TI.Season.Loci <- seq(1, length(CU.Temp.CHS.TI.Season.Loci$x), len = pmin(500, length(CU.Temp.CHS.TI.Season.Loci$x)))
isub.Season.CHS.TI.Season.Loci <- seq(1, length(CU.Season.CHS.TI.Season.Loci$x), len = pmin(500, length(CU.Season.CHS.TI.Season.Loci$x)))
isub.Neutral.CHS.TI.Season.Loci <- seq(1, length(CU.Neutral.CHS.TI.Season.Loci$x), len = pmin(500, length(CU.Neutral.CHS.TI.Season.Loci$x)))

ymax=0.04 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.CHS.TI.pdf")
plot(CU.Neutral.CHS.TI.Season.Loci$x[isub.Neutral.CHS.TI.Season.Loci], CU.Neutral.CHS.TI.Season.Loci$y[isub.Neutral.CHS.TI.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.Season.Loci$x[isub.Temp.CHS.TI.Season.Loci], CU.Temp.CHS.TI.Season.Loci$y[isub.Temp.CHS.TI.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.Season.Loci$x[isub.Season.CHS.TI.Season.Loci], CU.Season.CHS.TI.Season.Loci$y[isub.Season.CHS.TI.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.CHS.TI.dist <- cumimp(gf.CHS.TI.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.CHS.TI.dist <- cumimp(gf.CHS.TI.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.TI.dist <- cumimp(gf.CHS.TI.Neutral.SNPs, "dist")

isub.Temp.CHS.TI.dist <- seq(1, length(CU.Temp.CHS.TI.dist$x), len = pmin(500, length(CU.Temp.CHS.TI.dist$x)))
isub.Season.CHS.TI.dist <- seq(1, length(CU.Season.CHS.TI.dist$x), len = pmin(500, length(CU.Season.CHS.TI.dist$x)))
isub.Neutral.CHS.TI.dist <- seq(1, length(CU.Neutral.CHS.TI.dist$x), len = pmin(500, length(CU.Neutral.CHS.TI.dist$x)))

ymax=0.1 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.CHS.TI.pdf")
plot(CU.Neutral.CHS.TI.dist$x[isub.Neutral.CHS.TI.dist], CU.Neutral.CHS.TI.dist$y[isub.Neutral.CHS.TI.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.dist$x[isub.Temp.CHS.TI.dist], CU.Temp.CHS.TI.dist$y[isub.Temp.CHS.TI.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.dist$x[isub.Season.CHS.TI.dist], CU.Season.CHS.TI.dist$y[isub.Season.CHS.TI.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```


##### CHS.VS Cumulative plot
```
##Temp.Loci

CU.Temp.CHS.VS.Temp.Loci <- cumimp(gf.CHS.VS.TempLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Season.CHS.VS.Temp.Loci <- cumimp(gf.CHS.VS.SeasonLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.VS.Temp.Loci <- cumimp(gf.CHS.VS.Neutral.SNPs, "mean.temp.60d")

isub.Temp.CHS.VS.Temp.Loci <- seq(1, length(CU.Temp.CHS.VS.Temp.Loci$x), len = pmin(500, length(CU.Temp.CHS.VS.Temp.Loci$x)))
isub.Season.CHS.VS.Temp.Loci <- seq(1, length(CU.Season.CHS.VS.Temp.Loci$x), len = pmin(500, length(CU.Season.CHS.VS.Temp.Loci$x)))
isub.Neutral.CHS.VS.Temp.Loci <- seq(1, length(CU.Neutral.CHS.VS.Temp.Loci$x), len = pmin(500, length(CU.Neutral.CHS.VS.Temp.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.CHS.VS.pdf")
plot(CU.Neutral.CHS.VS.Temp.Loci$x[isub.Neutral.CHS.VS.Temp.Loci], CU.Neutral.CHS.VS.Temp.Loci$y[isub.Neutral.CHS.VS.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.Temp.Loci$x[isub.Temp.CHS.VS.Temp.Loci], CU.Temp.CHS.VS.Temp.Loci$y[isub.Temp.CHS.VS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.Temp.Loci$x[isub.Season.CHS.VS.Temp.Loci], CU.Season.CHS.VS.Temp.Loci$y[isub.Season.CHS.VS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.CHS.VS.Season.Loci <- cumimp(gf.CHS.VS.TempLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Season.CHS.VS.Season.Loci <- cumimp(gf.CHS.VS.SeasonLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.VS.Season.Loci <- cumimp(gf.CHS.VS.Neutral.SNPs, "days.above.6")

isub.Temp.CHS.VS.Season.Loci <- seq(1, length(CU.Temp.CHS.VS.Season.Loci$x), len = pmin(500, length(CU.Temp.CHS.VS.Season.Loci$x)))
isub.Season.CHS.VS.Season.Loci <- seq(1, length(CU.Season.CHS.VS.Season.Loci$x), len = pmin(500, length(CU.Season.CHS.VS.Season.Loci$x)))
isub.Neutral.CHS.VS.Season.Loci <- seq(1, length(CU.Neutral.CHS.VS.Season.Loci$x), len = pmin(500, length(CU.Neutral.CHS.VS.Season.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.CHS.VS.pdf")
plot(CU.Neutral.CHS.VS.Season.Loci$x[isub.Neutral.CHS.VS.Season.Loci], CU.Neutral.CHS.VS.Season.Loci$y[isub.Neutral.CHS.VS.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.Season.Loci$x[isub.Temp.CHS.VS.Season.Loci], CU.Temp.CHS.VS.Season.Loci$y[isub.Temp.CHS.VS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.Season.Loci$x[isub.Season.CHS.VS.Season.Loci], CU.Season.CHS.VS.Season.Loci$y[isub.Season.CHS.VS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.CHS.VS.dist <- cumimp(gf.CHS.VS.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.CHS.VS.dist <- cumimp(gf.CHS.VS.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.VS.dist <- cumimp(gf.CHS.VS.Neutral.SNPs, "dist")

isub.Temp.CHS.VS.dist <- seq(1, length(CU.Temp.CHS.VS.dist$x), len = pmin(500, length(CU.Temp.CHS.VS.dist$x)))
isub.Season.CHS.VS.dist <- seq(1, length(CU.Season.CHS.VS.dist$x), len = pmin(500, length(CU.Season.CHS.VS.dist$x)))
isub.Neutral.CHS.VS.dist <- seq(1, length(CU.Neutral.CHS.VS.dist$x), len = pmin(500, length(CU.Neutral.CHS.VS.dist$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.CHS.VS.pdf")
plot(CU.Neutral.CHS.VS.dist$x[isub.Neutral.CHS.VS.dist], CU.Neutral.CHS.VS.dist$y[isub.Neutral.CHS.VS.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.dist$x[isub.Temp.CHS.VS.dist], CU.Temp.CHS.VS.dist$y[isub.Temp.CHS.VS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.dist$x[isub.Season.CHS.VS.dist], CU.Season.CHS.VS.dist$y[isub.Season.CHS.VS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```

##### CHN Cumulative plot
```
##Temp.Loci

CU.Temp.CHN.Temp.Loci <- cumimp(gf.CHN.TempLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Season.CHN.Temp.Loci <- cumimp(gf.CHN.SeasonLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Neutral.CHN.Temp.Loci <- cumimp(gf.CHN.Neutral.SNPs, "mean.temp.60d")

isub.Temp.CHN.Temp.Loci <- seq(1, length(CU.Temp.CHN.Temp.Loci$x), len = pmin(500, length(CU.Temp.CHN.Temp.Loci$x)))
isub.Season.CHN.Temp.Loci <- seq(1, length(CU.Season.CHN.Temp.Loci$x), len = pmin(500, length(CU.Season.CHN.Temp.Loci$x)))
isub.Neutral.CHN.Temp.Loci <- seq(1, length(CU.Neutral.CHN.Temp.Loci$x), len = pmin(500, length(CU.Neutral.CHN.Temp.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.CHN.pdf")
plot(CU.Neutral.CHN.Temp.Loci$x[isub.Neutral.CHN.Temp.Loci], CU.Neutral.CHN.Temp.Loci$y[isub.Neutral.CHN.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.Temp.Loci$x[isub.Temp.CHN.Temp.Loci], CU.Temp.CHN.Temp.Loci$y[isub.Temp.CHN.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.Temp.Loci$x[isub.Season.CHN.Temp.Loci], CU.Season.CHN.Temp.Loci$y[isub.Season.CHN.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.CHN.Season.Loci <- cumimp(gf.CHN.TempLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Season.CHN.Season.Loci <- cumimp(gf.CHN.SeasonLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Neutral.CHN.Season.Loci <- cumimp(gf.CHN.Neutral.SNPs, "days.above.6")

isub.Temp.CHN.Season.Loci <- seq(1, length(CU.Temp.CHN.Season.Loci$x), len = pmin(500, length(CU.Temp.CHN.Season.Loci$x)))
isub.Season.CHN.Season.Loci <- seq(1, length(CU.Season.CHN.Season.Loci$x), len = pmin(500, length(CU.Season.CHN.Season.Loci$x)))
isub.Neutral.CHN.Season.Loci <- seq(1, length(CU.Neutral.CHN.Season.Loci$x), len = pmin(500, length(CU.Neutral.CHN.Season.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.CHN.pdf")
plot(CU.Neutral.CHN.Season.Loci$x[isub.Neutral.CHN.Season.Loci], CU.Neutral.CHN.Season.Loci$y[isub.Neutral.CHN.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.Season.Loci$x[isub.Temp.CHN.Season.Loci], CU.Temp.CHN.Season.Loci$y[isub.Temp.CHN.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.Season.Loci$x[isub.Season.CHN.Season.Loci], CU.Season.CHN.Season.Loci$y[isub.Season.CHN.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.CHN.dist <- cumimp(gf.CHN.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.CHN.dist <- cumimp(gf.CHN.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.CHN.dist <- cumimp(gf.CHN.Neutral.SNPs, "dist")

isub.Temp.CHN.dist <- seq(1, length(CU.Temp.CHN.dist$x), len = pmin(500, length(CU.Temp.CHN.dist$x)))
isub.Season.CHN.dist <- seq(1, length(CU.Season.CHN.dist$x), len = pmin(500, length(CU.Season.CHN.dist$x)))
isub.Neutral.CHN.dist <- seq(1, length(CU.Neutral.CHN.dist$x), len = pmin(500, length(CU.Neutral.CHN.dist$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.CHN.pdf")
plot(CU.Neutral.CHN.dist$x[isub.Neutral.CHN.dist], CU.Neutral.CHN.dist$y[isub.Neutral.CHN.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.dist$x[isub.Temp.CHN.dist], CU.Temp.CHN.dist$y[isub.Temp.CHN.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.dist$x[isub.Season.CHN.dist], CU.Season.CHN.dist$y[isub.Season.CHN.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```


##### CZ Cumulative plot
```
##Temp.Loci

CU.Temp.CZ.Temp.Loci <- cumimp(gf.CZ.TempLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Season.CZ.Temp.Loci <- cumimp(gf.CZ.SeasonLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Neutral.CZ.Temp.Loci <- cumimp(gf.CZ.Neutral.SNPs, "mean.temp.60d")

isub.Temp.CZ.Temp.Loci <- seq(1, length(CU.Temp.CZ.Temp.Loci$x), len = pmin(500, length(CU.Temp.CZ.Temp.Loci$x)))
isub.Season.CZ.Temp.Loci <- seq(1, length(CU.Season.CZ.Temp.Loci$x), len = pmin(500, length(CU.Season.CZ.Temp.Loci$x)))
isub.Neutral.CZ.Temp.Loci <- seq(1, length(CU.Neutral.CZ.Temp.Loci$x), len = pmin(500, length(CU.Neutral.CZ.Temp.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.CZ.pdf")
plot(CU.Neutral.CZ.Temp.Loci$x[isub.Neutral.CZ.Temp.Loci], CU.Neutral.CZ.Temp.Loci$y[isub.Neutral.CZ.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.Temp.Loci$x[isub.Temp.CZ.Temp.Loci], CU.Temp.CZ.Temp.Loci$y[isub.Temp.CZ.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.Temp.Loci$x[isub.Season.CZ.Temp.Loci], CU.Season.CZ.Temp.Loci$y[isub.Season.CZ.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.CZ.Season.Loci <- cumimp(gf.CZ.TempLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Season.CZ.Season.Loci <- cumimp(gf.CZ.SeasonLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Neutral.CZ.Season.Loci <- cumimp(gf.CZ.Neutral.SNPs, "days.above.6")

isub.Temp.CZ.Season.Loci <- seq(1, length(CU.Temp.CZ.Season.Loci$x), len = pmin(500, length(CU.Temp.CZ.Season.Loci$x)))
isub.Season.CZ.Season.Loci <- seq(1, length(CU.Season.CZ.Season.Loci$x), len = pmin(500, length(CU.Season.CZ.Season.Loci$x)))
isub.Neutral.CZ.Season.Loci <- seq(1, length(CU.Neutral.CZ.Season.Loci$x), len = pmin(500, length(CU.Neutral.CZ.Season.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.CZ.pdf")
plot(CU.Neutral.CZ.Season.Loci$x[isub.Neutral.CZ.Season.Loci], CU.Neutral.CZ.Season.Loci$y[isub.Neutral.CZ.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.Season.Loci$x[isub.Temp.CZ.Season.Loci], CU.Temp.CZ.Season.Loci$y[isub.Temp.CZ.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.Season.Loci$x[isub.Season.CZ.Season.Loci], CU.Season.CZ.Season.Loci$y[isub.Season.CZ.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.CZ.dist <- cumimp(gf.CZ.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.CZ.dist <- cumimp(gf.CZ.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.CZ.dist <- cumimp(gf.CZ.Neutral.SNPs, "dist")

isub.Temp.CZ.dist <- seq(1, length(CU.Temp.CZ.dist$x), len = pmin(500, length(CU.Temp.CZ.dist$x)))
isub.Season.CZ.dist <- seq(1, length(CU.Season.CZ.dist$x), len = pmin(500, length(CU.Season.CZ.dist$x)))
isub.Neutral.CZ.dist <- seq(1, length(CU.Neutral.CZ.dist$x), len = pmin(500, length(CU.Neutral.CZ.dist$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.CZ.pdf")
plot(CU.Neutral.CZ.dist$x[isub.Neutral.CZ.dist], CU.Neutral.CZ.dist$y[isub.Neutral.CZ.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.dist$x[isub.Temp.CZ.dist], CU.Temp.CZ.dist$y[isub.Temp.CZ.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.dist$x[isub.Season.CZ.dist], CU.Season.CZ.dist$y[isub.Season.CZ.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```

##### CHS Cumulative plot
```
##Temp.Loci

CU.Temp.CHS.Temp.Loci <- cumimp(gf.CHS.TempLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Season.CHS.Temp.Loci <- cumimp(gf.CHS.SeasonLoci.SNPs, "mean.temp.60d") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.Temp.Loci <- cumimp(gf.CHS.Neutral.SNPs, "mean.temp.60d")

isub.Temp.CHS.Temp.Loci <- seq(1, length(CU.Temp.CHS.Temp.Loci$x), len = pmin(500, length(CU.Temp.CHS.Temp.Loci$x)))
isub.Season.CHS.Temp.Loci <- seq(1, length(CU.Season.CHS.Temp.Loci$x), len = pmin(500, length(CU.Season.CHS.Temp.Loci$x)))
isub.Neutral.CHS.Temp.Loci <- seq(1, length(CU.Neutral.CHS.Temp.Loci$x), len = pmin(500, length(CU.Neutral.CHS.Temp.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.CHS.pdf")
plot(CU.Neutral.CHS.Temp.Loci$x[isub.Neutral.CHS.Temp.Loci], CU.Neutral.CHS.Temp.Loci$y[isub.Neutral.CHS.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.Temp.Loci$x[isub.Temp.CHS.Temp.Loci], CU.Temp.CHS.Temp.Loci$y[isub.Temp.CHS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.Temp.Loci$x[isub.Season.CHS.Temp.Loci], CU.Season.CHS.Temp.Loci$y[isub.Season.CHS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.CHS.Season.Loci <- cumimp(gf.CHS.TempLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Season.CHS.Season.Loci <- cumimp(gf.CHS.SeasonLoci.SNPs, "days.above.6") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.Season.Loci <- cumimp(gf.CHS.Neutral.SNPs, "days.above.6")

isub.Temp.CHS.Season.Loci <- seq(1, length(CU.Temp.CHS.Season.Loci$x), len = pmin(500, length(CU.Temp.CHS.Season.Loci$x)))
isub.Season.CHS.Season.Loci <- seq(1, length(CU.Season.CHS.Season.Loci$x), len = pmin(500, length(CU.Season.CHS.Season.Loci$x)))
isub.Neutral.CHS.Season.Loci <- seq(1, length(CU.Neutral.CHS.Season.Loci$x), len = pmin(500, length(CU.Neutral.CHS.Season.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.CHS.pdf")
plot(CU.Neutral.CHS.Season.Loci$x[isub.Neutral.CHS.Season.Loci], CU.Neutral.CHS.Season.Loci$y[isub.Neutral.CHS.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.Season.Loci$x[isub.Temp.CHS.Season.Loci], CU.Temp.CHS.Season.Loci$y[isub.Temp.CHS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.Season.Loci$x[isub.Season.CHS.Season.Loci], CU.Season.CHS.Season.Loci$y[isub.Season.CHS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.CHS.dist <- cumimp(gf.CHS.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.CHS.dist <- cumimp(gf.CHS.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.CHS.dist <- cumimp(gf.CHS.Neutral.SNPs, "dist")

isub.Temp.CHS.dist <- seq(1, length(CU.Temp.CHS.dist$x), len = pmin(500, length(CU.Temp.CHS.dist$x)))
isub.Season.CHS.dist <- seq(1, length(CU.Season.CHS.dist$x), len = pmin(500, length(CU.Season.CHS.dist$x)))
isub.Neutral.CHS.dist <- seq(1, length(CU.Neutral.CHS.dist$x), len = pmin(500, length(CU.Neutral.CHS.dist$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.CHS.pdf")
plot(CU.Neutral.CHS.dist$x[isub.Neutral.CHS.dist], CU.Neutral.CHS.dist$y[isub.Neutral.CHS.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.dist$x[isub.Temp.CHS.dist], CU.Temp.CHS.dist$y[isub.Temp.CHS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.dist$x[isub.Season.CHS.dist], CU.Season.CHS.dist$y[isub.Season.CHS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```

##### SE Cumulative plot
```
##Temp.Loci

CU.Temp.SE.Temp.Loci <- cumimp(gf.SE.TempLoci.SNPs, "mean.temp.60.days") ##find the cumulative importance for each gf.model output
CU.Season.SE.Temp.Loci <- cumimp(gf.SE.SeasonLoci.SNPs, "mean.temp.60.days") ##find the cumulative importance for each gf.model output
CU.Neutral.SE.Temp.Loci <- cumimp(gf.SE.Neutral.SNPs, "mean.temp.60.days")

isub.Temp.SE.Temp.Loci <- seq(1, length(CU.Temp.SE.Temp.Loci$x), len = pmin(500, length(CU.Temp.SE.Temp.Loci$x)))
isub.Season.SE.Temp.Loci <- seq(1, length(CU.Season.SE.Temp.Loci$x), len = pmin(500, length(CU.Season.SE.Temp.Loci$x)))
isub.Neutral.SE.Temp.Loci <- seq(1, length(CU.Neutral.SE.Temp.Loci$x), len = pmin(500, length(CU.Neutral.SE.Temp.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Temp.Loci.CumImp.SE.pdf")
plot(CU.Neutral.SE.Temp.Loci$x[isub.Neutral.SE.Temp.Loci], CU.Neutral.SE.Temp.Loci$y[isub.Neutral.SE.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.Temp.Loci$x[isub.Temp.SE.Temp.Loci], CU.Temp.SE.Temp.Loci$y[isub.Temp.SE.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.Temp.Loci$x[isub.Season.SE.Temp.Loci], CU.Season.SE.Temp.Loci$y[isub.Season.SE.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##Season.Loci

CU.Temp.SE.Season.Loci <- cumimp(gf.SE.TempLoci.SNPs, "days.6.degrees") ##find the cumulative importance for each gf.model output
CU.Season.SE.Season.Loci <- cumimp(gf.SE.SeasonLoci.SNPs, "days.6.degrees") ##find the cumulative importance for each gf.model output
CU.Neutral.SE.Season.Loci <- cumimp(gf.SE.Neutral.SNPs, "days.6.degrees")

isub.Temp.SE.Season.Loci <- seq(1, length(CU.Temp.SE.Season.Loci$x), len = pmin(500, length(CU.Temp.SE.Season.Loci$x)))
isub.Season.SE.Season.Loci <- seq(1, length(CU.Season.SE.Season.Loci$x), len = pmin(500, length(CU.Season.SE.Season.Loci$x)))
isub.Neutral.SE.Season.Loci <- seq(1, length(CU.Neutral.SE.Season.Loci$x), len = pmin(500, length(CU.Neutral.SE.Season.Loci$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.Season.Loci.CumImp.SE.pdf")
plot(CU.Neutral.SE.Season.Loci$x[isub.Neutral.SE.Season.Loci], CU.Neutral.SE.Season.Loci$y[isub.Neutral.SE.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.Season.Loci$x[isub.Temp.SE.Season.Loci], CU.Temp.SE.Season.Loci$y[isub.Temp.SE.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.Season.Loci$x[isub.Season.SE.Season.Loci], CU.Season.SE.Season.Loci$y[isub.Season.SE.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist

CU.Temp.SE.dist <- cumimp(gf.SE.TempLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Season.SE.dist <- cumimp(gf.SE.SeasonLoci.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral.SE.dist <- cumimp(gf.SE.Neutral.SNPs, "dist")

isub.Temp.SE.dist <- seq(1, length(CU.Temp.SE.dist$x), len = pmin(500, length(CU.Temp.SE.dist$x)))
isub.Season.SE.dist <- seq(1, length(CU.Season.SE.dist$x), len = pmin(500, length(CU.Season.SE.dist$x)))
isub.Neutral.SE.dist <- seq(1, length(CU.Neutral.SE.dist$x), len = pmin(500, length(CU.Neutral.SE.dist$x)))

ymax=0.08 #set the ymax so that it's the same for the the plots to be overlaid

pdf("CHP4.dist.CumImp.SE.pdf")
plot(CU.Neutral.SE.dist$x[isub.Neutral.SE.dist], CU.Neutral.SE.dist$y[isub.Neutral.SE.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.dist$x[isub.Temp.SE.dist], CU.Temp.SE.dist$y[isub.Temp.SE.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.dist$x[isub.Season.SE.dist], CU.Season.SE.dist$y[isub.Season.SE.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```





## CumImP: Combined plots
ymax=0.1

```
##TEMP

#limits including SE in plot
xmin=7
xmax=11.5

ymax = 0.15

#limits with Swiss transects only
xmin=7
xmax=11.5

ymax = 0.05


#pdf("CHP4.Temp.Loci.CumImp.AllTransects.pdf")
pdf("CHP4.Temp.Loci.CumImp.SwissTransects.pdf")
plot(CU.Neutral.CHS.VS.Temp.Loci$x[isub.Neutral.CHS.VS.Temp.Loci], CU.Neutral.CHS.VS.Temp.Loci$y[isub.Neutral.CHS.VS.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.Temp.Loci$x[isub.Temp.CHS.VS.Temp.Loci], CU.Temp.CHS.VS.Temp.Loci$y[isub.Temp.CHS.VS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.Temp.Loci$x[isub.Season.CHS.VS.Temp.Loci], CU.Season.CHS.VS.Temp.Loci$y[isub.Season.CHS.VS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHS.TI.Temp.Loci$x[isub.Neutral.CHS.TI.Temp.Loci], CU.Neutral.CHS.TI.Temp.Loci$y[isub.Neutral.CHS.TI.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.Temp.Loci$x[isub.Temp.CHS.TI.Temp.Loci], CU.Temp.CHS.TI.Temp.Loci$y[isub.Temp.CHS.TI.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.Temp.Loci$x[isub.Season.CHS.TI.Temp.Loci], CU.Season.CHS.TI.Temp.Loci$y[isub.Season.CHS.TI.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHN.Temp.Loci$x[isub.Neutral.CHN.Temp.Loci], CU.Neutral.CHN.Temp.Loci$y[isub.Neutral.CHN.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.Temp.Loci$x[isub.Temp.CHN.Temp.Loci], CU.Temp.CHN.Temp.Loci$y[isub.Temp.CHN.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.Temp.Loci$x[isub.Season.CHN.Temp.Loci], CU.Season.CHN.Temp.Loci$y[isub.Season.CHN.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CZ.Temp.Loci$x[isub.Neutral.CZ.Temp.Loci], CU.Neutral.CZ.Temp.Loci$y[isub.Neutral.CZ.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.Temp.Loci$x[isub.Temp.CZ.Temp.Loci], CU.Temp.CZ.Temp.Loci$y[isub.Temp.CZ.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.Temp.Loci$x[isub.Season.CZ.Temp.Loci], CU.Season.CZ.Temp.Loci$y[isub.Season.CZ.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHS.Temp.Loci$x[isub.Neutral.CHS.Temp.Loci], CU.Neutral.CHS.Temp.Loci$y[isub.Neutral.CHS.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.Temp.Loci$x[isub.Temp.CHS.Temp.Loci], CU.Temp.CHS.Temp.Loci$y[isub.Temp.CHS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.Temp.Loci$x[isub.Season.CHS.Temp.Loci], CU.Season.CHS.Temp.Loci$y[isub.Season.CHS.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.SE.Temp.Loci$x[isub.Neutral.SE.Temp.Loci], CU.Neutral.SE.Temp.Loci$y[isub.Neutral.SE.Temp.Loci], type = "s", ylab ="Cumulative Importance", xlab="Temp.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.Temp.Loci$x[isub.Temp.SE.Temp.Loci], CU.Temp.SE.Temp.Loci$y[isub.Temp.SE.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.Temp.Loci$x[isub.Season.SE.Temp.Loci], CU.Season.SE.Temp.Loci$y[isub.Season.SE.Temp.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))

legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##SEASON
#limits including SE in plot
xmin=50
xmax=280

ymax=0.2

#limits for CHonly

xmin=50
xmax=280

ymax=0.05


#pdf("CHP4.Season.Loci.CumImp.AllTransects.pdf")
pdf("CHP4.Season.Loci.CumImp.SwissTransects.pdf")

plot(CU.Neutral.CHS.VS.Season.Loci$x[isub.Neutral.CHS.VS.Season.Loci], CU.Neutral.CHS.VS.Season.Loci$y[isub.Neutral.CHS.VS.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.Season.Loci$x[isub.Temp.CHS.VS.Season.Loci], CU.Temp.CHS.VS.Season.Loci$y[isub.Temp.CHS.VS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.Season.Loci$x[isub.Season.CHS.VS.Season.Loci], CU.Season.CHS.VS.Season.Loci$y[isub.Season.CHS.VS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border

par(new=T)
plot(CU.Neutral.CHS.TI.Season.Loci$x[isub.Neutral.CHS.TI.Season.Loci], CU.Neutral.CHS.TI.Season.Loci$y[isub.Neutral.CHS.TI.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.Season.Loci$x[isub.Temp.CHS.TI.Season.Loci], CU.Temp.CHS.TI.Season.Loci$y[isub.Temp.CHS.TI.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.Season.Loci$x[isub.Season.CHS.TI.Season.Loci], CU.Season.CHS.TI.Season.Loci$y[isub.Season.CHS.TI.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHN.Season.Loci$x[isub.Neutral.CHN.Season.Loci], CU.Neutral.CHN.Season.Loci$y[isub.Neutral.CHN.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.Season.Loci$x[isub.Temp.CHN.Season.Loci], CU.Temp.CHN.Season.Loci$y[isub.Temp.CHN.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.Season.Loci$x[isub.Season.CHN.Season.Loci], CU.Season.CHN.Season.Loci$y[isub.Season.CHN.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CZ.Season.Loci$x[isub.Neutral.CZ.Season.Loci], CU.Neutral.CZ.Season.Loci$y[isub.Neutral.CZ.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.Season.Loci$x[isub.Temp.CZ.Season.Loci], CU.Temp.CZ.Season.Loci$y[isub.Temp.CZ.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.Season.Loci$x[isub.Season.CZ.Season.Loci], CU.Season.CZ.Season.Loci$y[isub.Season.CZ.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHS.Season.Loci$x[isub.Neutral.CHS.Season.Loci], CU.Neutral.CHS.Season.Loci$y[isub.Neutral.CHS.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.Season.Loci$x[isub.Temp.CHS.Season.Loci], CU.Temp.CHS.Season.Loci$y[isub.Temp.CHS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.Season.Loci$x[isub.Season.CHS.Season.Loci], CU.Season.CHS.Season.Loci$y[isub.Season.CHS.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.SE.Season.Loci$x[isub.Neutral.SE.Season.Loci], CU.Neutral.SE.Season.Loci$y[isub.Neutral.SE.Season.Loci], type = "s", ylab ="Cumulative Importance", xlab="Season.Loci", ylim=c(0,ymax), lty=1, lwd=1.5, col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.Season.Loci$x[isub.Temp.SE.Season.Loci], CU.Temp.SE.Season.Loci$y[isub.Temp.SE.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.Season.Loci$x[isub.Season.SE.Season.Loci], CU.Season.SE.Season.Loci$y[isub.Season.SE.Season.Loci], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))

dev.off()


## DISTANCE
#limits including SE in plot
xmin=0
xmax=600

ymax=0.15

#limits for CHonly

xmin=0
xmax=100

ymax=0.12

#pdf("CHP4.dist.CumImp.AllTransects.pdf")
pdf("CHP4.dist.CumImp.SwissTransects.pdf")
plot(CU.Neutral.CHS.VS.dist$x[isub.Neutral.CHS.VS.dist], CU.Neutral.CHS.VS.dist$y[isub.Neutral.CHS.VS.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.VS.dist$x[isub.Temp.CHS.VS.dist], CU.Temp.CHS.VS.dist$y[isub.Temp.CHS.VS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.VS.dist$x[isub.Season.CHS.VS.dist], CU.Season.CHS.VS.dist$y[isub.Season.CHS.VS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHS.TI.dist$x[isub.Neutral.CHS.TI.dist], CU.Neutral.CHS.TI.dist$y[isub.Neutral.CHS.TI.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.TI.dist$x[isub.Temp.CHS.TI.dist], CU.Temp.CHS.TI.dist$y[isub.Temp.CHS.TI.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.TI.dist$x[isub.Season.CHS.TI.dist], CU.Season.CHS.TI.dist$y[isub.Season.CHS.TI.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="darkorchid1", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHN.dist$x[isub.Neutral.CHN.dist], CU.Neutral.CHN.dist$y[isub.Neutral.CHN.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHN.dist$x[isub.Temp.CHN.dist], CU.Temp.CHN.dist$y[isub.Temp.CHN.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHN.dist$x[isub.Season.CHN.dist], CU.Season.CHN.dist$y[isub.Season.CHN.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="springgreen4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CZ.dist$x[isub.Neutral.CZ.dist], CU.Neutral.CZ.dist$y[isub.Neutral.CZ.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CZ.dist$x[isub.Temp.CZ.dist], CU.Temp.CZ.dist$y[isub.Temp.CZ.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CZ.dist$x[isub.Season.CZ.dist], CU.Season.CZ.dist$y[isub.Season.CZ.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="violetred", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.CHS.dist$x[isub.Neutral.CHS.dist], CU.Neutral.CHS.dist$y[isub.Neutral.CHS.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.CHS.dist$x[isub.Temp.CHS.dist], CU.Temp.CHS.dist$y[isub.Temp.CHS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.CHS.dist$x[isub.Season.CHS.dist], CU.Season.CHS.dist$y[isub.Season.CHS.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="slategray4", xlim=c(xmin,xmax))

par(new=T)
plot(CU.Neutral.SE.dist$x[isub.Neutral.SE.dist], CU.Neutral.SE.dist$y[isub.Neutral.SE.dist], type = "s", ylab ="Cumulative Importance", xlab="dist", ylim=c(0,ymax), lty=1, lwd=1.5, col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Temp.SE.dist$x[isub.Temp.SE.dist], CU.Temp.SE.dist$y[isub.Temp.SE.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))
par(new=T)  ##allows you to overlay the plots
plot(CU.Season.SE.dist$x[isub.Season.SE.dist], CU.Season.SE.dist$y[isub.Season.SE.dist], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=4, lwd=1.5, xaxt='n', yaxt='n', col="steelblue3", xlim=c(xmin,xmax))


#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()
```








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
        imp.vars <- imp.var.names <- c("mean.temp.60.days", "days.6.degrees", "dist")
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
            ylim <- c(0,0.5)   ###I'm editing this so that I can plot on one plot
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
pdf("CHN.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.TempLoci.SNPs)
dev.off()

pdf("CHN.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.SeasonLoci.SNPs)
dev.off()

pdf("CHN.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.Neutral.SNPs)
dev.off()




### CZ
pdf("CZ.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.TempLoci.SNPs)
dev.off()

pdf("CZ.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.SeasonLoci.SNPs)
dev.off()

pdf("CZ.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.Neutral.SNPs)
dev.off()





### CHS
pdf("CHS.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TempLoci.SNPs)
dev.off()

pdf("CHS.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.SeasonLoci.SNPs)
dev.off()

pdf("CHS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.Neutral.SNPs)
dev.off()





### CHS.VS
pdf("CHS.VS.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.TempLoci.SNPs)
dev.off()

pdf("CHS.VS.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.SeasonLoci.SNPs)
dev.off()

pdf("CHS.VS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.Neutral.SNPs)
dev.off()




### CHS.TI
pdf("CHS.TI.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.TempLoci.SNPs)
dev.off()

pdf("CHS.TI.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.SeasonLoci.SNPs)
dev.off()

pdf("CHS.TI.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.Neutral.SNPs)
dev.off()



##redefine the function with the new column names
### SE
pdf("SE.individual.TempSNPs.GF.pdf")
species.cumulative.plot(gf.SE.TempLoci.SNPs)
dev.off()

pdf("SE.individual.SeasonSNPs.GF.pdf")
species.cumulative.plot(gf.SE.SeasonLoci.SNPs)
dev.off()

pdf("SE.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.SE.Neutral.SNPs)
dev.off()




```
![alt_txt][Fig3]

[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png


# Part2: Find thresholds

For this part I will run GF with all the loci and identify the loci that are most associated with 

1: temp

2: season. 

Questions: 

Is there a threshold response in any of the transects? 

Is the threshold the same between transects? 


Set up gf for each trasect using all of the loci: 


### CHN

```

vcftools --vcf CHN.229.5265.recode.vcf --plink --out CHN.229.5265.plink

plink --file CHN.229.5265.plink --recode --recodeA --noweb --out CHN.229.5265.plink

##calculate MAF within each population

plink --file CHN.229.5265.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.5265.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.FullModel.MAF <- read.table("CHN.229.5265.plink.frq.strat", header=T)
head(CHN.FullModel.MAF)

CHN.FullModel.MAF <- CHN.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.FullModel.MAF2 <- melt(CHN.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.FullModel.MAF2)
head(CHN.FullModel.MAF2)


CHN.FullModel.MAF3 <- dcast(CHN.FullModel.MAF2, formula= CLST ~ SNP)
head(CHN.FullModel.MAF3)
colnames(CHN.FullModel.MAF3) <- paste("X", colnames(CHN.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.FullModel.MAF3, file="CHN.229.5265.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.CHN.FullModel <- read.csv("CHN.229.5265.FullModel.MAF.csv", header=T)
envGF.CHN.FullModel <- gf.CHN.FullModel[,2:8]
colnames(envGF.CHN.Adaptive)

CHN.FullModelSNPS <- CHN.FullModel.MAF3[,grep("X.", colnames(CHN.FullModel.MAF3))]
CHN.FullModelSNPS <- CHN.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHN.FullModel)/2)
maxLevel

gf.CHN.FullModel.SNPs <- gradientForest(cbind(envGF.CHN.FullModel, CHN.FullModelSNPS), predictor.vars=colnames(envGF.CHN.FullModel), response.vars=colnames(CHN.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##25 warnings (loci variable in less than 5 populations)


```


### CHS

```
###CHS

vcftools --vcf CHS.275.6339.recode.vcf --plink --out CHS.275.6339.plink

plink --file CHS.275.6339.plink --recode --recodeA --noweb --out CHS.275.6339.plink

##calculate MAF within each population

plink --file CHS.275.6339.plink --within CHS.PlinkCluster --freq --noweb --out CHS.275.6339.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.FullModel.MAF <- read.table("CHS.275.6339.plink.frq.strat", header=T)
head(CHS.FullModel.MAF)

CHS.FullModel.MAF <- CHS.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.FullModel.MAF2 <- melt(CHS.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.FullModel.MAF2)
head(CHS.FullModel.MAF2)


CHS.FullModel.MAF3 <- dcast(CHS.FullModel.MAF2, formula= CLST ~ SNP)
head(CHS.FullModel.MAF3)
colnames(CHS.FullModel.MAF3) <- paste("X", colnames(CHS.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.FullModel.MAF3, file="CHS.275.6339.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.CHS.FullModel <- read.csv("CHS.275.6339.FullModel.MAF.csv", header=T)
envGF.CHS.FullModel <- gf.CHS.FullModel[,2:9]
colnames(envGF.CHS.FullModel)

CHS.FullModelSNPS <- gf.CHS.FullModel[,grep("X.", colnames(gf.CHS.FullModel))]
CHS.FullModelSNPS <- CHS.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.FullModel)/2)
maxLevel

gf.CHS.FullModel.SNPs <- gradientForest(cbind(envGF.CHS.FullModel, CHS.FullModelSNPS), predictor.vars=colnames(envGF.CHS.FullModel), response.vars=colnames(CHS.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##16 warnings (loci variable in less than 5 populations)


```


### CZ

```

vcftools --vcf CZ.404.7288.recode.vcf --plink --out CZ.404.7288.plink

plink --file CZ.404.7288.plink --recode --recodeA --noweb --out CZ.404.7288.plink

##calculate MAF within each population

plink --file CZ.404.7288.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.7288.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.FullModel.MAF <- read.table("CZ.404.7288.plink.frq.strat", header=T)
head(CZ.FullModel.MAF)

CZ.FullModel.MAF <- CZ.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.FullModel.MAF2 <- melt(CZ.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.FullModel.MAF2)
head(CZ.FullModel.MAF2)


CZ.FullModel.MAF3 <- dcast(CZ.FullModel.MAF2, formula= CLST ~ SNP)
head(CZ.FullModel.MAF3)
colnames(CZ.FullModel.MAF3) <- paste("X", colnames(CZ.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.FullModel.MAF3, file="CZ.404.7288.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.CZ.FullModel <- read.csv("CZ.404.7288.FullModel.MAF.csv", header=T)
envGF.CZ.FullModel <- gf.CZ.FullModel[,2:12]
colnames(envGF.CZ.FullModel)

CZ.FullModelSNPS <- gf.CZ.FullModel[,grep("X.", colnames(gf.CZ.FullModel))]
CZ.FullModelSNPS <- CZ.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CZ.FullModel)/2)
maxLevel

gf.CZ.FullModel.SNPs <- gradientForest(cbind(envGF.CZ.FullModel, CZ.FullModelSNPS), predictor.vars=colnames(envGF.CZ.FullModel), response.vars=colnames(CZ.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##25 warnings (loci variable in less than 5 populations)


```


### CHS.VS

```
###CHS.VS

vcftools --vcf CHS.VS.135.5835.recode.vcf --plink --out CHS.VS.135.5835.plink

plink --file CHS.VS.135.5835.plink --recode --recodeA --noweb --out CHS.VS.135.5835.plink

##calculate MAF within each population

plink --file CHS.VS.135.5835.plink --within CHS.VS.PlinkCluster --freq --noweb --out CHS.VS.135.5835.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.FullModel.MAF <- read.table("CHS.VS.135.5835.plink.frq.strat", header=T)
head(CHS.VS.FullModel.MAF)

CHS.VS.FullModel.MAF <- CHS.VS.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.FullModel.MAF2 <- melt(CHS.VS.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.FullModel.MAF2)
head(CHS.VS.FullModel.MAF2)


CHS.VS.FullModel.MAF3 <- dcast(CHS.VS.FullModel.MAF2, formula= CLST ~ SNP)
head(CHS.VS.FullModel.MAF3)
colnames(CHS.VS.FullModel.MAF3) <- paste("X", colnames(CHS.VS.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.FullModel.MAF3, file="CHS.VS.135.5835.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.CHS.VS.FullModel <- read.csv("CHS.VS.135.5835.FullModel.MAF.csv", header=T)
envGF.CHS.VS.FullModel <- gf.CHS.VS.FullModel[,2:8]
colnames(envGF.CHS.VS.FullModel)

CHS.VS.FullModelSNPS <- gf.CHS.VS.FullModel[,grep("X.", colnames(gf.CHS.VS.FullModel))]
CHS.VS.FullModelSNPS <- CHS.VS.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.FullModel)/2)
maxLevel

gf.CHS.VS.FullModel.SNPs <- gradientForest(cbind(envGF.CHS.VS.FullModel, CHS.VS.FullModelSNPS), predictor.vars=colnames(envGF.CHS.VS.FullModel), response.vars=colnames(CHS.VS.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##50+ warnings (loci variable in less than 5 populations)


```



### CHS.TI
```

vcftools --vcf CHS.TI.140.5692.recode.vcf --plink --out CHS.TI.140.5692.plink

plink --file CHS.TI.140.5692.plink --recode --recodeA --noweb --out CHS.TI.140.5692.plink

##calculate MAF within each population

plink --file CHS.TI.140.5692.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.5692.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.FullModel.MAF <- read.table("CHS.TI.140.5692.plink.frq.strat", header=T)
head(CHS.TI.FullModel.MAF)

CHS.TI.FullModel.MAF <- CHS.TI.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.FullModel.MAF2 <- melt(CHS.TI.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.FullModel.MAF2)
head(CHS.TI.FullModel.MAF2)


CHS.TI.FullModel.MAF3 <- dcast(CHS.TI.FullModel.MAF2, formula= CLST ~ SNP)
head(CHS.TI.FullModel.MAF3)
colnames(CHS.TI.FullModel.MAF3) <- paste("X", colnames(CHS.TI.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.FullModel.MAF3, file="CHS.TI.140.5692.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.CHS.TI.FullModel <- read.csv("CHS.TI.140.5692.FullModel.MAF.csv", header=T)
envGF.CHS.TI.FullModel <- gf.CHS.TI.FullModel[,2:7]
colnames(envGF.CHS.TI.FullModel)

CHS.TI.FullModelSNPS <- gf.CHS.TI.FullModel[,grep("X.", colnames(gf.CHS.TI.FullModel))]
CHS.TI.FullModelSNPS <- CHS.TI.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.FullModel)/2)
maxLevel

gf.CHS.TI.FullModel.SNPs <- gradientForest(cbind(envGF.CHS.TI.FullModel, CHS.TI.FullModelSNPS), predictor.vars=colnames(envGF.CHS.TI.FullModel), response.vars=colnames(CHS.TI.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##25 warnings (loci variable in less than 5 populations)


```



### SE
```

vcftools --vcf SE132.2027.recode.vcf --plink --out SE.132.2027.plink

plink --file SE.132.2027.plink --recode --recodeA --noweb --out SE.132.2027.plink

##calculate MAF within each population

plink --file SE.132.2027.plink --within SE.PlinkCluster --freq --noweb --out SE.132.2027.plink 


######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.FullModel.MAF <- read.table("SE.132.2027.plink.frq.strat", header=T)
head(SE.FullModel.MAF)

SE.FullModel.MAF <- SE.FullModel.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.FullModel.MAF2 <- melt(SE.FullModel.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.FullModel.MAF2)
head(SE.FullModel.MAF2)


SE.FullModel.MAF3 <- dcast(SE.FullModel.MAF2, formula= CLST ~ SNP)
head(SE.FullModel.MAF3)
colnames(SE.FullModel.MAF3) <- paste("X", colnames(SE.FullModel.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.FullModel.MAF3, file="SE.132.2027.FullModel.MAF.csv")

```

Run GF

```
library(gradientForest)

gf.SE.FullModel <- read.csv("SE.132.2027.FullModel.MAF.csv", header=T)
envGF.SE.FullModel <- gf.SE.FullModel[,2:7]
colnames(envGF.SE.FullModel)

SE.FullModelSNPS <- gf.SE.FullModel[,grep("X.", colnames(gf.SE.FullModel))]
SE.FullModelSNPS <- SE.FullModelSNPS[,-1]
maxLevel <- log2(0.368*nrow(envGF.SE.FullModel)/2)
maxLevel

gf.SE.FullModel.SNPs <- gradientForest(cbind(envGF.SE.FullModel, SE.FullModelSNPS), predictor.vars=colnames(envGF.SE.FullModel), response.vars=colnames(SE.FullModelSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##50+ warnings (loci variable in less than 5 populations)


```


## Figure 1: R2 importance

Which of the variables is the most important across all the gradients. 

I'm comparing only dist, season and temp (ie. excluding MEM. And I'm rescaling R2 between 0 and 1 


```
R.sq.CHN.FullModel <- (rowMeans(gf.CHN.FullModel.SNPs$imp.rsq, na.rm=T))
R.sq.CHN.FullModel <- as.data.frame(R.sq.CHN.FullModel)
colnames(R.sq.CHN.FullModel) <- "CHN"

R.sq.CHS.FullModel <- (rowMeans(gf.CHS.FullModel.SNPs$imp.rsq, na.rm=T))
R.sq.CHS.FullModel <- as.data.frame(R.sq.CHS.FullModel)
colnames(R.sq.CHS.FullModel) <- "CHS"

R.sq.CHS.VS.FullModel <- rowMeans(gf.CHS.VS.FullModel.SNPs$imp.rsq, na.rm=T)
R.sq.CHS.VS.FullModel <- as.data.frame(R.sq.CHS.VS.FullModel)
colnames(R.sq.CHS.VS.FullModel) <- "CHS.VS"

R.sq.CHS.TI.FullModel <- rowMeans(gf.CHS.TI.FullModel.SNPs$imp.rsq, na.rm=T)
R.sq.CHS.TI.FullModel <- as.data.frame(R.sq.CHS.TI.FullModel)
colnames(R.sq.CHS.TI.FullModel) <- "CHS.TI"

R.sq.SE.FullModel <- rowMeans(gf.SE.FullModel.SNPs$imp.rsq, na.rm=T)
R.sq.SE.FullModel <- as.data.frame(R.sq.SE.FullModel)
colnames(R.sq.SE.FullModel) <- "SE"

R.sq.CZ.FullModel <- rowMeans(gf.CZ.FullModel.SNPs$imp.rsq, na.rm=T)
R.sq.CZ.FullModel <- as.data.frame(R.sq.CZ.FullModel)
colnames(R.sq.CZ.FullModel) <- "CZ"


R.sq.ALL <- R.sq.CHN.FullModel[1:5,]
R.sq.ALL <- as.data.frame(R.sq.ALL)
colnames(R.sq.ALL) <- "CHN"
row.names(R.sq.ALL) <- c("dist", "mean.temp.60d", "days.above.6", "MEM1", "MEM2")
R.sq.ALL$CHS <- R.sq.CHS.FullModel[1:5,]
R.sq.ALL$CHS.VS <- R.sq.CHS.VS.FullModel[1:5,]
R.sq.ALL$CHS.TI <- R.sq.CHS.TI.FullModel[1:5,]
R.sq.ALL$CZ <- R.sq.CZ.FullModel[1:5,]
R.sq.ALL$SE <- R.sq.SE.FullModel[1:5,]


library(reshape)
library(ggplot2)

R.sq.ALL <- as.matrix(R.sq.ALL)

R.sq.ALL.scale <- R.sq.ALL.scale[1:3,]   ##select only the variables of interest
R.sq.ALL.scale <- apply(R.sq.ALL.scale, MARGIN=2, FUN=function(X) (X - min(X))/diff(range(X)))   ##rescale between 0 and 1
R.sq.ALL.scale.melt <- melt(R.sq.ALL.scale)
R.sq.ALL.melt
R.sq.ALL.scale.melt <- R.sq.ALL.scale.melt[order(R.sq.ALL.scale.melt$X1),]
R.sq.ALL.scale.melt

pdf("GF.ALL.Figure1.R2_20171013.pdf")
ggplot(R.sq.ALL.scale.melt, aes(x=X2, y=X1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()
```




