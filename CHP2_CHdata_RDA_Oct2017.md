# RDA analysis of CH dataset

RDA: Effect of Climate & Geog on genetic variation

### Part1: RDA of full dataset for each region
      
      aim: partition variance into neutral, adaptive and residual

### Part 2: Effect of spatial scale

Compare the proportion of variance explained at different spatial scales


### Part 3: Identify the most important environmental variables

Of the 5 environmental variables, which is the most important in the different geographic regions, and globally. 


### Final: RDA of adaptive loci for each region  
 
 This was something I did before - a potential additional analysis with adaptive loci.. 
 
 
 

I can use RDA to partition the effects of climate and geog on genetic variation. 

Explanation and example: http://www.davidzeleny.net/anadat-r/doku.php/en:rda_examples

And Victoria's paper: http://www.amjbot.org/content/103/1/33.full.pdf+html

Paul Gugger tutorial: http://pgugger.al.umces.edu/assets/redundancy-analysis-for-landscape-genetics.pdf


Also, see these papers: 

http://onlinelibrary.wiley.com/doi/10.1111/mec.13409/full

http://onlinelibrary.wiley.com/doi/10.1111/mec.13811/full

http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2012.05709.x/full


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA


## Part 1

### Data:

I'm using the full CH dataset. I have to run the analysis seperately for all the different partitions of the data to investigate adaptation
at different geographic scales. 

#### 1. CHall.n932.loci7744

#### 2. CHN.n229.loci5265

#### 3. CHS.n275.loci6339

#### 4. CHS.VS.n135.loci5835

#### 5. CHS.TI.n140.loci5692

#### 6. CZ.n404.loci7288


### Input files:

Plink files need to be generated on the server. The data are too large to subset on my Mac. Gives the fault: 

```
After filtering, kept 940 out of 940 Individuals
Writing PLINK PED and MAP files ... 
Error:  Could not open temporary file.
```

####1. MAF of all loci

CHall
```
#Calculate MAF for the full dataset within region using PLINK

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets

###Use the *nosex file to create a file for subsetting the data within pop

plink --file CHall.932.7744.plink --freq --within CHall.932.cluster.pop --out CH932.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CH932.MAF <- read.table("CH932.pop.frq.strat", header=T)
CH932.MAF2 <- CH932.MAF[,c(3,2,6)]
summary(CH932.MAF2)
      CLST                SNP              MAF        
 abnd   :  7744   1000065:55:    81   Min.   :0.0000  
 agra   :  7744   1000448:60:    81   1st Qu.:0.0000  
 alpl   :  7744   1000622:42:    81   Median :0.1429  
 apla   :  7744   1000723:85:    81   Mean   :0.2310  
 arce   :  7744   1000806:91:    81   3rd Qu.:0.3750  
 bach   :  7744   1000841:7 :    81   Max.   :1.0000  
 (Other):580800   (Other)   :626778            
        

library("ggplot2")
library("reshape2")

CH932.MAF3 <- melt(CH932.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CH932.MAF3)
head(CH932.MAF3)


CH932.MAF4 <- dcast(CH932.MAF3, formula= CLST ~ SNP)
head(CH932.MAF4)


##Add X infront of all locusnames. 
colnames(CH932.MAF4) <- paste("X", colnames(CH932.MAF4), sep=".")
write.csv(CH932.MAF4, file="CH932.MAF.csv")
```


CHN
```
plink --file CHN.229.5265.plink --freq --within CHN229.cluster.pop --out CHN229.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHN229.MAF <- read.table("CHN229.pop.frq.strat", header=T)
CHN229.MAF2 <- CHN229.MAF[,c(3,2,6)]
summary(CHN229.MAF2)
     
library("ggplot2")
library("reshape2")

CHN229.MAF3 <- melt(CHN229.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN229.MAF3)
head(CHN229.MAF3)


CHN229.MAF4 <- dcast(CHN229.MAF3, formula= CLST ~ SNP)
head(CHN229.MAF4)


##Add X infront of all locusnames. 
colnames(CHN229.MAF4) <- paste("X", colnames(CHN229.MAF4), sep=".")
write.csv(CHN229.MAF4, file="CHN229.MAF.csv")

```

CHS
```
plink --file CHS.275.6339.plink --freq --within CHS283.cluster.pop --out CHS283.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS275.MAF <- read.table("CHS275.pop.frq.strat", header=T)
CHS275.MAF2 <- CHS275.MAF[,c(3,2,6)]
summary(CHS275.MAF2)
     
library("ggplot2")
library("reshape2")

CHS275.MAF3 <- melt(CHS275.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS275.MAF3)
head(CHS275.MAF3)


CHS275.MAF4 <- dcast(CHS275.MAF3, formula= CLST ~ SNP)
head(CHS275.MAF4)


##Add X infront of all locusnames. 
colnames(CHS275.MAF4) <- paste("X", colnames(CHS275.MAF4), sep=".")
write.csv(CHS275.MAF4, file="CHS275.MAF.csv")

```


CHS.VS
```
plink --file CHS.VS.135.5835.plink --freq --within CHS.VS.135.cluster.pop --out CHS.VS135.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS.VS135.MAF <- read.table("CHS.VS135.pop.frq.strat", header=T)
CHS.VS135.MAF2 <- CHS.VS135.MAF[,c(3,2,6)]
summary(CHS.VS135.MAF2)
     
library("ggplot2")
library("reshape2")

CHS.VS135.MAF3 <- melt(CHS.VS135.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS135.MAF3)
head(CHS.VS135.MAF3)


CHS.VS135.MAF4 <- dcast(CHS.VS135.MAF3, formula= CLST ~ SNP)
head(CHS.VS135.MAF4)


##Add X infront of all locusnames. 
colnames(CHS.VS135.MAF4) <- paste("X", colnames(CHS.VS135.MAF4), sep=".")
write.csv(CHS.VS135.MAF4, file="CHS.VS135.MAF.csv")

```


CHS.TI
```
plink --file CHS.TI.140.5692.plink --freq --within CHS.TI.140.cluster.pop --out CHS.TI140.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS.TI140.MAF <- read.table("CHS.TI140.pop.frq.strat", header=T)
CHS.TI140.MAF2 <- CHS.TI140.MAF[,c(3,2,6)]
summary(CHS.TI140.MAF2)
     
library("ggplot2")
library("reshape2")

CHS.TI140.MAF3 <- melt(CHS.TI140.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI140.MAF3)
head(CHS.TI140.MAF3)


CHS.TI140.MAF4 <- dcast(CHS.TI140.MAF3, formula= CLST ~ SNP)
head(CHS.TI140.MAF4)


##Add X infront of all locusnames. 
colnames(CHS.TI140.MAF4) <- paste("X", colnames(CHS.TI140.MAF4), sep=".")
write.csv(CHS.TI140.MAF4, file="CHS.TI140.MAF.csv")

```


CZ  
```
plink --file CZ.404.7288.plink --freq --within CZ404.cluster.pop --out CZ404.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CZ404.MAF <- read.table("CZ404.pop.frq.strat", header=T)
CZ404.MAF2 <- CZ404.MAF[,c(3,2,6)]
summary(CZ404.MAF2)
     
library("ggplot2")
library("reshape2")

CZ404.MAF3 <- melt(CZ404.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ404.MAF3)
head(CZ404.MAF3)


CZ404.MAF4 <- dcast(CZ404.MAF3, formula= CLST ~ SNP)
head(CZ404.MAF4)


##Add X infront of all locusnames. 
colnames(CZ404.MAF4) <- paste("X", colnames(CZ404.MAF4), sep=".")
write.csv(CZ404.MAF4, file="CZ404.MAF.csv")

```




#### 2. Geographic component: pcnm 

I used only geographic coordinates before. Now I've decided to use a PCNM which should describe geographic structure (akin to MEM)


I've already calculated MEMs for the gradient forest for all the transects: 

https://github.com/alexjvr1/Manuscripts/blob/90bb82edd50be25e4d467ca64b151798926bee9a/CHP2_GF.subsets.MACfiltered.md

So I'm changing the end of the code to write pcnms. All the code is given here, but it was run here: 

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017

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
 
library(vegan)
CHN.pcnm.scores <- pcnm(CHN.dxy, CHN.th)  ##calculate pcnm

test.scores(CHN.pcnm.scores, CHN.listw, nsim=999)  ##and select only positive pcnm

          stat  pval
1   0.84559173 0.001
2   0.90374681 0.001
3   0.68531162 0.001
4   0.43768794 0.005
5  -0.47001032 0.960
6  -0.07610158 0.560
7  -0.16121905 0.724
8  -0.33847617 0.885
9  -0.29365084 0.851
10 -0.14700643 0.691
11 -0.53066144 0.985
12 -0.77952568 0.999

#4 significant PCNM eigenfunctions with positive correlations.

write.table (CHN.pcnm.scores$vectors[,1], "CHN.scores_PCNM1.txt") 
write.table (CHN.pcnm.scores$vectors[,2], "CHN.scores_PCNM2.txt") 
write.table (CHN.pcnm.scores$vectors[,3], "CHN.scores_PCNM3.txt") 
write.table (CHN.pcnm.scores$vectors[,4], "CHN.scores_PCNM4.txt") 
```


#### CZ
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CZ <- read.csv("AllEnv.Data_CZ_20161010.csv", header=T)
head(env.data.CZ)

env.data.CZ <- env.data.CZ[-18,]##remove meie. This population was added by mistake

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
Number of regions: 37 
Number of nonzero links: 72 
Percentage nonzero weights: 5.259313 
Average number of links: 1.945946 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CZ.listw=nb2listw(CZ.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CZ.listw)
#[1] TRUE
 
library(vegan)
CZ.pcnm.scores <- pcnm(CZ.dxy, CZ.th)  ##calculate pcnm

test.scores(CZ.pcnm.scores, CZ.listw, nsim=999)  ##and select only positive pcnm

          stat  pval
1   0.95580479 0.001*
2   0.93050370 0.001*
3   0.91257634 0.001*
4   0.64966475 0.001*
5   0.25237401 0.041*
6   0.34391370 0.011*
7  -0.14166137 0.803
8   0.08943395 0.248
9  -0.11566399 0.709
10  0.11102298 0.218
11 -0.15455107 0.769
12 -0.32732604 0.977
13  0.16821893 0.104
14 -0.03065197 0.508
15 -0.63822650 0.999
16  0.34012456 0.005*
17 -0.52370184 0.998
18 -0.60240388 0.999
19 -0.67917757 0.999
20 -0.53936405 0.998
21 -0.50854789 0.997

#7 significant PCNM eigenfunctions with positive correlations.

write.table (CZ.pcnm.scores$vectors[,1], "CZ.scores_PCNM1.txt") 
write.table (CZ.pcnm.scores$vectors[,2], "CZ.scores_PCNM2.txt") 
write.table (CZ.pcnm.scores$vectors[,3], "CZ.scores_PCNM3.txt") 
write.table (CZ.pcnm.scores$vectors[,4], "CZ.scores_PCNM4.txt") 
write.table (CZ.pcnm.scores$vectors[,5], "CZ.scores_PCNM5.txt") 
write.table (CZ.pcnm.scores$vectors[,6], "CZ.scores_PCNM6.txt") 
write.table (CZ.pcnm.scores$vectors[,16], "CZ.scores_PCNM7.txt") 

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
 
library(vegan)
CHS.pcnm.scores <- pcnm(CHS.dxy, CHS.th)  ##calculate pcnm

test.scores(CHS.pcnm.scores, CHS.listw, nsim=999)  ##and select only positive pcnm

         stat  pval
1   0.9658591 0.001
2   0.9300624 0.001
3   0.4679170 0.001
4  -0.3328825 0.913
5  -0.0461216 0.587
6   0.4680426 0.006
7   0.3722407 0.024
8  -0.1211875 0.661
9  -0.2357535 0.840
10 -0.3741807 0.946

#3 significant PCNM eigenfunctions with positive correlations.

write.table (CHS.pcnm.scores$vectors[,1], "CHS.scores_PCNM1.txt") 
write.table (CHS.pcnm.scores$vectors[,2], "CHS.scores_PCNM2.txt") 
write.table (CHS.pcnm.scores$vectors[,3], "CHS.scores_PCNM3.txt") 
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
 
library(vegan)
CHS.TI.pcnm.scores <- pcnm(CHS.TI.dxy, CHS.TI.th)  ##calculate pcnm

test.scores(CHS.TI.pcnm.scores, CHS.TI.listw, nsim=999)  ##and select only positive pcnm

         stat  pval
1   0.92968201 0.001
2   0.80836905 0.001
3   0.73502128 0.001
4   0.41725039 0.005
5   0.27276369 0.090
6  -0.06803219 0.514
7  -0.27148032 0.771
8  -0.27717479 0.798
9  -0.74824234 0.987
10 -0.40149812 0.915

#4 significant PCNM eigenfunctions with positive correlations.

write.table (CHS.TI.pcnm.scores$vectors[,1], "CHS.TI.scores_PCNM1.txt") 
write.table (CHS.TI.pcnm.scores$vectors[,2], "CHS.TI.scores_PCNM2.txt") 
write.table (CHS.TI.pcnm.scores$vectors[,3], "CHS.TI.scores_PCNM3.txt")  
write.table (CHS.TI.pcnm.scores$vectors[,4], "CHS.TI.scores_PCNM4.txt")  
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
 
library(vegan)
CHS.VS.pcnm.scores <- pcnm(CHS.VS.dxy, CHS.VS.th)  ##calculate pcnm

test.scores(CHS.VS.pcnm.scores, CHS.VS.listw, nsim=999)  ##and select only positive pcnm

         stat  pval
1 0.395931977 0.037
2 0.037207343 0.038
3 0.009497033 0.300
4 0.468849442 0.019

#3 significant PCNM eigenfunctions with positive correlations.

write.table (CHS.VS.pcnm.scores$vectors[,1], "CHS.VS.scores_PCNM1.txt") 
write.table (CHS.VS.pcnm.scores$vectors[,2], "CHS.VS.scores_PCNM2.txt") 
write.table (CHS.VS.pcnm.scores$vectors[,4], "CHS.VS.scores_PCNM4.txt") 
```


#### CHall
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHall <- read.csv("Env.CHall.data.csv", header=T)
head(env.data.CHall)
#extract x and y
CHall.xy <- env.data.CHall[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
CHall.dxy <- distm(CHall.xy)
CHall.dxy <- as.dist(CHall.dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
CHall.th <- give.thresh(CHall.dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
CHall.nb1 <- mst.nb(CHall.dxy)
CHall.wh1 <- which(as.matrix(CHall.dxy)==CHall.th,arr.ind=TRUE)
plot(CHall.nb1,CHall.xy,pch=20,cex=2,lty=3)
lines(CHall.xy[CHall.wh1[1,],1],CHall.xy[CHall.wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
CHall.th 
#[1] 37752.95
CHall.nb1
Neighbour list object:
Number of regions: 81 
Number of nonzero links: 160 
Percentage nonzero weights: 2.438653 
Average number of links: 1.975309 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
CHall.listw=nb2listw(CHall.nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(CHall.listw)
#[1] TRUE
 
library(vegan)
CHall.pcnm.scores <- pcnm(CHall.dxy, CHall.th)  ##calculate pcnm

test.scores(CHall.pcnm.scores, CHall.listw, nsim=999)  ##and select only positive pcnm

           stat  pval
1   0.974630699 0.001*
2   0.983738764 0.001*
3   0.949606366 0.001*
4   0.935105256 0.001*
5   0.906896859 0.001*
6   0.894147516 0.001*
7   0.846225903 0.001*
8   0.891306648 0.001*
9   0.882729874 0.001*
10  0.533860219 0.001*
11  0.538148576 0.001*
12  0.220223091 0.022*
13  0.423930713 0.001*
14  0.265709642 0.008*
15 -0.190681120 0.936
16  0.045389272 0.310
17  0.119311054 0.110
18 -0.127488357 0.846
19 -0.173360967 0.938
20  0.038654599 0.336
21 -0.203487906 0.944
22 -0.158747679 0.897
23 -0.152646178 0.897
24 -0.266431324 0.998
25 -0.233618576 0.971
26 -0.420386320 0.999
27 -0.577584506 0.999
28 -0.123987718 0.837
29 -0.002309396 0.475
30 -0.234212862 0.965
31  0.058511897 0.245
32 -0.086798543 0.766
33 -0.321461715 0.995
34 -0.014731474 0.548
35 -0.059465574 0.949
36  0.096276663 0.056
37  0.026761631 0.225
38  0.034995430 0.227
39  0.218481211 0.025*
40 -0.495926269 0.982
41 -0.790173630 0.999
42 -0.771856198 0.999
43 -0.115110901 0.835
44 -0.377614344 0.988
45 -0.230270723 0.971
46 -0.414777341 0.979
47 -0.489151801 0.999

#4 significant PCNM eigenfunctions with positive correlations.

write.table (CHall.pcnm.scores$vectors[,1], "CHall.scores_PCNM1.txt") 
write.table (CHall.pcnm.scores$vectors[,2], "CHall.scores_PCNM2.txt") 
write.table (CHall.pcnm.scores$vectors[,3], "CHall.scores_PCNM3.txt") 
write.table (CHall.pcnm.scores$vectors[,4], "CHall.scores_PCNM4.txt") 
write.table (CHall.pcnm.scores$vectors[,5], "CHall.scores_PCNM5.txt") 
write.table (CHall.pcnm.scores$vectors[,6], "CHall.scores_PCNM6.txt") 
write.table (CHall.pcnm.scores$vectors[,7], "CHall.scores_PCNM7.txt") 
write.table (CHall.pcnm.scores$vectors[,8], "CHall.scores_PCNM8.txt") 
write.table (CHall.pcnm.scores$vectors[,9], "CHall.scores_PCNM9.txt") 
write.table (CHall.pcnm.scores$vectors[,10], "CHall.scores_PCNM10.txt") 
write.table (CHall.pcnm.scores$vectors[,11], "CHall.scores_PCNM11.txt") 
write.table (CHall.pcnm.scores$vectors[,12], "CHall.scores_PCNM12.txt") 
write.table (CHall.pcnm.scores$vectors[,13], "CHall.scores_PCNM13.txt") 
write.table (CHall.pcnm.scores$vectors[,14], "CHall.scores_PCNM14.txt") 
write.table (CHall.pcnm.scores$vectors[,39], "CHall.scores_PCNM15.txt") 
```


###OLD way: Geographic coordinates
paste coordinates into the MAF.csv file from Env.Data.all_20161025.csv

located: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2


i.e. lat, long




#### 3. Climate variables

```
1. solar.rad.60d (amt of solar radiation 60 days after egg laying)

2. pcpt.60d (precipitation in the 60days after egg laying)

3. day10cm (calender day on which only 10cm of snow remains)

4. temp.laying.date

5. shadow.days
```

paste the 5 variables from /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/Env.Data.all_20161025.csv




### Run RDA

#### 1. CHall

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHallData <- read.csv("CH932.MAF.csv", header=T)
CHallData <- CHallData[complete.cases(CHallData),]
GenData <- CHallData[,25:ncol(CHallData)]

GenData <- decostand(GenData, "hellinger") #detrend the genetic data 
Climate.Data <- CHallData[,3:24]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CH932full <- rda(GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 +PCNM10 + PCNM11 + PCNM12 +PCNM13 + PCNM14 +PCNM15 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s.

RDA.CH932full

RDA.CH932full
Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 +
PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 +
PCNM14 + PCNM15 + shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total          0.3387     1.0000     
Constrained    0.1941     0.5731   20
Unconstrained  0.1446     0.4269   60
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7    RDA8    RDA9   RDA10 
0.10352 0.02716 0.01238 0.00688 0.00502 0.00429 0.00371 0.00333 0.00312 0.00297 
  RDA11   RDA12   RDA13   RDA14   RDA15   RDA16   RDA17   RDA18   RDA19   RDA20 
0.00285 0.00262 0.00248 0.00236 0.00223 0.00215 0.00184 0.00180 0.00178 0.00162 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.020688 0.006219 0.005449 0.004814 0.004716 0.004201 0.003947 0.003725 
(Showed only 8 of all 60 unconstrained eigenvalues)


##calculate proportion explained (R2)

RsquareAdj(RDA.CH932full)

$r.squared
[1] 0.5731014

$adj.r.squared
[1] 0.4308019



anova(RDA.CH932full)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 + PCNM14 + PCNM15 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model    20  0.19410 4.0274  0.001 ***
Residual 60  0.14458                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CH932full.pdf")
plot(RDA.CH932full)
dev.off()
```


Best model overall
```
library(vegan)

RDA.CH932full.optModel <- ordistep(RDA.CH932full)


Step: GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM8 +      PCNM9 + PCNM10 + sol.rad.60d + pcpt.60d 

              Df    AIC        F Pr(>F)   
- PCNM10       1 456.57   2.6503  0.030 * 
- PCNM6        1 456.60   2.6771  0.025 * 
- PCNM3        1 456.70   2.7652  0.020 * 
- sol.rad.60d  1 457.45   3.4292  0.015 * 
- pcpt.60d     1 457.60   3.5644  0.005 **
- PCNM9        1 457.97   3.8999  0.005 **
- PCNM2        1 466.65  12.1461  0.005 **
- PCNM8        1 466.78  12.2729  0.005 **
- PCNM4        1 474.54  20.4435  0.005 **
- PCNM1        1 521.64  90.9973  0.005 **
- PCNM5        1 548.88 154.9492  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 +
PCNM6 + PCNM8 + PCNM9 + PCNM10 + sol.rad.60d + pcpt.60d, data =
Climate.Data)

               Inertia Proportion Rank
Total         430.9929     1.0000     
Constrained   222.5162     0.5163   11
Unconstrained 208.4767     0.4837   69
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8   RDA9  RDA10  RDA11 
137.92  30.70  16.39   9.43   5.99   4.75   4.41   3.75   3.31   3.14   2.72 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
30.806  7.662  6.838  6.675  6.126  5.736  5.162  4.763 
(Showed only 8 of all 69 unconstrained eigenvalues)

```



Partial out Geog and test which variable is the most important
```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CH932.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 +PCNM13 +PCNM14 + PCNM15), Climate.Data)

#head(summary(pRDA.CH940.geog))

pRDA.CH932.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 +
PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 +
PCNM13 + PCNM14 + PCNM15), data = Climate.Data)

              Inertia Proportion Rank
Total         0.33869    1.00000     
Conditional   0.17749    0.52405   15
Constrained   0.01661    0.04905    5
Unconstrained 0.14458    0.42690   60
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5 
0.007017 0.003147 0.002428 0.002020 0.002002 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.020688 0.006219 0.005449 0.004814 0.004716 0.004201 0.003947 0.003725 
(Showed only 8 of all 60 unconstrained eigenvalues)




RsquareAdj(pRDA.CH932.geog)

$r.squared
[1] 0.04905403

$adj.r.squared
[1] 0.01658972

anova(pRDA.CH932.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 + PCNM14 + PCNM15), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     5 0.016614 1.3789  0.006 **
Residual 60 0.144585                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CH940 GeneticData


pdf(file="pRDA.CH932.geog.pdf")
plot(pRDA.CH932.geog, main="pRDA CH932 (geog partialled out)")
dev.off()
```



Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CH932.geog)


Biplot scores for constraining variables

                     RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days      -0.45832 -0.1355 -0.37446  0.3868  0.50612   0
sol.rad.60d       0.34831  0.7203  0.25210 -0.1824  0.30492   0
pcpt.60d         -0.27648  0.6932  0.07362 -0.0963  0.01866   0
day10cm          -0.14166  0.5715  0.14609  0.4306  0.13571   0
temp.laying.date  0.06186  0.3177 -0.70016  0.1267 -0.23785   0


#eigenvalue (i.e. weight) of each RDA
pRDA.CH932.geog$CCA$eig

       RDA1        RDA2        RDA3        RDA4        RDA5 
0.007017060 0.003146749 0.002427694 0.002020352 0.002002091 
```



Partial out climate
```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CH932.climate <- rda(GenData~PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 + PCNM14 + PCNM15 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CH932.climate

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 +
PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 +
PCNM14 + PCNM15 + Condition(shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date), data = Climate.Data)

              Inertia Proportion Rank
Total         0.33869    1.00000     
Conditional   0.02951    0.08713    5
Constrained   0.16459    0.48597   15
Unconstrained 0.14458    0.42690   60
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7    RDA8    RDA9   RDA10 
0.09564 0.02335 0.01054 0.00486 0.00442 0.00350 0.00310 0.00293 0.00287 0.00271 
  RDA11   RDA12   RDA13   RDA14   RDA15 
0.00251 0.00230 0.00213 0.00202 0.00172 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.020688 0.006219 0.005449 0.004814 0.004716 0.004201 0.003947 0.003725 
(Showed only 8 of all 60 unconstrained eigenvalues)



#head(summary(pRDA.CH932.climate))

RsquareAdj(pRDA.CH932.climate)

$r.squared
[1] 0.4859667

$adj.r.squared
[1] 0.4045249

anova(pRDA.CH932.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 + PCNM14 + PCNM15 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model    15  0.16459 4.5535  0.001 ***
Residual 60  0.14458                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pdf(file="pRDA.CH932.climate.pdf")
plot(pRDA.CH932.climate, main="pRDA CH932 (climate partialled out)")
dev.off()
```

Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CHall.Climate.contribution <- (20.4271+18.04)  ##get these values from the RDA table. From the full and partial models run above
Venn.CHall.Geog.contribution <- (213.90415+18.04)
Venn.CHall.Geog.Climate.contribution <- 18.04

pdf("CHall.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CHall.Geog.contribution, area2=Venn.CHall.Climate.contribution, cross.area=Venn.CHall.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```


Find the most important loci associated with Climate
```
 summary(pRDA.CH932.geog)

Biplot scores for constraining variables

                     RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days      -0.45832 -0.1355 -0.37446  0.3868  0.50612   0
sol.rad.60d       0.34831  0.7203  0.25210 -0.1824  0.30492   0
pcpt.60d         -0.27648  0.6932  0.07362 -0.0963  0.01866   0
day10cm          -0.14166  0.5715  0.14609  0.4306  0.13571   0
temp.laying.date  0.06186  0.3177 -0.70016  0.1267 -0.23785   0


RDA1: pcpt.60d, day10cm, shadow.days

RDA2: day10cm, temp.laying.date
```


Best model when geography is partialled out
```
step.res.CHall.pgeog <- ordistep(pRDA.CH932.geog)

Step: GenData ~ sol.rad.60d + pcpt.60d + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 +
PCNM8 + PCNM9 + PCNM10 + PCNM11 + PCNM12 + PCNM13 + PCNM14 + PCNM15) 

                                                                                                                                           
                      Df    AIC    F      Pr(>F)
- sol.rad.60d         1   -116-690 2.2363 0.005**
- pcpt.60d            1   -116.611 2.3001 0.005**
- Condition(PCNM1-15) 15  -85.282
```

Best model when climate is partialled out
```
step.res.CHall.pclimate <- ordistep(pRDA.CH932.climate)
```


#### 2. CHN

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHNData <- read.csv("CHN229.MAF.csv", header=T)
GenData <- CHNData[,14:ncol(CHNData)]
GenData <- decostand(GenData, "hellinger") ##detrend the data

Climate.Data <- CHNData[,3:13]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHNfull <- rda(GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHNfull

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +
shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data
= Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Constrained   0.13205    0.61459    9
Unconstrained 0.08281    0.38541    9
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5     RDA6     RDA7     RDA8 
0.030828 0.024473 0.018272 0.014203 0.011790 0.010550 0.007933 0.007175 
    RDA9 
0.006824 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.014934 0.011844 0.010747 0.010134 0.008548 0.008055 0.007047 0.006180 
     PC9 
0.005317 
 


RsquareAdj(RDA.CHNfull)

$r.squared
[1] 0.6145933

$adj.r.squared
[1] 0.2291866



anova(RDA.CHNfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     9 0.132047 1.5947  0.001 ***
Residual  9 0.082806                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHNfull.pdf")
plot(RDA.CHNfull, main="RDA CHNfull")
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHNgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4), Climate.Data)

pRDA.CHNgeog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4),
data = Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Conditional   0.07850    0.36537    4
Constrained   0.05355    0.24922    5
Unconstrained 0.08281    0.38541    9
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5 
0.014513 0.011845 0.009541 0.009002 0.008644 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.014934 0.011844 0.010747 0.010134 0.008548 0.008055 0.007047 0.006180 
     PC9 
0.005317  

RsquareAdj(pRDA.CHNgeog)

$r.squared
[1] 0.2492243

$adj.r.squared
[1] 0.04514076


anova(pRDA.CHNgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4), data = Climate.Data)
         Df Variance     F Pr(>F)
Model     5 0.053547 1.164  0.129
Residual  9 0.082806


H0 not rejected: Climate does not explain a significant amount of CHN GeneticData





pdf(file="pRDA.CHNgeog.pdf")
plot(pRDA.CHNgeog, main="pRDA CHN(geog partialled out)")
dev.off()
```


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHNclimate <- rda(GenData~PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHNclimate

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +
Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date), data = Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Conditional   0.07712    0.35892    5
Constrained   0.05493    0.25567    4
Unconstrained 0.08281    0.38541    9
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4 
0.020120 0.014598 0.011917 0.008297 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.014934 0.011844 0.010747 0.010134 0.008548 0.008055 0.007047 0.006180 
     PC9 
0.005317  



RsquareAdj(pRDA.CHNclimate)

$r.squared
[1] 0.2556717

$adj.r.squared
[1] 0.1168337

anova(pRDA.CHNclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     4 0.054932 1.4926  0.004 **
Residual  9 0.082806                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: not rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHNclimate.pdf")
plot(pRDA.CHNclimate, main="pRDA (climate partialled out)")
dev.off()
```



Find the most important loci associated with Climate


```
summary(pRDA.CHNgeog)

Biplot scores for constraining variables

                    RDA1     RDA2    RDA3     RDA4     RDA5 PC1
shadow.days      -0.2301 -0.41578  0.1777 -0.01312  0.47282   0
sol.rad.60d      -0.6434 -0.23237 -0.4678  0.05344 -0.19547   0
pcpt.60d         -0.7412 -0.09528 -0.2354  0.36022 -0.05286   0
day10cm          -0.5316 -0.36279  0.1897  0.15015 -0.14827   0
temp.laying.date -0.6775  0.27792  0.2336 -0.48509 -0.01842   0


RDA1: temp.laying.date, shadow.days

RDA2: sol.rad.60d, pcpt.60d
```


Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CHNgeog)


Biplot scores for constraining variables

                    RDA1     RDA2    RDA3     RDA4     RDA5 PC1
shadow.days      -0.2301 -0.41578  0.1777 -0.01312  0.47282   0
sol.rad.60d      -0.6434 -0.23237 -0.4678  0.05344 -0.19547   0
pcpt.60d         -0.7412 -0.09528 -0.2354  0.36022 -0.05286   0
day10cm          -0.5316 -0.36279  0.1897  0.15015 -0.14827   0
temp.laying.date -0.6775  0.27792  0.2336 -0.48509 -0.01842   0


#eigenvalue (i.e. weight) of each RDA
pRDA.CHNgeog$CCA$eig

     RDA1        RDA2        RDA3        RDA4        RDA5 
0.014513463 0.011845253 0.009541467 0.009002205 0.008644273 
```



Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CHN.Climate.contribution <- (48.8465+19.42)  ##get these values from the RDA table. From the full and partial models run above
Venn.CHN.Geog.contribution <- (50.1661+19.42)
Venn.CHN.Geog.Climate.contribution <- 19.42

pdf("CHN.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CHN.Geog.contribution, area2=Venn.CHN.Climate.contribution, cross.area=Venn.CHN.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```



Best model overall


```
cd /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA/Sept2017/AnalysisFullDataset/Nov17.findOptimalModels

##R

RDA.CHNfull.optModel <- ordistep(RDA.CHNfull)

Start: GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + shadow.days + sol.rad.60d +      pcpt.60d + day10cm + temp.laying.date 

                   Df     AIC      F Pr(>F)   
- PCNM2             1 -28.424 0.9659  0.520   
- shadow.days       1 -28.365 0.9971  0.430   
- day10cm           1 -28.257 1.0540  0.360   
- sol.rad.60d       1 -28.223 1.0722  0.320   
- pcpt.60d          1 -28.147 1.1126  0.310   
- PCNM4             1 -27.930 1.2285  0.200   
- temp.laying.date  1 -27.979 1.2021  0.190   
- PCNM1             1 -27.197 1.6309  0.015 * 
- PCNM3             1 -26.860 1.8213  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: GenData ~ PCNM1 + PCNM3 + PCNM4 + shadow.days + sol.rad.60d +      pcpt.60d + day10cm + temp.laying.date 

                   Df     AIC      F Pr(>F)   
- shadow.days       1 -28.632 0.9889  0.515   
- pcpt.60d          1 -28.422 1.1113  0.305   
- temp.laying.date  1 -28.305 1.1801  0.235   
- sol.rad.60d       1 -28.270 1.2009  0.135   
- PCNM4             1 -28.090 1.3068  0.105   
- day10cm           1 -28.066 1.3213  0.105   
- PCNM3             1 -27.136 1.8893  0.010 **
- PCNM1             1 -27.477 1.6780  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: GenData ~ PCNM1 + PCNM3 + PCNM4 + sol.rad.60d + pcpt.60d + day10cm +      temp.laying.date 

                   Df     AIC      F Pr(>F)   
- pcpt.60d          1 -28.786 1.1224  0.230   
- sol.rad.60d       1 -28.672 1.1955  0.195   
- temp.laying.date  1 -28.691 1.1832  0.185   
- day10cm           1 -28.258 1.4642  0.030 * 
- PCNM4             1 -28.037 1.6101  0.020 * 
- PCNM1             1 -27.878 1.7158  0.015 * 
- PCNM3             1 -27.244 2.1476  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: GenData ~ PCNM1 + PCNM3 + PCNM4 + sol.rad.60d + day10cm + temp.laying.date 

                   Df     AIC      F Pr(>F)   
- temp.laying.date  1 -29.009 1.1769  0.215   
- sol.rad.60d       1 -28.842 1.2932  0.095 . 
- PCNM4             1 -28.151 1.7853  0.005 **
- PCNM1             1 -28.117 1.8100  0.005 **
- PCNM3             1 -27.669 2.1400  0.005 **
- day10cm           1 -27.267 2.4416  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: GenData ~ PCNM1 + PCNM3 + PCNM4 + sol.rad.60d + day10cm 

              Df     AIC      F Pr(>F)   
- sol.rad.60d  1 -29.199 1.2994  0.105   
- PCNM4        1 -28.590 1.7647  0.010 **
- PCNM1        1 -28.393 1.9184  0.005 **
- PCNM3        1 -28.183 2.0849  0.005 **
- day10cm      1 -27.275 2.8233  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Step: GenData ~ PCNM1 + PCNM3 + PCNM4 + day10cm 

          Df     AIC      F Pr(>F)   
- PCNM4    1 -29.042 1.6827  0.005 **
- PCNM3    1 -28.640 2.0178  0.005 **
- PCNM1    1 -28.500 2.1363  0.005 **
- day10cm  1 -27.859 2.6902  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> RDA.CHNfull.optModel
Call: rda(formula = GenData ~ PCNM1 + PCNM3 + PCNM4 + day10cm, data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Constrained   0.08073    0.37576    4
Unconstrained 0.13412    0.62424   14
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4 
0.029041 0.022018 0.017873 0.011802 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.016785 0.014188 0.013689 0.012529 0.010099 0.009617 0.009006 0.008232 
     PC9     PC10     PC11     PC12     PC13     PC14 
0.008164 0.007935 0.006904 0.006291 0.005855 0.004825 

```




Best model when space is partialled out
```

pRDA.CHNgeog <- rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4),
data = Climate.Data) 

pRDA.CHNgeog.optModel <- ordistep(pRDA.CHNgeog)


Step: GenData ~ pcpt.60d + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4) 

                                           Df     AIC      F Pr(>F)  
- pcpt.60d                                  1 -28.885 1.3855  0.065 .
- Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4)  4 -28.238                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


best Model when climate is partialled out
```
pRDA.CHNclimate <- rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date),
data = Climate.Data) 

pRDA.CHNclimate.optModel <- ordistep(pRDA.CHNclimate)


Step: GenData ~ PCNM1 + PCNM3 + Condition(shadow.days + sol.rad.60d +      pcpt.60d + day10cm + temp.laying.date) 

                                                                               Df AIC     F      Pr(>F)
- PCNM1                                                                         1 -27.492 1.6120 0.020*
- PCNM3                                                                         1 -27.002 1.9418 0.005**
- Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date)  5 -28.108

---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



#### 3. CHS

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHSData <- read.csv("CHS275.MAF.csv", header=T)
CHSData <- CHSData[complete.cases(CHSData),]  ##stba is removed here due to missing data
GenData <- CHSData[,13:ncol(CHSData)]
GenData <- decostand(GenData, "hellinger") ##detrend data

Climate.Data <- CHSData[,3:12]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHSfull <- rda(GenData ~ PCNM1 + PCNM2 + PCNM3 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHSfull

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

              Inertia Proportion Rank
Total          0.3083     1.0000     
Constrained    0.1570     0.5093    8
Unconstrained  0.1513     0.4907   15
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7    RDA8 
0.06531 0.02335 0.01610 0.01398 0.01228 0.01019 0.00839 0.00739 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025153 0.021440 0.015275 0.012001 0.009863 0.009675 0.009224 0.008483 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15 
0.007934 0.006660 0.006440 0.005588 0.004965 0.004347 0.004217  


RsquareAdj(RDA.CHSfull)

$r.squared
[1] 0.5092838

$adj.r.squared
[1] 0.2475686



anova(RDA.CHSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     8  0.15699 1.9459  0.001 ***
Residual 15  0.15127                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHSfull.pdf")
plot(RDA.CHSfull, main="RDA CHSfull")
dev.off()
```



```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHSgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(PCNM1 + PCNM2 + PCNM3), Climate.Data)

RsquareAdj(pRDA.CHSgeog)

$r.squared
[1] 0.2023947

$adj.r.squared
[1] 0.044646


anova(pRDA.CHSgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5 0.062389 1.2373  0.057 .
Residual 15 0.151266                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CHS GeneticData

pRDA.CHSgeog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.30826    1.00000     
Conditional   0.09460    0.30689    3
Constrained   0.06239    0.20239    5
Unconstrained 0.15127    0.49072   15
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5 
0.017872 0.016230 0.012429 0.008427 0.007432 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025153 0.021440 0.015275 0.012001 0.009863 0.009675 0.009224 0.008483 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15 
0.007934 0.006660 0.006440 0.005588 0.004965 0.004347 0.004217  


pdf(file="pRDA.CHSgeog.pdf")
plot(pRDA.CHSgeog, main="pRDA CHS(geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHSclimate <- rda(GenData~PCNM1 + PCNM2 + PCNM3 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHSclimate

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 +
Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date), data = Climate.Data)

              Inertia Proportion Rank
Total         0.30826    1.00000     
Conditional   0.10548    0.34219    5
Constrained   0.05151    0.16709    3
Unconstrained 0.15127    0.49072   15
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3 
0.025113 0.015476 0.010918 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025153 0.021440 0.015275 0.012001 0.009863 0.009675 0.009224 0.008483 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15 
0.007934 0.006660 0.006440 0.005588 0.004965 0.004347 0.004217  


RsquareAdj(pRDA.CHSclimate)

$r.squared
[1] 0.1670919

$adj.r.squared
[1] 0.08810106

anova(pRDA.CHSclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     3 0.051507 1.7025  0.002 **
Residual 15 0.151266                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: rejected -> Geography explains a significant portion of GenData

pdf(file="pRDA.CHSclimate.pdf")
plot(pRDA.CHSclimate, main="pRDA (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate


```
head(summary(pRDA.CHSgeog))

Biplot scores for constraining variables

                     RDA1      RDA2      RDA3     RDA4     RDA5 PC1
shadow.days       0.55907 -0.436678 -0.007770  0.30985  0.34680   0
sol.rad.60d      -0.51434  0.002098  0.612895 -0.16668  0.17302   0
pcpt.60d          0.06542 -0.143825  0.424469 -0.01164 -0.16270   0
day10cm          -0.12714  0.061116  0.725570  0.36524  0.08033   0
temp.laying.date -0.09510 -0.632217 -0.004846  0.55377 -0.41632   0

RDA1: temp.laying.date, shadow.days

RDA2: sol.rad.60d, day10cm
```


Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CHSgeog)


Biplot scores for constraining variables

                     RDA1      RDA2      RDA3     RDA4     RDA5 PC1
shadow.days       0.55907 -0.436678 -0.007770  0.30985  0.34680   0
sol.rad.60d      -0.51434  0.002098  0.612895 -0.16668  0.17302   0
pcpt.60d          0.06542 -0.143825  0.424469 -0.01164 -0.16270   0
day10cm          -0.12714  0.061116  0.725570  0.36524  0.08033   0
temp.laying.date -0.09510 -0.632217 -0.004846  0.55377 -0.41632   0


#eigenvalue (i.e. weight) of each RDA
pRDA.CHSgeog$CCA$eig

      RDA1        RDA2        RDA3        RDA4        RDA5 
0.017872167 0.016229721 0.012428778 0.008427029 0.007431517 
```



Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CHS.Climate.contribution <- (64.5578+44.88)  ##get these values from the RDA table. From the full and partial models run above
Venn.CHS.Geog.contribution <- (50.3074+44.88)
Venn.CHS.Geog.Climate.contribution <- 44.88

pdf("CHS.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CHS.Geog.contribution, area2=Venn.CHS.Climate.contribution, cross.area=Venn.CHS.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```



Best model overall
```
RDA.CHSfull.bestModel <- ordistep(RDA.CHSfull)


Step: GenData ~ PCNM1 + PCNM2 

        Df     AIC      F Pr(>F)   
- PCNM2  1 -30.808 1.8040  0.055 . 
- PCNM1  1 -26.822 5.9246  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


Best model when space is partialled out
```
step.res.CHS.pgeog <- ordistep(pRDA.CHSgeog)

Step: GenData ~ sol.rad.60d + pcpt.60d + temp.laying.date + Condition(PCNM1 +      PCNM2 + PCNM3) 

                                   Df     AIC      F Pr(>F)  
- temp.laying.date                  1 -29.317 1.3158  0.095 .
- pcpt.60d                          1 -29.225 1.3865  0.070 .
- sol.rad.60d                       1 -28.929 1.6142  0.030 *
- Condition(PCNM1 + PCNM2 + PCNM3)  3 -28.561                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```


Best model when climate is partialled out
```
step.res.CHS.pclimate <- ordistep(pRDA.CHSclimate)

Step: GenData ~ PCNM1 + Condition(shadow.days + sol.rad.60d + pcpt.60d +      day10cm + temp.laying.date) 

                                                                               Df AIC     F      Pr(>F)
- PCNM1                                                                         1 -27.317 2.3807 0.005**
- Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date)  5 -30.808

---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


```


#### 4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,12:ncol(CHS.VS.Data)]
GenData <- decostand(GenData, "hellinger") ##detrend data
Climate.Data <- CHS.VS.Data[,3:11]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.VSfull <- rda(GenData ~ PCNM1 + PCNM2 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.VSfull

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + shadow.days + sol.rad.60d
+ pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Constrained   0.18469    0.85575    7
Unconstrained 0.03113    0.14425    2
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7 
0.05298 0.03092 0.02816 0.02354 0.01883 0.01544 0.01481 

Eigenvalues for unconstrained axes:
     PC1      PC2 
0.016505 0.014626  




RsquareAdj(RDA.CHS.VSfull)
$r.squared
[1] 0.8557526

$adj.r.squared
[1] 0.3508868


anova(RDA.CHS.VSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     7 0.184687 1.695  0.008 **
Residual  2 0.031131                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.VSfull.noelev.pdf")
plot(RDA.CHS.VSfull)
dev.off()
```



```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.VS.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2), Climate.Data)

pRDA.CHS.VS.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Conditional   0.05274    0.24436    2
Constrained   0.13195    0.61139    5
Unconstrained 0.03113    0.14425    2
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5 
0.05043 0.02553 0.02232 0.01755 0.01612 

Eigenvalues for unconstrained axes:
     PC1      PC2 
0.016505 0.014626 


RsquareAdj(pRDA.CHS.VS.geog)
$r.squared
[1] 0.6113917

$adj.r.squared
[1] 0.3224228


anova(pRDA.CHS.VS.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5 0.131950 1.6954  0.088 .
Residual  2 0.031131                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 not rejected: Climate alone does not explain CHS.VS GeneticData


pdf(file="pRDA.CHS.VS.geog.pdf")
plot(pRDA.CHS.VS.geog, main="pRDA CHS.VS (geog partialled out)")
dev.off()
```




```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.VS.climate <- rda(GenData~PCNM1 + PCNM2+ Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.VS.climate

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Conditional   0.12718    0.58930    5
Constrained   0.05750    0.26645    2
Unconstrained 0.03113    0.14425    2
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2 
0.04081 0.01670 

Eigenvalues for unconstrained axes:
     PC1      PC2 
0.016505 0.014626 


#head(summary(pRDA.CHS.VS.climate))

RsquareAdj(pRDA.CHS.VS.climate)

$r.squared
[1] 0.2664485

$adj.r.squared
[1] 0.2749526


anova(pRDA.CHS.VS.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2 0.057505 1.8472  0.187
Residual  2 0.031131      

H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.pdf")
plot(pRDA.CHS.VS.climate, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
head(summary(pRDA.CHS.VS.geog))

Biplot scores for constraining variables

                    RDA1    RDA2    RDA3     RDA4     RDA5 PC1
shadow.days       0.2906  0.4324  0.3755  0.09537  0.08563   0
sol.rad.60d      -0.5865 -0.1228 -0.1705 -0.47359 -0.26945   0
pcpt.60d         -0.5241  0.5999 -0.5789 -0.08392  0.02238   0
day10cm          -0.4778  0.3572 -0.4011 -0.29523  0.14967   0
temp.laying.date -0.4201  0.4259 -0.4296  0.50117  0.44500   0

RDA1: day10cm, pcpt.60d, temp.laying.date

RDA2: pcpt.60d, temp.laying.date, day10cm
```


Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CHS.VS.geog)


Biplot scores for constraining variables

                    RDA1    RDA2    RDA3     RDA4     RDA5 PC1
shadow.days       0.2906  0.4324  0.3755  0.09537  0.08563   0
sol.rad.60d      -0.5865 -0.1228 -0.1705 -0.47359 -0.26945   0
pcpt.60d         -0.5241  0.5999 -0.5789 -0.08392  0.02238   0
day10cm          -0.4778  0.3572 -0.4011 -0.29523  0.14967   0
temp.laying.date -0.4201  0.4259 -0.4296  0.50117  0.44500   0


#eigenvalue (i.e. weight) of each RDA
pRDA.CHS.VS.geog$CCA$eig

     RDA1       RDA2       RDA3       RDA4       RDA5 
0.05043003 0.02553040 0.02231961 0.01755057 0.01611916 
```



Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CHS.VS.Climate.contribution <- (141.9559+-6.87)  ##get these values from the RDA table. From the full and partial models run above
Venn.CHS.VS.Geog.contribution <- (61.2333+-6.87)
Venn.CHS.VS.Geog.Climate.contribution <- -6.87

pdf("CHS.VS.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CHS.VS.Geog.contribution, area2=Venn.CHS.VS.Climate.contribution, cross.area=Venn.CHS.VS.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```


Best model overall
```
RDA.CHS.VSfull.bestModel <- ordistep(RDA.CHS.VSfull)

Step: GenData ~ PCNM1 + shadow.days + sol.rad.60d + pcpt.60d + temp.laying.date 

                   Df     AIC      F Pr(>F)   
- pcpt.60d          1 -14.373 1.4018  0.075 . 
- temp.laying.date  1 -14.108 1.5468  0.050 * 
- shadow.days       1 -13.963 1.6278  0.020 * 
- sol.rad.60d       1 -13.516 1.8850  0.010 **
- PCNM1             1 -13.203 2.0723  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```


Best model when space is partialled out
```
step.res.CHS.VS.pgeog <- ordistep(pRDA.CHS.VS.geog)

Step: GenData ~ 1 

```


Best model when climate is partialled out
```
step.res.CHS.VS.pclimate <- ordistep(pRDA.CHS.VS.climate)

Step: GenData ~ 1 
```




#### 5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI140.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,14:ncol(CHS.TI.Data)]
GenData <- decostand(GenData, "hellinger") ##detrend data
Climate.Data <- CHS.TI.Data[,3:13]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.TIfull <- rda(GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.TIfull

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +
shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data
= Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Constrained   0.20151    0.69687    9
Unconstrained 0.08765    0.30313    4
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7    RDA8    RDA9 
0.04670 0.03670 0.02742 0.02449 0.01578 0.01473 0.01367 0.01127 0.01075 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4 
0.03672 0.01934 0.01838 0.01321 


RsquareAdj(RDA.CHS.TIfull)
$r.squared
[1] 0.6968709

$adj.r.squared
[1] 0.01483055

anova(RDA.CHS.TIfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)
Model     9 0.201506 1.0217  0.471
Residual  4 0.087652              


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.TIfull.pdf")
plot(RDA.CHS.TIfull)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.TI.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4), Climate.Data)

pRDA.CHS.TI.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4),
data = Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Conditional   0.11263    0.38951    4
Constrained   0.08888    0.30737    5
Unconstrained 0.08765    0.30313    4
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5 
0.027112 0.018597 0.015661 0.014847 0.012661 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4 
0.03672 0.01934 0.01838 0.01321 


RsquareAdj(pRDA.CHS.TI.geog)
$r.squared
[1] 0.3073656

$adj.r.squared
[1] -0.1033438


anova(pRDA.CHS.TI.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5 0.088877 0.8112  0.803
Residual  4 0.087652                        


H0 not rejected: Climate alone does not explain CHS.TI GeneticData



pdf(file="pRDA.CHS.TI.geog.pdf")
plot(pRDA.CHS.TI.geog, main="pRDA CHS.TI (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.TI.climate <- rda(GenData~PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.TI.climate

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +
Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date), data = Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Conditional   0.12148    0.42011    5
Constrained   0.08003    0.27676    4
Unconstrained 0.08765    0.30313    4
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4 
0.030250 0.023329 0.014874 0.011575 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4 
0.03672 0.01934 0.01838 0.01321 

RsquareAdj(pRDA.CHS.TI.climate)
$r.squared
[1] 0.2767607

$adj.r.squared
[1] -0.04284855

#head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance     F Pr(>F)
Model     4 0.080028 0.913  0.628
Residual  4 0.087652 


H0: Rejected -> Geography explains a sig proportion of GenData

pdf(file="pRDA.CHS.TI.climate.pdf")
plot(pRDA.CHS.TI.climate, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```





Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.TI.geog))

Biplot scores for constraining variables

                   RDA1    RDA2     RDA3     RDA4     RDA5 PC1
shadow.days      0.2011 -0.1344 -0.16694  0.58450  0.18647   0
sol.rad.60d      0.1701  0.7226  0.17161  0.13981 -0.18513   0
pcpt.60d         0.4172  0.4533 -0.16464 -0.03591 -0.31509   0
day10cm          0.2936  0.5101  0.12119  0.28509  0.03155   0
temp.laying.date 0.2292 -0.3434  0.06574  0.41886 -0.04221   0





RDA1: shadow.days

RDA2: sol.rad.60d, day10cm, temp.laying.date
```

Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CHS.TI.geog)


Biplot scores for constraining variables

                   RDA1    RDA2     RDA3     RDA4     RDA5 PC1
shadow.days      0.2011 -0.1344 -0.16694  0.58450  0.18647   0
sol.rad.60d      0.1701  0.7226  0.17161  0.13981 -0.18513   0
pcpt.60d         0.4172  0.4533 -0.16464 -0.03591 -0.31509   0
day10cm          0.2936  0.5101  0.12119  0.28509  0.03155   0
temp.laying.date 0.2292 -0.3434  0.06574  0.41886 -0.04221   0



#eigenvalue (i.e. weight) of each RDA
pRDA.CHS.TI.geog$CCA$eig

      RDA1       RDA2       RDA3       RDA4       RDA5 
0.02711151 0.01859665 0.01566078 0.01484691 0.01266147 

```



Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CHS.TI.Climate.contribution <- (81.9553+28.72)  ##get these values from the RDA table. From the full and partial models run above
Venn.CHS.TI.Geog.contribution <- (76.6237+28.72)
Venn.CHS.TI.Geog.Climate.contribution <- 28.72

pdf("CHS.TI.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CHS.TI.Geog.contribution, area2=Venn.CHS.TI.Climate.contribution, cross.area=Venn.CHS.TI.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```

Best model overall
```
RDA.CHS.TIfull.bestModel <- ordistep(RDA.CHS.TIfull)


Step: GenData ~ sol.rad.60d + pcpt.60d 

              Df     AIC      F Pr(>F)  
- sol.rad.60d  1 -15.840 1.3467  0.095 .
- pcpt.60d     1 -15.675 1.4926  0.040 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


Best model when space is partialled out
```
step.res.CHS.TI.pgeog <- ordistep(pRDA.CHS.TI.geog)

Step: GenData ~ 1
```


Best model when climate is partialled out
```
step.res.CHS.TI.pclimate <- ordistep(pRDA.CHS.TIclimate)

Step: GenData ~ 1 
```


#### 6. CZ

```
library(vegan)

CZData <- read.csv("CZ404.MAF.csv", header=T)
CZData <- CZData[complete.cases(CZData),]  ##stba is removed here due to missing data
GenData <- CZData[,17:ncol(CZData)]
GenData <- decostand(GenData, "hellinger")
Climate.Data <- CZData[,3:16]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CZfull <- rda(GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. 

RDA.CZfull

Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 +
PCNM6 + PCNM7 + shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total          0.2568     1.0000     
Constrained    0.1249     0.4864   12
Unconstrained  0.1319     0.5136   24
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7    RDA8    RDA9   RDA10 
0.05082 0.01727 0.00949 0.00772 0.00704 0.00610 0.00586 0.00482 0.00438 0.00414 
  RDA11   RDA12 
0.00374 0.00353 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.016272 0.011501 0.009654 0.008946 0.007339 0.006615 0.005810 0.005438 
(Showed only 8 of all 24 unconstrained eigenvalues)


##calculate proportion explained (R2)

RsquareAdj(RDA.CZfull)

$r.squared
[1] 0.4864176

$adj.r.squared
[1] 0.2296264



anova(RDA.CZfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model    12  0.12492 1.8942  0.001 ***
Residual 24  0.13190                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CZfull.pdf")
plot(RDA.CZfull, main="CZ RDAfull")
dev.off()
```



```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CZ.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7), Climate.Data)

#head(summary(pRDA.CZ.geog.noelev))

pRDA.CZ.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 +
PCNM5 + PCNM6 + PCNM7), data = Climate.Data)

              Inertia Proportion Rank
Total         0.25682    1.00000     
Conditional   0.09584    0.37318    7
Constrained   0.02908    0.11323    5
Unconstrained 0.13190    0.51358   24
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4     RDA5 
0.008670 0.007009 0.005294 0.004295 0.003812 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.016272 0.011501 0.009654 0.008946 0.007339 0.006615 0.005810 0.005438 
(Showed only 8 of all 24 unconstrained eigenvalues)


RsquareAdj(pRDA.CZ.geog)

$r.squared
[1] 0.1132332

$adj.r.squared
[1] 0.007742354

anova(pRDA.CZ.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5  0.02908 1.0583  0.261
Residual 24  0.13190     


H0 rejected: Climate explains a significant amount of CZ GeneticData


pdf(file="pRDA.CZ.geog.pdf")
plot(pRDA.CZ.geog, main="pRDA CZ (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CZ.climate <- rda(GenData~PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CZ.climate
Call: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 +
Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm +
temp.laying.date), data = Climate.Data)

              Inertia Proportion Rank
Total         0.25682    1.00000     
Conditional   0.06405    0.24940    5
Constrained   0.04524    0.17617    4
Unconstrained 0.14752    0.57443   27
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2     RDA3     RDA4 
0.017897 0.015855 0.006479 0.005014 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.016724 0.012008 0.009774 0.009395 0.007985 0.006913 0.006228 0.005661 
(Showed only 8 of all 27 unconstrained eigenvalues)



#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CZ.climate)

$r.squared
[1] 0.1761724

$adj.r.squared
[1] 0.1057607


anova(pRDA.CZ.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     4 0.045244 2.0702  0.001 ***
Residual 27 0.147524                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CZ.climate.pdf")
plot(pRDA.CZ.climate, main="pRDA CZ (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
(head(summary(pRDA.CZ.geog)))

Biplot scores for constraining variables

                     RDA1     RDA2      RDA3    RDA4     RDA5 PC1
shadow.days      -0.07594 -0.04237 -0.422042 -0.7352  0.13552   0
sol.rad.60d      -0.64232 -0.32530  0.357004  0.3834  0.33305   0
pcpt.60d         -0.84298  0.20267  0.253135  0.1461  0.19477   0
day10cm          -0.68834 -0.07558 -0.007803 -0.1804  0.03830   0
temp.laying.date  0.08341  0.49105  0.558663 -0.2284 -0.06844   0


RDA1: pcpt.60d, day10cm

RDA2: sol.rad.60d
```


Calculate the proportion of variance explained by each environmental variable when space is partialled out
```

#This is calculated as
#sum(|env.var correlation with RDA|*eigenvalue of RDA) for all RDA / total variance  
#Total variance=calcualted in Full model

#Correlation of Env.variable with each RDA axis: 

summary(pRDA.CZ.geog)


Biplot scores for constraining variables

                     RDA1     RDA2      RDA3    RDA4     RDA5 PC1
shadow.days      -0.07594 -0.04237 -0.422042 -0.7352  0.13552   0
sol.rad.60d      -0.64232 -0.32530  0.357004  0.3834  0.33305   0
pcpt.60d         -0.84298  0.20267  0.253135  0.1461  0.19477   0
day10cm          -0.68834 -0.07558 -0.007803 -0.1804  0.03830   0
temp.laying.date  0.08341  0.49105  0.558663 -0.2284 -0.06844   0


#eigenvalue (i.e. weight) of each RDA
pRDA.CZ.geog$CCA$eig

       RDA1        RDA2        RDA3        RDA4        RDA5 
0.008670188 0.007008761 0.005294453 0.004294882 0.003812053 

```



Figure of the partitioning of variance
```
library(VennDiagram)

Venn.CZ.Climate.contribution <- (32.8959+63.06)  ##get these values from the RDA table. From the full and partial models run above
Venn.CZ.Geog.contribution <- (56.6319+63.06)
Venn.CZ.Geog.Climate.contribution <- 63.06

pdf("CZ.Venn.RDA.pdf")
draw.pairwise.venn(area1=Venn.CZ.Geog.contribution, area2=Venn.CZ.Climate.contribution, cross.area=Venn.CZ.Geog.Climate.contribution, category=c("Space", "Climate"), fill=c("blue", "green"))
dev.off()
```

Best model overall
```
RDA.CZfull.bestModel <- ordistep(RDA.CZfull)



```


Best model when space is partialled out
```
step.res.CZ.pgeog <- ordistep(pRDA.CZ.geog)

Step: GenData ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + pcpt.60d 

           Df     AIC      F Pr(>F)   
- PCNM4     1 -55.159 1.4761  0.070 . 
- pcpt.60d  1 -55.114 1.5159  0.055 . 
- PCNM2     1 -53.396 3.0612  0.015 * 
- PCNM1     1 -52.036 4.3364  0.005 **
- PCNM3     1 -50.942 5.3966  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


Best model when climate is partialled out
```
step.res.CZ.pclimate <- ordistep(pRDA.CZclimate)


```





## Part 2

I want to know whether the same environmental variables underlie adaptation across the different geographic regions. 

I will test this using RDA of the adaptive loci identified before (see Chp2_VennDiagrams.md)

I'm using only the loci identified by 2 or more methods (lfmm, bayenv2, XtX, or pcadapt) for this analysis. 

1. First I need to filter the vcf file to get the adaptive loci 

2. Convert to input files

3. Run full RDA

4. Run partial RDA 

This will be done for each region. ie. 6 analyses. 

### 1. Adpative loci

Rename the list of duplicated loci so that vcftools can recognise them

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/Venn/CHP2

```
CHN.duplicated.outliers <- read.table("CHN.duplicated.outliers", header=F)
CHN.duplicated.outliers <- as.data.frame(CHN.duplicated.outliers)
CHN.duplicated.outliers.names <- gsub("X.", "", CHN.duplicated.outliers$V1)
CHN.duplicated.outliers.names <- gsub("\\.", ":", CHN.duplicated.outliers.names)
write.table(CHN.duplicated.outliers.names, "CHN.duplicated.outliers.names", col.names=F, quote=F, row.names=F, sep=" ")

CHS.duplicated.outliers <- read.table("CHS.duplicated.outliers", header=F)
CHS.duplicated.outliers <- as.data.frame(CHS.duplicated.outliers)
CHS.duplicated.outliers.names <- gsub("X.", "", CHS.duplicated.outliers$V1)
CHS.duplicated.outliers.names <- gsub("\\.", ":", CHS.duplicated.outliers.names)
write.table(CHS.duplicated.outliers.names, "CHS.duplicated.outliers.names", col.names=F, quote=F, row.names=F, sep=" ")

CZ.duplicated.outliers <- read.table("CZ.duplicated.outliers", header=F)
CZ.duplicated.outliers <- as.data.frame(CZ.duplicated.outliers)
CZ.duplicated.outliers.names <- gsub("X.", "", CZ.duplicated.outliers$V1)
CZ.duplicated.outliers.names <- gsub("\\.", ":", CZ.duplicated.outliers.names)
write.table(CZ.duplicated.outliers.names, "CZ.duplicated.outliers.names", col.names=F, quote=F, row.names=F, sep=" ")


CHS.TI.duplicated.outliers <- read.table("CHS.TI.duplicated.outliers", header=F)
CHS.TI.duplicated.outliers <- as.data.frame(CHS.TI.duplicated.outliers)
CHS.TI.duplicated.outliers.names <- gsub("X.", "", CHS.TI.duplicated.outliers$V1)
CHS.TI.duplicated.outliers.names <- gsub("\\.", ":", CHS.TI.duplicated.outliers.names)
write.table(CHS.TI.duplicated.outliers.names, "CHS.TI.duplicated.outliers.names", col.names=F, quote=F, row.names=F, sep=" ")


CHS.VS.duplicated.outliers <- read.table("CHS.VS.duplicated.outliers", header=F)
CHS.VS.duplicated.outliers <- as.data.frame(CHS.VS.duplicated.outliers)
CHS.VS.duplicated.outliers.names <- gsub("X.", "", CHS.VS.duplicated.outliers$V1)
CHS.VS.duplicated.outliers.names <- gsub("\\.", ":", CHS.VS.duplicated.outliers.names)
write.table(CHS.VS.duplicated.outliers.names, "CHS.VS.duplicated.outliers.names", col.names=F, quote=F, row.names=F, sep=" ")
```


Filter the vcf file

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA/AdaptiveLociOnly
```
#cp /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/Venn/CHP2/*outliers.names .

vcftools --vcf CHN.229.9608.recode.vcf --snps CHN.duplicated.outliers.names --recode --recode-INFO-all --out CHN.229.25AdaptiveLoci

vcftools --vcf CHS.283.9608.recode.vcf --snps CHS.duplicated.outliers.names --recode --recode-INFO-all --out CHS.283.228AdaptiveLoci

vcftools --vcf CZ.404.9608.recode.vcf --snps CZ.duplicated.outliers.names --recode --recode-INFO-all --out CZ.404.46AdaptiveLoci

vcftools --vcf CHS.VS.135.9608.recode.vcf --snps CHS.VS.duplicated.outliers.names --recode --recode-INFO-all --out CHS.VS.135.356AdaptiveLoci

vcftools --vcf CHS.TI.148.9608.recode.vcf --snps CHS.TI.duplicated.outliers.names --recode --recode-INFO-all --out CHS.TI.148.285AdaptiveLoci

```

### 2. Input files:

convert all files to plink 

```
vcftools --vcf CHN.229.25AdaptiveLoci.recode.vcf --plink --out CHN.229.25AdaptiveLoci
vcftools --vcf CHS.283.228AdaptiveLoci.recode.vcf --plink --out CHS.283.228AdaptiveLoci
vcftools --vcf CZ.404.46AdaptiveLoci.recode.vcf --plink --out CZ.404.46AdaptiveLoci
vcftools --vcf CHS.TI.148.285AdaptiveLoci.recode.vcf --plink --out CHS.TI.148.285AdaptiveLoci
vcftools --vcf CHS.VS.135.356AdaptiveLoci.recode.vcf --plink --out CHS.VS.135.356AdaptiveLoci
```

Calculate the MAF of all loci

```
#Calculate MAF for the full dataset within region using PLINK

This was all done previously, so we can just copy the cluster.pop files from the input folder to here

cp ~/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/*cluster.pop .
```


########Still have to do this and correct the code for CHall
CHall
```
###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CH940.MAF <- read.table("CH940.pop.frq.strat", header=T)
CH940.MAF2 <- CH940.MAF[,c(3,2,6)]
summary(CH940.MAF2)
       CLST                SNP              MAF         
 abnd   :  9608   1000065:55:    82   Min.   :0.00000  
 agra   :  9608   1000356:65:    82   1st Qu.:0.00000  
 alpl   :  9608   1000622:42:    82   Median :0.09375  
 apla   :  9608   1000723:85:    82   Mean   :0.19825  
 arce   :  9608   1000806:91:    82   3rd Qu.:0.31250  
 bach   :  9608   1000841:7 :    82   Max.   :1.00000  
 (Other):730208   (Other)   :787364    

library("ggplot2")
library("reshape2")

CH940.MAF3 <- melt(CH940.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CH940.MAF3)
head(CH940.MAF3)


CH940.MAF4 <- dcast(CH940.MAF3, formula= CLST ~ SNP)
head(CH940.MAF4)


##Add X infront of all locusnames. 
colnames(CH940.MAF4) <- paste("X", colnames(CH940.MAF4), sep=".")
write.csv(CH940.MAF4, file="CH940.MAF.csv")
```


CHN
```
plink --file CHN.229.25AdaptiveLoci --freq --within CHN229.cluster.pop --out CHN229.25Adaptive.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHN229.MAF <- read.table("CHN229.25Adaptive.pop.frq.strat", header=T)
CHN229.MAF2 <- CHN229.MAF[,c(3,2,6)]
summary(CHN229.MAF2)
     
library("ggplot2")
library("reshape2")

CHN229.MAF3 <- melt(CHN229.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN229.MAF3)
head(CHN229.MAF3)


CHN229.MAF4 <- dcast(CHN229.MAF3, formula= CLST ~ SNP)
head(CHN229.MAF4)


##Add X infront of all locusnames. 
colnames(CHN229.MAF4) <- paste("X", colnames(CHN229.MAF4), sep=".")
write.csv(CHN229.MAF4, file="CHN229.MAF.csv")

```

CHS
```
plink --file CHS.283.228AdaptiveLoci --freq --within CHS283.cluster.pop --out CHS283.228Adaptive.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS283.MAF <- read.table("CHS283.228Adaptive.pop.frq.strat", header=T)
CHS283.MAF2 <- CHS283.MAF[,c(3,2,6)]
summary(CHS283.MAF2)
     
library("ggplot2")
library("reshape2")

CHS283.MAF3 <- melt(CHS283.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS283.MAF3)
head(CHS283.MAF3)


CHS283.MAF4 <- dcast(CHS283.MAF3, formula= CLST ~ SNP)
head(CHS283.MAF4)


##Add X infront of all locusnames. 
colnames(CHS283.MAF4) <- paste("X", colnames(CHS283.MAF4), sep=".")
write.csv(CHS283.MAF4, file="CHS283.MAF.csv")

```


CHS.VS
```
plink --file CHS.VS.135.356AdaptiveLoci --freq --within CHS.VS.cluster.pop --out CHS.VS.135.365Adaptive.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS.VS135.MAF <- read.table("CHS.VS.135.365Adaptive.pop.frq.strat", header=T)
CHS.VS135.MAF2 <- CHS.VS135.MAF[,c(3,2,6)]
summary(CHS.VS135.MAF2)
     
library("ggplot2")
library("reshape2")

CHS.VS135.MAF3 <- melt(CHS.VS135.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS135.MAF3)
head(CHS.VS135.MAF3)


CHS.VS135.MAF4 <- dcast(CHS.VS135.MAF3, formula= CLST ~ SNP)
head(CHS.VS135.MAF4)


##Add X infront of all locusnames. 
colnames(CHS.VS135.MAF4) <- paste("X", colnames(CHS.VS135.MAF4), sep=".")
write.csv(CHS.VS135.MAF4, file="CHS.VS135.MAF.csv")

```


CHS.TI
```
plink --file CHS.TI.148.285AdaptiveLoci --freq --within CHS.TI.148.cluster.pop --out CHS.TI.148.285Adaptive.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS.TI148.MAF <- read.table("CHS.TI.148.285Adaptive.pop.frq.strat", header=T)
CHS.TI148.MAF2 <- CHS.TI148.MAF[,c(3,2,6)]
summary(CHS.TI148.MAF2)
     
library("ggplot2")
library("reshape2")

CHS.TI148.MAF3 <- melt(CHS.TI148.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI148.MAF3)
head(CHS.TI148.MAF3)


CHS.TI148.MAF4 <- dcast(CHS.TI148.MAF3, formula= CLST ~ SNP)
head(CHS.TI148.MAF4)


##Add X infront of all locusnames. 
colnames(CHS.TI148.MAF4) <- paste("X", colnames(CHS.TI148.MAF4), sep=".")
write.csv(CHS.TI148.MAF4, file="CHS.TI148.MAF.csv")

```


CZ  
```
plink --file CZ.404.46AdaptiveLoci --freq --within CZ404.cluster.pop --out CZ.404.46Adaptive.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CZ404.MAF <- read.table("CZ.404.46Adaptive.pop.frq.strat", header=T)
CZ404.MAF2 <- CZ404.MAF[,c(3,2,6)]
summary(CZ404.MAF2)
     
library("ggplot2")
library("reshape2")

CZ404.MAF3 <- melt(CZ404.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ404.MAF3)
head(CZ404.MAF3)


CZ404.MAF4 <- dcast(CZ404.MAF3, formula= CLST ~ SNP)
head(CZ404.MAF4)


##Add X infront of all locusnames. 
colnames(CZ404.MAF4) <- paste("X", colnames(CZ404.MAF4), sep=".")
write.csv(CZ404.MAF4, file="CZ404.MAF.csv")

```




#### 3. Geographic coordinates

paste coordinates into the MAF.csv file from Env.Data.all_20161025.csv

located: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2

I will include elevation for geographic information for the CH data. 

i.e. lat, long


#### 4. Climate variables

```
1. solar.rad.60d (amt of solar radiation 60 days after egg laying)

2. pcpt.60d (precipitation in the 60days after egg laying)

3. day10cm (calender day on which only 10cm of snow remains)

4. temp.laying.date

5. shadow.days
```

paste the 5 variables from /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/Env.Data.all_20161025.csv

I removed stba from these analyses since there's missing environmental data. 


### 5. Run RDA: Adaptive loci

#####STILL NEED TO RUN CHall: waiting for LFMM data. 

#### 1. CHall

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHallData <- read.csv("CH940.MAF.csv", header=T)
CHallData <- CHallData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CHallData[,11:9618]
Climate.Data <- CHallData[,3:10]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CH940full.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. I got an error (Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd') when elev was after Long in the formula

RDA.CH940full.noelev

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Constrained   121.9575     0.3656    7
Unconstrained 211.6300     0.6344   30
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
67.10 16.52 13.46  7.97  6.88  5.34  4.69 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.582 19.869 12.236 11.725 10.567  9.725  9.000  8.742 
(Showed only 8 of all 30 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CH940full)

$r.squared
[1] 0.3655938

$adj.r.squared
[1] 0.2175657



anova(RDA.CH940full)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   121.96 2.4698  0.001 ***
Residual 30   211.63                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CH940full.pdf")
plot(RDA.CH940full.noelev)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CH940.geog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

#head(summary(pRDA.CH940.geog.noelev))

pRDA.CH940.geog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Conditional    75.6364     0.2267    2
Constrained    46.3212     0.1389    5
Unconstrained 211.6300     0.6344   30
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
19.893  8.720  7.426  5.415  4.869 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.582 19.869 12.236 11.725 10.567  9.725  9.000  8.742 
(Showed only 8 of all 30 unconstrained eigenvalues)



RsquareAdj(pRDA.CH940.geog.noelev)

$r.squared
[1] 0.1388576

$adj.r.squared
[1] 0.03501604

anova(pRDA.CH940.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     5   46.321 1.3133  0.005 **
Residual 30  211.630                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CH940 GeneticData


pdf(file="pRDA.CH940.geog.pdf")
plot(pRDA.CH940.geog.noelev, main="pRDA CH940 (geog partialled out)")
dev.off()
```


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CH940.climate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CH940.climate.noelev

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Conditional    75.1625     0.2253    5
Constrained    46.7950     0.1403    2
Unconstrained 211.6300     0.6344   30
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
37.78  9.01 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.582 19.869 12.236 11.725 10.567  9.725  9.000  8.742 
(Showed only 8 of all 30 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate.noelev))

RsquareAdj(pRDA.CH940.climate.noelev)

$r.squared
[1] 0.140278

$adj.r.squared
[1] 0.1132943

anova(pRDA.CH940.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     2   46.795 3.3168  0.001 ***
Residual 30  211.630                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1        

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CH940.climate.pdf")
plot(pRDA.CH940.climate.noelev, main="pRDA CH940 (climate partialled out)")
dev.off()
```


Find the most important loci associated with Climate
```
 summary(pRDA.CH940.geog.noelev)

Biplot scores for constraining variables

                       RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days       3.189e-01 -0.1027  0.03954  0.7566  0.08027   0
sol.rad.60d      -8.736e-05  0.7303 -0.32012 -0.2854 -0.51741   0
pcpt.60d          7.321e-01  0.3476 -0.33359 -0.2673 -0.26957   0
day10cm           5.950e-01  0.4360 -0.10138  0.1803 -0.05493   0
temp.laying.date  3.224e-01 -0.4331 -0.68384 -0.1886  0.38036   0

RDA1: pcpt.60d, sol.rad.60d

RDA2: sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog.noelev <- summary(pRDA.CH940.geog.noelev)

CH940.RDAloadings.noelev <- sum.pRDA.CH940.geog.noelev$species  ##write the loadings to data.frame
CH940.RDAloadings.noelev <- as.data.frame(CH940.RDAloadings.noelev) #change to df
names(CH940.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings.noelev$RDA1)
hist(CH940.RDAloadings.noelev$RDA2)
hist(CH940.RDAloadings.noelev$RDA3)
hist(CH940.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted.noelev <- CH940.RDAloadings.noelev[order(-abs(CH940.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent.noelev <- CH940.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted.noelev <- CH940.RDAloadings.noelev[order(-abs(CH940.RDAloadings.noelev$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent.noelev <- CH940.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent.noelev, "CH940.noelev.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent.noelev, "CH940.noelev.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```


#### 2. CHN

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHNData <- read.csv("CHN229.MAF.csv", header=T)
GenData <- CHNData[,11:ncol(CHNData)]
Climate.Data <- CHNData[,3:10]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHNfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHNfull.noelev

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total          0.7673     1.0000     
Constrained    0.3767     0.4910    7
Unconstrained  0.3906     0.5090   11
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7 
0.13736 0.10591 0.05750 0.03397 0.02288 0.01209 0.00701 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.10820 0.08713 0.06936 0.04217 0.02564 0.01950 0.01332 0.00900 0.00785 0.00482 
   PC11 
0.00362 


RsquareAdj(RDA.CHNfull.noelev)

$r.squared
[1] 0.4909507

$adj.r.squared
[1] 0.1670102


anova(RDA.CHNfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     7  0.37673 1.5156  0.021 *
Residual 11  0.39061                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHNfull.pdf")
plot(RDA.CHNfull.noelev)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHNgeog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(lat + long), Climate.Data)


RsquareAdj(pRDA.CHNgeog.noelev)

$r.squared
[1] 0.3248671

$adj.r.squared
[1] 0.1051662


anova(pRDA.CHNgeog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance     F Pr(>F)  
Model     5  0.24928 1.404  0.083 .
Residual 11  0.39061               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


H0 not rejected: Climate does not explain a significant amount of CHN GeneticData

pRDA.CHNgeog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

              Inertia Proportion Rank
Total          0.7673     1.0000     
Conditional    0.1274     0.1661    2
Constrained    0.2493     0.3249    5
Unconstrained  0.3906     0.5090   11
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5 
0.09661 0.07988 0.03573 0.02712 0.00994 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.10820 0.08713 0.06936 0.04217 0.02564 0.01950 0.01332 0.00900 0.00785 0.00482 
   PC11 
0.00362 



pdf(file="pRDA.CHNgeog.pdf")
plot(pRDA.CHNgeog.noelev, main="pRDA CHN(geog partialled out)")
dev.off()
```


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHNclimate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHNclimate.noelev

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total          0.7673     1.0000     
Conditional    0.2780     0.3623    5
Constrained    0.0987     0.1286    2
Unconstrained  0.3906     0.5091   11
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2 
0.07546 0.02324 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.10820 0.08713 0.06936 0.04217 0.02564 0.01950 0.01332 0.00900 0.00785 0.00482 
   PC11 
0.00362



RsquareAdj(pRDA.CHNclimate.noelev)

$r.squared
[1] 0.128625

$adj.r.squared
[1] 0.0499439

anova(pRDA.CHNclimate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2  0.09870 1.3897  0.164
Residual 11  0.39061  
          
          
H0: not rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHNclimate.pdf")
plot(pRDA.CHNclimate.noelev, main="pRDA (climate partialled out)")
dev.off()
```



Find the most important loci associated with Climate


```
summary(pRDA.CHNgeog.noelev)

Biplot scores for constraining variables

                   RDA1     RDA2     RDA3     RDA4    RDA5 PC1
shadow.days      0.6060 -0.27422  0.30041  0.33894  0.2015   0
sol.rad.60d      0.2882 -0.41274 -0.05727 -0.64998 -0.3496   0
pcpt.60d         0.5283 -0.10237  0.02474 -0.45830 -0.3773   0
day10cm          0.4645 -0.29550 -0.19325  0.05193 -0.2803   0
temp.laying.date 0.6619  0.07423 -0.56561 -0.16970  0.3641   0

```




#### 3. CHS

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHSData <- read.csv("CHS283.MAF.csv", header=T)
CHSData <- CHSData[complete.cases(CHSData),]  ##stba is removed here due to missing data
GenData <- CHSData[,11:ncol(CHSData)]
Climate.Data <- CHSData[,3:10]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHSfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHSfull.noelev
Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total         14.4611     1.0000     
Constrained    7.2578     0.5019    7
Unconstrained  7.2032     0.4981   16
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
3.413 1.268 0.689 0.628 0.497 0.446 0.315 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
1.2412 0.9521 0.8603 0.6993 0.5383 0.4349 0.3713 0.3389 0.2933 0.2925 0.2472 
  PC12   PC13   PC14   PC15   PC16 
0.2311 0.2122 0.1846 0.1673 0.1387 




RsquareAdj(RDA.CHSfull.noelev)

$r.squared
[1] 0.5018892

$adj.r.squared
[1] 0.2839657



anova(RDA.CHSfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   7.2578 2.3031  0.001 ***
Residual 16   7.2032                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHSfull.pdf")
plot(RDA.CHSfull.noelev)
dev.off()
```

![alt_txt][CHS.Fig1]

[CHS.Fig1]:https://cloud.githubusercontent.com/assets/12142475/23801930/0bf1bc5a-05b2-11e7-8954-fd562989513f.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHSgeog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(lat + long), Climate.Data)

RsquareAdj(pRDA.CHSgeog.noelev)

$r.squared
[1] 0.2089448

$adj.r.squared
[1] 0.05835999


anova(pRDA.CHSgeog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   3.0216 1.3423  0.019 *
Residual 16   7.2032                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CHS GeneticData

pRDA.CHSgeog.noelev
Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

              Inertia Proportion Rank
Total         14.4611     1.0000     
Conditional    4.2363     0.2929    2
Constrained    3.0216     0.2089    5
Unconstrained  7.2032     0.4981   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
1.0765 0.6715 0.5014 0.4511 0.3211 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
1.2412 0.9521 0.8603 0.6993 0.5383 0.4349 0.3713 0.3389 0.2933 0.2925 0.2472 
  PC12   PC13   PC14   PC15   PC16 
0.2311 0.2122 0.1846 0.1673 0.1387 




pdf(file="pRDA.CHSgeog.pdf")
plot(pRDA.CHSgeog.noelev, main="pRDA CHS(geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHSclimate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHSclimate.noelev

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total         14.4611     1.0000     
Conditional    5.4991     0.3803    5
Constrained    1.7587     0.1216    2
Unconstrained  7.2032     0.4981   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
1.1365 0.6222 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
1.2412 0.9521 0.8603 0.6993 0.5383 0.4349 0.3713 0.3389 0.2933 0.2925 0.2472 
  PC12   PC13   PC14   PC15   PC16 
0.2311 0.2122 0.1846 0.1673 0.1387  

RsquareAdj(pRDA.CHSclimate.noelev)

$r.squared
[1] 0.1216166

$adj.r.squared
[1] 0.07583967

anova(pRDA.CHSclimate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     2   1.7587 1.9532  0.003 **
Residual 16   7.2032                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: rejected -> Geography explains a significant portion of GenData

pdf(file="pRDA.CHSclimate.pdf")
plot(pRDA.CHSclimate.noelev, main="pRDA (climate partialled out)")
dev.off()
```


![alt_txt][pRDA.CHSclimate]

[pRDA.CHSclimate]:https://cloud.githubusercontent.com/assets/12142475/23802225/f27c7a7a-05b2-11e7-86b8-6c12e8b3f767.png


Find the most important loci associated with Climate


```
head(summary(pRDA.CHSgeog.noelev))

Biplot scores for constraining variables

                    RDA1    RDA2     RDA3     RDA4    RDA5 PC1
shadow.days      -0.3799 -0.1817  0.71558 -0.27610  0.2090   0
sol.rad.60d      -0.2149  0.2470 -0.60767  0.19128  0.5792   0
pcpt.60d         -0.4253  0.3217 -0.02614  0.04376 -0.1302   0
day10cm          -0.3185  0.4211 -0.24987 -0.44247  0.4243   0
temp.laying.date -0.7107 -0.3337 -0.08768 -0.40981 -0.3839   0
```




#### 4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,11:ncol(CHS.VS.Data)]
Climate.Data <- CHS.VS.Data[,3:10]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.VSfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.VSfull.noelev

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total         19.1804     1.0000     
Constrained   16.4681     0.8586    7
Unconstrained  2.7123     0.1414    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
5.355 3.188 2.503 1.738 1.609 1.141 0.934 

Eigenvalues for unconstrained axes:
   PC1    PC2 
1.6150 1.0973 



RsquareAdj(RDA.CHS.VSfull.noelev)
$r.squared
[1] 0.85859

$adj.r.squared
[1] 0.3636548

anova(RDA.CHS.VSfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     7  16.4681 1.7348  0.011 *
Residual  2   2.7123                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.VSfull.noelev.pdf")
plot(RDA.CHS.VSfull.noelev)
dev.off()
```



```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.VS.geog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

pRDA.CHS.VS.geog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

              Inertia Proportion Rank
Total         19.1804     1.0000     
Conditional    4.7242     0.2463    2
Constrained   11.7439     0.6123    5
Unconstrained  2.7123     0.1414    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
4.887 2.746 1.692 1.411 1.008 

Eigenvalues for unconstrained axes:
   PC1    PC2 
1.6150 1.0973 


RsquareAdj(pRDA.CHS.VS.geog.noelev)
$r.squared
[1] 0.6122848

$adj.r.squared
[1] 0.3326911


anova(pRDA.CHS.VS.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5  11.7439 1.7319  0.083 .
Residual  2   2.7123                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 not rejected: Climate alone does not explain CHS.VS GeneticData


pdf(file="pRDA.CHS.VS.geog.noelev.pdf")
plot(pRDA.CHS.VS.geog.noelev, main="pRDA CHS.VS (geog partialled out)")
dev.off()
```




```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.VS.climate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.VS.climate.noelev

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total         19.1804     1.0000     
Conditional   11.7555     0.6129    5
Constrained    4.7126     0.2457    2
Unconstrained  2.7123     0.1414    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
3.483 1.230 

Eigenvalues for unconstrained axes:
   PC1    PC2 
1.6150 1.0973 


#head(summary(pRDA.CHS.VS.climate))

RsquareAdj(pRDA.CHS.VS.climate.noelev)

$r.squared
[1] 0.2457

$adj.r.squared
[1] 0.2346524

anova(pRDA.CHS.VS.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2   4.7126 1.7375  0.216
Residual  2   2.7123   


H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.noelev.pdf")
plot(pRDA.CHS.VS.climate.noelev, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.VS.geog.noelev))

Biplot scores for constraining variables

                   RDA1     RDA2     RDA3     RDA4     RDA5 PC1
shadow.days      0.2786 -0.06524 -0.18470 -0.23073 -0.14493   0
sol.rad.60d      0.2056  0.24342 -0.24195  0.14293  0.49760   0
pcpt.60d         0.5070  0.75321 -0.19807  0.35447 -0.02805   0
day10cm          0.6165  0.63458 -0.06774  0.05259  0.24590   0
temp.laying.date 0.5920  0.42961  0.28835  0.31030 -0.51300   0
```




#### 5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI148.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,11:ncol(CHS.TI.Data)]
Climate.Data <- CHS.TI.Data[,3:10]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.TIfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.TIfull.noelev

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total         10.5378     1.0000     
Constrained    6.5350     0.6202    7
Unconstrained  4.0028     0.3798    6
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
1.7255 1.4654 1.1525 0.7753 0.5209 0.4703 0.4252 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
1.0375 1.0220 0.6755 0.5293 0.4146 0.3237 


RsquareAdj(RDA.CHS.TIfull.noelev)
$r.squared
[1] 0.6201516

$adj.r.squared
[1] 0.1769952

anova(RDA.CHS.TIfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     7   6.5350 1.3994  0.011 *
Residual  6   4.0028                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.TIfull.noelev.pdf")
plot(RDA.CHS.TIfull.noelev)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.TI.geog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

pRDA.CHS.TI.geog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

              Inertia Proportion Rank
Total         10.5378     1.0000     
Conditional    2.2451     0.2131    2
Constrained    4.2899     0.4071    5
Unconstrained  4.0028     0.3798    6
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
1.6975 1.1002 0.5380 0.4778 0.4764 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
1.0375 1.0220 0.6755 0.5293 0.4146 0.3237 


RsquareAdj(pRDA.CHS.TI.geog.noelev)
$r.squared
[1] 0.4070993

$adj.r.squared
[1] 0.1070243



anova(pRDA.CHS.TI.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   4.2899 1.2861  0.095 .
Residual  6   4.0028                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



H0 not rejected: Climate alone does not explain CHS.TI GeneticData



pdf(file="pRDA.CHS.TI.geog.noelev.pdf")
plot(pRDA.CHS.TI.geog.noelev, main="pRDA CHS.TI (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.TI.climate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.TI.climate.noelev

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total         10.5378     1.0000     
Conditional    4.4047     0.4180    5
Constrained    2.1303     0.2022    2
Unconstrained  4.0028     0.3798    6
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
1.2685 0.8618 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
1.0375 1.0220 0.6755 0.5293 0.4146 0.3237  

RsquareAdj(pRDA.CHS.TI.climate.noelev)
$r.squared
[1] 0.2021603

$adj.r.squared
[1] 0.1227593

#head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   2.1303 1.5966  0.065 .
Residual  6   4.0028                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

H0: Accepted -> Geography does not explains a sig proportion of GenData

pdf(file="pRDA.CHS.TI.climate.noelev.pdf")
plot(pRDA.CHS.TI.climate.noelev, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```





Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.TI.geog.noelev))

Biplot scores for constraining variables

                     RDA1     RDA2      RDA3     RDA4     RDA5 PC1
shadow.days       0.34968  0.77051 -0.007261 -0.42893 -0.07589   0
sol.rad.60d       0.14349 -0.67641  0.397620 -0.01601 -0.58221   0
pcpt.60d          0.39277  0.05651  0.696743  0.05771  0.03435   0
day10cm           0.35301 -0.24598  0.125802 -0.13424 -0.52315   0
temp.laying.date -0.08754  0.72435  0.095812  0.34925 -0.39001   0


```



#### 6. CZ

```
library(vegan)

CZData <- read.csv("CZ404.MAF.csv", header=T)
CZData <- CZData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CZData[,11:ncol(CZData)]
Climate.Data <- CZData[,3:10]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CZfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. 

RDA.CZfull.noelev

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

              Inertia Proportion Rank
Total          2.8016     1.0000     
Constrained    1.0613     0.3788    7
Unconstrained  1.7403     0.6212   30
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
0.6411 0.1349 0.1185 0.0921 0.0342 0.0236 0.0171 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
0.4198 0.1815 0.1446 0.1386 0.1213 0.0955 0.0853 0.0764 
(Showed only 8 of all 30 unconstrained eigenvalues)


##calculate proportion explained (R2)

RsquareAdj(RDA.CZfull.noelev)

$r.squared
[1] 0.3788244

$adj.r.squared
[1] 0.2338834



anova(RDA.CZfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   1.0613 2.6136  0.001 ***
Residual 30   1.7403                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CZfull.elev.noelev.pdf")
plot(RDA.CZfull.noelev)
dev.off()
```



```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CZ.geog.noelev <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

#head(summary(pRDA.CZ.geog.noelev))

pRDA.CZ.geog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

              Inertia Proportion Rank
Total          2.8016     1.0000     
Conditional    0.7000     0.2499    2
Constrained    0.3613     0.1290    5
Unconstrained  1.7403     0.6212   30
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2    RDA3    RDA4    RDA5 
0.14642 0.12356 0.04478 0.02888 0.01767 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
0.4198 0.1815 0.1446 0.1386 0.1213 0.0955 0.0853 0.0764 
(Showed only 8 of all 30 unconstrained eigenvalues)




RsquareAdj(pRDA.CZ.geog.noelev)

$r.squared
[1] 0.1289615

$adj.r.squared
[1] 0.02688544

anova(pRDA.CZ.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   0.3613 1.2457  0.099 .
Residual 30   1.7403                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 accepted: Climate alone does not explains a significant amount of CZ GeneticData


pdf(file="pRDA.CZ.geog.elev.noelev.pdf")
plot(pRDA.CZ.geog.noelev, main="pRDA CZ (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CZ.climate.noelev <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CZ.climate.noelev
Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

              Inertia Proportion Rank
Total          2.8016     1.0000     
Conditional    0.6999     0.2498    5
Constrained    0.3615     0.1290    2
Unconstrained  1.7403     0.6212   30
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2 
0.30843 0.05304 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
0.4198 0.1815 0.1446 0.1386 0.1213 0.0955 0.0853 0.0764 
(Showed only 8 of all 30 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CZ.climate.noelev)

$r.squared
[1] 0.129022

$adj.r.squared
[1] 0.1012994

anova(pRDA.CZ.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     2  0.36147 3.1156  0.001 ***
Residual 30  1.74031                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CZ.climate.noelev.pdf")
plot(pRDA.CZ.climate.noelev, main="pRDA CZ (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
(head(summary(pRDA.CZ.geog.noelev)))

Biplot scores for constraining variables

                    RDA1    RDA2     RDA3     RDA4    RDA5 PC1
shadow.days      -0.2035 -0.1860  0.08984  0.25995  0.7356   0
sol.rad.60d       0.4581 -0.6103 -0.42518 -0.34187 -0.3241   0
pcpt.60d         -0.3390 -0.7854 -0.31009 -0.09085 -0.2759   0
day10cm          -0.1142 -0.7304 -0.06059  0.17115  0.1017   0
temp.laying.date -0.6427  0.1863 -0.55871  0.38077 -0.1860   0

```

