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
1   0.95580479 0.001
2   0.93050370 0.001
3   0.91257634 0.001
4   0.64966475 0.001
5   0.25237401 0.041
6   0.34391370 0.011
7  -0.14166137 0.803
8   0.08943395 0.248
9  -0.11566399 0.709
10  0.11102298 0.218
11 -0.15455107 0.769
12 -0.32732604 0.977
13  0.16821893 0.104
14 -0.03065197 0.508
15 -0.63822650 0.999
16  0.34012456 0.005
17 -0.52370184 0.998
18 -0.60240388 0.999
19 -0.67917757 0.999
20 -0.53936405 0.998
21 -0.50854789 0.997

#6 significant PCNM eigenfunctions with positive correlations.

write.table (CHN.pcnm.scores$vectors[,1], "CHN.scores_PCNM1.txt") 
write.table (CHN.pcnm.scores$vectors[,2], "CHN.scores_PCNM2.txt") 
write.table (CHN.pcnm.scores$vectors[,3], "CHN.scores_PCNM3.txt") 
write.table (CHN.pcnm.scores$vectors[,4], "CHN.scores_PCNM4.txt") 
write.table (CHN.pcnm.scores$vectors[,5], "CHN.scores_PCNM5.txt") 
write.table (CHN.pcnm.scores$vectors[,6], "CHN.scores_PCNM6.txt") 

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


env.data.CHall <- rbind(env.data.CHN, env.data.CHS, env.data.CZ)
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
#[1] 42437.1
CHall.nb1
Neighbour list object:
Number of regions: 80 
Number of nonzero links: 158 
Percentage nonzero weights: 2.46875 
Average number of links: 1.975 

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
1   0.978870662 0.001*
2   0.979201633 0.001*
3   0.933694805 0.001*
4   0.917753544 0.001*
5   0.852246427 0.001*
6   0.900881326 0.001*
7   0.834207732 0.001*
8   0.856667583 0.001*
9   0.848587547 0.001*
10  0.645406423 0.001*
11  0.525168074 0.001*
12  0.166315768 0.043*
13  0.223339042 0.026*
14  0.014088782 0.406
15  0.036135417 0.355
16  0.009042856 0.410
17 -0.167323668 0.916
18 -0.279342795 0.994
19 -0.351123899 0.999
20 -0.394090816 0.999
21 -0.132939121 0.896
22 -0.007220108 0.502
23  0.061296277 0.253
24 -0.406238925 0.999
25 -0.069283508 0.731
26 -0.097663284 0.787
27 -0.269447987 0.992
28 -0.272271122 0.989
29 -0.456094441 0.999
30  0.176203280 0.029*
31  0.009238961 0.252
32 -0.060031799 0.948
33 -0.059769523 0.763
34  0.384102665 0.001*
35 -0.075602831 0.825
36 -0.001538547 0.621
37 -0.019890990 0.740
38 -0.831325686 0.999
39 -0.774231556 0.999
40 -0.215617762 0.968
41 -0.676525471 0.999
42 -0.265577634 0.982
43 -0.717997084 0.999
44 -0.428119796 0.997
45 -0.086623617 0.931

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
write.table (CHall.pcnm.scores$vectors[,30], "CHall.scores_PCNM14.txt") 
write.table (CHall.pcnm.scores$vectors[,34], "CHall.scores_PCNM15.txt") 
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
GenData <- CHallData[,10:ncol(CHallData)]
Climate.Data <- CHallData[,3:9]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CH932full <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s.

RDA.CH932full

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         430.9929     1.0000     
Constrained   172.3605     0.3999    7
Unconstrained 258.6324     0.6001   73
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
119.98  29.47   9.03   4.53   3.89   2.82   2.64 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
53.74 14.48  9.28  6.87  6.52   6.28  5.71  5.26 
(Showed only 8 of all 73 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CH932full)

$r.squared
[1] 0.3999149

$adj.r.squared
[1] 0.3423725



anova(RDA.CH932full)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   172.36 6.9499  0.001 ***
Residual 73   258.63                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CH932full.pdf")
plot(RDA.CH932full)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CH932.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

#head(summary(pRDA.CH940.geog))

pRDA.CH932.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

                Inertia Proportion Rank
Total         430.99291    1.00000     
Conditional   131.88102    0.30599    2
Constrained    40.47949    0.09392    5
Unconstrained 258.63241    0.60009   73
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
21.903  8.759  4.304  2.830  2.684 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
53.74 14.48  9.28  6.87  6.52  6.28  5.71  5.26 
(Showed only 8 of all 73 unconstrained eigenvalues)



RsquareAdj(pRDA.CH932.geog)

$r.squared
[1] 0.09392147

$adj.r.squared
[1] 0.0541741

anova(pRDA.CH932.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     5   40.479 2.2851  0.001 ***
Residual 73  258.632                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CH940 GeneticData


pdf(file="pRDA.CH932.geog.pdf")
plot(pRDA.CH932.geog, main="pRDA CH932 (geog partialled out)")
dev.off()
```


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CH932.climate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CH932.climate

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

                Inertia Proportion Rank
Total         430.99291    1.00000     
Conditional    38.46745    0.08925    5
Constrained   133.89305    0.31066    2
Unconstrained 258.63241    0.60009   73
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
109.09  24.80 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
53.74 14.48  9.28  6.87  6.52  6.28  5.71  5.26 
(Showed only 8 of all 73 unconstrained eigenvalues)


#head(summary(pRDA.CH932.climate))

RsquareAdj(pRDA.CH932.climate)

$r.squared
[1] 0.3106618

$adj.r.squared
[1] 0.3138359

anova(pRDA.CH932.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     2   133.89 18.896  0.001 ***
Residual 73   258.63                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CH932.climate.pdf")
plot(pRDA.CH932.climate, main="pRDA CH932 (climate partialled out)")
dev.off()
```


Find the most important loci associated with Climate
```
 summary(pRDA.CH932.geog)

Biplot scores for constraining variables

                    RDA1    RDA2    RDA3    RDA4    RDA5 PC1
shadow.days       0.3453 -0.1496 -0.3566  0.1332 -0.8445   0
sol.rad.60d      -0.2580  0.3047 -0.3726 -0.8022  0.2101   0
pcpt.60d          0.4786 -0.2806 -0.6136 -0.3718  0.2047   0
day10cm           0.3646  0.6714 -0.5720 -0.1345 -0.1135   0
temp.laying.date -0.1816 -0.3905 -0.7047  0.4556  0.1166   0


RDA1: pcpt.60d, day10cm, shadow.days

RDA2: day10cm, temp.laying.date
```



#### 2. CHN

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHNData <- read.csv("CHN229.MAF.csv", header=T)
GenData <- CHNData[,10:ncol(CHNData)]
Climate.Data <- CHNData[,3:9]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHNfull <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHNfull

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         191.7262     1.0000     
Constrained    92.3663     0.4818    7
Unconstrained  99.3599     0.5182   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
24.398 20.721 12.255 11.222  9.851  7.447  6.472 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
18.548 11.959 11.511  9.700  8.719  8.298  7.476  6.876  6.011  5.465  4.797  


RsquareAdj(RDA.CHNfull)

$r.squared
[1] 0.4817615

$adj.r.squared
[1] 0.1519733



anova(RDA.CHNfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   92.366 1.4608  0.001 ***
Residual 11   99.360                  
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

pRDA.CHNgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(lat + long), Climate.Data)

pRDA.CHNgeog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         191.7262     1.0000     
Conditional    39.8884     0.2080    2
Constrained    52.4779     0.2737    5
Unconstrained  99.3599     0.5182   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
14.300 11.402 11.117  8.817  6.842 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
18.548 11.959 11.511  9.700  8.719  8.298  7.476  6.876  6.011  5.465  4.797

RsquareAdj(pRDA.CHNgeog)

$r.squared
[1] 0.2737126

$adj.r.squared
[1] 0.0429183


anova(pRDA.CHNgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance     F Pr(>F)  
Model     5   52.478 1.162  0.099 .
Residual 11   99.360               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 not rejected: Climate does not explain a significant amount of CHN GeneticData





pdf(file="pRDA.CHNgeog.pdf")
plot(pRDA.CHNgeog, main="pRDA CHN(geog partialled out)")
dev.off()
```


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHNclimate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHNclimate

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         191.7262     1.0000     
Conditional    68.2651     0.3561    5
Constrained    24.1012     0.1257    2
Unconstrained  99.3599     0.5182   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
14.941  9.160 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
18.548 11.959 11.511  9.700  8.719  8.298  7.476  6.876  6.011  5.465  4.797  



RsquareAdj(pRDA.CHNclimate)

$r.squared
[1] 0.1257064

$adj.r.squared
[1] 0.04358942

anova(pRDA.CHNclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   24.101 1.3341  0.076 .
Residual 11   99.360                
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

                    RDA1     RDA2     RDA3      RDA4      RDA5 PC1
shadow.days      -0.6372  0.33047  0.33471  0.246755 -0.008364   0
sol.rad.60d      -0.2536 -0.76583  0.25466 -0.046557  0.289053   0
pcpt.60d         -0.4723 -0.51210 -0.01391  0.002262  0.396309   0
day10cm          -0.4603 -0.09573  0.23986 -0.301758  0.228056   0
temp.laying.date -0.7094 -0.22601 -0.22636 -0.319271 -0.466171   0


RDA1: temp.laying.date, shadow.days

RDA2: sol.rad.60d, pcpt.60d
```




#### 3. CHS

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHSData <- read.csv("CHS275.MAF.csv", header=T)
CHSData <- CHSData[complete.cases(CHSData),]  ##stba is removed here due to missing data
GenData <- CHSData[,10:ncol(CHSData)]
Climate.Data <- CHSData[,3:9]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHSfull <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHSfull

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         316.3975     1.0000     
Constrained   146.7810     0.4639    7
Unconstrained 169.6165     0.5361   16
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
67.07 18.33 16.49 13.72 12.46  9.73  8.98 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
23.457 21.232 17.695 14.813 10.934 10.191  9.838  9.364  8.827  7.752  7.579 
  PC12   PC13   PC14   PC15   PC16 
 6.942  6.050  5.689  4.708  4.545 


RsquareAdj(RDA.CHSfull)

$r.squared
[1] 0.4639132

$adj.r.squared
[1] 0.2293752



anova(RDA.CHSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)    
Model     7   146.78 1.978  0.001 ***
Residual 16   169.62                 
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

pRDA.CHSgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(lat + long), Climate.Data)

RsquareAdj(pRDA.CHSgeog)

$r.squared
[1] 0.1986423

$adj.r.squared
[1] 0.0340785


anova(pRDA.CHSgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5    62.85 1.1857  0.046 *
Residual 16   169.62                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CHS GeneticData

pRDA.CHSgeog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         316.3975     1.0000     
Conditional    83.9311     0.2653    2
Constrained    62.8499     0.1986    5
Unconstrained 169.6165     0.5361   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
16.495 13.851 13.555  9.934  9.015 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
23.457 21.232 17.695 14.813 10.934 10.191  9.838  9.364  8.827  7.752  7.579 
  PC12   PC13   PC14   PC15   PC16 
 6.942  6.050  5.689  4.708  4.545




pdf(file="pRDA.CHSgeog.pdf")
plot(pRDA.CHSgeog, main="pRDA CHS(geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHSclimate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHSclimate

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         316.3975     1.0000     
Conditional   109.4408     0.3459    5
Constrained    37.3402     0.1180    2
Unconstrained 169.6165     0.5361   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
21.927 15.413 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
23.457 21.232 17.695 14.813 10.934 10.191  9.838  9.364  8.827  7.752  7.579 
  PC12   PC13   PC14   PC15   PC16 
 6.942  6.050  5.689  4.708  4.545 

RsquareAdj(pRDA.CHSclimate)

$r.squared
[1] 0.1180168

$adj.r.squared
[1] 0.06517423

anova(pRDA.CHSclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     2    37.34 1.7612  0.001 ***
Residual 16   169.62                  
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

                     RDA1     RDA2     RDA3    RDA4      RDA5 PC1
shadow.days      -0.57138 -0.36197  0.55148  0.2124 -0.051183   0
sol.rad.60d      -0.02309  0.64415 -0.36900  0.3378 -0.427501   0
pcpt.60d         -0.32004  0.37876  0.08764 -0.2157  0.061851   0
day10cm          -0.12807  0.64213  0.24394  0.4782  0.006326   0
temp.laying.date -0.67551  0.01769 -0.12396  0.1134  0.674231   0

RDA1: temp.laying.date, shadow.days

RDA2: sol.rad.60d, day10cm
```





#### 4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,10:ncol(CHS.VS.Data)]
Climate.Data <- CHS.VS.Data[,3:9]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.VSfull <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.VSfull

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         228.9971     1.0000     
Constrained   196.2912     0.8572    7
Unconstrained  32.7059     0.1428    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
56.77 33.82 27.43 25.30 19.07 18.10 15.81 

Eigenvalues for unconstrained axes:
   PC1    PC2 
18.203 14.503 




RsquareAdj(RDA.CHS.VSfull)
$r.squared
[1] 0.8571775

$adj.r.squared
[1] 0.3572988


anova(RDA.CHS.VSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7  196.291 1.7148  0.001 ***
Residual  2   32.706                  
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

pRDA.CHS.VS.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

pRDA.CHS.VS.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         228.9971     1.0000     
Conditional    65.7180     0.2870    2
Constrained   130.5732     0.5702    5
Unconstrained  32.7059     0.1428    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
47.97 26.50 20.58 19.04 16.49 

Eigenvalues for unconstrained axes:
   PC1    PC2 
18.203 14.503 

RsquareAdj(pRDA.CHS.VS.geog)
$r.squared
[1] 0.5701957

$adj.r.squared
[1] 0.2740365


anova(pRDA.CHS.VS.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5  130.573 1.5969  0.087 .
Residual  2   32.706                
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

pRDA.CHS.VS.climate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.VS.climate

Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         228.9971     1.0000     
Conditional   135.0835     0.5899    5
Constrained    61.2076     0.2673    2
Unconstrained  32.7059     0.1428    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
42.77 18.44 

Eigenvalues for unconstrained axes:
   PC1    PC2 
18.203 14.503 


#head(summary(pRDA.CHS.VS.climate))

RsquareAdj(pRDA.CHS.VS.climate)

$r.squared
[1] 0.2672856

$adj.r.squared
[1] 0.2800421

anova(pRDA.CHS.VS.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2   61.208 1.8715  0.151
Residual  2   32.706   

H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.pdf")
plot(pRDA.CHS.VS.climate, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.VS.geog))

Biplot scores for constraining variables

                   RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days      0.2814 -0.1945  0.18520 0.16543 -0.10742   0
sol.rad.60d      0.2770  0.2697 -0.02120 0.09208  0.51938   0
pcpt.60d         0.5170  0.8100  0.13979 0.21681 -0.01227   0
day10cm          0.5726  0.5415 -0.09465 0.44351  0.15408   0
temp.laying.date 0.4772  0.5430 -0.13052 0.06771 -0.65840   0

RDA1: day10cm, pcpt.60d, temp.laying.date

RDA2: pcpt.60d, temp.laying.date, day10cm
```




#### 5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI140.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,10:ncol(CHS.TI.Data)]
Climate.Data <- CHS.TI.Data[,3:9]

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.TIfull <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.TIfull

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         266.1176     1.0000     
Constrained   164.1611     0.6169    7
Unconstrained 101.9564     0.3831    6
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
39.28 35.46 26.50 24.45 14.05 13.03 11.39 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
27.919 21.777 17.034 15.974 10.499  8.754 


RsquareAdj(RDA.CHS.TIfull)
$r.squared
[1] 0.6168745

$adj.r.squared
[1] 0.1698948

anova(RDA.CHS.TIfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     7   164.16 1.3801  0.005 **
Residual  6   101.96                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.TIfull.pdf")
plot(RDA.CHS.TIfull)
dev.off()
```


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.TI.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

pRDA.CHS.TI.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         266.1176     1.0000     
Conditional    64.9251     0.2440    2
Constrained    99.2361     0.3729    5
Unconstrained 101.9564     0.3831    6
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
33.31 25.74 15.03 13.50 11.66 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
27.919 21.777 17.034 15.974 10.499  8.754 


RsquareAdj(pRDA.CHS.TI.geog)
$r.squared
[1] 0.3729032

$adj.r.squared
[1] 0.06338319



anova(pRDA.CHS.TI.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance     F Pr(>F)
Model     5   99.236 1.168  0.181
Residual  6  101.956               


H0 not rejected: Climate alone does not explain CHS.TI GeneticData



pdf(file="pRDA.CHS.TI.geog.pdf")
plot(pRDA.CHS.TI.geog, main="pRDA CHS.TI (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.TI.climate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.TI.climate

CCall: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         266.1176     1.0000     
Conditional   110.6767     0.4159    5
Constrained    53.4844     0.2010    2
Unconstrained 101.9564     0.3831    6
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
29.004 24.480 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6 
27.919 21.777 17.034 15.974 10.499  8.754 

RsquareAdj(pRDA.CHS.TI.climate)
$r.squared
[1] 0.2009804

$adj.r.squared
[1] 0.1190668

#head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   53.484 1.5737  0.048 *
Residual  6  101.956                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

H0: Rejected -> Geography explains a sig proportion of GenData

pdf(file="pRDA.CHS.TI.climate.pdf")
plot(pRDA.CHS.TI.climate, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```





Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.TI.geog))

Biplot scores for constraining variables

                    RDA1    RDA2   RDA3    RDA4     RDA5 PC1
shadow.days       0.6023  0.4185 0.1044  0.3742 -0.46578   0
sol.rad.60d      -0.2516 -0.8185 0.4245  0.2192  0.11839   0
pcpt.60d          0.3800 -0.1492 0.6673 -0.1839  0.04035   0
day10cm           0.1219 -0.5413 0.1967  0.3797 -0.04461   0
temp.laying.date  0.2237  0.5442 0.2352  0.5215  0.37699   0


RDA1: shadow.days

RDA2: sol.rad.60d, day10cm, temp.laying.date
```



#### 6. CZ

```
library(vegan)

CZData <- read.csv("CZ404.MAF.csv", header=T)
CZData <- CZData[complete.cases(CZData),]  ##stba is removed here due to missing data
GenData <- CZData[,10:ncol(CZData)]
Climate.Data <- CZData[,3:9]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CZfull <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. 

RDA.CZfull

Call: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d +
pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)

               Inertia Proportion Rank
Total         310.3385     1.0000     
Constrained   114.6941     0.3696    7
Unconstrained 195.6444     0.6304   30
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
64.16 15.48 12.01  7.43  6.44  4.91  4.26 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
24.782 18.471 11.336 10.941  9.803  9.155  8.394  7.905 
(Showed only 8 of all 30 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CZfull)

$r.squared
[1] 0.3695774

$adj.r.squared
[1] 0.2224788



anova(RDA.CZfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   114.69 2.5124  0.001 ***
Residual 30   195.64                  
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

pRDA.CZ.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), Climate.Data)

#head(summary(pRDA.CZ.geog.noelev))

pRDA.CZ.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         310.3385     1.0000     
Conditional    71.5225     0.2305    2
Constrained    43.1716     0.1391    5
Unconstrained 195.6444     0.6304   30
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
18.747  8.108  6.875  4.996  4.445 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
24.782 18.471 11.336 10.941  9.803  9.155  8.394  7.905 
(Showed only 8 of all 30 unconstrained eigenvalues)


RsquareAdj(pRDA.CZ.geog)

$r.squared
[1] 0.1391112

$adj.r.squared
[1] 0.03598595

anova(pRDA.CZ.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     5   43.172 1.324  0.004 **
Residual 30  195.644                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CZ GeneticData


pdf(file="pRDA.CZ.geog.pdf")
plot(pRDA.CZ.geog, main="pRDA CZ (geog partialled out)")
dev.off()
```



```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CZ.climate <- rda(GenData~lat+long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CZ.climate
Call: rda(formula = GenData ~ lat + long + Condition(shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         310.3385     1.0000     
Conditional    70.3460     0.2267    5
Constrained    44.3481     0.1429    2
Unconstrained 195.6444     0.6304   30
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
36.13  8.22 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
24.782 18.471 11.336 10.941  9.803  9.155  8.394  7.905 
(Showed only 8 of all 30 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CZ.climate)

$r.squared
[1] 0.1429024

$adj.r.squared
[1] 0.1166358

anova(pRDA.CZ.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     2   44.348 3.4002  0.001 ***
Residual 30  195.644                  
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

                    RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days      0.31508 -0.1141  0.04581  0.7586  0.05346   0
sol.rad.60d      0.01513  0.7493 -0.30059 -0.2955 -0.49557   0
pcpt.60d         0.73838  0.3458 -0.33076 -0.2685 -0.25699   0
day10cm          0.60361  0.4257 -0.09102  0.1824 -0.05277   0
temp.laying.date 0.30882 -0.4304 -0.70349 -0.1659  0.36929   0

RDA1: pcpt.60d, day10cm

RDA2: sol.rad.60d
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

