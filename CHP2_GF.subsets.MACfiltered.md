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

#### CHS.VS
```
#Calculating MEM variables
#install.packages("tripack")
#install.packages("spacemakeR", repos="http://R-Forge.R-project.org")

library(spacemakeR)


env.data.CHS.VS <- read.csv("AllEnv.Data_CHS.VS_20161010.csv", header=T)
head(env.data.CHS.VS)
#extract x and y
xy <- env.data.CHS.VS[, c("long","lat")]

#install.packages("geosphere")
library(geosphere) #calculate a matrix of geographic distances
dxy <- distm(xy)
dxy <- as.dist(dxy)

#Function that returns the maximum distance of the minimum spanning tree based on a distance matrix.
th <- give.thresh(dxy)
#Function to compute neighborhood based on the minimum spanning tree. Returns an object of the class nb (see spdep package).
nb1 <- mst.nb(dxy)
wh1 <- which(as.matrix(dxy)==th,arr.ind=TRUE)
plot(nb1,xy,pch=20,cex=2,lty=3)
lines(xy[wh1[1,],1],xy[wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")
#thershold distance
th 
#[1] 41004.58
nb1
Neighbour list object:
Number of regions: 10 
Number of nonzero links: 18 
Percentage nonzero weights: 18 
Average number of links: 1.8 

#install.packages("spdep")
library(spdep)
#transform nb to listw (spdep package)
listw=nb2listw(nb1, glist=NULL, style="W", zero.policy=NULL)
#The can.be.simmed helper function checks whether a spatial weights object is similar to
#symmetric and can be so transformed to yield real eigenvalues or for Cholesky decomposition.
can.be.simmed(listw)
#[1] TRUE
 
#Function to compute Moran's eigenvectors of a listw object
#This functions compute eigenvector's of a doubly centered spatial weighting matrix. 
#Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
#scores=scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
#MEM.autocor: A string indicating if all MEMs must be returned or only those corresponding to positive or negative autocorrelation.
#Only positive correlations:
scores=scores.listw(listw, echo = FALSE, MEM.autocor = "positive")
	#listw not symmetric, (w+t(w)) used in the place of w 

#Function to compute and test Moran's I for eigenvectors of spatial weighting matrices. 
#This function tests Moran's I for each eigenvector of a spatial weighting matrix
test.scores(scores,listw,nsim=999)
         stat   pval
1 0.998929647 0.0010
2 0.750000000 0.0015
3 0.523833554 0.0300
4 0.008291494 0.3580


#3 significant MEM eigenfunctions with positive correlations. OBS. use first half (2 in our case)

write.table (scores$vectors[,1], "scores_MEM1.txt") 
write.table (scores$vectors[,2], "scores_MEM2.txt") 
```

## Great circle distance between all sampling points

```
##Geographic distance between all CH populations
##using package Rdist

library(fields)
#Import .csv with coordinates (done above)

setwd("/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017/CHS.VS")

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
bli <- as.data.frame(distance.matrix.CHS.VS)
head(bli)
colnames(bli) <- env.data.CHS.VS$site
rownames(bli) <- env.data.CHS.VS$site

bli[lower.tri(bli,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
bli=as.data.frame(as.table(as.matrix(bli)))  #Turn into a 3-column table
bli
bli=na.omit(bli)  #Get rid of the junk we flagged above
bli
colnames(bli)<-c("site1", "site2", "dist(km)")
head(bli)

bli2 <- bli[sort(bli$site2),]

head(bli2)


##write to csv
write.csv(bli, file="distance.CHS.VS.csv",row.names=F)
```

## Create a SNP input file 

I'm creating two input files.

1. Adaptive Loci (found by at least two methods)

2. Neutral loci (All - adaptive)


##### 1. Adaptive loci
```

#Filter the Adaptive loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.AdaptiveLoci.names --recode --recode-INFO-all --out CHS.VS.135.358

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.358.recode.vcf --plink --out CHS.VS.135.358.plink

plink --file CHS.VS.135.358.plink --recode --recodeA --out CHS.VS.135.358.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHall932.9608.PlinkCluster

And calculate MAF with Plink

```

plink --file CHS.VS.135.358.plink --within CHS.VS.PlinkCluster --freq --out CHS.VS.135.358

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.Adaptive.MAF <- read.table("CHS.VS.135.358.frq.strat", header=T)
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
write.csv(CHS.VS.Adaptive.MAF3, file="CHS.VS.135.358.Adaptive.MAF.csv")
```


##### 2. Neutral loci
```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --exclude CHS.VS.NeutralLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.Neutral.recode.vcf --plink --out CHS.VS.135.Neutral.plink

plink --file CHS.VS.135.Neutral.plink --recode --recodeA --out CHS.VS.135.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
Neutral.loci.names.map <- read.table("CHS.VS.135.Neutral.plink.map", header=F)
Neutral.loci.names <-  Neutral.loci.names.map$V2
Neutral.loci.names <- as.data.frame(Neutral.loci.names)
CHS.VS.1000.loci.names <- Neutral.loci.names[sample(nrow(Neutral.loci.names),1000),]
CHS.VS.1000.loci.names <- as.data.frame(CHS.VS.1000.loci.names)
summary(CHS.VS.1000.loci.names)
write.table(CHS.VS.1000.loci.names, "CHS.VS.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.VS.135.Neutral.recode.vcf --snps CHS.VS.1000.loci.names --recode --recode-INFO-all --out CHS.VS.1000NeutralLoci

vcftools --vcf CHS.VS.1000NeutralLoci.recode.vcf --plink --out CHS.VS.1000NeutralLoci.plink
plink --file CHS.VS.1000NeutralLoci.plink --recode --recodeA --out CHS.VS.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster

And calculate MAF with Plink

```

plink --file CHS.VS.1000NeutralLoci.plink --within CHS.VS.PlinkCluster --freq --out CHS.VS.1000NeutralLoci

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
write.csv(CHS.VS.Neutral.MAF3, file="CHS.VS.135.Neutral.Neutral.MAF.csv")
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



