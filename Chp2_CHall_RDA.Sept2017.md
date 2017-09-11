# RDA analysis of CH dataset

RDA: Effect of Climate & Geog on genetic variation


I can use RDA to partition the effects of climate and geog on genetic variation. 

Explanation and example: http://www.davidzeleny.net/anadat-r/doku.php/en:rda_examples

And Victoria's paper: http://www.amjbot.org/content/103/1/33.full.pdf+html

Paul Gugger tutorial: http://pgugger.al.umces.edu/assets/redundancy-analysis-for-landscape-genetics.pdf


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA

### Data:

I'm using the full CH dataset. I have to run the analysis seperately for all the different partitions of the data to investigate adaptation
at different geographic scales. 

#### 1. CHall.n940.loci9608

#### 2. CHN.n229.loci9608

#### 3. CHS.n283.loci9608

#### 4. CHS.VS.n135.loci9608

#### 5. CHS.TI.n148.loci9608

#### 6. CZ.n404.loci9608


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

plink --file CHall.940.9608.plink --freq --within CH940.cluster.pop --out CH940.pop


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
plink --file CHN.229.9608.plink --freq --within CHN229.cluster.pop --out CHN229.pop


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
plink --file CHS.283.9608.plink --freq --within CHS283.cluster.pop --out CHS283.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS283.MAF <- read.table("CHS283.pop.frq.strat", header=T)
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
plink --file CHS.VS.135.plink --freq --within CHS.VS.cluster.pop --out CHS.VS135.pop


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
plink --file CHS.TI.148.plink --freq --within CHS.TI.148.cluster.pop --out CHS.TI148.pop


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

CHS.TI148.MAF <- read.table("CHS.TI148.pop.frq.strat", header=T)
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
plink --file CZ.404.9608.plink --freq --within CZ404.cluster.pop --out CZ404.pop


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




#### 2. Geographic coordinates

paste coordinates into the MAF.csv file from Env.Data.all_20161025.csv

located: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2

I will include elevation for geographic information for the CH data. 

i.e. lat, long, elev


#### 3. Climate variables

```
1. solar.rad.60d (amt of solar radiation 60 days after egg laying)

2. pcpt.60d (precipitation in the 60days after egg laying)

3. day10cm (calender day on which only 10cm of snow remains)

4. temp.laying.date

5. shadow.days
```

paste the 5 variables from /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/Env.Data.all_20161025.csv

Note that stba has missing data here. I could remove this population from the rest of these analyses. 


### Run RDA

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
GenData <- CHNData[,11:9618]
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
Total         218.3893     1.0000     
Constrained   104.8478     0.4801    7
Unconstrained 113.5415     0.5199   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
27.471 23.112 14.411 12.965 10.979  8.698  7.212 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
21.112 13.473 13.220 11.204  9.805  9.420  8.675  7.677  6.944  6.368  5.643 


RsquareAdj(RDA.CHNfull.noelev)

$r.squared
[1] 0.4800959

$adj.r.squared
[1] 0.1492478



anova(RDA.CHNfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   104.85 1.4511  0.001 ***
Residual 11   113.54                  
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
[1] 0.2751348

$adj.r.squared
[1] 0.0436666


anova(pRDA.CHNgeog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   60.086 1.1642  0.077 .
Residual 11  113.541                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1  


H0 not rejected: Climate does not explain a significant amount of CHN GeneticData

pRDA.CHNgeog.noelev

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         218.3893     1.0000     
Conditional    44.7613     0.2050    2
Constrained    60.0865     0.2751    5
Unconstrained 113.5415     0.5199   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
16.393 13.400 12.578 10.030  7.685 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
21.112 13.473 13.220 11.204  9.805  9.420  8.675  7.677  6.944  6.368  5.643 



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
Total         218.3893     1.0000     
Conditional    77.4433     0.3546    5
Constrained    27.4045     0.1255    2
Unconstrained 113.5415     0.5199   11
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
17.053 10.351 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
21.112 13.473 13.220 11.204  9.805  9.420  8.675  7.677  6.944  6.368  5.643 



RsquareAdj(pRDA.CHNclimate.noelev)

$r.squared
[1] 0.1254847

$adj.r.squared
[1] 0.0428631

anova(pRDA.CHNclimate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   27.405 1.3275  0.076 .
Residual 11  113.541                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: not rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHNclimate.pdf")
plot(pRDA.CHNclimate.noelev, main="pRDA (climate partialled out)")
dev.off()
```



Find the most important loci associated with Climate


```
summary(pRDA.CHNgeog.noelev)

Biplot scores for constraining variables

                    RDA1     RDA2     RDA3      RDA4      RDA5 PC1
shadow.days      -0.6486 -0.41034 0.151712  0.275719 -0.002203   0
sol.rad.60d      -0.1709  0.51585 0.641601 -0.062784 -0.300858   0
pcpt.60d         -0.4192  0.46447 0.294863 -0.009642 -0.405071   0
day10cm          -0.4323 -0.04717 0.288758 -0.284975 -0.260702   0
temp.laying.date -0.7060  0.31165 0.002971 -0.373020  0.436136   0

RDA1: shadow.days, temp.laying.date

RDA2: shadow.days, sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHNgeog.elev <- summary(pRDA.CHNgeog.noelev)

CHN.RDAloadings.noelev <- sum.pRDA.CHNgeog.noelev$species  ##write the loadings to data.frame
CHN.RDAloadings.noelev <- as.data.frame(CHN.RDAloadings.noelev) #change to df
names(CHN.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHN.RDAloadings.noelev$RDA1)
hist(CHN.RDAloadings.noelev$RDA2)
hist(CHN.RDAloadings.noelev$RDA3)
hist(CHN.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CHN.RDAloadings.sorted.noelev <- CHN.RDAloadings.noelev[order(-abs(CHN.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHN.RDA1.outliers.1percent.noelev <- CHN.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CHN.RDAloadings.sorted.noelev <- CHN.RDAloadings.noelev[order(-abs(CHN.RDAloadings.noelev$RDA2)),]  #same for RDA2
CHN.RDA2.outliers.1percent.noelev <- CHN.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CHN.RDA1.outliers.1percent.noelev, "CHN.noelev.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHN.RDA2.outliers.1percent.noelev, "CHN.noelev.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```


#### 3. CHS

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHSData <- read.csv("CHS283.MAF.csv", header=T)
CHSData <- CHSData[complete.cases(CHSData),]  ##stba is removed here due to missing data
GenData <- CHSData[,11:9618]
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
Total         345.1880     1.0000     
Constrained   157.5261     0.4563    7
Unconstrained 187.6619     0.5437   16
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
70.21 19.86 17.75 15.34 13.59 10.80  9.98 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
26.236 23.414 19.410 16.216 12.249 11.444 10.847 10.146  9.667  8.944  8.349 
  PC12   PC13   PC14   PC15   PC16 
 7.614  6.633  6.187  5.277  5.029 


RsquareAdj(RDA.CHSfull.noelev)

$r.squared
[1] 0.4563487

$adj.r.squared
[1] 0.2185012



anova(RDA.CHSfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   157.53 1.9187  0.001 ***
Residual 16   187.66                  
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
[1] 0.2001038

$adj.r.squared
[1] 0.03309021


anova(pRDA.CHSgeog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   69.073 1.1778  0.038 *
Residual 16  187.662                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CHS GeneticData

pRDA.CHSgeog.noelev
Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long), data =
Climate.Data)

               Inertia Proportion Rank
Total         345.1880     1.0000     
Conditional    88.4526     0.2562    2
Constrained    69.0734     0.2001    5
Unconstrained 187.6619     0.5437   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
17.766 15.449 14.844 11.016  9.999 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
26.236 23.414 19.410 16.216 12.249 11.444 10.847 10.146  9.667  8.944  8.349 
  PC12   PC13   PC14   PC15   PC16 
 7.614  6.633  6.187  5.277  5.029 




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
Total         345.1880     1.0000     
Conditional   117.0711     0.3392    5
Constrained    40.4550     0.1172    2
Unconstrained 187.6619     0.5437   16
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2 
23.725 16.730 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
26.236 23.414 19.410 16.216 12.249 11.444 10.847 10.146  9.667  8.944  8.349 
  PC12   PC13   PC14   PC15   PC16 
 7.614  6.633  6.187  5.277  5.029 

RsquareAdj(pRDA.CHSclimate.noelev)

$r.squared
[1] 0.1171971

$adj.r.squared
[1] 0.06291862

anova(pRDA.CHSclimate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     2   40.455 1.7246  0.002 **
Residual 16  187.662                 
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

                     RDA1    RDA2     RDA3     RDA4     RDA5 PC1
shadow.days      -0.52558 -0.3297  0.60984 -0.22506 -0.04276   0
sol.rad.60d      -0.08959  0.6056 -0.45177 -0.36674 -0.36762   0
pcpt.60d         -0.34804  0.3543  0.05319  0.22785  0.05012   0
day10cm          -0.17525  0.6628  0.16698 -0.46269  0.07053   0
temp.laying.date -0.67169 -0.0683 -0.09052 -0.04383  0.68811   0

RDA1: shadow.days, temp.laying.date

RDA2: sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHSgeog.noelev <- summary(pRDA.CHSgeog.noelev)

CHS.RDAloadings.noelev <- sum.pRDA.CHSgeog.noelev$species  ##write the loadings to data.frame
CHS.RDAloadings.noelev <- as.data.frame(CHS.RDAloadings.noelev) #change to df
names(CHS.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.RDAloadings.noelev$RDA1)
hist(CHS.RDAloadings.noelev$RDA2)
hist(CHS.RDAloadings.noelev$RDA3)
hist(CHS.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.RDAloadings.sorted.noelev <- CHS.RDAloadings.noelev[order(-abs(CHS.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.RDA1.outliers.1percent.noelev <- CHS.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CHS.RDAloadings.sorted.noelev <- CHS.RDAloadings.noelev[order(-abs(CHS.RDAloadings.noelev$RDA2)),]  #same for RDA2
CHS.RDA2.outliers.1percent.noelev <- CHS.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.RDA1.outliers.1percent.noelev, "CHS.noelev.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.RDA2.outliers.1percent.noelev, "CHS.noelev.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```



#### 4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,11:9618]
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
Total         252.5933     1.0000     
Constrained   216.1574     0.8558    7
Unconstrained  36.4359     0.1442    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
59.66 37.72 30.56 28.42 22.23 19.96 17.62 

Eigenvalues for unconstrained axes:
   PC1    PC2 
20.389 16.047 




RsquareAdj(RDA.CHS.VSfull.noelev)
$r.squared
[1] 0.8557528

$adj.r.squared
[1] 0.3508874


anova(RDA.CHS.VSfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     7  216.157 1.695  0.003 **
Residual  2   36.436                
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
Total         252.5933     1.0000     
Conditional    73.3697     0.2905    2
Constrained   142.7877     0.5653    5
Unconstrained  36.4359     0.1442    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
50.79 29.32 23.72 20.94 18.02 

Eigenvalues for unconstrained axes:
   PC1    PC2 
20.389 16.047 

RsquareAdj(pRDA.CHS.VS.geog.noelev)
$r.squared
[1] 0.565287

$adj.r.squared
[1] 0.2631457


anova(pRDA.CHS.VS.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5  142.788 1.5676  0.076 .
Residual  2   36.436                
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
Total         252.5933     1.0000     
Conditional   149.6384     0.5924    5
Constrained    66.5190     0.2633    2
Unconstrained  36.4359     0.1442    2
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
45.85 20.67 

Eigenvalues for unconstrained axes:
   PC1    PC2 
20.389 16.047


#head(summary(pRDA.CHS.VS.climate))

RsquareAdj(pRDA.CHS.VS.climate.noelev)

$r.squared
[1] 0.2633442

$adj.r.squared
[1] 0.2679682

anova(pRDA.CHS.VS.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2   66.519 1.8256  0.173
Residual  2   36.436              

H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.noelev.pdf")
plot(pRDA.CHS.VS.climate.noelev, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.VS.geog.noelev))

Biplot scores for constraining variables

                   RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days      0.2818 -0.2034  0.17278 0.17170 -0.10070   0
sol.rad.60d      0.2795  0.2831 -0.05314 0.11498  0.50381   0
pcpt.60d         0.5049  0.8066  0.18660 0.22142 -0.02531   0
day10cm          0.5698  0.5592 -0.09371 0.43904  0.10819   0
temp.laying.date 0.4639  0.5463 -0.03048 0.02092 -0.68020   0

RDA1: day10cm, pcpt.60d

RDA2: day10cm, pcpt.60d, temp.laying.date
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHS.VS.geog.noelev <- summary(pRDA.CHS.VS.geog.noelev)

CHS.VS.RDAloadings.noelev <- sum.pRDA.CHS.VS.geog.noelev$species  ##write the loadings to data.frame
CHS.VS.RDAloadings.noelev <- as.data.frame(CHS.VS.RDAloadings.noelev) #change to df
names(CHS.VS.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.VS.RDAloadings.noelev$RDA1)
hist(CHS.VS.RDAloadings.noelev$RDA2)
hist(CHS.VS.RDAloadings.noelev$RDA3)
hist(CHS.VS.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.VS.RDAloadings.sorted.noelev <- CHS.VS.RDAloadings.noelev[order(-abs(CHS.VS.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.VS.RDA1.outliers.1percent.noelev <- CHS.VS.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CHS.VS.RDAloadings.sorted.noelev <- CHS.VS.RDAloadings.noelev[order(-abs(CHS.VS.RDAloadings.noelev$RDA2)),]  #same for RDA2
CHS.VS.RDA2.outliers.1percent.noelev <- CHS.VS.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.VS.RDA1.outliers.1percent.noelev, "CHS.noelev.VS.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.VS.RDA2.outliers.1percent.noelev, "CHS.noelev.VS.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```


#### 5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI148.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,11:9618]
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
Total         303.4368     1.0000     
Constrained   187.6470     0.6184    7
Unconstrained 115.7899     0.3816    6
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7 
45.98 40.71 30.82 27.58 15.63 14.07 12.86 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6 
33.07 25.14 18.53 17.92 11.66  9.48 


RsquareAdj(RDA.CHS.TIfull.noelev)
$r.squared
[1] 0.6184054

$adj.r.squared
[1] 0.1732117

anova(RDA.CHS.TIfull.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     7   187.65 1.3891  0.003 **
Residual  6   115.79                 
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
Total         303.4368     1.0000     
Conditional    72.1685     0.2378    2
Constrained   115.4785     0.3806    5
Unconstrained 115.7899     0.3816    6
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
40.03 30.90 16.81 14.62 13.12 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6 
33.07 25.14 18.53 17.92 11.66  9.48 


RsquareAdj(pRDA.CHS.TI.geog.noelev)
$r.squared
[1] 0.3805684

$adj.r.squared
[1] 0.07394976



anova(pRDA.CHS.TI.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5   115.48 1.1968  0.189
Residual  6   115.79              


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
Total         303.4368     1.0000     
Conditional   124.9119     0.4117    5
Constrained    62.7350     0.2067    2
Unconstrained 115.7899     0.3816    6
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
35.70 27.03 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6 
33.07 25.14 18.53 17.92 11.66  9.48 

RsquareAdj(pRDA.CHS.TI.climate.noelev)
$r.squared
[1] 0.2067483

$adj.r.squared
[1] 0.1292689

#head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   62.735 1.6254  0.049 *
Residual  6  115.790                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

H0: Rejected -> Geography explains a sig proportion of GenData

pdf(file="pRDA.CHS.TI.climate.noelev.pdf")
plot(pRDA.CHS.TI.climate.noelev, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```





Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.TI.geog.noelev))

Biplot scores for constraining variables

                    RDA1     RDA2    RDA3    RDA4     RDA5 PC1
shadow.days      0.43153 -0.60200 -0.1351 -0.3388  0.47339   0
sol.rad.60d      0.01517  0.86277 -0.4238 -0.2111 -0.08190   0
pcpt.60d         0.40842  0.04162 -0.6537  0.2248 -0.03140   0
day10cm          0.27926  0.46930 -0.2147 -0.3778  0.07448   0
temp.laying.date 0.02548 -0.58652 -0.3001 -0.5038 -0.35643   0


RDA1: shadow.days, pcpt.60d

RDA2: sol.rad.60d, shadow.days, temp.laying.date
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHS.TI.geog.noelev <- summary(pRDA.CHS.TI.geog.noelev)

CHS.TI.RDAloadings.noelev <- sum.pRDA.CHS.TI.geog.noelev$species  ##write the loadings to data.frame
CHS.TI.RDAloadings.noelev <- as.data.frame(CHS.TI.RDAloadings.noelev) #change to df
names(CHS.TI.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.TI.RDAloadings.noelev$RDA1)
hist(CHS.TI.RDAloadings.noelev$RDA2)
hist(CHS.TI.RDAloadings.noelev$RDA3)
hist(CHS.TI.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.TI.RDAloadings.sorted.noelev <- CHS.TI.RDAloadings.noelev[order(-abs(CHS.TI.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.TI.RDA1.outliers.1percent.noelev <- CHS.TI.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CHS.TI.RDAloadings.sorted.noelev <- CHS.TI.RDAloadings.noelev[order(-abs(CHS.TI.RDAloadings.noelev$RDA2)),]  #same for RDA2
CHS.TI.RDA2.outliers.1percent.noelev <- CHS.TI.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.TI.RDA1.outliers.1percent.noelev, "CHS.noelev.TI.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.TI.RDA2.outliers.1percent.noelev, "CHS.noelev.TI.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

#### 6. CZ

```
library(vegan)

CZData <- read.csv("CZ404.MAF.csv", header=T)
CZData <- CZData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CZData[,11:9618]
Climate.Data <- CZData[,3:10]
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CZfull.noelev <- rda(GenData ~ lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. 

RDA.CZfull

Call: rda(formula = GenData ~ elev + lat + long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Constrained   131.7473     0.3949    8
Unconstrained 201.8402     0.6051   29
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7  RDA8 
69.93 16.77 13.70  8.03  7.53  6.23  4.88  4.69 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.022 17.900 11.864 11.718 10.560  9.379  8.851  7.973 
(Showed only 8 of all 29 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CZfull.noelev)

$r.squared
[1] 0.3655938

$adj.r.squared
[1] 0.2175657



anova(RDA.CZfull.noelev)

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


RsquareAdj(pRDA.CZ.geog.noelev)

$r.squared
[1] 0.1388576

$adj.r.squared
[1] 0.03501604

anova(pRDA.CZ.geog.noelev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     5   46.321 1.3133  0.005 **
Residual 30  211.630                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CZ GeneticData


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


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CZ.climate.noelev)

$r.squared
[1] 0.140278

$adj.r.squared
[1] 0.1132943

anova(pRDA.CZ.climate.noelev)

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

pdf(file="pRDA.CZ.climate.noelev.pdf")
plot(pRDA.CZ.climate.noelev, main="pRDA CZ (climate partialled out)")
dev.off()
```




Find the most important loci associated with Climate
```
(head(summary(pRDA.CZ.geog.noelev)))

Biplot scores for constraining variables

                       RDA1    RDA2     RDA3    RDA4     RDA5 PC1
shadow.days       3.189e-01 -0.1027  0.03954  0.7566  0.08027   0
sol.rad.60d      -8.736e-05  0.7303 -0.32012 -0.2854 -0.51741   0
pcpt.60d          7.321e-01  0.3476 -0.33359 -0.2673 -0.26957   0
day10cm           5.950e-01  0.4360 -0.10138  0.1803 -0.05493   0
temp.laying.date  3.224e-01 -0.4331 -0.68384 -0.1886  0.38036   0

RDA1: sol.rad.60d, pcpt.60d

RDA2: sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CZ.geog.noelev <- summary(pRDA.CZ.geog.noelev)

CZ.RDAloadings.noelev <- sum.pRDA.CZ.geog.noelev$species  ##write the loadings to data.frame
CZ.RDAloadings.noelev <- as.data.frame(CZ.RDAloadings.noelev) #change to df
names(CZ.RDAloadings.noelev)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CZ.RDAloadings.noelev$RDA1)
hist(CZ.RDAloadings.noelev$RDA2)
hist(CZ.RDAloadings.noelev$RDA3)
hist(CZ.RDAloadings.noelev$RDA4)

##select the top 100 loci (~1%) for each RDA

CZ.RDAloadings.sorted.noelev <- CZ.RDAloadings.noelev[order(-abs(CZ.RDAloadings.noelev$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CZ.RDA1.outliers.1percent.noelev <- CZ.RDAloadings.sorted.noelev[1:100, 1:2]  ##write to df

CZ.RDAloadings.sorted.noelev <- CZ.RDAloadings.noelev[order(-abs(CZ.RDAloadings.noelev$RDA2)),]  #same for RDA2
CZ.RDA2.outliers.1percent.noelev <- CZ.RDAloadings.sorted.noelev[1:100, c(1:2)]  

#write both to file. 
write.table(CZ.RDA1.outliers.1percent.noelev, "CZ.noelev.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CZ.RDA2.outliers.1percent.noelev, "CZ.noelev.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```
