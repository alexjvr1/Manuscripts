#RDA analysis of CH dataset

RDA: Effect of Climate & Geog on genetic variation


I can use RDA to partition the effects of climate and geog on genetic variation. 

Explanation and example: http://www.davidzeleny.net/anadat-r/doku.php/en:rda_examples

And Victoria's paper: http://www.amjbot.org/content/103/1/33.full.pdf+html

Paul Gugger tutorial: http://pgugger.al.umces.edu/assets/redundancy-analysis-for-landscape-genetics.pdf

###Data:

I'm using the full CH dataset. I have to run the analysis seperately for all the different partitions of the data to investigate adaptation
at different geographic scales. 

####1. CHall.n940.loci9608

####2. CHN.n229.loci9608

####3. CHS.n283.loci9608

####4. CHS.VS.n135.loci9608

####5. CHS.TI.n148.loci9608

####6. CZ.n404.loci9608


###Input files:

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




####2. Geographic coordinates

paste coordinates into the MAF.csv file from Env.Data.all_20161025.csv

located: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2

I will include elevation for geographic information for the CH data. 

i.e. lat, long, elev


####3. Climate variables

```
1. solar.rad.60d (amt of solar radiation 60 days after egg laying)

2. pcpt.60d (precipitation in the 60days after egg laying)

3. day10cm (calender day on which only 10cm of snow remains)

4. temp.laying.date

5. shadow.days
```

paste the 5 variables from /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/Env.Data.all_20161025.csv

Note that stba has missing data here. I could remove this population from the rest of these analyses. 


###Run RDA

####1. CHall

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHallData <- read.csv("CH940.MAF.csv", header=T)
CHallData <- CHallData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CHallData[,11:9618]
Climate.Data <- CHallData[,3:7]
Climate.Data$Lat <- CHallData$lat
Climate.Data$Long <- CHallData$long
Climate.Data$elev <- CHallData$elev
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CH940full <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. I got an error (Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd') when elev was after Long in the formula

RDA.CH940full

Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         461.1563     1.0000     
Constrained   184.0588     0.3991    8
Unconstrained 601.2492     1.3038   70
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8 
123.58  31.44  10.17   5.02   4.71   3.36   3.09   2.70 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 20.48 19.92 19.28 18.48 17.89 17.29 17.19 
(Showed only 8 of all 70 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CH940full)

$r.squared
[1] 0.3991245

$adj.r.squared
[1] 0.3323605



anova(RDA.CH940full)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     8   184.06 2.7551  0.001 ***
Residual 72   601.25                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CH940full.elev.pdf")
plot(RDA.CH940full)
dev.off()
```

![alt_txt][CH940.full.fig1]
[CH940.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23798623/c56ded96-05a5-11e7-8085-70fb6aad317d.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CH940.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long +elev), Climate.Data)

#head(summary(pRDA.CH940.geog))

pRDA.CH940.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

                Inertia Proportion Rank
Total         461.15630    1.00000     
Conditional   144.57238    0.31350    3
Constrained    39.48639    0.08562    5
Unconstrained 277.09753    0.60088   72
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
22.377  6.884  4.001  3.336  2.888 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 15.51 10.10  7.53  7.03  6.97  6.33  5.72 
(Showed only 8 of all 72 unconstrained eigenvalues)


RsquareAdj(pRDA.CH940.geog)

$r.squared
[1] 0.08562474

$adj.r.squared
[1] 0.04560755

anova(pRDA.CH940.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance     F Pr(>F)    
Model     5   39.486 2.052  0.001 ***
Residual 72  277.098                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CH940 GeneticData


pdf(file="pRDA.CH940.geog.elev.pdf")
plot(pRDA.CH940.geog, main="pRDA CH940 (geog partialled out)")
dev.off()
```


![alt_txt][CH940.Fig2]
[CH940.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23798791/6cee94bc-05a6-11e7-91e8-587383b3428c.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CH940.climate <- rda(GenData~Lat+Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CH940.climate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

                Inertia Proportion Rank
Total         461.15630    1.00000     
Conditional    41.31089    0.08958    5
Constrained   142.74788    0.30954    3
Unconstrained 277.09753    0.60088   72
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3 
112.42  26.51   3.82 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 15.51 10.10  7.53  7.03  6.97  6.33  5.72 
(Showed only 8 of all 72 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CH940.climate)

$r.squared
[1] 0.3095434

$adj.r.squared
[1] 0.303474

anova(pRDA.CH940.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     3   142.75 12.364  0.001 ***
Residual 72   277.10                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1         

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CH940.climate.elev.pdf")
plot(pRDA.CH940.climate, main="pRDA CH940 (climate partialled out)")
dev.off()
```

![alt_txt][CH940.Fig3]
[CH940.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23798964/19945ab2-05a7-11e7-8465-30d221431940.png



Find the most important loci associated with Climate
```
 summary(pRDA.CH940.geog)

Importance of components:
                          RDA1    RDA2    RDA3    RDA4    RDA5     PC1     PC2
Eigenvalue            22.37711 6.88381 4.00051 3.33649 2.88847 55.6741 15.5130
Proportion Explained   0.07068 0.02174 0.01264 0.01054 0.00912  0.1759  0.0490
Cumulative Proportion  0.07068 0.09243 0.10506 0.11560 0.12473  0.3006  0.3496

                    RDA1     RDA2       RDA3    RDA4     RDA5 PC1
shadow.days       0.3497 -0.21845  0.0544813 -0.6602 -0.62512   0
sol.rad.60d      -0.3521 -0.25972 -0.0752741  0.5668 -0.26445   0
pcpt.60d          0.4492 -0.70416  0.0004636  0.1741 -0.01842   0
day10cm           0.2255  0.04426  0.3289740 -0.1241 -0.23575   0
temp.laying.date -0.1623 -0.67753  0.3294060 -0.4541  0.32005   0

RDA1: pcpt.60d, sol.rad.60d, shadow.days

RDA2: pcpt.60d, temp.laying.date
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```



####2. CHN

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHNData <- read.csv("CHN229.MAF.csv", header=T)
GenData <- CHNData[,10:9617]
Climate.Data <- CHNData[,3:7]
Climate.Data$Lat <- CHNData$lat
Climate.Data$Long <- CHNData$long
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHNfull <- rda(GenData ~ Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHNfull
Call: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d +
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


RsquareAdj(RDA.CHNfull)


anova(RDA.CHNfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   104.85 1.4511  0.001 ***
Residual 11   113.54                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHNfull.pdf")
plot(RDA.CHNfull)
dev.off()
```

![alt_txt][CHN.Fig1]
[CHN.Fig1]:https://cloud.githubusercontent.com/assets/12142475/23724999/2bffd0ea-044f-11e7-9924-2cedcc7c2b08.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHNgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(Lat + Long), Climate.Data)

RsquareAdj(RDA.CHNgeog)

anova(pRDA.CHNgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5   60.086 1.1642  0.081 .
Residual 11  113.541                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 not rejected: Climate does not explain CHN GeneticData

pRDA.CHNgeog
Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long), data =
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
plot(pRDA.CHNgeog, main="pRDA CHN(geog partialled out)")
dev.off()
```

![alt_txt][pRDA.CHNgeog]
[pRDA.CHNgeog]:https://cloud.githubusercontent.com/assets/12142475/23725154/bd04708c-044f-11e7-93c7-ffd31113b2cc.png

```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHNclimate <- rda(GenData~Lat+Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)
pRDA.CHNclimate
Call: rda(formula = GenData ~ Lat + Long + Condition(shadow.days +
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

RsquareAdj(RDA.CHNclimate)


anova(pRDA.CHNclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   27.405 1.3275  0.073 .
Residual 11  113.541                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: not rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHNclimate.pdf")
plot(pRDA.CHNclimate, main="pRDA (climate partialled out)")
dev.off()
```


![alt_txt][pRDA.CHNclimate]
[pRDA.CHNclimate]:https://cloud.githubusercontent.com/assets/12142475/23725291/5545050a-0450-11e7-9868-a1fe3ab1707c.png


Find the most important loci associated with Climate

##NBNB Still need to do this for the CHNdata
```
 summary(pRDA.geog)

Accumulated constrained eigenvalues
Importance of components:
                       RDA1  RDA2   RDA3    RDA4    RDA5
Eigenvalue            15.55 9.094 4.4877 2.71164 1.96009
Proportion Explained   0.46 0.269 0.1328 0.08022 0.05799
Cumulative Proportion  0.46 0.729 0.8618 0.94201 1.00000

                 RDA1     RDA2    RDA3     RDA4     RDA5 PC1
bio5.scaled   0.10875 -0.16439  0.2365 -0.07146  0.10935   0
bio15.scaled  0.51115  0.23431  0.3214 -0.04714 -0.05866   0
bio13.scaled  0.47356 -0.23902 -0.3016 -0.17095 -0.53324   0
bio18.scaled  0.47439 -0.06841  0.2092  0.02536 -0.26091   0
bio2.scaled  -0.05063 -0.13750  0.2565 -0.06603 -0.02417   0

RDA1: BIO13, BIO15, BIO18

RDA2: BIO13, BIO15
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```


####3. CHS

###Don't have this input file yet

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

GenData <- read.csv("SE.132.MAF.csv", header=T)
GenData <- GenData[,15:2055]

Climate.Data <- read.csv("SE.145.MAF.csv", header=T)
names(Climate.Data)
Climate.variables <- Climate.Data[,1:14]
names(Climate.variables)
Climate.Data <- Climate.variables[,10:14]
Climate.Data$Lat <- Climate.variables$Lat
Climate.Data$Long <- Climate.variables$Long
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.SEfull <- rda(GenData ~ Lat + Long +bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled, Climate.Data) 
RDA.SEfull
Call: rda(formula = GenData ~ Lat + Long + bio5.scaled + bio15.scaled +
bio13.scaled + bio18.scaled + bio2.scaled, data = Climate.Data)

              Inertia Proportion Rank
Total         80.7640     1.0000     
Constrained   61.8140     0.7654    7
Unconstrained 18.9499     0.2346    7
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7 
25.969 13.331  7.957  4.612  4.274  3.625  2.046 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7 
5.153 3.555 2.913 2.256 1.959 1.732 1.381 

anova(RDA.SEfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     7   68.333 3.2732  0.001 ***
Residual  7   20.876                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

plot(RDA.SEfull)
```

alt_txt

```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.geog <- rda(GenData~bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled+ Condition(Lat + Long), Climate.Data)
head(summary(pRDA.geog))

anova(pRDA.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled + Condition(Lat + Long), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     5   33.803 2.2668  0.003 **
Residual  7   20.876                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate does explain GeneticData

pRDA.geog
Call: rda(formula = GenData ~ bio5.scaled + bio15.scaled + bio13.scaled
+ bio18.scaled + bio2.scaled + Condition(Lat + Long), data =
Climate.Data)

              Inertia Proportion Rank
Total         89.2091     1.0000     
Conditional   34.5300     0.3871    2
Constrained   33.8026     0.3789    5
Unconstrained 20.8765     0.2340    7
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
15.549  9.094  4.488  2.712  1.960 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7 
6.190 3.872 3.460 2.250 2.179 1.662 1.264 


plot(pRDA.geog, main="pRDA (geog partialled out)")
alt_txt

##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.climate <- rda(GenData~Lat+Long + Condition(bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled), Climate.Data)
pRDA.climate
Call: rda(formula = GenData ~ Lat + Long + Condition(bio5.scaled +
bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled), data =
Climate.Data)

               Inertia Proportion Rank
Total         89.20911    1.00000     
Conditional   60.53395    0.67856    5
Constrained    7.79869    0.08742    2
Unconstrained 20.87648    0.23402    7
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2 
5.364 2.435 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7 
6.190 3.872 3.460 2.250 2.179 1.662 1.264 


head(summary(pRDA.climate))

anova(pRDA.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + Condition(bio5.scaled + bio15.scaled + bio13.scaled + bio18.scaled + bio2.scaled), data = Climate.Data)


         Df Variance      F Pr(>F)
Model     2   7.7987 1.3075  0.173
Residual  7  20.8765          

H0: not rejected -> Geography alone does not explain GenData

plot(pRDA.climate, main="pRDA (climate partialled out)")
alt_txt

Find the most important loci associated with Climate

 summary(pRDA.geog)

Accumulated constrained eigenvalues
Importance of components:
                       RDA1  RDA2   RDA3    RDA4    RDA5
Eigenvalue            15.55 9.094 4.4877 2.71164 1.96009
Proportion Explained   0.46 0.269 0.1328 0.08022 0.05799
Cumulative Proportion  0.46 0.729 0.8618 0.94201 1.00000

                 RDA1     RDA2    RDA3     RDA4     RDA5 PC1
bio5.scaled   0.10875 -0.16439  0.2365 -0.07146  0.10935   0
bio15.scaled  0.51115  0.23431  0.3214 -0.04714 -0.05866   0
bio13.scaled  0.47356 -0.23902 -0.3016 -0.17095 -0.53324   0
bio18.scaled  0.47439 -0.06841  0.2092  0.02536 -0.26091   0
bio2.scaled  -0.05063 -0.13750  0.2565 -0.06603 -0.02417   0

RDA1: BIO13, BIO15, BIO18

RDA2: BIO13, BIO15
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

####4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,10:9617]
Climate.Data <- CHS.VS.Data[,3:7]
Climate.Data$Lat <- CHS.VS.Data$lat
Climate.Data$Long <- CHS.VS.Data$long
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.VSfull <- rda(GenData ~ Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.VSfull

Call: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d +
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

anova(RDA.CHS.VSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     7  216.157 1.695  0.002 **
Residual  2   36.436                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.VSfull.pdf")
plot(RDA.CHS.VSfull)
dev.off()
```

![alt_txt][CHS.VS.full.fig1]
[CHS.VS.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23725606/82fcd602-0451-11e7-908d-5636aa5ea84e.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.VS.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long), Climate.Data)

pRDA.CHS.VS.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long), data =
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
20.389 16.04

anova(pRDA.CHS.VS.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     5  142.788 1.5676  0.089 .
Residual  2   36.436                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 not rejected: Climate alone does not explain CHS.VS GeneticData


pdf(file="pRDA.CHS.VS.geog.pdf")
plot(pRDA.CHS.VS.geog, main="pRDA CHS.VS (geog partialled out)")
dev.off()
```


![alt_txt][CHS.VS.Fig2]
[CHS.VS.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23725723/ef6053c8-0451-11e7-92ea-7ab92a39c62f.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.VS.climate <- rda(GenData~Lat+Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.VS.climate
Call: rda(formula = GenData ~ Lat + Long + Condition(shadow.days +
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


head(summary(pRDA.CHS.VS.climate))

anova(pRDA.CHS.VS.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2   66.519 1.8256  0.156
Residual  2   36.436          

H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.pdf")
plot(pRDA.CHS.VS.climate, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```

![alt_txt][CHS.VS.Fig3]
[CHS.VS.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23725815/4ec75ce4-0452-11e7-9835-ebea040a6a2b.png



Find the most important loci associated with Climate

###STILL NEED TO DO THIS FOR CHS.VS
```
 summary(pRDA.geog)

Accumulated constrained eigenvalues
Importance of components:
                       RDA1  RDA2   RDA3    RDA4    RDA5
Eigenvalue            15.55 9.094 4.4877 2.71164 1.96009
Proportion Explained   0.46 0.269 0.1328 0.08022 0.05799
Cumulative Proportion  0.46 0.729 0.8618 0.94201 1.00000

                 RDA1     RDA2    RDA3     RDA4     RDA5 PC1
bio5.scaled   0.10875 -0.16439  0.2365 -0.07146  0.10935   0
bio15.scaled  0.51115  0.23431  0.3214 -0.04714 -0.05866   0
bio13.scaled  0.47356 -0.23902 -0.3016 -0.17095 -0.53324   0
bio18.scaled  0.47439 -0.06841  0.2092  0.02536 -0.26091   0
bio2.scaled  -0.05063 -0.13750  0.2565 -0.06603 -0.02417   0

RDA1: BIO13, BIO15, BIO18

RDA2: BIO13, BIO15
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

####5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI148.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,10:9617]
Climate.Data <- CHS.TI.Data[,3:7]
Climate.Data$Lat <- CHS.TI.Data$lat
Climate.Data$Long <- CHS.TI.Data$long
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.TIfull <- rda(GenData ~ Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.TIfull

Call: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d +
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

anova(RDA.CHS.TIfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     7   187.65 1.3891  0.005 **
Residual  6   115.79                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.TIfull.pdf")
plot(RDA.CHS.TIfull)
dev.off()
```

![alt_txt][CHS.TI.full.fig1]
[CHS.TI.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23725994/edd26400-0452-11e7-809a-4d013669b188.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.TI.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long), Climate.Data)

pRDA.CHS.TI.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long), data =
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


anova(pRDA.CHS.TI.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5   115.48 1.1968  0.152
Residual  6   115.79


H0 not rejected: Climate alone does not explain CHS.TI GeneticData



pdf(file="pRDA.CHS.TI.geog.pdf")
plot(pRDA.CHS.TI.geog, main="pRDA CHS.TI (geog partialled out)")
dev.off()
```


![alt_txt][CHS.TI.Fig2]
[CHS.TI.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23726094/5156dcc2-0453-11e7-9108-268e0a98b4d6.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.TI.climate <- rda(GenData~Lat+Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.TI.climate
Call: rda(formula = GenData ~ Lat + Long + Condition(shadow.days +
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



head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     2   62.735 1.6254  0.057 .
Residual  6  115.790                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       

H0: Not Rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHS.TI.climate.pdf")
plot(pRDA.CHS.TI.climate, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```

![alt_txt][CHS.TI.Fig3]
[CHS.TI.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23726201/bb81b0a4-0453-11e7-8289-c68e78425571.png



Find the most important loci associated with Climate

###STILL NEED TO DO THIS FOR CH940
```
 summary(pRDA.geog)

Accumulated constrained eigenvalues
Importance of components:
                       RDA1  RDA2   RDA3    RDA4    RDA5
Eigenvalue            15.55 9.094 4.4877 2.71164 1.96009
Proportion Explained   0.46 0.269 0.1328 0.08022 0.05799
Cumulative Proportion  0.46 0.729 0.8618 0.94201 1.00000

                 RDA1     RDA2    RDA3     RDA4     RDA5 PC1
bio5.scaled   0.10875 -0.16439  0.2365 -0.07146  0.10935   0
bio15.scaled  0.51115  0.23431  0.3214 -0.04714 -0.05866   0
bio13.scaled  0.47356 -0.23902 -0.3016 -0.17095 -0.53324   0
bio18.scaled  0.47439 -0.06841  0.2092  0.02536 -0.26091   0
bio2.scaled  -0.05063 -0.13750  0.2565 -0.06603 -0.02417   0

RDA1: BIO13, BIO15, BIO18

RDA2: BIO13, BIO15
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

####6. CZ

#####Still need this input file

```
library(vegan)

CHallData <- read.csv("CH940.MAF.csv", header=T)
CHallData <- CHallData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CHallData[,11:9618]
Climate.Data <- CHallData[,3:7]
Climate.Data$Lat <- CHallData$lat
Climate.Data$Long <- CHallData$long
Climate.Data$elev <- CHallData$elev
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CH940full <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. I got an error (Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd') when elev was after Long in the formula

RDA.CH940full

Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         461.1563     1.0000     
Constrained   184.0588     0.3991    8
Unconstrained 601.2492     1.3038   70
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8 
123.58  31.44  10.17   5.02   4.71   3.36   3.09   2.70 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 20.48 19.92 19.28 18.48 17.89 17.29 17.19 
(Showed only 8 of all 70 unconstrained eigenvalues)

##calculate proportion explained (R2)

RsquareAdj(RDA.CH940full)

$r.squared
[1] 0.3991245

$adj.r.squared
[1] 0.3323605



anova(RDA.CH940full)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     8   184.06 2.7551  0.001 ***
Residual 72   601.25                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CH940full.elev.pdf")
plot(RDA.CH940full)
dev.off()
```

![alt_txt][CH940.full.fig1]
[CH940.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23798623/c56ded96-05a5-11e7-8085-70fb6aad317d.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CH940.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long +elev), Climate.Data)

#head(summary(pRDA.CH940.geog))

pRDA.CH940.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

                Inertia Proportion Rank
Total         461.15630    1.00000     
Conditional   144.57238    0.31350    3
Constrained    39.48639    0.08562    5
Unconstrained 277.09753    0.60088   72
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
22.377  6.884  4.001  3.336  2.888 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 15.51 10.10  7.53  7.03  6.97  6.33  5.72 
(Showed only 8 of all 72 unconstrained eigenvalues)


RsquareAdj(pRDA.CH940.geog)

$r.squared
[1] 0.08562474

$adj.r.squared
[1] 0.04560755

anova(pRDA.CH940.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance     F Pr(>F)    
Model     5   39.486 2.052  0.001 ***
Residual 72  277.098                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CH940 GeneticData


pdf(file="pRDA.CH940.geog.elev.pdf")
plot(pRDA.CH940.geog, main="pRDA CH940 (geog partialled out)")
dev.off()
```


![alt_txt][CH940.Fig2]
[CH940.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23798791/6cee94bc-05a6-11e7-91e8-587383b3428c.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CH940.climate <- rda(GenData~Lat+Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CH940.climate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

                Inertia Proportion Rank
Total         461.15630    1.00000     
Conditional    41.31089    0.08958    5
Constrained   142.74788    0.30954    3
Unconstrained 277.09753    0.60088   72
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3 
112.42  26.51   3.82 

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
55.67 15.51 10.10  7.53  7.03  6.97  6.33  5.72 
(Showed only 8 of all 72 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CH940.climate)

$r.squared
[1] 0.3095434

$adj.r.squared
[1] 0.303474

anova(pRDA.CH940.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     3   142.75 12.364  0.001 ***
Residual 72   277.10                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1         

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CH940.climate.elev.pdf")
plot(pRDA.CH940.climate, main="pRDA CH940 (climate partialled out)")
dev.off()
```

![alt_txt][CH940.Fig3]
[CH940.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23798964/19945ab2-05a7-11e7-8465-30d221431940.png



Find the most important loci associated with Climate
```
 summary(pRDA.CH940.geog)

Importance of components:
                          RDA1    RDA2    RDA3    RDA4    RDA5     PC1     PC2
Eigenvalue            22.37711 6.88381 4.00051 3.33649 2.88847 55.6741 15.5130
Proportion Explained   0.07068 0.02174 0.01264 0.01054 0.00912  0.1759  0.0490
Cumulative Proportion  0.07068 0.09243 0.10506 0.11560 0.12473  0.3006  0.3496

                    RDA1     RDA2       RDA3    RDA4     RDA5 PC1
shadow.days       0.3497 -0.21845  0.0544813 -0.6602 -0.62512   0
sol.rad.60d      -0.3521 -0.25972 -0.0752741  0.5668 -0.26445   0
pcpt.60d          0.4492 -0.70416  0.0004636  0.1741 -0.01842   0
day10cm           0.2255  0.04426  0.3289740 -0.1241 -0.23575   0
temp.laying.date -0.1623 -0.67753  0.3294060 -0.4541  0.32005   0

RDA1: pcpt.60d, sol.rad.60d, shadow.days

RDA2: pcpt.60d, temp.laying.date
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CH940.geog <- summary(pRDA.CH940.geog)

CH940.RDAloadings <- sum.pRDA.CH940.geog$species  ##write the loadings to data.frame
CH940.RDAloadings <- as.data.frame(CH940.RDAloadings) #change to df
names(CH940.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CH940.RDAloadings$RDA1)
hist(CH940.RDAloadings$RDA2)
hist(CH940.RDAloadings$RDA3)
hist(CH940.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CH940.RDA1.outliers.1percent <- CH940.RDAloadings.sorted[1:100, 1:2]  ##write to df

CH940.RDAloadings.sorted <- CH940.RDAloadings[order(-abs(CH940.RDAloadings$RDA2)),]  #same for RDA2
CH940.RDA2.outliers.1percent <- CH940.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CH940.RDA1.outliers.1percent, "CH940.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CH940.RDA2.outliers.1percent, "CH940.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```
