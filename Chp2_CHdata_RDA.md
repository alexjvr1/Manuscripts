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

Accumulated constrained eigenvalues
Importance of components:
                         RDA1   RDA2   RDA3   RDA4    RDA5
Eigenvalue            22.3771 6.8838 4.0005 3.3365 2.88847
Proportion Explained   0.5667 0.1743 0.1013 0.0845 0.07315
Cumulative Proportion  0.5667 0.7410 0.8424 0.9268 1.00000

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
GenData <- CHNData[,11:9618]
Climate.Data <- CHNData[,3:7]
Climate.Data$Lat <- CHNData$lat
Climate.Data$Long <- CHNData$long
Climate.Data$elev <- CHNData$elev

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHNfull <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHNfull
Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         218.3893     1.0000     
Constrained   116.3646     0.5328    8
Unconstrained 102.0247     0.4672   10
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5   RDA6   RDA7   RDA8 
27.579 23.436 16.344 13.696 11.119  8.764  8.334  7.094 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10 
20.592 13.234 11.357 10.658  9.426  8.675  8.514  7.397  6.460  5.711 


RsquareAdj(RDA.CHNfull)

$r.squared
[1] 0.5328311

$adj.r.squared
[1] 0.1590959



anova(RDA.CHNfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     8   116.36 1.4257  0.001 ***
Residual 10   102.03                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHNfull.elev.pdf")
plot(RDA.CHNfull)
dev.off()
```

![alt_txt][CHN.Fig1]
[CHN.Fig1]:https://cloud.githubusercontent.com/assets/12142475/23800453/c1e573cc-05ac-11e7-8b04-e0c335d1edcd.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHNgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(Lat + Long + elev), Climate.Data)

RsquareAdj(pRDA.CHNgeog)

$r.squared
[1] 0.2689343

$adj.r.squared
[1] 0.04241981


anova(pRDA.CHNgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5   58.732 1.1513   0.13
Residual 10  102.025       


H0 not rejected: Climate does not explain a significant amoutn of CHN GeneticData

pRDA.CHNgeog
Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

               Inertia Proportion Rank
Total         218.3893     1.0000     
Conditional    57.6322     0.2639    3
Constrained    58.7324     0.2689    5
Unconstrained 102.0247     0.4672   10
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
18.579 12.879 10.688  8.817  7.769 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10 
20.592 13.234 11.357 10.658  9.426  8.675  8.514  7.397  6.460  5.711 




pdf(file="pRDA.CHNgeog.elev.pdf")
plot(pRDA.CHNgeog, main="pRDA CHN(geog partialled out)")
dev.off()
```

![alt_txt][pRDA.CHNgeog]
[pRDA.CHNgeog]:https://cloud.githubusercontent.com/assets/12142475/23800567/323bcdce-05ad-11e7-9ea7-11de9c4df627.png

```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHNclimate <- rda(GenData~Lat+Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHNclimate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         218.3893     1.0000     
Conditional    77.4433     0.3546    5
Constrained    38.9213     0.1782    3
Unconstrained 102.0247     0.4672   10
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3 
17.216 11.984  9.722 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10 
20.592 13.234 11.357 10.658  9.426  8.675  8.514  7.397  6.460  5.711 

RsquareAdj(pRDA.CHNclimate)

$r.squared
[1] 0.1782199

$adj.r.squared
[1] 0.05271127

anova(pRDA.CHNclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     3   38.921 1.2716  0.067 .
Residual 10  102.025                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: not rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHNclimate.elev.pdf")
plot(pRDA.CHNclimate, main="pRDA (climate partialled out)")
dev.off()
```


![alt_txt][pRDA.CHNclimate]
[pRDA.CHNclimate]:https://cloud.githubusercontent.com/assets/12142475/23801293/eef1e2b2-05af-11e7-9318-7a09fc6317fd.png


Find the most important loci associated with Climate

##NBNB Still need to do this for the CHNdata
```
summary(pRDA.CHNgeog)

Accumulated constrained eigenvalues
Importance of components:
                         RDA1    RDA2    RDA3   RDA4   RDA5
Eigenvalue            18.5794 12.8795 10.6880 8.8168 7.7687
Proportion Explained   0.3163  0.2193  0.1820 0.1501 0.1323
Cumulative Proportion  0.3163  0.5356  0.7176 0.8677 1.0000

                     RDA1    RDA2     RDA3      RDA4    RDA5 PC1
shadow.days      -0.61073  0.3666 -0.05435  0.347371 -0.2364   0
sol.rad.60d       0.04329  0.1876 -0.32866 -0.155540  0.2277   0
pcpt.60d         -0.21112 -0.0454 -0.17428  0.074965  0.2094   0
day10cm          -0.30440  0.2813  0.21105 -0.007648  0.1068   0
temp.laying.date -0.51582 -0.1804  0.16413 -0.414999 -0.4870   0

RDA1: shadow.days, temp.laying.date

RDA2: shadow.days, day10cm
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHNgeog <- summary(pRDA.CHNgeog)

CHN.RDAloadings <- sum.pRDA.CHNgeog$species  ##write the loadings to data.frame
CHN.RDAloadings <- as.data.frame(CHN.RDAloadings) #change to df
names(CHN.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHN.RDAloadings$RDA1)
hist(CHN.RDAloadings$RDA2)
hist(CHN.RDAloadings$RDA3)
hist(CHN.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CHN.RDAloadings.sorted <- CHN.RDAloadings[order(-abs(CHN.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHN.RDA1.outliers.1percent <- CHN.RDAloadings.sorted[1:100, 1:2]  ##write to df

CHN.RDAloadings.sorted <- CHN.RDAloadings[order(-abs(CHN.RDAloadings$RDA2)),]  #same for RDA2
CHN.RDA2.outliers.1percent <- CHN.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CHN.RDA1.outliers.1percent, "CHN.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHN.RDA2.outliers.1percent, "CHN.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```


####3. CHS

See this tutorial for the interpretation: REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger redundancy-analysis-for-landscape-genetics.pdf on mac

```
library(vegan)

CHSData <- read.csv("CHS283.MAF.csv", header=T)
CHSData <- CHSData[complete.cases(CHSData),]  ##stba is removed here due to missing data
GenData <- CHSData[,11:9618]
Climate.Data <- CHSData[,3:7]
Climate.Data$Lat <- CHSData$lat
Climate.Data$Long <- CHSData$long
Climate.Data$elev <- CHSData$elev

names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHSfull <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) 

RDA.CHSfull
Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         345.1880     1.0000     
Constrained   169.2622     0.4903    8
Unconstrained 175.9258     0.5097   15
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7  RDA8 
70.28 19.89 17.86 15.38 15.25 12.58  9.98  8.03 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
25.085 23.143 19.370 16.162 11.444 10.913 10.201  9.823  9.002  8.942  7.785 
  PC12   PC13   PC14   PC15 
 6.646  6.252  5.925  5.235 


RsquareAdj(RDA.CHSfull)

$r.squared
[1] 0.4903478

$adj.r.squared
[1] 0.2185333



anova(RDA.CHSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)    
Model     8   169.26 1.804  0.001 ***
Residual 15   175.93                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHSfull.elev.pdf")
plot(RDA.CHSfull)
dev.off()
```

![alt_txt][CHS.Fig1]
[CHS.Fig1]:https://cloud.githubusercontent.com/assets/12142475/23801930/0bf1bc5a-05b2-11e7-8954-fd562989513f.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHSgeog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date+ Condition(Lat + Long + elev), Climate.Data)

RsquareAdj(pRDA.CHSgeog)

$r.squared
[1] 0.191701

$adj.r.squared
[1] 0.02508945


anova(pRDA.CHSgeog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5   66.173 1.1284   0.16
Residual 15  175.926     


H0 not rejected: Climate does not explain a significant amoutn of CHN GeneticData

pRDA.CHSgeog
Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

               Inertia Proportion Rank
Total         345.1880     1.0000     
Conditional   103.0893     0.2986    3
Constrained    66.1729     0.1917    5
Unconstrained 175.9258     0.5097   15
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
17.649 15.714 13.663 10.309  8.838 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
25.085 23.143 19.370 16.162 11.444 10.913 10.201  9.823  9.002  8.942  7.785 
  PC12   PC13   PC14   PC15 
 6.646  6.252  5.925  5.235 




pdf(file="pRDA.CHSgeog.elev.pdf")
plot(pRDA.CHSgeog, main="pRDA CHS(geog partialled out)")
dev.off()
```

![alt_txt][pRDA.CHSgeog]
[pRDA.CHSgeog]:https://cloud.githubusercontent.com/assets/12142475/23802059/75764916-05b2-11e7-8dc8-71d3e8c19f74.png

```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHSclimate <- rda(GenData~Lat+Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHSclimate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         345.1880     1.0000     
Conditional   117.0711     0.3392    5
Constrained    52.1911     0.1512    3
Unconstrained 175.9258     0.5097   15
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3 
23.895 17.096 11.201 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
25.085 23.143 19.370 16.162 11.444 10.913 10.201  9.823  9.002  8.942  7.785 
  PC12   PC13   PC14   PC15 
 6.646  6.252  5.925  5.235 

RsquareAdj(pRDA.CHSclimate)

$r.squared
[1] 0.1511962

$adj.r.squared
[1] 0.06295074

anova(pRDA.CHSclimate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     3   52.191 1.4833  0.014 *
Residual 15  175.926                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          
          
H0: rejected -> Geography explains a significant portion of GenData

pdf(file="pRDA.CHSclimate.elev.pdf")
plot(pRDA.CHSclimate, main="pRDA (climate partialled out)")
dev.off()
```


![alt_txt][pRDA.CHSclimate]
[pRDA.CHSclimate]:https://cloud.githubusercontent.com/assets/12142475/23802225/f27c7a7a-05b2-11e7-86b8-6c12e8b3f767.png


Find the most important loci associated with Climate

##NBNB Still need to do this for the CHNdata
```
head(summary(pRDA.CHSgeog))

Accumulated constrained eigenvalues
Importance of components:
                         RDA1    RDA2    RDA3   RDA4   RDA5
Eigenvalue            18.5794 12.8795 10.6880 8.8168 7.7687
Proportion Explained   0.3163  0.2193  0.1820 0.1501 0.1323
Cumulative Proportion  0.3163  0.5356  0.7176 0.8677 1.0000

 Biplot scores for constraining variables

                     RDA1     RDA2     RDA3     RDA4    RDA5 PC1
shadow.days      -0.58258  0.58674 -0.15416 -0.06295 -0.3059   0
sol.rad.60d       0.07458 -0.52244  0.19003  0.47917 -0.1448   0
pcpt.60d         -0.25403 -0.19847 -0.35009  0.02455  0.1738   0
day10cm           0.01894 -0.06537 -0.13778  0.03732 -0.1734   0
temp.laying.date -0.64893 -0.12578  0.06734 -0.64539  0.1613   0

RDA1: shadow.days, temp.laying.date

RDA2: shadow.days, sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHSgeog <- summary(pRDA.CHSgeog)

CHS.RDAloadings <- sum.pRDA.CHSgeog$species  ##write the loadings to data.frame
CHS.RDAloadings <- as.data.frame(CHS.RDAloadings) #change to df
names(CHS.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.RDAloadings$RDA1)
hist(CHS.RDAloadings$RDA2)
hist(CHS.RDAloadings$RDA3)
hist(CHS.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.RDAloadings.sorted <- CHS.RDAloadings[order(-abs(CHS.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.RDA1.outliers.1percent <- CHS.RDAloadings.sorted[1:100, 1:2]  ##write to df

CHS.RDAloadings.sorted <- CHS.RDAloadings[order(-abs(CHS.RDAloadings$RDA2)),]  #same for RDA2
CHS.RDA2.outliers.1percent <- CHS.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.RDA1.outliers.1percent, "CHS.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.RDA2.outliers.1percent, "CHS.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```



####4. CHS.VS

```
library(vegan)

CHS.VS.Data <- read.csv("CHS.VS135.MAF.csv", header=T)
CHS.VS.Data <- CHS.VS.Data[complete.cases(CHS.VS.Data),]  ##stba is removed here due to missing data
GenData <- CHS.VS.Data[,11:9618]
Climate.Data <- CHS.VS.Data[,3:7]
Climate.Data$Lat <- CHS.VS.Data$lat
Climate.Data$Long <- CHS.VS.Data$long
Climate.Data$elev <- CHS.VS.Data$elev
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.VSfull <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.VSfull

Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

                Inertia Proportion Rank
Total         252.59331    1.00000     
Constrained   235.62747    0.93283    8
Unconstrained  16.96584    0.06717    1
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7  RDA8 
61.39 37.74 30.66 29.41 24.79 19.96 17.62 14.05 

Eigenvalues for unconstrained axes:
   PC1 
16.966 


RsquareAdj(RDA.CHS.VSfull)
$r.squared
[1] 0.9328334

$adj.r.squared
[1] 0.3955005


anova(RDA.CHS.VSfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     8  235.627 1.736  0.008 **
Residual  1   16.966                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.VSfull.pdf")
plot(RDA.CHS.VSfull)
dev.off()
```

![alt_txt][CHS.VS.full.fig1]
[CHS.VS.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23802427/97841474-05b3-11e7-8d44-cd92159a3fee.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.VS.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long +elev), Climate.Data)

pRDA.CHS.VS.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

                Inertia Proportion Rank
Total         252.59331    1.00000     
Conditional   112.97328    0.44725    3
Constrained   122.65419    0.48558    5
Unconstrained  16.96584    0.06717    1
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
40.47 27.44 22.55 18.03 14.16 

Eigenvalues for unconstrained axes:
   PC1 
16.966 

$r.squared
[1] 0.4855797

$adj.r.squared
[1] 0.22462


anova(pRDA.CHS.VS.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5  122.654 1.4459  0.204
Residual  1   16.966 


H0 not rejected: Climate alone does not explain CHS.VS GeneticData


pdf(file="pRDA.CHS.VS.geog.pdf")
plot(pRDA.CHS.VS.geog, main="pRDA CHS.VS (geog partialled out)")
dev.off()
```


![alt_txt][CHS.VS.Fig2]
[CHS.VS.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23802573/0b9f49a0-05b4-11e7-886b-01d9b7b4defb.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.VS.climate <- rda(GenData~Lat+Long +elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.VS.climate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

                Inertia Proportion Rank
Total         252.59331    1.00000     
Conditional   149.63843    0.59241    5
Constrained    85.98905    0.34042    3
Unconstrained  16.96584    0.06717    1
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3 
47.53 21.53 16.93 

Eigenvalues for unconstrained axes:
   PC1 
16.966 


#head(summary(pRDA.CHS.VS.climate))

RsquareAdj(pRDA.CHS.VS.climate)

$r.squared
[1] 0.3404249

$adj.r.squared
[1] 0.3125813

anova(pRDA.CHS.VS.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     3   85.989 1.6895  0.264
Residual  1   16.966              

H0: Not Rejected -> Geography alone doesn't explain GenData

pdf(file="pRDA.CHS.VS.climate.pdf")
plot(pRDA.CHS.VS.climate, main="pRDA CHS.VS (climate partialled out)")
dev.off()
```

![alt_txt][CHS.VS.Fig3]
[CHS.VS.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23802708/67e21b8e-05b4-11e7-936a-9170ff6b281e.png



Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.VS.geog))

Accumulated constrained eigenvalues
Importance of components:
                         RDA1    RDA2    RDA3    RDA4    RDA5
Eigenvalue            40.4675 27.4375 22.5537 18.0311 14.1643
Proportion Explained   0.3299  0.2237  0.1839  0.1470  0.1155
Cumulative Proportion  0.3299  0.5536  0.7375  0.8845  1.0000

                     RDA1      RDA2     RDA3      RDA4     RDA5 PC1
shadow.days       0.27453  0.196130  0.21966  0.098751  0.02212   0
sol.rad.60d      -0.02280 -0.014690 -0.05157 -0.515855 -0.03210   0
pcpt.60d         -0.27810  0.176723 -0.15135  0.007086  0.07282   0
day10cm          -0.08859 -0.043333  0.11277 -0.141309 -0.10641   0
temp.laying.date -0.07741  0.009801 -0.22783  0.664696 -0.07906   0

RDA1: shadow.days, pcpt.60d

RDA2: shadow.days, pcpt.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHS.VS.geog <- summary(pRDA.CHS.VS.geog)

CHS.VS.RDAloadings <- sum.pRDA.CHS.VS.geog$species  ##write the loadings to data.frame
CHS.VS.RDAloadings <- as.data.frame(CHS.VS.RDAloadings) #change to df
names(CHS.VS.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.VS.RDAloadings$RDA1)
hist(CHS.VS.RDAloadings$RDA2)
hist(CHS.VS.RDAloadings$RDA3)
hist(CHS.VS.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.VS.RDAloadings.sorted <- CHS.VS.RDAloadings[order(-abs(CHS.VS.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.VS.RDA1.outliers.1percent <- CHS.VS.RDAloadings.sorted[1:100, 1:2]  ##write to df

CHS.VS.RDAloadings.sorted <- CHS.VS.RDAloadings[order(-abs(CHS.VS.RDAloadings$RDA2)),]  #same for RDA2
CHS.VS.RDA2.outliers.1percent <- CHS.VS.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.VS.RDA1.outliers.1percent, "CHS.VS.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.VS.RDA2.outliers.1percent, "CHS.VS.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

####5. CHS.TI

```
library(vegan)

CHS.TI.Data <- read.csv("CHS.TI148.MAF.csv", header=T)
CHS.TI.Data <- CHS.TI.Data[complete.cases(CHS.TI.Data),]  ##stba is removed here due to missing data
GenData <- CHS.TI.Data[,11:9618]
Climate.Data <- CHS.TI.Data[,3:7]
Climate.Data$Lat <- CHS.TI.Data$lat
Climate.Data$Long <- CHS.TI.Data$long
Climate.Data$elev <- CHS.TI.Data$elev
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CHS.TIfull <- rda(GenData ~ elev + Lat + Long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s 
RDA.CHS.TIfull

Call: rda(formula = GenData ~ elev + Lat + Long + shadow.days +
sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data =
Climate.Data)

               Inertia Proportion Rank
Total         303.4368     1.0000     
Constrained   210.9920     0.6953    8
Unconstrained  92.4448     0.3047    5
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7  RDA8 
46.22 40.74 35.99 27.61 19.21 14.53 13.92 12.77 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5 
27.568 23.412 17.927 14.014  9.523 

RsquareAdj(RDA.CHS.TIfull)
$r.squared
[1] 0.6953409

$adj.r.squared
[1] 0.2078862

anova(RDA.CHS.TIfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Lat + Long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     8  210.992 1.4265  0.003 **
Residual  5   92.445                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CHS.TIfull.pdf")
plot(RDA.CHS.TIfull)
dev.off()
```

![alt_txt][CHS.TI.full.fig1]
[CHS.TI.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23803001/6abec66c-05b5-11e7-8b04-71fca5e9b5b3.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CHS.TI.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long +elev), Climate.Data)

pRDA.CHS.TI.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(Lat + Long + elev), data =
Climate.Data)

               Inertia Proportion Rank
Total         303.4368     1.0000     
Conditional    98.4304     0.3244    3
Constrained   112.5617     0.3710    5
Unconstrained  92.4448     0.3047    5
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3  RDA4  RDA5 
38.66 26.71 19.36 14.66 13.17 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5 
27.568 23.412 17.927 14.014  9.523 


RsquareAdj(pRDA.CHS.TI.geog)
$r.squared
[1] 0.3709558

$adj.r.squared
[1] 0.08618564


anova(pRDA.CHS.TI.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(Lat + Long + elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     5  112.562 1.2176  0.171
Residual  5   92.445 


H0 not rejected: Climate alone does not explain CHS.TI GeneticData



pdf(file="pRDA.CHS.TI.geog.pdf")
plot(pRDA.CHS.TI.geog, main="pRDA CHS.TI (geog partialled out)")
dev.off()
```


![alt_txt][CHS.TI.Fig2]
[CHS.TI.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23803098/bcd8240c-05b5-11e7-89ed-c659e6d6a900.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CHS.TI.climate <- rda(GenData~Lat+Long +elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CHS.TI.climate
Call: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         303.4368     1.0000     
Conditional   124.9119     0.4117    5
Constrained    86.0801     0.2837    3
Unconstrained  92.4448     0.3047    5
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3 
36.00 27.03 23.05 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5 
27.568 23.412 17.927 14.014  9.523 

RsquareAdj(pRDA.CHS.TI.climate)
$r.squared
[1] 0.2836837

$adj.r.squared
[1] 0.1639434

#head(summary(pRDA.CHS.TI.climate))

anova(pRDA.CHS.TI.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Lat + Long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     3   86.080 1.5519  0.063 .
Residual  5   92.445                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       

H0: Not Rejected -> Geography alone does not explain GenData

pdf(file="pRDA.CHS.TI.climate.pdf")
plot(pRDA.CHS.TI.climate, main="pRDA CHS.TI (climate partialled out)")
dev.off()
```

![alt_txt][CHS.TI.Fig3]
[CHS.TI.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23803199/12397ee6-05b6-11e7-9214-24205ceae95d.png



Find the most important loci associated with Climate
```
 head(summary(pRDA.CHS.TI.geog))

Accumulated constrained eigenvalues
Importance of components:
                         RDA1    RDA2    RDA3    RDA4   RDA5
Eigenvalue            38.6573 26.7132 19.3624 14.6589 13.170
Proportion Explained   0.3434  0.2373  0.1720  0.1302  0.117
Cumulative Proportion  0.3434  0.5807  0.7528  0.8830  1.000


                     RDA1     RDA2     RDA3     RDA4     RDA5 PC1
shadow.days       0.60147  0.30551 -0.33814 -0.19323  0.54025   0
sol.rad.60d      -0.35527 -0.11283  0.08379  0.31205 -0.02584   0
pcpt.60d          0.33034 -0.06300 -0.39163  0.54204 -0.06361   0
day10cm           0.04061  0.00895  0.06302  0.02395  0.15971   0
temp.laying.date  0.19051  0.57288 -0.53335 -0.33446 -0.23098   0


RDA1: shadow.days

RDA2: temp.laying.date
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CHS.TI.geog <- summary(pRDA.CHS.TI.geog)

CHS.TI.RDAloadings <- sum.pRDA.CHS.TI.geog$species  ##write the loadings to data.frame
CHS.TI.RDAloadings <- as.data.frame(CHS.TI.RDAloadings) #change to df
names(CHS.TI.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CHS.TI.RDAloadings$RDA1)
hist(CHS.TI.RDAloadings$RDA2)
hist(CHS.TI.RDAloadings$RDA3)
hist(CHS.TI.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CHS.TI.RDAloadings.sorted <- CHS.TI.RDAloadings[order(-abs(CHS.TI.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CHS.TI.RDA1.outliers.1percent <- CHS.TI.RDAloadings.sorted[1:100, 1:2]  ##write to df

CHS.TI.RDAloadings.sorted <- CHS.TI.RDAloadings[order(-abs(CHS.TI.RDAloadings$RDA2)),]  #same for RDA2
CHS.TI.RDA2.outliers.1percent <- CHS.TI.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CHS.TI.RDA1.outliers.1percent, "CHS.TI.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CHS.TI.RDA2.outliers.1percent, "CHS.TI.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```

####6. CZ

```
library(vegan)

CHallData <- read.csv("CZ404.MAF.csv", header=T)
CHallData <- CHallData[complete.cases(CHallData),]  ##stba is removed here due to missing data
GenData <- CHallData[,11:9618]
Climate.Data <- CHallData[,3:10]
#Climate.Data$Lat <- CHallData$lat
#Climate.Data$Long <- CHallData$long
#Climate.Data$elev <- CHallData$elev
names(Climate.Data)
```


```
##1. Run Full RDA model to determine how much of the variation is explainable by the expanatory variables we have
##H0: climate data does not affect genotype

RDA.CZfull <- rda(GenData ~ elev + lat + long +shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, Climate.Data) ##runs in 20s. 

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

RsquareAdj(RDA.CZfull)

$r.squared
[1] 0.3949408

$adj.r.squared
[1] 0.2280279



anova(RDA.CZfull)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + lat + long + shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date, data = Climate.Data)
         Df Variance      F Pr(>F)    
Model     8   131.75 2.3661  0.001 ***
Residual 29   201.84                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##to see which variables are most important, we can plot the results in a biplot

pdf(file="RDA.CZfull.elev.pdf")
plot(RDA.CZfull)
dev.off()
```

![alt_txt][CZ.full.fig1]
[CZ.full.fig1]:https://cloud.githubusercontent.com/assets/12142475/23803552/338e3dd8-05b7-11e7-9d71-d147c99a23db.png


```
##Partial out geog
H0: Climate does not explain genetic data

pRDA.CZ.geog <- rda(GenData~shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long +elev), Climate.Data)

#head(summary(pRDA.CZ.geog))

pRDA.CZ.geog

Call: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d +
day10cm + temp.laying.date + Condition(lat + long + elev), data =
Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Conditional    84.7235     0.2540    3
Constrained    47.0238     0.1410    5
Unconstrained 201.8402     0.6051   29
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1   RDA2   RDA3   RDA4   RDA5 
20.994  8.460  7.550  5.069  4.950 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.022 17.900 11.864 11.718 10.560  9.379  8.851  7.973 
(Showed only 8 of all 29 unconstrained eigenvalues)


RsquareAdj(pRDA.CZ.geog)

$r.squared
[1] 0.140964

$adj.r.squared
[1] 0.03987665

anova(pRDA.CZ.geog)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date + Condition(lat + long + elev), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     5   47.024 1.3513  0.005 **
Residual 29  201.840                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


H0 rejected: Climate explains a significant amount of CZ GeneticData


pdf(file="pRDA.CZ.geog.elev.pdf")
plot(pRDA.CZ.geog, main="pRDA CZ (geog partialled out)")
dev.off()
```


![alt_txt][CZ.Fig2]
[CZ.Fig2]:https://cloud.githubusercontent.com/assets/12142475/23803843/1dc60d18-05b8-11e7-97b3-fdf1c32ff38d.png


```
##Partial out climate

H0: Geog alone does not explain Genetic data

pRDA.CZ.climate <- rda(GenData~lat+long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), Climate.Data)

pRDA.CZ.climate
Call: rda(formula = GenData ~ lat + long + elev + Condition(shadow.days
+ sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data =
Climate.Data)

               Inertia Proportion Rank
Total         333.5875     1.0000     
Conditional    75.1625     0.2253    5
Constrained    56.5848     0.1696    3
Unconstrained 201.8402     0.6051   29
Inertia is variance 

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3 
40.85  9.34  6.39 

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
26.022 17.900 11.864 11.718 10.560  9.379  8.851  7.973 
(Showed only 8 of all 29 unconstrained eigenvalues)


#head(summary(pRDA.CH940.climate))

RsquareAdj(pRDA.CZ.climate)

$r.squared
[1] 0.169625

$adj.r.squared
[1] 0.1237565

anova(pRDA.CZ.climate)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ lat + long + elev + Condition(shadow.days + sol.rad.60d + pcpt.60d + day10cm + temp.laying.date), data = Climate.Data)
         Df Variance    F Pr(>F)    
Model     3   56.585 2.71  0.001 ***
Residual 29  201.840                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

H0: Rejected -> Geography explains a significant amount of GenData

pdf(file="pRDA.CZ.climate.elev.pdf")
plot(pRDA.CZ.climate, main="pRDA CZ (climate partialled out)")
dev.off()
```

![alt_txt][CZ.Fig3]
[CZ.Fig3]:https://cloud.githubusercontent.com/assets/12142475/23803964/937a0118-05b8-11e7-89c1-720883a40a06.png



Find the most important loci associated with Climate
```
(head(summary(pRDA.CZ.geog)))

Accumulated constrained eigenvalues
Importance of components:
                         RDA1   RDA2   RDA3   RDA4   RDA5
Eigenvalue            20.9944 8.4603 7.5504 5.0686 4.9501
Proportion Explained   0.4465 0.1799 0.1606 0.1078 0.1053
Cumulative Proportion  0.4465 0.6264 0.7869 0.8947 1.0000

                    RDA1     RDA2     RDA3     RDA4     RDA5 PC1
shadow.days       0.3785  0.04372 -0.07447 -0.39420 -0.59888   0
sol.rad.60d      -0.2597  0.41389  0.12271 -0.16903  0.32231   0
pcpt.60d          0.4291 -0.01140  0.21903 -0.01499  0.08971   0
day10cm           0.3856  0.17103 -0.06955 -0.10843 -0.27725   0
temp.laying.date  0.3328 -0.23521  0.76408  0.39983 -0.16824   0

RDA1: pcpt.60d

RDA2: sol.rad.60d
```

Select the top 100 loci (~1%) associated with RDA1 and RDA2
```
#Find the loadings of of loci on each RDA

sum.pRDA.CZ.geog <- summary(pRDA.CZ.geog)

CZ.RDAloadings <- sum.pRDA.CZ.geog$species  ##write the loadings to data.frame
CZ.RDAloadings <- as.data.frame(CZ.RDAloadings) #change to df
names(CZ.RDAloadings)
par(mfrow=c(2,2))  #check that loadings on each RDA is normally distributed
hist(CZ.RDAloadings$RDA1)
hist(CZ.RDAloadings$RDA2)
hist(CZ.RDAloadings$RDA3)
hist(CZ.RDAloadings$RDA4)

##select the top 100 loci (~1%) for each RDA

CZ.RDAloadings.sorted <- CZ.RDAloadings[order(-abs(CZ.RDAloadings$RDA1)),] #sort by absolute value of RDA1. be sure to add "-" in front of "abs" to order from largest to smallest
CZ.RDA1.outliers.1percent <- CZ.RDAloadings.sorted[1:100, 1:2]  ##write to df

CZ.RDAloadings.sorted <- CZ.RDAloadings[order(-abs(CZ.RDAloadings$RDA2)),]  #same for RDA2
CZ.RDA2.outliers.1percent <- CZ.RDAloadings.sorted[1:100, c(1:2)]  

#write both to file. 
write.table(CZ.RDA1.outliers.1percent, "CZ.RDA1.outliers", sep=" ", col.names=T, row.names=T, quote=F)
write.table(CZ.RDA2.outliers.1percent, "CZ.RDA2.outliers", sep=" ", col.names=T, row.names=T, quote=F)

```
