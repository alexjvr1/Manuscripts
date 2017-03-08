#RDA analysis of CH dataset

RDA: Effect of Climate & Geog on genetic variation


I can use RDA to partition the effects of climate and geog on genetic variation. 

Explanation and example: http://www.davidzeleny.net/anadat-r/doku.php/en:rda_examples

And Victoria's paper: http://www.amjbot.org/content/103/1/33.full.pdf+html

###Data:

I'm using the full CH dataset. I have to run the analysis seperately for all the different partitions of the data to investigate adaptation
at different geographic scales. 

####1. CHall

940 indivs

9608 loci

####2. CHN

229 indivs

9608 loci

####3. CHS

283 indivs

9608 loci

####4. CHS.VS

135 indivs

9608 loci

####5. CHS.TI

148 indivs

9608 loci

####6. CZ

404 indivs

9608 loci



###Input files:

####1. MAF of all loci

```
#Calculate MAF for the full dataset within region using PLINK

/Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/RDA/

###Use the *nosex file to create a file for subsetting the data. Here I've made 1=within region, 2=within pop

plink --file SE132.FINAL.plink --freq --within SE132.region.cluster --out SE132.cluster

plink --file SE132.FINAL.plink --freq --within SE132.pop.cluster --out SE132.pop
```


```
###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

SE.MAF <- read.table("SE132.pop.frq.strat", header=T)
SE.MAF2 <- SE.MAF[,c(3,2,6)]
summary(SE.MAF2)
         CLST               SNP             MAF        
 X1.Sk.Ho  : 2081   100865:23 :   15   Min.   :0.0000  
 X1.Sk.SF  : 2081   101108:100:   15   1st Qu.:0.0000  
 X1.Sk.SL  : 2081   101142:72 :   15   Median :0.1111  
 X2.Upp.Gra: 2081   101270:3  :   15   Mean   :0.2028  
 X2.Upp.K  : 2081   101367:84 :   15   3rd Qu.:0.3333  
 X2.Upp.O  : 2081   101609:30 :   15   Max.   :1.0000  
 (Other)   :18729   (Other)   :31125     

library("ggplot2")
library("reshape2")

SE.MAF3 <- melt(SE.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.MAF3)
head(SE.MAF3)


SE132.MAF4 <- dcast(SE.MAF3, formula= CLST ~ SNP)
head(SE145.MAF4)
```


```
##Add X infront of all locusnames. 
colnames(SE132.MAF4) <- paste("X", colnames(SE132.MAF4), sep=".")
write.csv(SE132.MAF4, file="SE.132.MAF.csv")
```


####2. Geographic coordinates




####3. Climate variables

```
1. solar.rad.60d (amt of solar radiation 60 days after egg laying)

2. pcpt.60d (precipitation in the 60days after egg laying)

3. day10cm (calender day on which only 10cm of snow remains)

4. temp.laying.date

5. shadow.days
```





###Run RDA

####1. CHall

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


####2. CHN

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


####3. CHS

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

####4. CHS.VS

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

####5. CHS.TI

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

####6. CZ

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

