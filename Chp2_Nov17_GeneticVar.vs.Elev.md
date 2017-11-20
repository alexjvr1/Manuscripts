# Genetic variation across elevation

## Aim: 

1. Determine how genetic variation is distributed across elevation 

2. Proportion of genetic variation explained by geographic dist vs elevation 

I'm comparing results from RAD and mtDNA for all analyses to get an impression of historic vs contemporary geneflow

## Methods

### 1. Linear model of diversity measures vs elevation and elevation^2 and Transect

  - Are they significantly different between transects?
  
  - Is diversity distributed linearly across elevation? 
  
  - Figures: fit lm (as with the environmental variables before)
  
  
### 2. RDA to partition variance

  - RDA of detrended Fst vs distance + elevation 
  
  - pRDA of detrended Fst vs elevation
  
  - pRDA of detrended Fst vs distance
  
compare mtDNA vs RAD (contemporary vs historic geneflow)


## Analyses

One of the measures of diversity is estimated populations size (Ne). I'll compare this to the estimated census size based on nr of 
clutches counted at each site. 

To calculate Ne from RAD data: 




### 1. Linear models

I will use a linear model to test the whether region or elevation have a significant effect on each measure of genetic diversity.

I'm running both a straighforward linear model, and a model including the squared term of elevation for possible quadratic 
(i.e. non-linear) relationships between diversity and elevation.


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GeneticVar.vs.Elev


#### RAD: 

Measures of diversity: Fis, Hs, Ho, Ne, clutch size

First I have to calculate the Ne for each population: 


I'm using NeEstimator for this: http://www.molecularfisherieslaboratory.com.au/neestimator-software/

This estimates genetic effective population size from allele frequency data.

I converted to Genepop format using pgdspider. I first tried with only the 1000 Neutral loci, and then with the full 9608 loci. 
I used the LD and excess Het estimate with random mating. All analyses gave me Inf estimates of Ne. So I'm leaving this out. 


First add the centered and square term for elevation
```
RAD.data <- read.table("RAD.diversityMeasures", header=T)

RAD.data$elev.c <- RAD.data$Elev - mean(RAD.data$Elev) ##center the data
RAD.data$elev.c <- RAD.data$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

RAD.data$elev.c.squared <- RAD.data$elev.c^2
```

Run the models

```
m1.RAD.Fis <- lm(RAD.data$FIS~RAD.data$elev.c + RAD.data$Region + RAD.data$Region:RAD.data$elev.c)

m1.RAD.Hs <- lm(RAD.data$Hs~RAD.data$elev.c + RAD.data$Region + RAD.data$Region:RAD.data$elev.c)

m1.RAD.Ho <- lm(RAD.data$Ho~RAD.data$elev.c + RAD.data$Region + RAD.data$Region:RAD.data$elev.c)

m2.RAD.Fis <- lm(RAD.data$FIS~RAD.data$elev.c +RAD.data$Region + RAD.data$elev.c.squared + RAD.data$Region:RAD.data$elev.c + RAD.data$Region:RAD.data$elev.c.squared)

m2.RAD.Hs <- lm(RAD.data$Hs~RAD.data$elev.c + RAD.data$Region + RAD.data$elev.c.squared + RAD.data$Region:RAD.data$elev.c + RAD.data$Region:RAD.data$elev.c.squared)

m2.RAD.Ho <- lm(RAD.data$Ho~RAD.data$elev.c + RAD.data$Region + RAD.data$elev.c.squared + RAD.data$Region:RAD.data$elev.c + RAD.data$Region:RAD.data$elev.c.squared)

```



Select model
```
AIC(m1.RAD.Fis)
[1] -188.5216
AIC(m2.RAD.Fis)
[1] -182.84

AIC(m1.RAD.Hs)
[1] -395.9505
AIC(m2.RAD.Hs)
[1] -389.7367

AIC(m1.RAD.Ho)
[1] -369.4961
AIC(m2.RAD.Ho)
[1] -366.6641
```

All m1 models

Assess model
```
library(car)

Anova(m1.RAD.Fis)

Anova Table (Type II tests)

Response: RAD.data$FIS
                                 Sum Sq Df F value  Pr(>F)   
RAD.data$elev.c                 0.04561  1  9.5613 0.00284 **
RAD.data$Region                 0.02148  3  1.5010 0.22172   
RAD.data$elev.c:RAD.data$Region 0.00630  3  0.4400 0.72509   
Residuals                       0.33870 71                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Anova(m1.RAD.Hs)

Anova Table (Type II tests)

Response: RAD.data$Hs
                                   Sum Sq Df F value   Pr(>F)    
RAD.data$elev.c                 0.0003684  1  1.0666  0.30521    
RAD.data$Region                 0.0183735  3 17.7348 1.08e-08 ***
RAD.data$elev.c:RAD.data$Region 0.0029904  3  2.8865  0.04154 *  
Residuals                       0.0245189 71                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Anova(m1.RAD.Ho)
Anova Table (Type II tests)

Response: RAD.data$Ho
                                  Sum Sq Df F value    Pr(>F)    
RAD.data$elev.c                 0.001356  1  2.8095   0.09811 .  
RAD.data$Region                 0.020559  3 14.1975 2.418e-07 ***
RAD.data$elev.c:RAD.data$Region 0.001513  3  1.0451   0.37803    
Residuals                       0.034271 71                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

And plot
```
library(broom)
library(ggplot2)



##RAD.Fis
Fis.fit.RAD <- lm(formula = RAD.data$FIS ~ RAD.data$elev.c + 
    RAD.data$Region + 
    RAD.data$Region:RAD.data$elev.c)


fit.Fis.RAD <- augment(Fis.fit.RAD)
fit.Fis.RAD$colour <- gsub("CHS.TI", "thistle3", fit.Fis.RAD$RAD.data.Region)
fit.Fis.RAD$colour <- gsub("CHS.VS", "mediumpurple4", fit.Fis.RAD$colour)
fit.Fis.RAD$colour <- gsub("CHN", "palegreen3", fit.Fis.RAD$colour)
fit.Fis.RAD$colour <- gsub("CZ", "skyblue3", fit.Fis.RAD$colour)

fit.Fis.RAD$elev <- RAD.data$Elev

p1.Fis.RAD <- ggplot(fit.Fis.RAD, aes(x=elev, y=RAD.data.FIS , group=RAD.data.Region)) + geom_point(color=fit.Fis.RAD$colour) + geom_line(aes(x=elev, y=.fitted, colour=RAD.data.Region, group=RAD.data.Region, color=fit.Fis.RAD$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CZ"="skyblue3", "CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))



##RAD.Hs

Hs.fit.RAD <- lm(formula = RAD.data$Hs ~ RAD.data$elev.c + 
    RAD.data$Region + 
    RAD.data$Region:RAD.data$elev.c)


fit.Hs.RAD <- augment(Hs.fit.RAD)
fit.Hs.RAD$colour <- gsub("CHS.TI", "thistle3", fit.Hs.RAD$RAD.data.Region)
fit.Hs.RAD$colour <- gsub("CHS.VS", "mediumpurple4", fit.Hs.RAD$colour)
fit.Hs.RAD$colour <- gsub("CHN", "palegreen3", fit.Hs.RAD$colour)
fit.Hs.RAD$colour <- gsub("CZ", "skyblue3", fit.Hs.RAD$colour)

fit.Hs.RAD$elev <- RAD.data$Elev

p2.Hs.RAD <- ggplot(fit.Hs.RAD, aes(x=elev, y=RAD.data.Hs , group=RAD.data.Region)) + geom_point(color=fit.Hs.RAD$colour) + geom_line(aes(x=elev, y=.fitted, colour=RAD.data.Region, group=RAD.data.Region, color=fit.Hs.RAD$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CZ"="skyblue3", "CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))



##RAD.Ho

Ho.fit.RAD <- lm(formula = RAD.data$Ho ~ RAD.data$elev.c + 
    RAD.data$Region + 
    RAD.data$Region:RAD.data$elev.c)


fit.Ho.RAD <- augment(Ho.fit.RAD)
fit.Ho.RAD$colour <- gsub("CHS.TI", "thistle3", fit.Ho.RAD$RAD.data.Region)
fit.Ho.RAD$colour <- gsub("CHS.VS", "mediumpurple4", fit.Ho.RAD$colour)
fit.Ho.RAD$colour <- gsub("CHN", "palegreen3", fit.Ho.RAD$colour)
fit.Ho.RAD$colour <- gsub("CZ", "skyblue3", fit.Ho.RAD$colour)

fit.Ho.RAD$elev <- RAD.data$Elev

p3.Ho.RAD <- ggplot(fit.Ho.RAD, aes(x=elev, y=RAD.data.Ho , group=RAD.data.Region)) + geom_point(color=fit.Ho.RAD$colour) + geom_line(aes(x=elev, y=.fitted, colour=RAD.data.Region, group=RAD.data.Region, color=fit.Ho.RAD$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CZ"="skyblue3", "CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))
```

multiplot
```
pdf("RADvsElev.pdf")
multiplot(p1.Fis.RAD, p2.Hs.RAD, cols=1)
dev.off()


pdf("RADvsElev.pdf")
multiplot(p1.Fis.RAD, p3.Ho.RAD, cols=1)
dev.off()

```



#### mtDNA

Haplotype diversity, nucleotide diversity

First add the centered and square term for elevation
```
mtDNA.data <- read.table("mtDNA.diversityMeasures.67pops", header=T)

mtDNA.data$elev.c <- mtDNA.data$Elev - mean(mtDNA.data$Elev) ##center the data
mtDNA.data$elev.c <- mtDNA.data$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

mtDNA.data$elev.c.squared <- mtDNA.data$elev.c^2
```

Run the models

```
m1.mtDNA.Hd <- lm(mtDNA.data$Hd~mtDNA.data$elev.c + mtDNA.data$Region + mtDNA.data$Region:mtDNA.data$elev.c)

m1.mtDNA.nd <- lm(mtDNA.data$nd~mtDNA.data$elev.c + mtDNA.data$Region + mtDNA.data$Region:mtDNA.data$elev.c)

m2.mtDNA.Hd <- lm(mtDNA.data$Hd~mtDNA.data$elev.c + mtDNA.data$Region + mtDNA.data$elev.c.squared + mtDNA.data$Region:mtDNA.data$elev.c + mtDNA.data$Region:mtDNA.data$elev.c.squared)

m2.mtDNA.nd <- lm(mtDNA.data$nd~mtDNA.data$elev.c + mtDNA.data$Region + mtDNA.data$elev.c.squared + mtDNA.data$Region:mtDNA.data$elev.c + mtDNA.data$Region:mtDNA.data$elev.c.squared)
```



Select model
```
AIC(m1.mtDNA.Hd)
[1] 41.18423
AIC(m1.mtDNA.nd)
[1] -557.0955
AIC(m2.mtDNA.Hd)
[1] 48.56382
AIC(m2.mtDNA.nd)
[1] -550.5264

```

Assess model
```
library(car)

Anova(m1.mtDNA.Hd)

Anova Table (Type II tests)

Response: mtDNA.data$Hd
                                    Sum Sq Df F value  Pr(>F)  
mtDNA.data$elev.c                   0.1906  1  2.0280 0.15969  
mtDNA.data$Region                   0.7280  3  2.5823 0.06183 .
mtDNA.data$elev.c:mtDNA.data$Region 0.7238  3  2.5673 0.06294 .
Residuals                           5.5447 59                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Anova(m1.mtDNA.nd)

Anova Table (Type II tests)

Response: mtDNA.data$nd
                                        Sum Sq Df F value Pr(>F)  
mtDNA.data$elev.c                   0.00000411  1  0.3301 0.5678  
mtDNA.data$Region                   0.00011323  3  3.0329 0.0362 *
mtDNA.data$elev.c:mtDNA.data$Region 0.00001785  3  0.4782 0.6987  
Residuals                           0.00073422 59                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

And plot
```
library(broom)
library(ggplot2)



##mtDNA.Hd
fit.Hd.mtDNA$elev <- mtDNA.data$Elev

Hd.fit.mtDNA <- lm(formula = mtDNA.data$Hd ~ mtDNA.data$elev.c + 
    mtDNA.data$Region + 
    mtDNA.data$Region:mtDNA.data$elev.c)


fit.Hd.mtDNA <- augment(Hd.fit.mtDNA)
fit.Hd.mtDNA$colour <- gsub("CHS.TI", "thistle3", fit.Hd.mtDNA$mtDNA.data.Region)
fit.Hd.mtDNA$colour <- gsub("CHS.VS", "mediumpurple4", fit.Hd.mtDNA$colour)
fit.Hd.mtDNA$colour <- gsub("CHN", "palegreen3", fit.Hd.mtDNA$colour)
fit.Hd.mtDNA$colour <- gsub("CZ", "skyblue3", fit.Hd.mtDNA$colour)


p1.Hd.mtDNA <- ggplot(fit.Hd.mtDNA, aes(x=elev, y=mtDNA.data.Hd , group=mtDNA.data.Region)) + geom_point(color=fit.Hd.mtDNA$colour) + geom_line(aes(x=elev, y=.fitted, colour=mtDNA.data.Region, group=mtDNA.data.Region, color=fit.Hd.mtDNA$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CZ"="skyblue3", "CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))


##mtDNA.nd

fit.nd.mtDNA$elev <- mtDNA.data$Elev

nd.fit.mtDNA <- lm(formula = mtDNA.data$nd ~ mtDNA.data$elev.c + 
    mtDNA.data$Region + 
    mtDNA.data$Region:mtDNA.data$elev.c)


fit.nd.mtDNA <- augment(nd.fit.mtDNA)
fit.nd.mtDNA$colour <- gsub("CHS.TI", "thistle3", fit.nd.mtDNA$mtDNA.data.Region)
fit.nd.mtDNA$colour <- gsub("CHS.VS", "mediumpurple4", fit.nd.mtDNA$colour)
fit.nd.mtDNA$colour <- gsub("CHN", "palegreen3", fit.nd.mtDNA$colour)
fit.nd.mtDNA$colour <- gsub("CZ", "skyblue3", fit.nd.mtDNA$colour)


p2.nd.mtDNA <- ggplot(fit.nd.mtDNA, aes(x=elev, y=mtDNA.data.nd , group=mtDNA.data.Region)) + geom_point(color=fit.nd.mtDNA$colour) + geom_line(aes(x=elev, y=.fitted, colour=mtDNA.data.Region, group=mtDNA.data.Region, color=fit.nd.mtDNA$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CZ"="skyblue3", "CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))


```

multiplot
```
pdf("mtDNAvsElev.pdf")
multiplot(p1.Hd.mtDNA, p2.nd.mtDNA, cols=1)
dev.off()
```

### 2. RDA

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA/Sept2017/AnalysisFullDataset

#### RADdata

MAF input files have been prepared before: https://github.com/alexjvr1/Manuscripts/blob/master/CHP2_CHdata_RDA_Oct2017.md

I have to add elevation and geographic distance. Since I want to understand how genetic diversity is distributed up elevation,
I will measure all distances from the lowest elevation population

CHN: rade

CHS.VS: orge

CHS.TI: scai

CZ: ente


First I need to calculate the geographic distance between all populations. And then take the distance from the lowest elevation and most north/most southern population. (see above)
```
library(fields)
library(reshape2)

CHNData <- read.csv("CHN229.MAF.csv", header=T)

CHN_long.lat <- cbind(CHNData$long, CHNData$lat)

CHN.dist.matrix <- rdist.earth(CHN_long.lat, miles=F)
CHN.m.dist <- as.matrix(CHN.dist.matrix)

colnames(CHN.m.dist) <- (CHNData$X.CLST)
rownames(CHN.m.dist) <- (CHNData$X.CLST)

head(CHN.m.dist)
CHN.m2.dist <- melt(CHN.m.dist)[melt(upper.tri(CHN.m.dist))$value,]

names(CHN.m2.dist) <- c("c1", "c2", "distance")
CHN.m2.dist$log.km <- log(CHN.m2.dist$distance)

CHN.dist.rade <- subset(CHN.m2.dist, c1==c("rade"), select=c1:log.km)

new.row <- c("rade", "rade", "0", "0")
CHN.dist.rade.new <- rbind(new.row,CHN.dist.rade)
CHNData$X.CLST
CHN.dist.rade.new$c2
CHNData$log.km <- CHN.dist.rade.new$log.km
```


CHS.TI
```
CHS.TIData <- read.csv("CHS.TI140.MAF.csv", header=T)

CHS.TI_long.lat <- cbind(CHS.TIData$long, CHS.TIData$lat)

CHS.TI.dist.matrix <- rdist.earth(CHS.TI_long.lat, miles=F)
CHS.TI.m.dist <- as.matrix(CHS.TI.dist.matrix)

colnames(CHS.TI.m.dist) <- (CHS.TIData$X.CLST)
rownames(CHS.TI.m.dist) <- (CHS.TIData$X.CLST)

head(CHS.TI.m.dist)
CHS.TI.m2.dist <- melt(CHS.TI.m.dist)[melt(upper.tri(CHS.TI.m.dist))$value,]

names(CHS.TI.m2.dist) <- c("c1", "c2", "distance")
CHS.TI.m2.dist$log.km <- log(CHS.TI.m2.dist$distance)

CHS.TI.dist.scai <- subset(CHS.TI.m2.dist, c1==c("scai"), select=c1:log.km)

new.row <- c("scai", "scai", "0", "0")
CHS.TI.dist.scai.new <- rbind(new.row,CHS.TI.dist.scai)
CHS.TIData$X.CLST
CHS.TI.dist.scai.new$c2
CHS.TIData$log.km <- CHS.TI.dist.scai.new$log.km
```


CHS.VS
```
CHS.VSData <- read.csv("CHS.VS135.MAF.csv", header=T)

CHS.VS_long.lat <- cbind(CHS.VSData$long, CHS.VSData$lat)

CHS.VS.dist.matrix <- rdist.earth(CHS.VS_long.lat, miles=F)
CHS.VS.m.dist <- as.matrix(CHS.VS.dist.matrix)

colnames(CHS.VS.m.dist) <- (CHS.VSData$X.CLST)
rownames(CHS.VS.m.dist) <- (CHS.VSData$X.CLST)

head(CHS.VS.m.dist)
CHS.VS.m2.dist <- melt(CHS.VS.m.dist)[melt(upper.tri(CHS.VS.m.dist))$value,]

names(CHS.VS.m2.dist) <- c("c1", "c2", "distance")
CHS.VS.m2.dist$log.km <- log(CHS.VS.m2.dist$distance)

CHS.VS.dist.orge <- subset(CHS.VS.m2.dist, c1==c("orge"), select=c1:log.km)

new.row <- c("orge", "orge", "0", "0")
CHS.VS.dist.orge.new <- rbind(new.row,CHS.VS.dist.orge)
CHS.VSData$X.CLST
CHS.VS.dist.orge.new$c2
CHS.VSData$log.km <- CHS.VS.dist.orge.new$log.km
```

CZ
```
CZData <- read.csv("CZ404.MAF.csv", header=T)

CZ_long.lat <- cbind(CZData$long, CZData$lat)

CZ.dist.matrix <- rdist.earth(CZ_long.lat, miles=F)
CZ.m.dist <- as.matrix(CZ.dist.matrix)

colnames(CZ.m.dist) <- (CZData$X.CLST)
rownames(CZ.m.dist) <- (CZData$X.CLST)

head(CZ.m.dist)
CZ.m2.dist <- melt(CZ.m.dist)[melt(upper.tri(CZ.m.dist))$value,]

names(CZ.m2.dist) <- c("c1", "c2", "distance")
CZ.m2.dist$log.km <- log(CZ.m2.dist$distance)

CZ.dist.ente <- subset(CZ.m2.dist, c1==c("ente"), select=c1:log.km)

new.row <- c("ente", "ente", "0", "0")
CZ.dist.ente.new <- rbind(new.row,CZ.dist.ente)
CZData$X.CLST
CZ.dist.ente.new$c2
CZData$log.km <- CZ.dist.ente.new$log.km
```


##### Run RDA models

###### CHN
```
library(vegan)

GenData <- CHNData[,15:5279]
GenData <- decostand(GenData, "hellinger")

Climate.Data <- CHNData$Elev
Climate.Data
Climate.Data <- as.data.frame(Climate.Data)
colnames(Climate.Data) <- "elev"
Climate.Data$log.km <- CHNData$log.km
Climate.Data <- as.data.frame(Climate.Data)
class(Climate.Data$log.km) <- "numeric"
Climate.Data <- as.data.frame(Climate.Data)
summary(Climate.Data)
```

CHN Full RDA
```
RDA.CHN.dist.elev <- rda(GenData ~ elev + log.km, Climate.Data)

RDA.CHN.dist.elev

Call: rda(formula = GenData ~ elev + log.km, data = Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Constrained   0.03425    0.15940    2
Unconstrained 0.18061    0.84060   16
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2 
0.022545 0.011702 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025321 0.023201 0.017360 0.015786 0.013331 0.011585 0.010067 0.009347 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16 
0.009033 0.008145 0.007880 0.006800 0.006685 0.006050 0.005194 0.004822 


RsquareAdj(RDA.CHN.dist.elev)

$r.squared
[1] 0.1593952

$adj.r.squared
[1] 0.05431955


anova(RDA.CHN.dist.elev)


Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + log.km, data = Climate.Data)
         Df Variance     F Pr(>F)   
Model     2 0.034247 1.517  0.004 **
Residual 16 0.180607                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

pCHN.distOnly
```
pCHN.distOnly <- rda(GenData ~ log.km + Condition(elev), Climate.Data) 

pCHN.distOnly

Call: rda(formula = GenData ~ log.km + Condition(elev), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Conditional   0.02028    0.09439    1
Constrained   0.01397    0.06501    1
Unconstrained 0.18061    0.84060   16
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.013967 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025321 0.023201 0.017360 0.015786 0.013331 0.011585 0.010067 0.009347 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16 
0.009033 0.008145 0.007880 0.006800 0.006685 0.006050 0.005194 0.004822

RsquareAdj(pCHN.distOnly)

$r.squared
[1] 0.06500662

$adj.r.squared
[1] 0.01320228

anova(pCHN.distOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ log.km + Condition(elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     1 0.013967 1.2373  0.155
Residual 16 0.180607     

```

pCHN.elevOnly
```
pCHN.elevOnly <- rda(GenData ~elev + Condition(log.km), Climate.Data)

pCHN.elevOnly

Call: rda(formula = GenData ~ elev + Condition(log.km), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21485    1.00000     
Conditional   0.01962    0.09131    1
Constrained   0.01463    0.06809    1
Unconstrained 0.18061    0.84060   16
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.014628 

Eigenvalues for unconstrained axes:
     PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8 
0.025321 0.023201 0.017360 0.015786 0.013331 0.011585 0.010067 0.009347 
     PC9     PC10     PC11     PC12     PC13     PC14     PC15     PC16 
0.009033 0.008145 0.007880 0.006800 0.006685 0.006050 0.005194 0.004822 




RsquareAdj(pCHN.elevOnly)

$r.squared
[1] 0.06808528

$adj.r.squared
[1] 0.01646204

anova(pCHN.elevOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Condition(log.km), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     1 0.014628 1.2959  0.077 .
Residual 16 0.180607                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
```


###### CHS.TI
```
library(vegan)

GenData <- CHS.TIData[,15:5706]
GenData <- decostand(GenData, "hellinger")

Climate.Data <- CHS.TIData$elev
Climate.Data
Climate.Data <- as.data.frame(Climate.Data)
colnames(Climate.Data) <- "elev"
Climate.Data$log.km <- CHS.TIData$log.km
Climate.Data <- as.data.frame(Climate.Data)
class(Climate.Data$log.km) <- "numeric"
Climate.Data <- as.data.frame(Climate.Data)
summary(Climate.Data)
```

CHS.TI Full RDA
```
RDA.CHS.TI.dist.elev <- rda(GenData ~ elev + log.km, Climate.Data)

RDA.CHS.TI.dist.elev

Call: rda(formula = GenData ~ elev + log.km, data = Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Constrained   0.05668    0.19602    2
Unconstrained 0.23248    0.80398   11
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1    RDA2 
0.03314 0.02354 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.05747 0.03864 0.02596 0.02037 0.01727 0.01603 0.01502 0.01339 0.01149 0.00921 
   PC11 
0.00763 




RsquareAdj(RDA.CHS.TI.dist.elev)

$r.squared
[1] 0.1960223

$adj.r.squared
[1] 0.04984454


anova(RDA.CHS.TI.dist.elev)


Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + log.km, data = Climate.Data)
         Df Variance     F Pr(>F)  
Model     2 0.056681 1.341  0.085 .
Residual 11 0.232477               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

pCHS.TI.distOnly
```
pCHS.TI.distOnly <- rda(GenData ~ log.km + Condition(elev), Climate.Data) 

pCHS.TI.distOnly

Call: rda(formula = GenData ~ log.km + Condition(elev), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Conditional   0.02835    0.09805    1
Constrained   0.02833    0.09797    1
Unconstrained 0.23248    0.80398   11
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1 
0.02833 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.05747 0.03864 0.02596 0.02037 0.01727 0.01603 0.01502 0.01339 0.01149 0.00921 
   PC11 
0.00763 



RsquareAdj(pCHS.TI.distOnly)

$r.squared
[1] 0.09797342

$adj.r.squared
[1] 0.02695825

anova(pCHS.TI.distOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ log.km + Condition(elev), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     1  0.02833 1.3405  0.144
Residual 11  0.23248     

```

pCHS.TI.elevOnly
```
pCHS.TI.elevOnly <- rda(GenData ~elev + Condition(log.km), Climate.Data)

pCHS.TI.elevOnly

Call: rda(formula = GenData ~ elev + Condition(log.km), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.28916    1.00000     
Conditional   0.03303    0.11422    1
Constrained   0.02365    0.08181    1
Unconstrained 0.23248    0.80398   11
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.023655 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
0.05747 0.03864 0.02596 0.02037 0.01727 0.01603 0.01502 0.01339 0.01149 0.00921 
   PC11 
0.00763 




RsquareAdj(pCHS.TI.elevOnly)

$r.squared
[1] 0.081806

$adj.r.squared
[1] 0.009443545


anova(pCHS.TI.elevOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Condition(log.km), data = Climate.Data)
         Df Variance      F Pr(>F)
Model     1 0.023655 1.1193  0.276
Residual 11 0.232477        
```

###### CHS.VS
```
library(vegan)

GenData <- CHS.VSData[,13:5847]
GenData <- decostand(GenData, "hellinger")

Climate.Data <- CHS.VSData$elev
Climate.Data
Climate.Data <- as.data.frame(Climate.Data)
colnames(Climate.Data) <- "elev"
Climate.Data$log.km <- CHS.VSData$log.km
Climate.Data <- as.data.frame(Climate.Data)
class(Climate.Data$log.km) <- "numeric"
Climate.Data <- as.data.frame(Climate.Data)
summary(Climate.Data)
```

CHS.VS Full RDA
```
RDA.CHS.VS.dist.elev <- rda(GenData ~ elev + log.km, Climate.Data)

RDA.CHS.VS.dist.elev

Call: rda(formula = GenData ~ elev + log.km, data = Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Constrained   0.05477    0.25379    2
Unconstrained 0.16105    0.74621    7
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2 
0.029482 0.025292 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7 
0.04447 0.02786 0.02539 0.02101 0.01720 0.01417 0.01095 


RsquareAdj(RDA.CHS.VS.dist.elev)

$r.squared
[1] 0.2537949

$adj.r.squared
[1] 0.0405935

anova(RDA.CHS.VS.dist.elev)


Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + log.km, data = Climate.Data)
         Df Variance      F Pr(>F)
Model     2 0.054774 1.1904  0.171
Residual  7 0.161045               

```

pCHS.VS.distOnly
```
pCHS.VS.distOnly <- rda(GenData ~ log.km + Condition(elev), Climate.Data) 

pCHS.VS.distOnly

Call: rda(formula = GenData ~ log.km + Condition(elev), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Conditional   0.02737    0.12684    1
Constrained   0.02740    0.12696    1
Unconstrained 0.16105    0.74621    7
Inertia is variance 

Eigenvalues for constrained axes:
  RDA1 
0.0274 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7 
0.04447 0.02786 0.02539 0.02101 0.01720 0.01417 0.01095 



RsquareAdj(pCHS.VS.distOnly)

$r.squared
[1] 0.1269586

$adj.r.squared
[1] 0.02290261

anova(pCHS.VS.distOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ log.km + Condition(elev), data = Climate.Data)
         Df Variance     F Pr(>F)
Model     1  0.02740 1.191  0.248
Residual  7  0.16104     

```

pCHS.VS.elevOnly
```
pCHS.VS.elevOnly <- rda(GenData ~elev + Condition(log.km), Climate.Data)

pCHS.VS.elevOnly

Call: rda(formula = GenData ~ elev + Condition(log.km), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.21582    1.00000     
Conditional   0.02581    0.11959    1
Constrained   0.02896    0.13421    1
Unconstrained 0.16105    0.74621    7
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.028965 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7 
0.04447 0.02786 0.02539 0.02101 0.01720 0.01417 0.01095 




RsquareAdj(pCHS.VS.elevOnly)

$r.squared
[1] 0.1342098

$adj.r.squared
[1] 0.03106016


anova(pCHS.VS.elevOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Condition(log.km), data = Climate.Data)
         Df Variance     F Pr(>F)
Model     1 0.028965 1.259  0.178
Residual  7 0.161045     
```


###### CZ
```
library(vegan)

GenData <- CZData[,18:7305]
GenData <- decostand(GenData, "hellinger")

Climate.Data <- CZData$elev
Climate.Data
Climate.Data <- as.data.frame(Climate.Data)
colnames(Climate.Data) <- "elev"
Climate.Data$log.km <- CZData$log.km
Climate.Data <- as.data.frame(Climate.Data)
class(Climate.Data$log.km) <- "numeric"
Climate.Data <- as.data.frame(Climate.Data)
summary(Climate.Data)
```

CZ Full RDA
```
RDA.CZ.dist.elev <- rda(GenData ~ elev + log.km, Climate.Data)

RDA.CZ.dist.elev

Call: rda(formula = GenData ~ elev + log.km, data = Climate.Data)

              Inertia Proportion Rank
Total         0.25682    1.00000     
Constrained   0.03628    0.14128    2
Unconstrained 0.22053    0.85872   34
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1     RDA2 
0.028764 0.007520 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
0.04093 0.01756 0.01277 0.01100 0.01034 0.00823 0.00777 0.00719 
(Showed only 8 of all 34 unconstrained eigenvalues)



RsquareAdj(RDA.CZ.dist.elev)

$r.squared
[1] 0.1412819

$adj.r.squared
[1] 0.09076911

anova(RDA.CZ.dist.elev)


Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + log.km, data = Climate.Data)
         Df Variance     F Pr(>F)    
Model     2 0.036284 2.797  0.001 ***
Residual 34 0.220534                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

pCZ.distOnly
```
pCZ.distOnly <- rda(GenData ~ log.km + Condition(elev), Climate.Data) 

pCZ.distOnly

Call: rda(formula = GenData ~ log.km + Condition(elev), data =
Climate.Data)

              Inertia Proportion Rank
Total         0.25682    1.00000     
Conditional   0.01903    0.07411    1
Constrained   0.01725    0.06717    1
Unconstrained 0.22053    0.85872   34
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.017251 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
0.04093 0.01756 0.01277 0.01100 0.01034 0.00823 0.00777 0.00719 
(Showed only 8 of all 34 unconstrained eigenvalues)

RsquareAdj(pCZ.distOnly)

$r.squared
[1] 0.06717223

$adj.r.squared
[1] 0.04311341

anova(pCZ.distOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ log.km + Condition(elev), data = Climate.Data)
         Df Variance      F Pr(>F)   
Model     1 0.017251 2.6596  0.002 **
Residual 34 0.220534                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

pCZ.elevOnly
```
pCZ.elevOnly <- rda(GenData ~elev + Condition(log.km), Climate.Data)

pCZ.elevOnly

Call: rda(formula = GenData ~ elev + Condition(log.km), data =
Climate.Data)

               Inertia Proportion Rank
Total         0.256818   1.000000     
Conditional   0.026948   0.104930    1
Constrained   0.009336   0.036352    1
Unconstrained 0.220534   0.858718   34
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.009336 

Eigenvalues for unconstrained axes:
    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
0.04093 0.01756 0.01277 0.01100 0.01034 0.00823 0.00777 0.00719 
(Showed only 8 of all 34 unconstrained eigenvalues)




RsquareAdj(pCZ.elevOnly)

$r.squared
[1] 0.03635195

$adj.r.squared
[1] 0.01141255


anova(pCZ.elevOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ elev + Condition(log.km), data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     1 0.009336 1.4393   0.06 .
Residual 34 0.220534                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

#### mtDNA






define multiplot
```
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
```
