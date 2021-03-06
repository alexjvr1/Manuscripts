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


### 3. Fst across elevation 

Isolation by elevation: 

Comparison vs distance
 
Comparison between mtDNA and RAD data

Mantel and partial mantel tests to determine if Fst is correlated with elevation or distance. 

### 4. Connectivity within each region

I'll use a network analysis to visualise how connected populations are across elevation. i.e. are high elevation populations isolated? 



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


###### Test correlation and p-value

I can get R-values for the correlation and test whether they are significant


CHN
```
library(Hmisc)

CHN.elev <- subset(RAD.data, Region=="CHN",)

cor.test(CHN.elev$elev.c, CHN.elev$FIS, type="pearson")

	Pearson's product-moment correlation

data:  CHN.elev$elev.c and CHN.elev$FIS
t = -0.84295, df = 17, p-value = 0.411
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.5999302  0.2793195
sample estimates:
      cor 
-0.200302 

cor.test(CHN.elev$elev.c, CHN.elev$Hs, type="pearson")

	Pearson's product-moment correlation

data:  CHN.elev$elev.c and CHN.elev$Hs
t = -0.25738, df = 17, p-value = 0.8
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.5022970  0.4033207
sample estimates:
        cor 
-0.06230177 

cor.test(CHN.elev$elev.c, CHN.elev$Ho, type="pearson")

	Pearson's product-moment correlation

data:  CHN.elev$elev.c and CHN.elev$Ho
t = 0.21188, df = 17, p-value = 0.8347
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.4125038  0.4940147
sample estimates:
       cor 
0.05132113 
```


CHS.VS
```
CHS.VS.elev <- subset(RAD.data, Region=="CHS.VS",)

cor.test(CHS.VS.elev$elev.c, CHS.VS.elev$FIS, type="pearson")

	Pearson's product-moment correlation

data:  CHS.VS.elev$elev.c and CHS.VS.elev$FIS
t = -1.7315, df = 8, p-value = 0.1216
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.8667912  0.1601732
sample estimates:
       cor 
-0.5221073 


cor.test(CHS.VS.elev$elev.c, CHS.VS.elev$Hs, type="pearson")

	Pearson's product-moment correlation

data:  CHS.VS.elev$elev.c and CHS.VS.elev$Hs
t = -4.1896, df = 8, p-value = 0.003039
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.9583366 -0.4165586
sample estimates:
       cor 
-0.8288082 

cor.test(CHS.VS.elev$elev.c, CHS.VS.elev$Ho, type="pearson")

	Pearson's product-moment correlation

data:  CHS.VS.elev$elev.c and CHS.VS.elev$Ho
t = -0.11241, df = 8, p-value = 0.9133
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6530097  0.6050437
sample estimates:
        cor 
-0.03971036 

```


CHS.TI
```
CHS.TI.elev <- subset(RAD.data, Region=="CHS.TI",)

cor.test(CHS.TI.elev$elev.c, CHS.TI.elev$FIS, type="pearson")

	Pearson's product-moment correlation

data:  CHS.TI.elev$elev.c and CHS.TI.elev$FIS
t = -0.99518, df = 13, p-value = 0.3378
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6849688  0.2850503
sample estimates:
       cor 
-0.2660633 


cor.test(CHS.TI.elev$elev.c, CHS.TI.elev$Hs, type="pearson")

	Pearson's product-moment correlation

data:  CHS.TI.elev$elev.c and CHS.TI.elev$Hs
t = 1.4697, df = 13, p-value = 0.1654
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.1671183  0.7455610
sample estimates:
      cor 
0.3774581 


cor.test(CHS.TI.elev$elev.c, CHS.TI.elev$Ho, type="pearson")

	Pearson's product-moment correlation

data:  CHS.TI.elev$elev.c and CHS.TI.elev$Ho
t = 2.5228, df = 13, p-value = 0.02547
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.08640765 0.83912519
sample estimates:
      cor 
0.5732945 

```


CZ
```
CZ.elev <- subset(RAD.data, Region=="CZ",)

cor.test(CZ.elev$elev.c, CZ.elev$FIS, type="pearson")

	Pearson's product-moment correlation

data:  CZ.elev$elev.c and CZ.elev$FIS
t = -2.3069, df = 33, p-value = 0.02748
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.6279237 -0.0449986
sample estimates:
       cor 
-0.3726569 



cor.test(CZ.elev$elev.c, CZ.elev$Hs, type="pearson")

	Pearson's product-moment correlation

data:  CZ.elev$elev.c and CZ.elev$Hs
t = -1.3039, df = 33, p-value = 0.2013
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.5164922  0.1208150
sample estimates:
       cor 
-0.2213431 


cor.test(CZ.elev$elev.c, CZ.elev$Ho, type="pearson")

	Pearson's product-moment correlation

data:  CZ.elev$elev.c and CZ.elev$Ho
t = 0.73267, df = 33, p-value = 0.4689
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.2158301  0.4411627
sample estimates:
      cor 
0.1265161 

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


###### test for correlation

CHN
```
CHN.elev <- subset(mtDNA.data, Region=="CHN",)

cor.test(CHN.elev$elev.c, CHN.elev$Hd, type="pearson")

	Pearson's product-moment correlation

data:  CHN.elev$elev.c and CHN.elev$Hd
t = 2.5595, df = 13, p-value = 0.02376
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.09466821 0.84157228
sample estimates:
      cor 
0.5788593 

cor.test(CHN.elev$elev.c, CHN.elev$nd, type="pearson")

	Pearson's product-moment correlation

data:  CHN.elev$elev.c and CHN.elev$nd
t = 2.5032, df = 13, p-value = 0.02643
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.08198632 0.83780278
sample estimates:
      cor 
0.5702975 


```


CHS.TI
```
CHS.TI.elev <- subset(mtDNA.data, Region=="CHS.TI",)

cor.test(CHS.TI.elev$elev.c, CHS.TI.elev$Hd, type="pearson")

	Pearson's product-moment correlation

data:  CHS.TI.elev$elev.c and CHS.TI.elev$Hd
t = -1.4386, df = 13, p-value = 0.1739
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.7419863  0.1748850
sample estimates:
       cor 
-0.3705767 

cor.test(CHS.TI.elev$elev.c, CHS.TI.elev$nd, type="pearson")

	Pearson's product-moment correlation

data:  CHS.TI.elev$elev.c and CHS.TI.elev$nd
t = 0.1089, df = 13, p-value = 0.9149
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.4896450  0.5341917
sample estimates:
       cor 
0.03019038 

```




CHS.VS
```
CHS.VS.elev <- subset(mtDNA.data, Region=="CHS.VS",)

cor.test(CHS.VS.elev$elev.c, CHS.VS.elev$Hd, type="pearson")

	Pearson's product-moment correlation

data:  CHS.VS.elev$elev.c and CHS.VS.elev$Hd
t = -1.0616, df = 5, p-value = 0.337
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.8933974  0.4788215
sample estimates:
       cor 
-0.4288837 

cor.test(CHS.VS.elev$elev.c, CHS.VS.elev$nd, type="pearson")

	Pearson's product-moment correlation

data:  CHS.VS.elev$elev.c and CHS.VS.elev$nd
t = -1.0576, df = 5, p-value = 0.3386
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.8930735  0.4800555
sample estimates:
       cor 
-0.4275752 

```





CZ
```
CZ.elev <- subset(mtDNA.data, Region=="CZ",)

cor.test(CZ.elev$elev.c, CZ.elev$Hd, type="pearson")

	Pearson's product-moment correlation

data:  CZ.elev$elev.c and CZ.elev$Hd
t = -1.3984, df = 28, p-value = 0.173
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.5638616  0.1153972
sample estimates:
       cor 
-0.2554939 


cor.test(CZ.elev$elev.c, CZ.elev$nd, type="pearson")

	Pearson's product-moment correlation

data:  CZ.elev$elev.c and CZ.elev$nd
t = -0.51593, df = 28, p-value = 0.61
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.4418630  0.2727635
sample estimates:
        cor 
-0.09704187 
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


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GeneticVar.vs.Elev


HD input files have been prepared before: 

R: mtDNA.data$Hd

I have to add elevation and geographic distance. Since I want to understand how genetic diversity is distributed up elevation, I will measure all distances from the lowest elevation population

CHN: berm
CHS.VS: grma
CHS.TI: mart
CZ: kebe


CHN
```
library(fields)
library(reshape2)

CHNData <- subset(mtDNA.data, Region=="CHN",)
CHN_long.lat <- cbind(CHNData$Long, CHNData$Lat)

CHN.dist.matrix <- rdist.earth(CHN_long.lat, miles=F)
CHN.m.dist <- as.matrix(CHN.dist.matrix)

colnames(CHN.m.dist) <- (CHNData$Pop)
rownames(CHN.m.dist) <- (CHNData$Pop)

head(CHN.m.dist)
CHN.m2.dist <- melt(CHN.m.dist)[melt(upper.tri(CHN.m.dist))$value,]

names(CHN.m2.dist) <- c("c1", "c2", "distance")
CHN.m2.dist$log.km <- log(CHN.m2.dist$distance)

CHN.dist.berm <- subset(CHN.m2.dist, c1==c("berm"), select=c1:log.km)

new.row <- c("berm", "berm", "0", "0")
CHN.dist.berm.new <- rbind(new.row,CHN.dist.berm)
CHNData$Pop
CHN.dist.berm.new$c2
CHNData$log.km <- CHN.dist.berm.new$log.km
CHNData$log.km <- as.numeric(CHNData$log.km)
```

CHN Full RDA
```
RDA.CHN.dist.elev <- rda(CHNData$Hd ~ Elev + log.km, CHNData)

RDA.CHN.dist.elev

Call: rda(formula = CHNData$Hd ~ Elev + log.km, data = CHNData)

              Inertia Proportion Rank
Total         0.08098    1.00000     
Constrained   0.02828    0.34920    1
Unconstrained 0.05270    0.65080    1
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1 
0.02828 

Eigenvalues for unconstrained axes:
   PC1 
0.0527 


RsquareAdj(RDA.CHN.dist.elev)

$r.squared
[1] 0.3491992

$adj.r.squared
[1] 0.2407324

anova(RDA.CHN.dist.elev)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHNData$Hd ~ Elev + log.km, data = CHNData)
         Df Variance      F Pr(>F)  
Model     2 0.028280 3.2194  0.064 .
Residual 12 0.052705                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

pCHN.distOnly
```
pCHN.distOnly <- rda(CHNData$Hd ~log.km + Condition(Elev), CHNData)

pCHN.distOnly

Call: rda(formula = CHNData$Hd ~ log.km + Condition(Elev), data =
CHNData)

               Inertia Proportion Rank
Total         0.080984   1.000000     
Conditional   0.027136   0.335078    1
Constrained   0.001144   0.014121    1
Unconstrained 0.052705   0.650801    1
Inertia is variance 

Eigenvalues for constrained axes:
     RDA1 
0.0011436 

Eigenvalues for unconstrained axes:
   PC1 
0.0527 



RsquareAdj(pCHN.distOnly)

$r.squared
[1] 0.01412113

$adj.r.squared
[1] -0.04319783

anova(pCHN.distOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHNData$Hd ~ log.km + Condition(Elev), data = CHNData)
         Df Variance      F Pr(>F)
Model     1 0.001144 0.2604  0.646
Residual 12 0.052705        
```

pCHN.elevOnly
```
pCHN.elevOnly <- rda(CHNData$Hd ~Elev + Condition(log.km), CHNData)

pCHN.elevOnly

Call: rda(formula = CHNData$Hd ~ Elev + Condition(log.km), data =
CHNData)

               Inertia Proportion Rank
Total         0.080984   1.000000     
Conditional   0.009229   0.113955    1
Constrained   0.019051   0.235245    1
Unconstrained 0.052705   0.650801    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.019051 

Eigenvalues for unconstrained axes:
   PC1 
0.0527 



RsquareAdj(pCHN.elevOnly)

$r.squared
[1] 0.2352446

$adj.r.squared
[1] 0.1949352

anova(pCHN.elevOnly)

Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHNData$Hd ~ Elev + Condition(log.km), data = CHNData)
         Df Variance      F Pr(>F)  
Model     1 0.019051 4.3376  0.053 .
Residual 12 0.052705                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


CHS.TI
```
CHS.TIData <- subset(mtDNA.data, Region=="CHS.TI",)
CHS.TI_long.lat <- cbind(CHS.TIData$Long, CHS.TIData$Lat)

CHS.TI.dist.matrix <- rdist.earth(CHS.TI_long.lat, miles=F)
CHS.TI.m.dist <- as.matrix(CHS.TI.dist.matrix)

colnames(CHS.TI.m.dist) <- (CHS.TIData$Pop)
rownames(CHS.TI.m.dist) <- (CHS.TIData$Pop)

head(CHS.TI.m.dist)
CHS.TI.m2.dist <- melt(CHS.TI.m.dist)[melt(upper.tri(CHS.TI.m.dist))$value,]

names(CHS.TI.m2.dist) <- c("c1", "c2", "distance")
CHS.TI.m2.dist$log.km <- log(CHS.TI.m2.dist$distance)

CHS.TI.dist.mart <- subset(CHS.TI.m2.dist, c1==c("mart"), select=c1:log.km)

new.row <- c("mart", "mart", "0", "0")
CHS.TI.dist.mart.new <- rbind(new.row,CHS.TI.dist.mart)
CHS.TIData$Pop
CHS.TI.dist.mart.new$c2
CHS.TIData$log.km <- CHS.TI.dist.mart.new$log.km
CHS.TIData$log.km <- as.numeric(CHS.TIData$log.km)

```

CHS.TI Full RDA
```
RDA.CHS.TI.dist.elev <- rda(CHS.TIData$Hd ~ Elev + log.km, CHS.TIData)
RDA.CHS.TI.dist.elev
Call: rda(formula = CHS.TIData$Hd ~ Elev + log.km, data = CHS.TIData)

              Inertia Proportion Rank
Total         0.10273    1.00000     
Constrained   0.01423    0.13851    1
Unconstrained 0.08850    0.86149    1
Inertia is variance 

Eigenvalues for constrained axes:
   RDA1 
0.01423 

Eigenvalues for unconstrained axes:
   PC1 
0.0885 

RsquareAdj(RDA.CHS.TI.dist.elev)
$r.squared
[1] 0.1385137

$adj.r.squared
[1] -0.005067377

anova(RDA.CHS.TI.dist.elev)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.TIData$Hd ~ Elev + log.km, data = CHS.TIData)
         Df Variance      F Pr(>F)
Model     2 0.014230 0.9647  0.396
Residual 12 0.088501    
```

pCHS.TI.distOnly
```
pCHS.TI.distOnly <- rda(CHS.TIData$Hd ~ log.km + Condition(Elev), CHS.TIData) 
pCHS.TI.distOnly

Call: rda(formula = CHS.TIData$Hd ~ log.km + Condition(Elev), data =
CHS.TIData)

                Inertia Proportion Rank
Total         0.1027304  1.0000000     
Conditional   0.0141077  0.1373271    1
Constrained   0.0001219  0.0011866    1
Unconstrained 0.0885008  0.8614863    1
Inertia is variance 

Eigenvalues for constrained axes:
     RDA1 
0.0001219 

Eigenvalues for unconstrained axes:
   PC1 
0.0885 

RsquareAdj(pCHS.TI.distOnly)
$r.squared
[1] 0.001186575

$adj.r.squared
[1] -0.07603503

anova(pCHS.TI.distOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.TIData$Hd ~ log.km + Condition(Elev), data = CHS.TIData)
         Df Variance      F Pr(>F)
Model     1 0.000122 0.0165  0.888
Residual 12 0.088501  
```

pCHS.TI.elevOnly
```
pCHS.TI.elevOnly <- rda(CHS.TIData$Hd ~Elev + Condition(log.km), CHS.TIData)
pCHS.TI.elevOnly
Call: rda(formula = CHS.TIData$Hd ~ Elev + Condition(log.km), data =
CHS.TIData)

               Inertia Proportion Rank
Total         0.102730   1.000000     
Conditional   0.006907   0.067233    1
Constrained   0.007323   0.071281    1
Unconstrained 0.088501   0.861486    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.007323 

Eigenvalues for unconstrained axes:
   PC1 
0.0885 

RsquareAdj(pCHS.TI.elevOnly)
$r.squared
[1] 0.07128101

$adj.r.squared
[1] -0.0005487082

anova(pCHS.TI.elevOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.TIData$Hd ~ Elev + Condition(log.km), data = CHS.TIData)
         Df Variance      F Pr(>F)
Model     1 0.007323 0.9929  0.344
Residual 12 0.088501    
```


CHS.VS
```
CHS.VSData <- subset(mtDNA.data, Region=="CHS.VS",)
CHS.VS_long.lat <- cbind(CHS.VSData$Long, CHS.VSData$Lat)

CHS.VS.dist.matrix <- rdist.earth(CHS.VS_long.lat, miles=F)
CHS.VS.m.dist <- as.matrix(CHS.VS.dist.matrix)

colnames(CHS.VS.m.dist) <- (CHS.VSData$Pop)
rownames(CHS.VS.m.dist) <- (CHS.VSData$Pop)

head(CHS.VS.m.dist)
CHS.VS.m2.dist <- melt(CHS.VS.m.dist)[melt(upper.tri(CHS.VS.m.dist))$value,]

names(CHS.VS.m2.dist) <- c("c1", "c2", "distance")
CHS.VS.m2.dist$log.km <- log(CHS.VS.m2.dist$distance)

CHS.VS.dist.grma <- subset(CHS.VS.m2.dist, c1==c("grma"), select=c1:log.km)

new.row <- c("grma", "grma", "0", "0")
CHS.VS.dist.grma.new <- rbind(new.row,CHS.VS.dist.grma)
CHS.VSData$Pop
CHS.VS.dist.grma.new$c2
CHS.VSData$log.km <- CHS.VS.dist.grma.new$log.km
CHS.VSData$log.km <- as.numeric(CHS.VSData$log.km)

```

CHS.VS Full RDA
```
RDA.CHS.VS.dist.elev <- rda(CHS.VSData$Hd ~ Elev + log.km, CHS.VSData)
RDA.CHS.VS.dist.elev
Call: rda(formula = CHS.VSData$Hd ~ Elev + log.km, data = CHS.VSData)

              Inertia Proportion Rank
Total         0.11694    1.00000     
Constrained   0.02258    0.19305    1
Unconstrained 0.09437    0.80695    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.022575 

Eigenvalues for unconstrained axes:
    PC1 
0.09437 

RsquareAdj(RDA.CHS.VS.dist.elev)
$r.squared
[1] 0.1930486

$adj.r.squared
[1] -0.210427

anova(RDA.CHS.VS.dist.elev)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.VSData$Hd ~ Elev + log.km, data = CHS.VSData)
         Df Variance      F Pr(>F)
Model     2 0.022575 0.4785  0.665
Residual  4 0.094366    
```

pCHS.VS.distOnly
```
pCHS.VS.distOnly <- rda(CHS.VSData$Hd ~ log.km + Condition(Elev), CHS.VSData) 
pCHS.VS.distOnly

Call: rda(formula = CHS.VSData$Hd ~ log.km + Condition(Elev), data =
CHS.VSData)

               Inertia Proportion Rank
Total         0.116941   1.000000     
Conditional   0.021510   0.183941    1
Constrained   0.001065   0.009107    1
Unconstrained 0.094366   0.806951    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.001065 

Eigenvalues for unconstrained axes:
    PC1 
0.09437 

RsquareAdj(pCHS.VS.distOnly)
$r.squared
[1] 0.009107431

$adj.r.squared
[1] -0.2311565

anova(pCHS.VS.distOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.VSData$Hd ~ log.km + Condition(Elev), data = CHS.VSData)
         Df Variance      F Pr(>F)
Model     1 0.001065 0.0451  0.921
Residual  4 0.094366   
```

pCHS.VS.elevOnly
```
pCHS.VS.elevOnly <- rda(CHS.VSData$Hd ~Elev + Condition(log.km), CHS.VSData)
pCHS.VS.elevOnly
Call: rda(formula = CHS.VSData$Hd ~ Elev + Condition(log.km), data =
CHS.VSData)

               Inertia Proportion Rank
Total         0.116941   1.000000     
Conditional   0.009428   0.080619    1
Constrained   0.013148   0.112429    1
Unconstrained 0.094366   0.806951    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.013148 

Eigenvalues for unconstrained axes:
    PC1 
0.09437 

RsquareAdj(pCHS.VS.elevOnly)
$r.squared
[1] 0.1124293

$adj.r.squared
[1] -0.1071702

anova(pCHS.VS.elevOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CHS.VSData$Hd ~ Elev + Condition(log.km), data = CHS.VSData)
         Df Variance      F Pr(>F)
Model     1 0.013148 0.5573  0.509
Residual  4 0.094366     
```


CZ
```
CZData <- subset(mtDNA.data, Region=="CZ",)
CZ_long.lat <- cbind(CZData$Long, CZData$Lat)

CZ.dist.matrix <- rdist.earth(CZ_long.lat, miles=F)
CZ.m.dist <- as.matrix(CZ.dist.matrix)

colnames(CZ.m.dist) <- (CZData$Pop)
rownames(CZ.m.dist) <- (CZData$Pop)

head(CZ.m.dist)
CZ.m2.dist <- melt(CZ.m.dist)[melt(upper.tri(CZ.m.dist))$value,]

names(CZ.m2.dist) <- c("c1", "c2", "distance")
CZ.m2.dist$log.km <- log(CZ.m2.dist$distance)

CZ.dist.kebe <- subset(CZ.m2.dist, c1==c("kebe"), select=c1:log.km)

new.row <- c("kebe", "kebe", "0", "0")
CZ.dist.kebe.new <- rbind(new.row,CZ.dist.kebe)
CZData$Pop
CZ.dist.kebe.new$c2
CZData$log.km <- CZ.dist.kebe.new$log.km
CZData$log.km <- as.numeric(CZData$log.km)

```

CZ Full RDA
```
RDA.CZ.dist.elev <- rda(CZData$Hd ~ Elev + log.km, CZData)
RDA.CZ.dist.elev
Call: rda(formula = CZData$Hd ~ Elev + log.km, data = CZData)

               Inertia Proportion Rank
Total         0.109845   1.000000     
Constrained   0.009354   0.085156    1
Unconstrained 0.100491   0.914844    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.009354 

Eigenvalues for unconstrained axes:
    PC1 
0.10049 

RsquareAdj(RDA.CZ.dist.elev)
$r.squared
[1] 0.08515555

$adj.r.squared
[1] 0.0173893

anova(RDA.CZ.dist.elev)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CZData$Hd ~ Elev + log.km, data = CZData)
         Df Variance      F Pr(>F)
Model     2 0.009354 1.2566   0.29
Residual 27 0.100491  
```

pCZ.distOnly
```
pCZ.distOnly <- rda(CZData$Hd ~ log.km + Condition(Elev), CZData) 
pCZ.distOnly

Call: rda(formula = CZData$Hd ~ log.km + Condition(Elev), data =
CZData)

               Inertia Proportion Rank
Total         0.109845   1.000000     
Conditional   0.007170   0.065277    1
Constrained   0.002184   0.019878    1
Unconstrained 0.100491   0.914844    1
Inertia is variance 

Eigenvalues for constrained axes:
     RDA1 
0.0021835 

Eigenvalues for unconstrained axes:
    PC1 
0.10049 

RsquareAdj(pCZ.distOnly)
$r.squared
[1] 0.01987841

$adj.r.squared
[1] -0.01450489

anova(pCZ.distOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CZData$Hd ~ log.km + Condition(Elev), data = CZData)
         Df Variance      F Pr(>F)
Model     1 0.002184 0.5867  0.495
Residual 27 0.100491   
```

pCZ.elevOnly
```
pCZ.elevOnly <- rda(CZData$Hd ~Elev + Condition(log.km), CZData)
pCZ.elevOnly
Call: rda(formula = CZData$Hd ~ Elev + Condition(log.km), data =
CZData)

               Inertia Proportion Rank
Total         0.116941   1.000000     
Conditional   0.009428   0.080619    1
Constrained   0.013148   0.112429    1
Unconstrained 0.094366   0.806951    1
Inertia is variance 

Eigenvalues for constrained axes:
    RDA1 
0.013148 

Eigenvalues for unconstrained axes:
    PC1 
0.09437 

RsquareAdj(pCZ.elevOnly)
$r.squared
[1] 0.1124293

$adj.r.squared
[1] -0.1071702

anova(pCZ.elevOnly)
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = CZData$Hd ~ Elev + Condition(log.km), data = CZData)
         Df Variance      F Pr(>F)
Model     1 0.013148 0.5573  0.509
Residual  4 0.094366     
```




## 3. Isolation by Dist vs Isolation by Elev

Within a transect. Is Fst higher across elevation than across distance? 

IBD has previously been calculated using the RAD data

I'll have to add in all the mtDNA data later... 

##### Mantel assesment

I originally just ran some linear models to test for IBD and isolation by elevation. (I've also partialled out the effects of geog and elev on Genetic distance using RDA). 
However, Josh was not happy about this. So I'm running some Mantel and partial mantel tests to see whether there's IBD and IBelevation while controlling for the other. 

Mantel test to test for IBD and IBElevation in each region

Two questions: 

1. To test for scale in IBD patterns, we use pRDA to determine if genetic distance increases with geog distance partialling out elev

2. To test for gene flow across elevation, pRDA partiallying out geographic distance. 


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats


CHN 
```
library(vegan)

CHN.fst.matrix <- as.matrix(CHN.fst) ##create an fst matrix
CHN.elev.dist.matrix <- as.matrix(CHN.elev.dist)  ##elevational differences as matrix

CHN.mantel.geog <- mantel(CHN.fst.matrix, CHN.m.dist, method="pearson", permutations=999) ##mantel fst x geog

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHN.fst.matrix, ydis = CHN.m.dist, method = "pearson",      permutations = 999) 

Mantel statistic r: 0.3584 
      Significance: 0.002 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.169 0.217 0.268 0.292 
Permutation: free
Number of permutations: 999


CHN.mantel.elev <- mantel(CHN.fst.matrix, CHN.elev.dist.matrix, method="pearson", permutations=999)  ##mantel fst x elev


Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHN.fst.matrix, ydis = CHN.elev.dist.matrix, method = "pearson",      permutations = 999) 

Mantel statistic r: 0.3215 
      Significance: 0.001 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.132 0.171 0.209 0.236 
Permutation: free
Number of permutations: 999

CHN.mantel.geog.pElev <- mantel.partial(CHN.fst.matrix, CHN.m.dist, CHN.elev.dist.matrix, method="pearson", permutations=999)  ##pMantel fst x geog, elev constant

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHN.fst.matrix, ydis = CHN.m.dist, zdis = CHN.elev.dist.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r: 0.3857 
      Significance: 0.002 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.154 0.205 0.252 0.284 
Permutation: free
Number of permutations: 999

CHN.mantel.elev.pGeog <- mantel.partial(CHN.fst.matrix, CHN.elev.dist.matrix, CHN.m.dist, method="pearson", permutations=999) ##pMantel fst x elev, geog constant

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHN.fst.matrix, ydis = CHN.elev.dist.matrix,      zdis = CHN.m.dist, method = "pearson", permutations = 999) 

Mantel statistic r: 0.3525 
      Significance: 0.001 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.145 0.186 0.226 0.263 
Permutation: free
Number of permutations: 999

```


CHS.TI 
```
#library(vegan)
library(fields)

First make sure that oalp is removed from all datasets (i.e. I need 15x15 datasets.)

CHS.TI.noOalp_lon.lat <- CHS.TI_lon.lat[-10,]
CHS.TI.dist.noOalp.matrix <- rdist.earth(CHS.TI.noOalp_lon.lat, miles=F)
CHS.TI.dist.noOalp.matrix <- as.matrix(CHS.TI.dist.noOalp.matrix) ##geographic distance


CHS.TI.elev.fst.matrix <- as.matrix(CHS.TI.elev.fst) ##create an fst matrix
CHS.TI.elev.dist.matrix <- as.matrix(CHS.TI.elev.dist)  ##elevational differences as matrix

CHS.TI.mantel.geog <- mantel(CHS.TI.elev.fst.matrix, CHS.TI.dist.noOalp.matrix, method="pearson", permutations=999) ##mantel fst x geog

CHS.TI.mantel.geog 

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHS.TI.elev.fst.matrix, ydis = CHS.TI.dist.noOalp.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r: 0.2319 
      Significance: 0.026 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.147 0.183 0.222 0.259 
Permutation: free
Number of permutations: 999

CHS.TI.mantel.elev <- mantel(CHS.TI.elev.fst.matrix, CHS.TI.elev.dist.matrix, method="pearson", permutations=999)  ##mantel fst x elev

CHS.TI.mantel.elev
Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHS.TI.elev.fst.matrix, ydis = CHS.TI.elev.dist.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r: 0.1429 
      Significance: 0.039 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0932 0.1274 0.1736 0.2022 
Permutation: free
Number of permutations: 999

CHS.TI.mantel.geog.pElev <- mantel.partial(CHS.TI.elev.fst.matrix, CHS.TI.dist.noOalp.matrix, CHS.TI.elev.dist.matrix, method="pearson", permutations=999)  ##pMantel fst x geog, elev constant

CHS.TI.mantel.geog.pElev

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHS.TI.elev.fst.matrix, ydis = CHS.TI.dist.noOalp.matrix,      zdis = CHS.TI.elev.dist.matrix, method = "pearson", permutations = 999) 

Mantel statistic r: 0.1927 
      Significance: 0.058 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.154 0.200 0.232 0.269 
Permutation: free
Number of permutations: 999

CHS.TI.mantel.elev.pGeog <- mantel.partial(CHS.TI.elev.fst.matrix, CHS.TI.elev.dist.matrix, CHS.TI.dist.noOalp.matrix, method="pearson", permutations=999) ##pMantel fst x elev, geog constant

CHS.TI.mantel.elev.pGeog

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHS.TI.elev.fst.matrix, ydis = CHS.TI.elev.dist.matrix,      zdis = CHS.TI.dist.noOalp.matrix, method = "pearson", permutations = 999) 

Mantel statistic r: 0.05626 
      Significance: 0.245 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.123 0.161 0.198 0.223 
Permutation: free
Number of permutations: 999
```



Mantel test to test for IBD and IBElevation in each region

Two questions: 

1. To test for scale in IBD patterns, we use pRDA to determine if genetic distance increases with geog distance partialling out elev

2. To test for gene flow across elevation, pRDA partiallying out geographic distance. 





I need to do this within all regions, and for CHall. 


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

CHS.VS 
```
#library(vegan)

CHS.VS.fst.matrix <- as.matrix(CHS.VS.fst) ##create an fst matrix 
CHS.VS.elev.dist.matrix <- as.matrix(CHS.VS.elev.dist)  ##elevational differences as matrix

CHS.VS.mantel.geog <- mantel(CHS.VS.fst.matrix, CHS.VS.m.dist, method="pearson", permutations=999) ##mantel fst x geog

CHS.VS.mantel.geog 

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHS.VS.fst.matrix, ydis = CHS.VS.m.dist, method = "pearson",      permutations = 999) 

Mantel statistic r: 0.04773 
      Significance: 0.378 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.407 0.478 0.534 0.591 
Permutation: free
Number of permutations: 999

CHS.VS.mantel.elev <- mantel(CHS.VS.fst.matrix, CHS.VS.elev.dist.matrix, method="pearson", permutations=999)  ##mantel fst x elev

CHS.VS.mantel.elev

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHS.VS.fst.matrix, ydis = CHS.VS.elev.dist.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r: 0.02675 
      Significance: 0.408 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.269 0.352 0.402 0.475 
Permutation: free
Number of permutations: 999

CHS.VS.mantel.geog.pElev <- mantel.partial(CHS.VS.fst.matrix, CHS.VS.m.dist, CHS.VS.elev.dist.matrix, method="pearson", permutations=999)  ##pMantel fst x geog, elev constant

CHS.VS.mantel.geog.pElev

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHS.VS.fst.matrix, ydis = CHS.VS.m.dist,      zdis = CHS.VS.elev.dist.matrix, method = "pearson", permutations = 999) 

Mantel statistic r: 0.04776 
      Significance: 0.401 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.410 0.496 0.556 0.634 
Permutation: free
Number of permutations: 999

CHS.VS.mantel.elev.pGeog <- mantel.partial(CHS.VS.fst.matrix, CHS.VS.elev.dist.matrix, CHS.VS.m.dist, method="pearson", permutations=999) ##pMantel fst x elev, geog constant

CHS.VS.mantel.elev.pGeog

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHS.VS.fst.matrix, ydis = CHS.VS.elev.dist.matrix,      zdis = CHS.VS.m.dist, method = "pearson", permutations = 999) 

Mantel statistic r: 0.02681 
      Significance: 0.376 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.274 0.351 0.396 0.445 
Permutation: free
Number of permutations: 999


```

Mantel test to test for IBD and IBElevation in each region

Two questions: 

1. To test for scale in IBD patterns, we use pRDA to determine if genetic distance increases with geog distance partialling out elev

2. To test for gene flow across elevation, pRDA partiallying out geographic distance. 





I need to do this within all regions, and for CHall. 


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

CZ 
```
#library(vegan)

CZ.fst.matrix <- as.matrix(CZ.fst) ##create an fst matrix 
CZ.elev.dist.matrix <- as.matrix(CZ.elev.dist)  ##elevational differences as matrix

CZ.mantel.geog <- mantel(CZ.fst.matrix, CZ.m.dist, method="pearson", permutations=999) ##mantel fst x geog

CZ.mantel.geog 

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CZ.fst.matrix, ydis = CZ.m.dist, method = "pearson",      permutations = 999) 

Mantel statistic r: 0.6381 
      Significance: 0.001 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.141 0.179 0.213 0.267 
Permutation: free
Number of permutations: 999

CZ.mantel.elev <- mantel(CZ.fst.matrix, CZ.elev.dist.matrix, method="pearson", permutations=999)  ##mantel fst x elev

CZ.mantel.elev

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CZ.fst.matrix, ydis = CZ.elev.dist.matrix, method = "pearson",      permutations = 999) 

Mantel statistic r: 0.1519 
      Significance: 0.002 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0503 0.0651 0.0767 0.1104 
Permutation: free
Number of permutations: 999

CZ.mantel.geog.pElev <- mantel.partial(CZ.fst.matrix, CZ.m.dist, CZ.elev.dist.matrix, method="pearson", permutations=999)  ##pMantel fst x geog, elev constant

CZ.mantel.geog.pElev

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CZ.fst.matrix, ydis = CZ.m.dist, zdis = CZ.elev.dist.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r: 0.6304 
      Significance: 0.001 

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.141 0.184 0.210 0.233 
Permutation: free
Number of permutations: 999

CZ.mantel.elev.pGeog <- mantel.partial(CZ.fst.matrix, CZ.elev.dist.matrix, CZ.m.dist, method="pearson", permutations=999) ##pMantel fst x elev, geog constant

CZ.mantel.elev.pGeog

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CZ.fst.matrix, ydis = CZ.elev.dist.matrix,      zdis = CZ.m.dist, method = "pearson", permutations = 999) 

Mantel statistic r: -0.08314 
      Significance: 0.978 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0559 0.0723 0.0876 0.1090 
Permutation: free
Number of permutations: 999



```


Mantel test to test for IBD and IBElevation in each region

Two questions: 

1. To test for scale in IBD patterns, we use pRDA to determine if genetic distance increases with geog distance partialling out elev

2. To test for gene flow across elevation, pRDA partiallying out geographic distance. 





I need to do this within all regions, and for CHall. 


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

CHall 
```
#library(vegan)

CHall.fst.matrix <- as.matrix(CHall.fst)
dim(CHall.fst.matrix)
[1] 82 82

CHall.elev.diff.matrix <- as.matrix(CHall.elev.diff)
dim(CHall.elev.diff.matrix)
[1] 82 82

CHall.geog.dist.matrix <- as.matrix(CHall.m.dist)
dim(CHall.geog.dist.matrix)


CHall.mantel.geog <- mantel(CHall.fst.matrix, CHall.geog.dist.matrix, method="pearson", permutations=999) ##mantel fst x geog

CHall.mantel.geog 

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHall.fst.matrix, ydis = CHall.geog.dist.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r:   0.3 
      Significance: 0.001 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0424 0.0547 0.0663 0.0841 
Permutation: free
Number of permutations: 999

CHall.mantel.elev <- mantel(CHall.fst.matrix, CHall.elev.diff.matrix, method="pearson", permutations=999)  ##mantel fst x elev

CHall.mantel.elev

Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel(xdis = CHall.fst.matrix, ydis = CHall.elev.diff.matrix,      method = "pearson", permutations = 999) 

Mantel statistic r:  0.05 
      Significance: 0.042 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0354 0.0474 0.0578 0.0719 
Permutation: free
Number of permutations: 999

CHall.mantel.geog.pElev <- mantel.partial(CHall.fst.matrix, CHall.geog.dist.matrix, CHall.elev.diff.matrix, method="pearson", permutations=999)  ##pMantel fst x geog, elev constant

CHall.mantel.geog.pElev

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHall.fst.matrix, ydis = CHall.geog.dist.matrix,      zdis = CHall.elev.diff.matrix, method = "pearson", permutations = 999) 

Mantel statistic r: 0.2987 
      Significance: 0.001 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0402 0.0527 0.0623 0.0778 
Permutation: free
Number of permutations: 999

CHall.mantel.elev.pGeog <- mantel.partial(CHall.fst.matrix, CHall.elev.diff.matrix, CHall.geog.dist.matrix, method="pearson", permutations=999) ##pMantel fst x elev, geog constant

CHall.mantel.elev.pGeog

Partial Mantel statistic based on Pearson's product-moment correlation 

Call:
mantel.partial(xdis = CHall.fst.matrix, ydis = CHall.elev.diff.matrix,      zdis = CHall.geog.dist.matrix, method = "pearson", permutations = 999) 

Mantel statistic r: 0.04018 
      Significance: 0.066 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0335 0.0440 0.0548 0.0657 
Permutation: free
Number of permutations: 999



```



### 3.2 

## Isolation by Elevation

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

CHN
```
library(reshape2)
CHall.elev <- read.table("CHall.82pops.coords.elev", header=T)


CHN.elev <- subset(CHall.elev, region=="CHN",)
CHN.elevonly <- CHN.elev$elev
CHN.elev.dist <- dist(CHN.elevonly, method="euclidean")

m.CHN.elev <- as.matrix(CHN.elev.dist)
m2.CHN.elev <- melt(m.CHN.elev)[melt(upper.tri(m.CHN.elev))$value,]
colnames(m2.CHN.elev) <- c("site1", "site2", "elev")
m2.CHN.elev$log.elev <- log(m2.CHN.elev$value)

pdf("CHN.Isol.Elev")
plot(m2.CHN$IBD~m2.CHN.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)", ylim=c(0,0.6))
abline(fit <- lm(m2.CHN$IBD~m2.CHN.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation plot - CHN")
dev.off()

##test significance
fit <- lm(m2.CHN$IBD~m2.CHN.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHN$IBD
                      Df   Sum Sq   Mean Sq F value    Pr(>F)    
m2.CHN.elev$log.elev   1 0.005876 0.0058756   13.21 0.0003692 ***
Residuals            169 0.075168 0.0004448                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


## Isolation by Elevation

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

CHS.VS
```
CHS.VS.elev <- subset(CHall.elev, region=="CHS.VS",)
CHS.VS.elev <- CHS.VS.elev[-9,]  ##remove oalp from analysis

CHS.VS.elevonly <- CHS.VS.elev$elev
CHS.VS.elev.dist <- dist(CHS.VS.elevonly, method="euclidean")

m.CHS.VS.elev <- as.matrix(CHS.VS.elev.dist)
m2.CHS.VS.elev <- melt(m.CHS.VS.elev)[melt(upper.tri(m.CHS.VS.elev))$value,]
colnames(m2.CHS.VS.elev) <- c("site1", "site2", "elev")
m2.CHS.VS.elev$log.elev <- log(m2.CHS.VS.elev$elev)

pdf("CHS.VS.Isol.Elev")
plot(m2.CHS.VS$IBD~m2.CHS.VS.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)", ylim=c(0,0.6))
abline(fit <- lm(m2.CHS.VS$IBD~m2.CHS.VS.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation plot - CHS.VS")
dev.off()


##test significance
fit <- lm(m2.CHS.VS$IBD~m2.CHS.VS.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.VS$IBD
                        Df   Sum Sq    Mean Sq F value Pr(>F)
m2.CHS.VS.elev$log.elev  1 0.000051 0.00005118  0.0679 0.7956
Residuals               43 0.032399 0.00075347 

```


CZ
```
CZ.elev <- subset(CHall.elev, region=="CZ",)

#one of the values (wise vs egel) = -Inf because these sites are at exactly the same elevation. 
I’m changing the elevation of egel to 446

CZ.elev.2$elev[CZ.elev.2$site=="egel" & CZ.elev.2$elev==445] <- 446


CZ.elevonly <- CZ.elev.2$elev
CZ.elev.dist <- dist(CZ.elevonly, method="euclidean")

m.CZ.elev <- as.matrix(CZ.elev.dist)
m2.CZ.elev <- melt(m.CZ.elev)[melt(upper.tri(m.CZ.elev))$value,]
colnames(m2.CZ.elev) <- c("site1", "site2", "elev")
m2.CZ.elev$log.elev <- log(m2.CZ.elev$elev)


pdf("CZ.Isol.Elev")
plot(m2.CZ$IBD~m2.CZ.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)", ylim=c(0,0.6))
abline(fit <- lm(m2.CZ$IBD~m2.CZ.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation plot - CZ")
dev.off()


And test the significance
fit <- lm(m2.CZ$IBD~m2.CZ.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CZ$IBD
                     Df  Sum Sq   Mean Sq F value   Pr(>F)   
m2.CZ.elev$log.elev   1 0.03005 0.0300522  8.8367 0.003059 **
Residuals           664 2.25816 0.0034008                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```


CHS.TI
```
##I had to remove oalp from the genind object and recalculate the Fst
library(hierfstat)

obj <- seppop(CHall.L3$CHS.TI)

CHS.TI.genind.oalpremoved <- repool(obj$agra, obj$arce, obj$cava, obj$forn, obj$gola, obj$gott, obj$lucm, obj$magn, obj$mart, obj$pozz, obj$sali, obj$scai, obj$star, obj$stba, obj$zeni)

CHS.TI.elev.fst <- pairwise.fst(CHS.TI.genind.oalpremoved, pop=NULL, res.type=c("dist", "matrix"))


m.CHS.TI.fst.elev <- as.matrix(CHS.TI.elev.fst)
m2.CHS.TI.fst.elev <- melt(m.CHS.TI.fst.elev)[melt(upper.tri(m.CHS.TI.fst.elev))$value,]
names(m2.CHS.TI.fst.elev) <- c("c1", "c2", "distance")
m2.CHS.TI.fst.elev$IBD <- (m2.CHS.TI.fst.elev$distance/(1-m2.CHS.TI.fst.elev$distance))



CHS.TI.elev <- subset(CHall.elev, region=="CHS.TI",)

CHS.TI.elevonly <- CHS.TI.elev$elev
CHS.TI.elev.dist <- dist(CHS.TI.elevonly, method="euclidean")

m.CHS.TI.elev <- as.matrix(CHS.TI.elev.dist)
m2.CHS.TI.elev <- melt(m.CHS.TI.elev)[melt(upper.tri(m.CHS.TI.elev))$value,]
colnames(m2.CHS.TI.elev) <- c("site1", "site2", "elev")
m2.CHS.TI.elev$log.elev <- log(m2.CHS.TI.elev$elev)

pdf("CHS.TI.Isol.Elev")
plot(m2.CHS.TI.fst.elev$IBD~m2.CHS.TI.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)", ylim=c(0,0.6))
abline(fit <- lm(m2.CHS.TI.fst.elev$IBD~m2.CHS.TI.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation plot - CHS.TI")
dev.off()

##test significance

fit <- lm(m2.CHS.TI.fst.elev$IBD~m2.CHS.TI.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.TI.fst.elev$IBD
                         Df   Sum Sq    Mean Sq F value Pr(>F)
m2.CHS.TI.elev$log.elev   1 0.000008 0.00000844  0.0075  0.931
Residuals               103 0.115268 0.00111910  

```


### Check significance of IBD calculated previously

https://github.com/alexjvr1/Manuscripts/blob/19f7dd1845d4fa852e6d729f0ee456e274b0182c/5.CHP2_CH_LandscapeGenomics.md

```
##significance of IBD

#CHall
fit <- lm(m2.CHall$IBD~CHall.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHall$IBD
                       Df Sum Sq Mean Sq F value    Pr(>F)    
CHall.m2.dist$log.km    1  4.079  4.0789  404.96 < 2.2e-16 ***
Residuals            3319 33.430  0.0101                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#CHN
fit <- lm(m2.CHN$IBD~CHN.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHN$IBD
                    Df   Sum Sq   Mean Sq F value   Pr(>F)    
CHN.m2.dist$log.km   1 0.010874 0.0108743   26.19 8.34e-07 ***
Residuals          169 0.070170 0.0004152                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#CHS
fit <- lm(m2.CHS$IBD~CHS.m2.dist$log.km)
anova(fit)
Analysis of Variance Table

Response: m2.CHS$IBD
                    Df  Sum Sq  Mean Sq F value    Pr(>F)    
CHS.m2.dist$log.km   1 0.10629 0.106290   115.3 < 2.2e-16 ***
Residuals          323 0.29777 0.000922                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#CHS.VS
fit <- lm(m2.CHS.VS$IBD~CHS.VS.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.VS$IBD
                      Df    Sum Sq    Mean Sq F value Pr(>F)
CHS.VS.m2.dist$log.km  1 0.0019346 0.00193458   2.726  0.106
Residuals             43 0.0305158 0.00070967     

#CHS.TI
fit <- lm(m2.CHS.TI$IBD~CHS.TI.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.TI$IBD
                       Df   Sum Sq   Mean Sq F value    Pr(>F)    
CHS.TI.m2.dist$log.km   1 0.013235 0.0132350  12.242 0.0006597 ***
Residuals             118 0.127566 0.0010811                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#CZ
fit <- lm(m2.CZ$IBD~CZ.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CZ$IBD
                   Df  Sum Sq Mean Sq F value    Pr(>F)    
CZ.m2.dist$log.km   1 0.69433 0.69433  289.25 < 2.2e-16 ***
Residuals         664 1.59388 0.00240                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```


#### mtDNA IBD vs IBelevation


## Isolation by Elevation mtDNA

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats


CHN
```
library(reshape2)
library(fields)

Fst.info <- read.table("Fst.data", header=T)

##geog dist
CHN.mtDNA.elev <- subset(Fst.info, Region=="CHN",)
CHN.mtDNA_lon.lat <- cbind(CHN.mtDNA.elev$Long, CHN.mtDNA.elev$Lat)
CHN.mtDNA.distance.matrix <- rdist.earth(CHN.mtDNA_lon.lat, miles=F)
CHN.mtDNA.m.dist <- as.matrix(CHN.mtDNA.distance.matrix)

CHN.mtDNA.m2.dist <- melt(CHN.mtDNA.m.dist)[melt(upper.tri(CHN.mtDNA.m.dist))$value,]
names(CHN.mtDNA.m2.dist) <- c("c1", "c2", "distance")
CHN.mtDNA.m2.dist$log.km <- log(CHN.mtDNA.m2.dist$distance)



##elev dist
CHN.mtDNA.elev <- subset(Fst.info, Region=="CHN",)
CHN.mtDNA.elevonly <- CHN.mtDNA.elev$Elev
CHN.mtDNA.elev.dist <- dist(CHN.mtDNA.elevonly, method="euclidean")

m.CHN.mtDNA.elev <- as.matrix(CHN.mtDNA.elev.dist)
m2.CHN.mtDNA.elev <- melt(m.CHN.mtDNA.elev)[melt(upper.tri(m.CHN.mtDNA.elev))$value,]
colnames(m2.CHN.mtDNA.elev) <- c("site1", "site2", "elev")
m2.CHN.mtDNA.elev$log.elev <- log(m2.CHN.mtDNA.elev$elev)


##PhiSt
CHN.phiST <- read.table("CHN.phiST", header=T)
CHN.phiST <- as.matrix(CHN.phiST)
CHN.mtDNA.names <- colnames(CHN.phiST)

m2.CHN.mtDNA.phist <- melt(CHN.phiST)[melt(lower.tri(CHN.phiST))$value,]
names(m2.CHN.mtDNA.phist) <- c("c1", "c2", "distance")
m2.CHN.mtDNA.phist$IBD <- (m2.CHN.mtDNA.phist$distance/(1-m2.CHN.mtDNA.phist$distance))


##IBD
pdf("CHN.mtDNA.IBD.pdf")
plot(m2.CHN.mtDNA.phist$IBD~CHN.mtDNA.m2.dist$log.km, pch=20, cex=0.5, xlab="distance (km)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHN.mtDNA.phist$IBD~CHN.mtDNA.m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance (mtDNA) plot - CHN")
dev.off()

##test significance IBD
fit <- lm(m2.CHN.mtDNA.phist$IBD~CHN.mtDNA.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHN.mtDNA.phist$IBD
                          Df  Sum Sq   Mean Sq F value Pr(>F)
CHN.mtDNA.m2.dist$log.km   1 0.00017 0.0001671  0.0129 0.9097
Residuals                103 1.33029 0.0129154 



##IB Elevation
pdf("CHN.mtDNA.Isol.Elev.pdf")
plot(m2.CHN.mtDNA.phist$IBD~m2.CHN.mtDNA.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHN.mtDNA.phist$IBD~m2.CHN.mtDNA.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation (mtDNA) plot - CHN")
dev.off()


##test significance IBelev
fit <- lm(m2.CHN.mtDNA.phist$IBD~m2.CHN.mtDNA.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHN.mtDNA.phist$IBD
                            Df  Sum Sq  Mean Sq F value Pr(>F)
m2.CHN.mtDNA.elev$log.elev   1 0.00396 0.003959  0.3074 0.5805
Residuals                  103 1.32650 0.012879
```



CHS.TI
```
##geog dist
CHS.TI.mtDNA.elev <- subset(Fst.info, Region=="CHS.TI",)
CHS.TI.mtDNA_lon.lat <- cbind(CHS.TI.mtDNA.elev$Long, CHS.TI.mtDNA.elev$Lat)
CHS.TI.mtDNA.distance.matrix <- rdist.earth(CHS.TI.mtDNA_lon.lat, miles=F)
CHS.TI.mtDNA.m.dist <- as.matrix(CHS.TI.mtDNA.distance.matrix)

CHS.TI.mtDNA.m2.dist <- melt(CHS.TI.mtDNA.m.dist)[melt(upper.tri(CHS.TI.mtDNA.m.dist))$value,]
names(CHS.TI.mtDNA.m2.dist) <- c("c1", "c2", "distance")
CHS.TI.mtDNA.m2.dist$log.km <- log(CHS.TI.mtDNA.m2.dist$distance)



##elev dist
CHS.TI.mtDNA.elev <- subset(Fst.info, Region=="CHS.TI",)
CHS.TI.mtDNA.elevonly <- CHS.TI.mtDNA.elev$Elev
CHS.TI.mtDNA.elev.dist <- dist(CHS.TI.mtDNA.elevonly, method="euclidean")

m.CHS.TI.mtDNA.elev <- as.matrix(CHS.TI.mtDNA.elev.dist)
m2.CHS.TI.mtDNA.elev <- melt(m.CHS.TI.mtDNA.elev)[melt(upper.tri(m.CHS.TI.mtDNA.elev))$value,]
colnames(m2.CHS.TI.mtDNA.elev) <- c("site1", "site2", "elev")
m2.CHS.TI.mtDNA.elev$log.elev <- log(m2.CHS.TI.mtDNA.elev$elev)


##PhiSt
CHS.TI.phiST <- read.table("CHS.TI.phiST", header=T)
CHS.TI.phiST <- as.matrix(CHS.TI.phiST)
CHS.TI.mtDNA.names <- colnames(CHS.TI.phiST)

m2.CHS.TI.mtDNA.phist <- melt(CHS.TI.phiST)[melt(lower.tri(CHS.TI.phiST))$value,]
names(m2.CHS.TI.mtDNA.phist) <- c("c1", "c2", "distance")
m2.CHS.TI.mtDNA.phist$IBD <- (m2.CHS.TI.mtDNA.phist$distance/(1-m2.CHS.TI.mtDNA.phist$distance))


##IBD
pdf("CHS.TI.mtDNA.IBD.pdf")
plot(m2.CHS.TI.mtDNA.phist$IBD~CHS.TI.mtDNA.m2.dist$log.km, pch=20, cex=0.5, xlab="distance (km)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHS.TI.mtDNA.phist$IBD~CHS.TI.mtDNA.m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance (mtDNA) plot - CHS.TI")
dev.off()

##test significance IBD
fit <- lm(m2.CHS.TI.mtDNA.phist$IBD~CHS.TI.mtDNA.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.TI.mtDNA.phist$IBD
                            Df  Sum Sq Mean Sq F value Pr(>F)
CHS.TI.mtDNA.m2.dist$log.km  1   0.104 0.10429   0.076 0.7835
Residuals                   89 122.206 1.37310



##IB Elevation
pdf("CHS.TI.mtDNA.Isol.Elev.pdf")
plot(m2.CHS.TI.mtDNA.phist$IBD~m2.CHS.TI.mtDNA.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHS.TI.mtDNA.phist$IBD~m2.CHS.TI.mtDNA.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation (mtDNA) plot - CHS.TI")
dev.off()


##test significance IBelev
fit <- lm(m2.CHS.TI.mtDNA.phist$IBD~m2.CHS.TI.mtDNA.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.TI.mtDNA.phist$IBD
                              Df  Sum Sq Mean Sq F value Pr(>F)
m2.CHS.TI.mtDNA.elev$log.elev  1   0.179 0.17859  0.1301 0.7191
Residuals                     89 122.132 1.37227    
```

CHS.VS
```
##geog dist
CHS.VS.mtDNA.elev <- subset(Fst.info, Region=="CHS.VS",)
CHS.VS.mtDNA_lon.lat <- cbind(CHS.VS.mtDNA.elev$Long, CHS.VS.mtDNA.elev$Lat)
CHS.VS.mtDNA.distance.matrix <- rdist.earth(CHS.VS.mtDNA_lon.lat, miles=F)
CHS.VS.mtDNA.m.dist <- as.matrix(CHS.VS.mtDNA.distance.matrix)

CHS.VS.mtDNA.m2.dist <- melt(CHS.VS.mtDNA.m.dist)[melt(upper.tri(CHS.VS.mtDNA.m.dist))$value,]
names(CHS.VS.mtDNA.m2.dist) <- c("c1", "c2", "distance")
CHS.VS.mtDNA.m2.dist$log.km <- log(CHS.VS.mtDNA.m2.dist$distance)



##elev dist
CHS.VS.mtDNA.elev <- subset(Fst.info, Region=="CHS.VS",)
CHS.VS.mtDNA.elevonly <- CHS.VS.mtDNA.elev$Elev
CHS.VS.mtDNA.elev.dist <- dist(CHS.VS.mtDNA.elevonly, method="euclidean")

m.CHS.VS.mtDNA.elev <- as.matrix(CHS.VS.mtDNA.elev.dist)
m2.CHS.VS.mtDNA.elev <- melt(m.CHS.VS.mtDNA.elev)[melt(upper.tri(m.CHS.VS.mtDNA.elev))$value,]
colnames(m2.CHS.VS.mtDNA.elev) <- c("site1", "site2", "elev")
m2.CHS.VS.mtDNA.elev$log.elev <- log(m2.CHS.VS.mtDNA.elev$elev)


##PhiSt
CHS.VS.phiST <- read.table("CHS.VS.phiST", header=T)
CHS.VS.phiST <- as.matrix(CHS.VS.phiST)
CHS.VS.mtDNA.names <- colnames(CHS.VS.phiST)

m2.CHS.VS.mtDNA.phist <- melt(CHS.VS.phiST)[melt(lower.tri(CHS.VS.phiST))$value,]
names(m2.CHS.VS.mtDNA.phist) <- c("c1", "c2", "distance")
m2.CHS.VS.mtDNA.phist$IBD <- (m2.CHS.VS.mtDNA.phist$distance/(1-m2.CHS.VS.mtDNA.phist$distance))


##IBD
pdf("CHS.VS.mtDNA.IBD.pdf")
plot(m2.CHS.VS.mtDNA.phist$IBD~CHS.VS.mtDNA.m2.dist$log.km, pch=20, cex=0.5, xlab="distance (km)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHS.VS.mtDNA.phist$IBD~CHS.VS.mtDNA.m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance (mtDNA) plot - CHS.VS")
dev.off()

##test significance IBD
fit <- lm(m2.CHS.VS.mtDNA.phist$IBD~CHS.VS.mtDNA.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.VS.mtDNA.phist$IBD
                            Df     Sum Sq    Mean Sq F value Pr(>F)
CHS.VS.mtDNA.m2.dist$log.km  1 1.1277e+10 1.1277e+10  0.1192 0.7337
Residuals                   19 1.7982e+12 9.4644e+10  



##IB Elevation
pdf("CHS.VS.mtDNA.Isol.Elev.pdf")
plot(m2.CHS.VS.mtDNA.phist$IBD~m2.CHS.VS.mtDNA.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CHS.VS.mtDNA.phist$IBD~m2.CHS.VS.mtDNA.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation (mtDNA) plot - CHS.VS")
dev.off()


##test significance IBelev
fit <- lm(m2.CHS.VS.mtDNA.phist$IBD~m2.CHS.VS.mtDNA.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CHS.VS.mtDNA.phist$IBD
                              Df     Sum Sq    Mean Sq F value Pr(>F)
m2.CHS.VS.mtDNA.elev$log.elev  1 9.8464e+09 9.8464e+09   0.104 0.7507
Residuals                     19 1.7997e+12 9.4720e+10      
```


CZ
```
##geog dist
CZ.mtDNA.elev <- subset(Fst.info, Region=="CZ",)
CZ.mtDNA_lon.lat <- cbind(CZ.mtDNA.elev$Long, CZ.mtDNA.elev$Lat)
CZ.mtDNA.distance.matrix <- rdist.earth(CZ.mtDNA_lon.lat, miles=F)
CZ.mtDNA.m.dist <- as.matrix(CZ.mtDNA.distance.matrix)

CZ.mtDNA.m2.dist <- melt(CZ.mtDNA.m.dist)[melt(upper.tri(CZ.mtDNA.m.dist))$value,]
names(CZ.mtDNA.m2.dist) <- c("c1", "c2", "distance")
CZ.mtDNA.m2.dist$log.km <- log(CZ.mtDNA.m2.dist$distance)



##elev dist

CZ.mtDNA.elev <- subset(Fst.info, Region=="CZ",)
CZ.mtDNA.elevonly <- CZ.mtDNA.elev$Elev
CZ.mtDNA.elev.dist <- dist(CZ.mtDNA.elevonly, method="euclidean")

m.CZ.mtDNA.elev <- as.matrix(CZ.mtDNA.elev.dist)
m2.CZ.mtDNA.elev <- melt(m.CZ.mtDNA.elev)[melt(upper.tri(m.CZ.mtDNA.elev))$value,]
colnames(m2.CZ.mtDNA.elev) <- c("site1", "site2", "elev")
m2.CZ.mtDNA.elev$log.elev <- log(m2.CZ.mtDNA.elev$elev)


##PhiSt
CZ.phiST <- read.table("CZ.phiST", header=T)
CZ.phiST <- as.matrix(CZ.phiST)
CZ.mtDNA.names <- colnames(CZ.phiST)

m2.CZ.mtDNA.phist <- melt(CZ.phiST)[melt(lower.tri(CZ.phiST))$value,]
names(m2.CZ.mtDNA.phist) <- c("c1", "c2", "distance")
m2.CZ.mtDNA.phist$IBD <- (m2.CZ.mtDNA.phist$distance/(1-m2.CZ.mtDNA.phist$distance))


##IBD
pdf("CZ.mtDNA.IBD.pdf")
plot(m2.CZ.mtDNA.phist$IBD~CZ.mtDNA.m2.dist$log.km, pch=20, cex=0.5, xlab="distance (km)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CZ.mtDNA.phist$IBD~CZ.mtDNA.m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance (mtDNA) plot - CZ")
dev.off()

##test significance IBD
fit <- lm(m2.CZ.mtDNA.phist$IBD~CZ.mtDNA.m2.dist$log.km)
anova(fit)

Analysis of Variance Table

Response: m2.CZ.mtDNA.phist$IBD
                         Df     Sum Sq    Mean Sq F value Pr(>F)
CZ.mtDNA.m2.dist$log.km   1 5.7630e+10 5.7630e+10  0.9666 0.3262
Residuals               376 2.2418e+13 5.9624e+10     



##IB Elevation
pdf("CZ.mtDNA.Isol.Elev.pdf")
plot(m2.CZ.mtDNA.phist$IBD~m2.CZ.mtDNA.elev$log.elev, pch=20, cex=0.5, xlab="elevation (m)", ylab="Fst/(1-Fst)")
abline(fit <- lm(m2.CZ.mtDNA.phist$IBD~m2.CZ.mtDNA.elev$log.elev))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by elevation (mtDNA) plot - CZ")
dev.off()


##test significance IBelev
fit <- lm(m2.CZ.mtDNA.phist$IBD~m2.CZ.mtDNA.elev$log.elev)
anova(fit)

Analysis of Variance Table

Response: m2.CZ.mtDNA.phist$IBD
                           Df     Sum Sq    Mean Sq F value Pr(>F)
m2.CZ.mtDNA.elev$log.elev   1 8.1185e+09 8.1185e+09  0.1359 0.7126
Residuals                 376 2.2468e+13 5.9755e+10  
```












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


### 4. Connectivity across elevation

*** I didn't finish this. Converting between packages and adding population and coordinate information was annoying. I'll add this into a next paper. 

I've already described the Fst. But now I will use a network analysis to determine how connected populations are across elevation. ie. are high elevation populations more isolated? 

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/Dataset1_forPCA/summstats

I'm working in this directory because I already have all the genind objects created in this R project. 

I'll use the package EcoGenetics. 


```
#install.packages("EcoGenetics")
library(EcoGenetics)




```




