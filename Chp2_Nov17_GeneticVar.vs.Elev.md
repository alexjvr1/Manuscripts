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
[1] -146.3384
AIC(m2.RAD.Fis)
[1] -139.51

AIC(m1.RAD.Hs)
[1] -384.6909
AIC(m2.RAD.Hs)
[1] -377.7738

AIC(m1.RAD.Ho)
[1] -298.8624
AIC(m2.RAD.Ho)
[1] -292.3717
```

All m1 models

Assess model
```
library(car)

Anova(m1.RAD.Fis)

Anova Table (Type II tests)

Response: RAD.data$FIS
                                 Sum Sq Df F value Pr(>F)
RAD.data$elev.c                 0.01293  1  1.5139 0.2225
RAD.data$Region                 0.04600  3  1.7951 0.1556
RAD.data$elev.c:RAD.data$Region 0.01545  3  0.6029 0.6152
Residuals                       0.62357 73  


Anova(m1.RAD.Hs)

Anova Table (Type II tests)

Response: RAD.data$Hs
                                  Sum Sq Df F value    Pr(>F)    
RAD.data$elev.c                 0.001313  1  2.9147   0.09203 .  
RAD.data$Region                 0.022724  3 16.8178 2.105e-08 ***
RAD.data$elev.c:RAD.data$Region 0.003732  3  2.7620   0.04809 *  
Residuals                       0.032879 73                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


Anova(m1.RAD.Ho)
Anova Table (Type II tests)

Response: RAD.data$Ho
                                  Sum Sq Df F value    Pr(>F)    
RAD.data$elev.c                 0.000039  1  0.0304    0.8622    
RAD.data$Region                 0.033237  3  8.5256 6.299e-05 ***
RAD.data$elev.c:RAD.data$Region 0.004098  3  1.0512    0.3752    
Residuals                       0.094864 73                      
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
multiplot(p1.Fis.RAD, p2.Hs.RAD, p3.Ho.RAD, cols=1)
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

CHall: scai

CHN: rade

CHS.VS: orge

CHS.TI: scai

CZ: ente

```



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
