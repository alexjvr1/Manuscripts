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


#### RAD: 

Measures of diversity: Fis, Hs, Ho, Ne, clutch size


```
CHN.CHS.CZ.EnvData <- read.table("CHall.CHN.CHS.CZ.EnvData", header=T) ##read in all the EnvData for CHS.VS and CHS.TI

CHN.CHS.CZ.EnvData$elev.c <- CHN.CHS.CZ.EnvData$elev - mean(CHN.CHS.CZ.EnvData$elev) ##center the data
CHN.CHS.CZ.EnvData$elev.c <- CHN.CHS.CZ.EnvData$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

CHN.CHS.CZonly.EnvData <- CHN.CHS.CZ.EnvData[!(CHN.CHS.CZ.EnvData$Transect %in% ("CHall")),] ##remove CHS, so that we include only VS and TI in the model

CHN.CHS.CZonly.EnvData$elev.c.squared <- CHN.CHS.CZonly.EnvData$elev.c^2
```



#### mtDNA

Haplotype diversity, nucleotide diversity



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

