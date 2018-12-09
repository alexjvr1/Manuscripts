# Plot MaxEnt results

Use results from Charlotte & Elsa's 3rd year project MaxEnt output to display MaxEnt output.

MaxEnt run for Current climate, 2050 and 2070. 
See [here](https://github.com/alexjvr1/Manuscripts/blob/master/ENM_CharlottesCode.md) for Charlotte's code. 


Some resources for MaxEnt in R

https://www.molecularecologist.com/2013/04/species-distribution-models-in-r/

Function to import maxent models. From [here](https://rdrr.io/github/johnbaums/rmaxent/src/R/import_maxent.R)
```
import_maxent <- function(dir, lambdas, html) {
  l <- ifelse(missing(lambdas), 
         list.files(dir, '\\.lambdas$', full.names=TRUE)[1], 
         file.path(dir, basename(lambdas)))
  h <- ifelse(missing(html), 
              list.files(dir, '\\.html$', full.names=TRUE)[1], 
              file.path(dir, basename(html)))
  m <- methods::new('MaxEnt',
           lambdas=readLines(l),
           results=t(utils::read.csv(file.path(dir, 'maxentResults.csv'))[, -1]),
           path=dir,
           html=h,
           hasabsence=file.exists(file.path(dir, 'absence')))
  m@absence <- read.csv(file.path(dir, 'absence'))[, -(1:3)]
  m@presence <- utils::read.csv(file.path(dir, 'presence'))[, -(1:3)]
  m
}

```

### Aims: 

1. Calculate niche overlap between regions within Switzerland. (using Schroeder's D)

2. Calculate change in suitable cells between current and future climates. Display as boxplots

3. Plot projections of combined models for a specific time period (present, 2050, 2070)



## 1.Niche Overlap


Use current predicted habitat to determine the overlap in niche between 1. CHall, 2. CHN, 3. CHSW, 4. CHSE. 

CHall = model calculated from sites across all of CH

CHN = model based on sites north of the Alps

CHSW = South West CH (Valais)

CHSE = South East CH (Ticino)




## 2. Change in suitable cells

Calculate the change in suitable habitat under future climate conditions. 

1. Import MaxEnt model for current distribution. 

2. Calculate number of suitable cells based on threshold

3. Import maxEnt model for future climates

4. Calculate number of suitable cells. 

5. Calculate % change for 2050 and 2070. 

6. Boxplot of the results. 



#### CH Full

working directory: 
/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland


load libraries and MaxEnt
```
library(raster)
library(dismo)
library(sp)
library(jsonlite)
library(rgdal)
library(spaa)

## MaxEnt needs to be downloaded and maxent.jar should be moved to the dismo package directory:  
## system.file("java", package="dismo")  ##use this to find the directory

MaxEnt <- system.file("java",package="dismo")  ##get MaxEnt
```

Import MaxEnt output into R
```
#install.packages("devtools")
library(devtools)
install_github("johnbaums/rmaxent")


##2050


CCSM.45MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland")
```

Get Future climate data for 2050 and 2070
```
CCSM.26_2050=getData('CMIP5', var='bio', res=2.5, rcp=26, model='CC', year=50)
CCSM.45_2050=getData('CMIP5', var='bio', res=2.5, rcp=45, model='CC', year=50)
CCSM.60_2050=getData('CMIP5', var='bio', res=2.5, rcp=60, model='CC', year=50)
CCSM.85_2050=getData('CMIP5', var='bio', res=2.5, rcp=85, model='CC', year=50)
CSM1.26_2050=getData('CMIP5', var='bio', res=2.5, rcp=26, model='BC', year=50)
CSM1.45_2050=getData('CMIP5', var='bio', res=2.5, rcp=45, model='BC', year=50)
CSM1.60_2050=getData('CMIP5', var='bio', res=2.5, rcp=60, model='BC', year=50)
CSM1.85_2050=getData('CMIP5', var='bio', res=2.5, rcp=85, model='BC', year=50)


CCSM.26_2070=getData('CMIP5', var='bio', res=2.5, rcp=26, model='CC', year=70)
CCSM.45_2070=getData('CMIP5', var='bio', res=2.5, rcp=45, model='CC', year=70)
CCSM.60_2070=getData('CMIP5', var='bio', res=2.5, rcp=60, model='CC', year=70)
CCSM.85_2070=getData('CMIP5', var='bio', res=2.5, rcp=85, model='CC', year=70)
CSM1.26_2070=getData('CMIP5', var='bio', res=2.5, rcp=26, model='BC', year=70)
CSM1.45_2070=getData('CMIP5', var='bio', res=2.5, rcp=45, model='BC', year=70)
CSM1.60_2070=getData('CMIP5', var='bio', res=2.5, rcp=60, model='BC', year=70)
CSM1.85_2070=getData('CMIP5', var='bio', res=2.5, rcp=85, model='BC', year=70)


#rename all the bioclim variables
names(CCSM.26_2050)=names(climate)
names(CCSM.45_2050)=names(climate)
names(CCSM.60_2050)=names(climate)
names(CCSM.85_2050)=names(climate)
names(CSM1.26_2050)=names(climate)
names(CSM1.45_2050)=names(climate)
names(CSM1.60_2050)=names(climate)
names(CSM1.85_2050)=names(climate)

names(CCSM.26_2070)=names(climate)
names(CCSM.45_2070)=names(climate)
names(CCSM.60_2070)=names(climate)
names(CCSM.85_2070)=names(climate)
names(CSM1.26_2070)=names(climate)
names(CSM1.45_2070)=names(climate)
names(CSM1.60_2070)=names(climate)
names(CSM1.85_2070)=names(climate)

```
