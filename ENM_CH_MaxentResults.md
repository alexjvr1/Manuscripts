# Plot MaxEnt results

Use results from Charlotte & Elsa's 3rd year project MaxEnt output to display MaxEnt output.

MaxEnt run for Current climate, 2050 and 2070. 
See [here](https://github.com/alexjvr1/Manuscripts/blob/master/ENM_CharlottesCode.md) for Charlotte's code. 


Some resources for MaxEnt in R

https://www.molecularecologist.com/2013/04/species-distribution-models-in-r/



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

1. Plot current distribution. 

2. Calculate number of suitable cells. 

3. Plot future climates

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




```

