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

Upload models for CHN and CHSW
```
CHN <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntNorthSubset/Present/")

CHSW <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntSouthWestSubset/Present/")
```

Subset the bioclim models
```
#N
Present.Bio.N <- crop(Subset.Bio, extent(6.1,10.1,46.7,47.5))
#SE
Present.Bio.SE <- crop(Subset.Bio, extent(8.0,10.1,46,47.5))
#SW
Present.Bio.SW <- crop(Subset.Bio, extent(6.1,8,46,47.5))
```



Calculate the niche overlap between the different regions using the R package [ecospat](https://cran.r-project.org/web/packages/ecospat/vignettes/vignette_ecospat_package.pdf)
```
#install.packages("ecospat")
library(ecospat)


##PCA of the BioClim variables for the CHN and CHSW data combined. This is for the entire region based on the extent of the bioclim data. 

df.Present.Bio.N <- as.data.frame(Present.Bio.N@data@values[, 1:5])
df.Present.Bio.SW <- as.data.frame(Present.Bio.SW@data@values[, 1:5])


pca.env <- dudi.pca(rbind(df.Present.Bio.N, df.Present.Bio.SW),scannf=F,nf=2)  #scannf=F: don't display screeplot; nf=2: keep two axes

Duality diagramm
class: pca dudi
$call: dudi.pca(df = rbind(df.Present.Bio.N, df.Present.Bio.SW), scannf = F, 
    nf = 2)

$nf: 2 axis-components saved
$rank: 5
eigen values: 2.747 1.653 0.3828 0.1936 0.02299
  vector length mode    content       
1 $cw    5      numeric column weights
2 $lw    3480   numeric row weights   
3 $eig   5      numeric eigen values  

  data.frame nrow ncol content             
1 $tab       3480 5    modified array      
2 $li        3480 2    row coordinates     
3 $l1        3480 2    row normed scores   
4 $co        5    2    column coordinates  
5 $c1        5    2    column normed scores
other elements: cent norm 




##plot the PCA

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
```

![alt_txt][PCA]

[PCA]:https://user-images.githubusercontent.com/12142475/49727061-51dbf200-fc67-11e8-9821-ece211399001.png



Predict the scores on the axes
```
# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for the CHN sites with known R.temp populations
CHN.presence.bio <- CHN@presence #write all the bioclim vars for known sites to a dataframe (make sure they're in the same order and format as with the full dataset (i.e. not transformed in this case). Here CHN is the model I imported in before (see above) 
pca.sp.N <- dudi.pca(CHN.presence.bio, scannf=F, nf=2)  ##pca on these variables 
scores.sp.N <- pca.sp.N$li

# PCA scores for the CHSW sites with known R.temp populations
SW.presence.bio <- CHSW@presence #write all the bioclim vars for known sites to a dataframe
pca.sp.SW <- dudi.pca(SW.presence.bio, scannf=F, nf=2)  ##pca on these variables 
scores.sp.SW <- pca.sp.SW$li


# PCA scores for the whole CHN study area
scores.clim.N <- suprow(pca.env,df.Present.Bio.N[,1:5])$li
dim(df.Present.Bio.N) #we added North first (rbind command above). Check how big the df and if it corresponds with the PCA rows
[1] 1824
dim(scores.clim.N)
[1] 1824    2


# PCA scores for the whole CHSW study area
scores.clim.SW <- suprow(pca.env,df.Present.Bio.SW[,1:5])$li
dim(df.Present.Bio.SW)
[1] 1656    5
dim(scores.clim.SW)
[1] 1656    2
```


Calculate the Occurrence Densities Grid with ecospat.grid.clim.dyn()
```
#For a species in CHN
# gridding the CHN niche
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim,
glob1=scores.clim.N,
sp=scores.sp.N, R=100,
th.sp=0.05)


For a species in the invaded range (Australia)
11
# gridding the invasive niche
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
glob1=scores.clim.inv,
sp=scores.sp.inv, R=100,
th.sp=0)
```



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
library(rJava)  #rJava package
library(maptools)


## MaxEnt needs to be downloaded and maxent.jar should be moved to the dismo package directory:  
## system.file("java", package="dismo")  ##use this to find the directory

MaxEnt <- system.file("java",package="dismo")  ##get MaxEnt
```

Import MaxEnt models run by Charlotte into R
```
#install.packages("devtools")
library(devtools)
install_github("johnbaums/rmaxent")
```


Currentmodel
```
SwitzerlandFull <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland")

```

ModelPredictions 
```
pred.SwitzerlandFull <- predict(SwitzerlandFull, Subset.Bio)

```


Evaluate model and choose suitability threshold
```
# withold 20% of the data for testing the model
Rtempocc=xy1
fold <- kfold(Rtempocc, k=5)
Rtemptest <- Rtempocc[fold == 1, ]
Rtemptrain <- Rtempocc[fold != 1, ]


#fit the maxent model
Rtemp.me <- maxent(Subset.Bio, Rtemptrain)

# response curves
response(Rtemp.me)
# predict to entire dataset
Rtemp.pred <- predict(Rtemp.me, Subset.Bio)
```

Evaluate the model
```
#testing the model
# background data
bg <- randomPoints(Subset.Bio, 1000) #background "pseudoabsences"

#simplest way to use 'evaluate'
e1 <- evaluate(Rtemp.me, p=Rtemptest, a=bg, x=Subset.Bio)

plot(e1, 'ROC')
pvtest <- data.frame(extract(Subset.Bio, Rtemptest))
avtest <- data.frame(extract(Subset.Bio, bg))
e2 <- evaluate(Rtemp.me, p=pvtest, a=avtest)
```




Upload all the models that Charlotte ran before
```
##2050

CCSM.26MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CCSM2.26")

CCSM.45MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CCSM2.45")

CCSM.60MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CCSM2.60")

CCSM.85MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CCSM2.80")

CSM1.26MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CSM1.26")

CSM1.45MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CSM1.45")

CSM1.60MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CSM1.60")

CSM1.85MaxEnt_2050.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2050_Switerland/CSM1.85")
```


```
##2070

CCSM.26MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CCSM2.26")

CCSM.45MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CCSM2.45")

CCSM.60MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CCSM2.60")

CCSM.85MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CCSM2.80")

CSM1.26MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CSM1.26")

CSM1.45MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CSM1.45")

CSM1.60MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CSM1.60")

CSM1.85MaxEnt_2070.FULL <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntFullSwitzerland/MaxEntCroppedBiosFullFuture_2070_Switerland/CSM1.85")

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


#names to use to compare with Charlotte's models
#rename all the bioclim variables

Charlotte.bioclimnames.cc2650 <- c("cc26bi508", "cc26bi509","cc26bi5013", "cc26bi5015", "cc26bi5018")
Charlotte.bioclimnames.cc4550 <- c("cc45bi508", "cc45bi509","cc45bi5013", "cc45bi5015", "cc45bi5018")
Charlotte.bioclimnames.cc6050 <- c("cc60bi508", "cc60bi509","cc60bi5013", "cc60bi5015", "cc60bi5018")
Charlotte.bioclimnames.cc8550 <- c("cc85bi508", "cc85bi509","cc85bi5013", "cc85bi5015", "cc85bi5018")
Charlotte.bioclimnames.bc2650 <- c("bc26bi508", "bc26bi509","bc26bi5013", "bc26bi5015", "bc26bi5018")
Charlotte.bioclimnames.bc4550 <- c("bc45bi508", "bc45bi509","bc45bi5013", "bc45bi5015", "bc45bi5018")
Charlotte.bioclimnames.bc6050 <- c("bc60bi508", "bc60bi509","bc60bi5013", "bc60bi5015", "bc60bi5018")
Charlotte.bioclimnames.bc8550 <- c("bc85bi508", "bc85bi509","bc85bi5013", "bc85bi5015", "bc85bi5018")


Charlotte.bioclimnames.cc2670 <- c("cc26bi708", "cc26bi709","cc26bi7013", "cc26bi7015", "cc26bi7018")
Charlotte.bioclimnames.cc4570 <- c("cc45bi708", "cc45bi709","cc45bi7013", "cc45bi7015", "cc45bi7018")
Charlotte.bioclimnames.cc6070 <- c("cc60bi708", "cc60bi709","cc60bi7013", "cc60bi7015", "cc60bi7018")
Charlotte.bioclimnames.cc8570 <- c("cc85bi708", "cc85bi709","cc85bi7013", "cc85bi7015", "cc85bi7018")
Charlotte.bioclimnames.bc2670 <- c("bc26bi708", "bc26bi709","bc26bi7013", "bc26bi7015", "bc26bi7018")
Charlotte.bioclimnames.bc4570 <- c("bc45bi708", "bc45bi709","bc45bi7013", "bc45bi7015", "bc45bi7018")
Charlotte.bioclimnames.bc6070 <- c("bc60bi708", "bc60bi709","bc60bi7013", "bc60bi7015", "bc60bi7018")
Charlotte.bioclimnames.bc8570 <- c("bc85bi708", "bc85bi709","bc85bi7013", "bc85bi7015", "bc85bi7018")

#names(CCSM.26_2050_CH5)=Charlotte.bioclimnames.cc2650
#names(CCSM.45_2050_CH5)=Charlotte.bioclimnames.cc4550
#names(CCSM.60_2050_CH5)=Charlotte.bioclimnames.cc6050
#names(CCSM.85_2050_CH5)=Charlotte.bioclimnames.cc8550
#names(CSM1.26_2050_CH5)=Charlotte.bioclimnames.bc2650
#names(CSM1.45_2050_CH5)=Charlotte.bioclimnames.bc4550
#names(CSM1.60_2050_CH5)=Charlotte.bioclimnames.bc6050
#names(CSM1.85_2050_CH5)=Charlotte.bioclimnames.bc8550

##I think these models might be swapped around. I've renamed the CCSM models with the bc variable names.. The same happened with CCSM.60MaxEnt_2050.FULL & CSM1.80MaxEnt_2050.FULL
#names(CCSM.26_2070_CH5)=Charlotte.bioclimnames.bc2670
#names(CCSM.45_2070_CH5)=Charlotte.bioclimnames.bc4570
#names(CCSM.60_2070_CH5)=Charlotte.bioclimnames.bc6070
#names(CCSM.85_2070_CH5)=Charlotte.bioclimnames.bc8570
#names(CSM1.26_2070_CH5)=Charlotte.bioclimnames.cc2670
#names(CSM1.45_2070_CH5)=Charlotte.bioclimnames.cc4570
#names(CSM1.60_2070_CH5)=Charlotte.bioclimnames.cc6070
#names(CSM1.85_2070_CH5)=Charlotte.bioclimnames.cc8570



##crop extent to CH
CCSM.26_2050_CH <- crop(CCSM.26_2050,extent(5.8,10.6,45.5,47.9))
CCSM.45_2050_CH <- crop(CCSM.45_2050,extent(5.8,10.6,45.5,47.9))
CCSM.60_2050_CH <- crop(CCSM.60_2050,extent(5.8,10.6,45.5,47.9))
CCSM.85_2050_CH <- crop(CCSM.85_2050,extent(5.8,10.6,45.5,47.9))
CSM1.26_2050_CH <- crop(CSM1.26_2050,extent(5.8,10.6,45.5,47.9))
CSM1.45_2050_CH <- crop(CSM1.45_2050,extent(5.8,10.6,45.5,47.9))
CSM1.60_2050_CH <- crop(CSM1.60_2050,extent(5.8,10.6,45.5,47.9))
CSM1.85_2050_CH <- crop(CSM1.85_2050,extent(5.8,10.6,45.5,47.9))

CCSM.26_2070_CH <- crop(CCSM.26_2070,extent(5.8,10.6,45.5,47.9))
CCSM.45_2070_CH <- crop(CCSM.45_2070,extent(5.8,10.6,45.5,47.9))
CCSM.60_2070_CH <- crop(CCSM.60_2070,extent(5.8,10.6,45.5,47.9))
CCSM.85_2070_CH <- crop(CCSM.85_2070,extent(5.8,10.6,45.5,47.9))
CSM1.26_2070_CH <- crop(CSM1.26_2070,extent(5.8,10.6,45.5,47.9))
CSM1.45_2070_CH <- crop(CSM1.45_2070,extent(5.8,10.6,45.5,47.9))
CSM1.60_2070_CH <- crop(CSM1.60_2070,extent(5.8,10.6,45.5,47.9))
CSM1.85_2070_CH <- crop(CSM1.85_2070,extent(5.8,10.6,45.5,47.9))

##keep only the 6 bioclim variables
important.bios <- c("bio8", "bio9", "bio13", "bio15", "bio18")

CCSM.26_2050_CH5 <- subset(CCSM.26_2050_CH,important.bios)
CCSM.45_2050_CH5 <- subset(CCSM.45_2050_CH,important.bios)
CCSM.60_2050_CH5 <- subset(CCSM.60_2050_CH,important.bios)
CCSM.85_2050_CH5 <- subset(CCSM.85_2050_CH,important.bios)
CSM1.26_2050_CH5 <- subset(CSM1.26_2050_CH,important.bios)
CSM1.45_2050_CH5 <- subset(CSM1.45_2050_CH,important.bios)
CSM1.60_2050_CH5 <- subset(CSM1.60_2050_CH,important.bios)
CSM1.85_2050_CH5 <- subset(CSM1.85_2050_CH,important.bios)

CCSM.26_2070_CH5 <- subset(CCSM.26_2070_CH,important.bios)
CCSM.45_2070_CH5 <- subset(CCSM.45_2070_CH,important.bios)
CCSM.60_2070_CH5 <- subset(CCSM.60_2070_CH,important.bios)
CCSM.85_2070_CH5 <- subset(CCSM.85_2070_CH,important.bios)
CSM1.26_2070_CH5 <- subset(CSM1.26_2070_CH,important.bios)
CSM1.45_2070_CH5 <- subset(CSM1.45_2070_CH,important.bios)
CSM1.60_2070_CH5 <- subset(CSM1.60_2070_CH,important.bios)
CSM1.85_2070_CH5 <- subset(CSM1.85_2070_CH,important.bios)

```


Run predictions for all the models: 
```
pred.CCSM.26_2050_CH <- predict(CCSM.26MaxEnt_2050.FULL, CCSM.26_2050_CH5)
pred.CCSM.45_2050_CH <- predict(CCSM.45MaxEnt_2050.FULL, CCSM.45_2050_CH5)
pred.CCSM.60_2050_CH <- predict(CCSM.60MaxEnt_2050.FULL, CCSM.60_2050_CH5) #wrong variables
pred.CCSM.85_2050_CH <- predict(CCSM.85MaxEnt_2050.FULL, CCSM.85_2050_CH5) #wrong variables
pred.CSM1.26_2050_CH <- predict(CSM1.26MaxEnt_2050.FULL, CSM1.26_2050_CH5)
pred.CSM1.45_2050_CH <- predict(CSM1.45MaxEnt_2050.FULL, CSM1.45_2050_CH5)
pred.CSM1.60_2050_CH <- predict(CSM1.60MaxEnt_2050.FULL, CSM1.60_2050_CH5)
pred.CSM1.85_2050_CH <- predict(CSM1.85MaxEnt_2050.FULL, CSM1.85_2050_CH5) 

pred.CCSM.26_2070_CH <- predict(CCSM.26MaxEnt_2070.FULL, CCSM.26_2070_CH5)
pred.CCSM.45_2070_CH <- predict(CCSM.45MaxEnt_2070.FULL, CCSM.45_2070_CH5) #wrong variables
pred.CCSM.60_2070_CH <- predict(CCSM.60MaxEnt_2070.FULL, CCSM.60_2070_CH5)
pred.CCSM.85_2070_CH <- predict(CCSM.85MaxEnt_2070.FULL, CCSM.85_2070_CH5)
pred.CSM1.26_2070_CH <- predict(CSM1.26MaxEnt_2070.FULL, CSM1.26_2070_CH5)
pred.CSM1.45_2070_CH <- predict(CSM1.45MaxEnt_2070.FULL, CSM1.45_2070_CH5)
pred.CSM1.60_2070_CH <- predict(CSM1.60MaxEnt_2070.FULL, CSM1.60_2070_CH5)
pred.CSM1.85_2070_CH <- predict(CSM1.85MaxEnt_2070.FULL, CSM1.85_2070_CH5) #wrong variables

```


Now we have all the models to plot the change in suitable cells

i.e. present climate predict, and all future climates predict. 

Determine the number of suitable cells in each model. We want to change the presence probability into binary suitability values. There are different ways to choose a thresholds for this. I'd suggest using 0.7 for this data. It's quite stringent, but let's see what it looks like. 

The predicted outputs all have an @data@values layer which gives you the raw suitability for each cell in the dataset. We can simply count all the cells above our threshold for this: 
```
#Switzerland Full

suitability.df <- ()

sum(pred.SwitzerlandFull@data@values > 0.7)
[1] 1104

sum(pred.CCSM.26_2050_CH@data@values > 0.7)
sum(pred.CCSM.45_2050_CH@data@values > 0.7)
sum(pred.CCSM.60_2050_CH@data@values > 0.7)
sum(pred.CCSM.85_2050_CH@data@values > 0.7)
sum(pred.CSM1.26_2050_CH@data@values > 0.7)
sum(pred.CSM1.45_2050_CH@data@values > 0.7)
sum(pred.CSM1.60_2050_CH@data@values > 0.7)
sum(pred.CSM1.85_2050_CH@data@values > 0.7)

sum(pred.CCSM.26_2070_CH@data@values > 0.7)
sum(pred.CCSM.45_2070_CH@data@values > 0.7)
sum(pred.CCSM.60_2070_CH@data@values > 0.7)
sum(pred.CCSM.85_2070_CH@data@values > 0.7)
sum(pred.CSM1.26_2070_CH@data@values > 0.7)
sum(pred.CSM1.45_2070_CH@data@values > 0.7)
sum(pred.CSM1.60_2070_CH@data@values > 0.7)
sum(pred.CSM1.85_2070_CH@data@values > 0.7)

```


Create a table in excel with all the info and import back into R

In the example here I've copied the CHfull results 4 times to represent the other regions, so the data is exactly the same. The last column is the proportional change relative to the present predictions. 
```
Cells.Comparison <- read.table("SuitableCellsComparison", header=T)

head(Cells.Comparison)
  region    time nr.cells proportion change.in.suitable.cells
1  Chall present     1104  0.1655172                 100.0000
2  Chall    2050     1194  0.1790105                 108.1522
3  Chall    2050     1130  0.1694153                 102.3551
4  Chall    2050     1127  0.1689655                 102.0833
5  Chall    2050     1257  0.1884558                 113.8587
6  Chall    2050     1207  0.1809595                 109.3297
```

Plot
```

```


## 3. Plot Niche

Aim: 

1. Plot niche based on current data

2. plot future niche (2050)

3. Plot future niche (2070)



We need to plot the predicted niche suitability onto a map. We'll decide on a suitability threshold to turn the predicted values into binary 1s and 0s (suitable, not suitable). This will also make combining the models easy for aims 2 and 3. 



Import all models and bioclim data
```
Upload models for CHN and CHSW

CHN <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntNorthSubset/Present/")

CHSW <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntSouthWestSubset/Present/")

CHSE <- import_maxent("/Users/alexjvr/2016RADAnalysis/2018StudentENMproject/MaxEntSouthEastSubset/Present/")


Subset Bioclim vars
#N
Present.Bio.N <- crop(Subset.Bio, extent(6.1,10.1,46.7,47.5))
#SE
Present.Bio.SE <- crop(Subset.Bio, extent(8.0,10.1,46,47.5))
#SW
Present.Bio.SW <- crop(Subset.Bio, extent(6.1,8,46,47.5))
```


predict maps for Rtemp Full, North, SoutWest, and SouthEast
```
library(dismo)

Rtemp.pred <- predict(Rtemp.me, Subset.Bio)
Rtemp.N.pred <- predict(CHN, Present.Bio.N)
Rtemp.SW.pred <- predict(CHSW, Present.Bio.SW)
Rtemp.SE.pred <- predict(CHSE, Present.Bio.SE)
```


Extract the predicted value and cell coordinates for each grid point
```
library(raster)
pAll <- rasterToPoints(Rtemp.pred) 

summary(pAll)  #see below. This gives you coordinates in cols 1&2, and the prediction value in col3 (layer). 

       x                y             layer        
 Min.   : 5.812   Min.   :45.52   Min.   :0.01539  
 1st Qu.: 6.979   1st Qu.:46.10   1st Qu.:0.18942  
 Median : 8.188   Median :46.71   Median :0.38617  
 Mean   : 8.188   Mean   :46.71   Mean   :0.41567  
 3rd Qu.: 9.396   3rd Qu.:47.31   3rd Qu.:0.60553  
 Max.   :10.562   Max.   :47.90   Max.   :0.96043 

```

Add a binary column to this with 1 for all cells > 0.7 (or whatever your threshold is). For the boxplot above we used 0.7, so remember to check that you're using the same threshold throughout. 
```
pAll$thresh <- pAll$layer
pAll[which(pAll$layer>=0.7),4] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pAll[which(pAll$layer<0.7),4] <- 0 #0 to everything else
```

Change the coordinates to spatial data points: 
```
library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)
# and their dependencies


pAll.coords <- pAll[which(pAll$thresh==1),1:2]
coordinates(pAll.coords) <- pAll.coords[,1:2]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #define the projection. You can confirm this by looking at your rasterfile
proj4string(pAll.coords) <- crs.geo  #change the coordinates to a spatial dataframe with the correct projection
```




Draw a map of CH: 
```
require(spatial.tools)
elevation<-getData("alt", country = "CH")
x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(x)
slope <- terrain(elevation, opt = "slope")
aspect <- terrain(elevation, opt = "aspect")
hill <- hillShade(slope, aspect, 40, 270)
plot(hill, col = grey(0:100/100), legend = FALSE, main = "Switzerland")

#add the spatial data points
plot(pAll.coords, pch = 20, cex = 1, add = TRUE) #you can add different colours like before when you plotted the CH map. 
```
Some of the points are outside of CH. This might be because the extent of the raster files were always in a rectangle that encompasses CH, rather than only CH. You can decide whether to keep these or not. 



#### Maps of the regions


```
pCHN <- as.data.frame(rasterToPoints(Rtemp.N.pred))
pCHSW <- as.data.frame(rasterToPoints(Rtemp.SW.pred))
pCHSE <- as.data.frame(rasterToPoints(Rtemp.SE.pred))


pCHN$thresh <- pCHN$layer
pCHN[which(pCHN$layer>=0.7),4] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHN[which(pCHN$layer<0.7),4] <- 0 #0 to everything else
pCHN.coords <- pCHN[which(pCHN$thresh==1),1:2]
coordinates(pCHN.coords) <- pCHN.coords[,1:2]  ##convert to spatial points
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #define the projection. You can confirm this by looking at your rasterfile
proj4string(pCHN.coords) <- crs.geo  #change the coordinates to a spatial dataframe with the correct projection


pCHSW$thresh <- pCHSW$layer
pCHSW[which(pCHSW$layer>=0.7),4] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHSW[which(pCHSW$layer<0.7),4] <- 0 #0 to everything else
pCHSW.coords <- pCHSW[which(pCHSW$thresh==1),1:2]
coordinates(pCHSW.coords) <- pCHSW.coords[,1:2]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #define the projection. You can confirm this by looking at your rasterfile
proj4string(pCHSW.coords) <- crs.geo  #change the coordinates to a spatial dataframe with the correct projection


pCHSE$thresh <- pCHSE$layer
pCHSE[which(pCHSE$layer>=0.7),4] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHSE[which(pCHSE$layer<0.7),4] <- 0 #0 to everything else
pCHSE.coords <- pCHSE[which(pCHSE$thresh==1),1:2]
coordinates(pCHSE.coords) <- pCHSE.coords[,1:2]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #define the projection. You can confirm this by looking at your rasterfile
proj4string(pCHSE.coords) <- crs.geo  #change the coordinates to a spatial dataframe with the correct projection

```



Draw a map of CH: 
```
require(spatial.tools)
elevation<-getData("alt", country = "CH")
x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(x)
slope <- terrain(elevation, opt = "slope")
aspect <- terrain(elevation, opt = "aspect")
hill <- hillShade(slope, aspect, 40, 270)
plot(hill, col = grey(0:100/100), legend = FALSE, main = "Switzerland")

#add the spatial data points
plot(pCHN.coords, pch = 20, cex = 1, col="blue", add = TRUE) #you can add different colours like before when you plotted the CH map. Maybe match the PCA and barplots?
plot(pCHSE.coords, pch = 20, cex = 1, col="purple", add = TRUE)
plot(pCHSW.coords, pch = 20, cex = 1, col="green", add = TRUE)
```
You can use Rcolours and pch to change this figure. 


Threshold of 0.5

![alt_txt][thresh0.5]

[thresh0.5]:https://user-images.githubusercontent.com/12142475/49741636-bd858580-fc8e-11e8-8518-df50d7a64533.png




Threshold of 0.7

![alt_txt][thresh0.7]

[thresh0.7]:https://user-images.githubusercontent.com/12142475/49741635-bd858580-fc8e-11e8-9dba-a5575f6fdebe.png



#### Plot future climates

We want to aggregate the future climates and change the colour on the plot according to the number of models that found a cell to be suitable (>0.7). So we'll add the binary suitability threshold, then add all the columns together, then create a final column called 0, 1, 2 or 3 depending on the number of models that supported that cell (up to 25%, 25-50%, 50-75%, 75-100%). We'll assign three colours to the rows we want to map (classes >25% suitability). 

SwitzerlandFull
```
#predicted maps made above
#stack all the raster layers

CHfull.2050 <- stack(pred.CCSM.26_2050_CH, pred.CCSM.45_2050_CH, pred.CSM1.26_2050_CH, pred.CSM1.45_2050_CH, pred.CSM1.60_2050_CH, pred.CSM1.85_2050_CH) #pred.CCSM.60_2050_CH and pred.CCSM.85_2050_CH removed from this code, because they need to be rerun

pCHfull.2050 <- as.data.frame(rasterToPoints(CHfull.2050)) #this gives you xy coordinates, and a layer with suitability values for each of the future climate scenarios. 

pCHfull.2050[,9:14] <- pCHfull.2050[,3:8]  # add an additional column for each layer. These will be our thresholds. 
pCHfull.2050[which(pCHfull.2050$layer.1.1>=0.7),9] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.1.1<0.7),9] <- 0 #0 to everything else
pCHfull.2050[which(pCHfull.2050$layer.2.1>=0.7),10] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.2.1<0.7),10] <- 0 #0 to everything else
pCHfull.2050[which(pCHfull.2050$layer.3.1>=0.7),11] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.3.1<0.7),11] <- 0 #0 to everything else
pCHfull.2050[which(pCHfull.2050$layer.4.1>=0.7),12] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.4.1<0.7),12] <- 0 #0 to everything else
pCHfull.2050[which(pCHfull.2050$layer.5.1>=0.7),13] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.5.1<0.7),13] <- 0 #0 to everything else
pCHfull.2050[which(pCHfull.2050$layer.6.1>=0.7),14] <- 1 #assign 1 to all cells above the threshold in a new column called thresh
pCHfull.2050[which(pCHfull.2050$layer.6.1<0.7),14] <- 0 #0 to everything else


#Add a column with aggregated presence values
pCHfull.2050$thresh <- rowSums(pCHfull.2050[9:14])


#extract all the cells where >25% of the models support presence
pCHfull.2050.coords <- pCHfull.2050[which(pCHfull.2050$thresh>1),1:2]
pCHfull.2050.names <- as.data.frame(pCHfull.2050[which(pCHfull.2050$thresh>1),15]) #get the value of these (presence val) so that we can colour accordingly
colnames(pCHfull.2050.names) <- "thresh"

pCHfull.2050.names$col <- pCHfull.2050.names$thresh
pCHfull.2050.names <- within(pCHfull.2050.names, col[thresh>5] <- "darkslategray")
pCHfull.2050.names <- within(pCHfull.2050.names, col[thresh==4] <- "darkslategray3")
pCHfull.2050.names <- within(pCHfull.2050.names, col[thresh==5] <- "darkslategray3")
pCHfull.2050.names <- within(pCHfull.2050.names, col[thresh==2] <- "darkslategray1")
pCHfull.2050.names <- within(pCHfull.2050.names, col[thresh==3] <- "darkslategray1")

coordinates(pCHfull.2050.coords) <- pCHfull.2050.coords[,1:2]
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #define the projection. You can confirm this by looking at your rasterfile
proj4string(pCHfull.2050.coords) <- crs.geo  #change the coordinates to a spatial dataframe with the correct projection

col<- as.character(pCHfull.2050.names$col) #define colours
```


Draw map
Draw a map of CH: 
```
require(spatial.tools)
elevation<-getData("alt", country = "CH")
x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(x)
slope <- terrain(elevation, opt = "slope")
aspect <- terrain(elevation, opt = "aspect")
hill <- hillShade(slope, aspect, 40, 270)
plot(hill, col = grey(0:100/100), legend = FALSE, main = "Switzerland")

#add the spatial data points
plot(pCHfull.2050.coords, pch = 20, cex = 1, col=col, add = TRUE)

```


![alt_txt][concat.map]

[concat.map]:https://user-images.githubusercontent.com/12142475/49748151-da28ba00-fc9c-11e8-841c-c9bd59f04fda.png
