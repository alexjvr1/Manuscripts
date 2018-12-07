# ENM in MaxEnt

Code by Charlotte (3rd year project, 2018). 

```
#Switzerland 

#Installing packages
install.packages(c('raster', 'rgdal', 'dismo', 'rJava', 'jsonlite'))
install.packages('sp') #Spacial Data Analysis 

#Libraries
library(raster)
library(dismo)
library(sp)
library(jsonlite)
library(rgdal)
library(spaa)

#Getting MaxEnt
MaxEnt <- system.file("java",package="dismo")
```


### Get data

Climate 
```
#Climate Data 
climate <- getData('worldclim', var='bio', res=2.5) #No Need to convert to Rasta- What set of Bioclim is this?
climate2 <- crop(climate, extent(5.8,10.6,45.5,47.9)) ##crop to map to Switzerland
plot(Climate2)

```

Occurrance data based on sampled sites and GBIF
```
# Alex Occurance data - Switzerland
#xy <- Envi[, c("long","lat")]  #extract coordinates
#xy.forbioclim <- xy[,c(1,2)] #rename xy for rest of code
#spdf <- SpatialPointsDataFrame(coords = xy.forbioclim, data=xy, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#spdf = switzerland, swespdf = Sweden


#Gibif Data Download 
RtempGeo = gbif("Rana", "temporaria", geo=F)    ##download all the data from GBIF - this took around 3 hours on my home wifi, should be faster at Uni
acgeo <- RtempGeo
dim(RtempGeo)    #dimensions of the data. Each row is a record, each column some info.. 

#Switzerland GIBIF
##Subset for Switzerlansw
acgeo.SW <- acgeo[which(acgeo$lat >=25),]
acgeo.SW <- acgeo.SW[which(acgeo.SW$lat <48),]
plot(wrld_simpl, xlim=c(6,11), ylim=c(45.5,48), axes=T)
points(acgeo.SW$lon, acgeo.SW$lat, col='red', pch=20, cex=0.75)

##Cleaning- cleaned down to Switzerland, after 1970 and by human observation as is Alex's data
acgeo.SW <- subset(acgeo.SW, basisOfRecord=="HUMAN_OBSERVATION")
acgeo.SW <- subset(acgeo.SW, country=="Switzerland")
acgeo.SW <- subset(acgeo.SW, year>=1970)

##Replot using this subset
plot(wrld_simpl, xlim=c(6,11), ylim=c(45.5,48), axes=T)
points(acgeo.SW$lon, acgeo.SW$lat, col='red', pch=20, cex=0.75)

#Create CSV
write.csv(acgeo.SW, "C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/acgeo_SW.csv")

#Merg With Alex's Data!
GIBIF.SW <- read.csv(file.choose())#Combined presence data
xyGIBSW <- GIBIF.SW[, c("long","lat")]  #extract coordinates
xy1 <- xyGIBSW[,c(1,2)] #rename xy for rest of code
SW <- SpatialPointsDataFrame(coords = xy1, data=xy1, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
```


Reduce the number of climate variables in the dateset by removing correlated variables. 
```
#Correlation of Climatatic Factors 
install.packages("corrr")  
library(corrr)

###Reduce number of variables using 
#Prep Data 
values.SW <- extract(climate2,SW)#Climate Files and Presence
df.SW <- cbind.data.frame(coordinates(SW),values.SW)
head(df.SW)

write.csv(df.SW, "C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/df_SW.csv")

#Library
install.packages('corrplot')
install.packages('caret')
library(corrplot)
library(caret)

datMy.SW <- read.csv("C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/df_SW.csv", header=TRUE)#df_SW.csv
datMy.SW <- datMy.SW[, 3:21]
datMy.scale.SW <- scale(datMy.SW[1:ncol(datMy.SW)], center=T, scale = T) #scale all the features (from feature 2 bacause feature 1 is the predictor output)
corMatMy.SW <- cor(datMy.scale.SW, method = "spearman") #compute the correlation matrix
corrplot(corMatMy.SW, order = "hclust", tl.cex = 1)
  #visualize the matrix, clustering features by correlation index.

highlyCor.SW <- findCorrelation(corMatMy.SW, 0.80) #After inspection, apply correlation filter at 0.70,
#then we remove all the variable correlated with more 0.8.
datMyFiltered.scale.SW <- datMy.scale.SW[,-highlyCor.SW]
corMatMy.SW <- cor(datMyFiltered.scale.SW)
corrplot(corMatMy.SW, order = "hclust", tl.cex=1)

##Subset by Genome 
Genomic <- read.csv(file.choose())#swiss.csv
SouthEast <- subset(Genomic, region=="SouthEast")
SouthWest <- subset(Genomic, region=="SouthWest")
North <- subset(Genomic, region=="North")

#Cutting down bioclim data 
#Subsetting Bioclims
bio8 <- file.choose()
bio9 <- file.choose()
bio13 <- file.choose()
bio15 <- file.choose()
bio18 <- file.choose()

Subset.Bio <- stack(bio8,bio9,bio13,bio15,bio18)
Subset.Bio <- crop(Subset.Bio, extent(5.8,10.6,45.5,47.9))
```


### Run MaxEnt Switzerland

First this was run before reducing the number of climate variables
```
meGIBIFSwi <- maxent(climate2, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output')
```


#### Full model for Switzerland. 
```
SwitzerlandFull <- maxent(Subset.Bio, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Full Switzerland/Full Switzerland- 5 Bios')

#Prediction
Pred.Pres.Full <- predict(SwitzerlandFull, Subset.Bio)#What goes here??? Why doesn't it work!?
plot(Pred.Pres.Full)
```


#### SUBSETS FOR N/SE/SW

Subset the climate data
```
#N
Present.Bio.N <- crop(Subset.Bio, extent(6.1,10.1,46.7,47.5))
#SE
Present.Bio.SE <- crop(Subset.Bio, extent(8.0,10.1,46,47.5))
#SW
Present.Bio.SW <- crop(Subset.Bio, extent(6.1,8,46,47.5))
```

### Future climates

Download all BioClim data for future climate scenarios. 

2050
```
Full Future 
#Library 

#Future 2050
#CCSM.26 2050
bio8CCSM.26 <- file.choose()
bio9CCSM.26 <- file.choose()
bio13CCSM.26 <- file.choose()
bio15CCSM.26 <- file.choose()
bio18CCSM.26 <- file.choose()

CCSM.26 <- stack(bio8CCSM.26,bio9CCSM.26,bio13CCSM.26,bio15CCSM.26,bio18CCSM.26) 
CCSM.26 <- crop(CCSM.26,extent(5.8,10.6,45.5,47.9))

CCSM.26MaxEnt <- maxent(CCSM.26, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.26')

#Predict
Pred.CCSM.26.Full <- predict(CCSM.26MaxEnt, CCSM.26)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.26.Full)

#CCSM.45 2050
bio8CCSM.45 <- file.choose()
bio9CCSM.45 <- file.choose()
bio13CCSM.45 <- file.choose()
bio15CCSM.45 <- file.choose()
bio18CCSM.45 <- file.choose()

CCSM.45 <- stack(bio8CCSM.45,bio9CCSM.45,bio13CCSM.45,bio15CCSM.45,bio18CCSM.45) 
CCSM.45 <- crop(CCSM.45,extent(5.8,10.6,45.5,47.9))

CCSM.45MaxEnt <- maxent(CCSM.45, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.45')

Pred.CCSM.45.Full <- predict(CCSM.45MaxEnt, CCSM.45)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.45.Full)
     
#CCSM.60 2050
bio8CCSM.60 <- file.choose()
bio9CCSM.60 <- file.choose()
bio13CCSM.60 <- file.choose()
bio15CCSM.60 <- file.choose()
bio18CCSM.60 <- file.choose()

CCSM.60 <- stack(bio8CCSM.60,bio9CCSM.60,bio13CCSM.60,bio15CCSM.60,bio18CCSM.60) 
CCSM.60 <- crop(CCSM.60,extent(5.8,10.6,45.5,47.9))

CCSM.60MaxEnt <- maxent(CCSM.60, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.60')

Pred.CCSM.60.Full <- predict(CCSM.60MaxEnt, CCSM.60)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.60.Full)

#CCSM.85 2050
bio8CCSM.85 <- file.choose()
bio9CCSM.85 <- file.choose()
bio13CCSM.85 <- file.choose()
bio15CCSM.85 <- file.choose()
bio18CCSM.85 <- file.choose()

CCSM.85 <- stack(bio8CCSM.85,bio9CCSM.85,bio13CCSM.85,bio15CCSM.85,bio18CCSM.85) 
CCSM.85 <- crop(CCSM.60,extent(5.8,10.6,45.5,47.9))

CCSM.85MaxEnt <- maxent(CCSM.85, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.80')

Pred.CCSM.85.Full <- predict(CCSM.85MaxEnt, CCSM.85)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.85.Full)

#CSM1.26 2050
bio8csm1.26 <- file.choose()
bio9csm1.26 <- file.choose()
bio13csm1.26 <- file.choose()
bio15csm1.26 <- file.choose()
bio18csm1.26 <- file.choose()

CSM1.26 <- stack(bio8csm1.26,bio9csm1.26,bio13csm1.26,bio15csm1.26,bio18csm1.26) 
CSM1.26 <- crop(CSM1.26,extent(5.8,10.6,45.5,47.9))

CSM1.26MaxEnt <- maxent(CSM1.26, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.26')

Pred.CSM1.26.Full <- predict(CSM1.26MaxEnt, CSM1.26)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.26.Full)

#CSM1.45 2050
bio8csm1.45 <- file.choose()
bio9csm1.45 <- file.choose()
bio13csm1.45 <- file.choose()
bio15csm1.45 <- file.choose()
bio18csm1.45 <- file.choose()

CSM1.45 <- stack(bio8csm1.45,bio9csm1.45,bio13csm1.45,bio15csm1.45,bio18csm1.45) 
CSM1.45 <- crop(CSM1.45,extent(5.8,10.6,45.5,47.9))

CSM1.45MaxEnt <- maxent(CSM1.45, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.45')

Pred.CSM1.45.Full <- predict(CSM1.45MaxEnt, CSM1.45)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.45.Full)

#CSM1.60
bio8csm1.60 <- file.choose()
bio9csm1.60 <- file.choose()
bio13csm1.60 <- file.choose()
bio15csm1.60 <- file.choose()
bio18csm1.60 <- file.choose()

CSM1.60 <- stack(bio8csm1.60,bio9csm1.60,bio13csm1.60,bio15csm1.60,bio18csm1.60) 
CSM1.60 <- crop(CSM1.60,extent(5.8,10.6,45.5,47.9))

CSM1.60MaxEnt <- maxent(CSM1.60, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.60')

Pred.CSM1.60.Full <- predict(CSM1.60MaxEnt, CSM1.60)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.60.Full)

#CSM1.85
bio8csm1.85 <- file.choose()
bio9csm1.85 <- file.choose()
bio13csm1.85 <- file.choose()
bio15csm1.85 <- file.choose()
bio18csm1.85 <- file.choose()

CSM1.85 <- stack(bio8csm1.85,bio9csm1.85,bio13csm1.85,bio15csm1.85,bio18csm1.85) 
CSM1.85 <- crop(CSM1.85,extent(5.8,10.6,45.5,47.9))

CSM1.85MaxEnt <- maxent(CSM1.85, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.85')

Pred.CSM1.85.Full <- predict(CSM1.85MaxEnt, CSM1.85)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.85.Full)
```

2070
```
##Future 2070
#CCSM.26 2070
bio8CCSM.267 <- file.choose()
bio9CCSM.267 <- file.choose()
bio13CCSM.267 <- file.choose()
bio15CCSM.267 <- file.choose()
bio18CCSM.267 <- file.choose()

CCSM.267 <- stack(bio8CCSM.267,bio9CCSM.267,bio13CCSM.267,bio15CCSM.267,bio18CCSM.267) 
CCSM.267 <- crop(CCSM.267,extent(5.8,10.6,45.5,47.9))

CCSM.267MaxEnt <- maxent(CCSM.267, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.26')

Pred.CCSM.267.Full <- predict(CCSM.267MaxEnt, CCSM.267)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.267.Full)

#CCSM.45 2070
bio8CCSM.457 <- file.choose()
bio9CCSM.457 <- file.choose()
bio13CCSM.457 <- file.choose()
bio15CCSM.457 <- file.choose()
bio18CCSM.457 <- file.choose()

CCSM.457 <- stack(bio8CCSM.457,bio9CCSM.457,bio13CCSM.457,bio15CCSM.457,bio18CCSM.457) 
CCSM.457 <- crop(CCSM.457,extent(5.8,10.6,45.5,47.9))

CCSM.457MaxEnt <- maxent(CCSM.457, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.45')

Pred.CCSM.457.Full <- predict(CCSM.457MaxEnt, CCSM.457)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.457.Full)

#CCSM.60
bio8CCSM.607 <- file.choose()
bio9CCSM.607 <- file.choose()
bio13CCSM.607 <- file.choose()
bio15CCSM.607 <- file.choose()
bio18CCSM.607 <- file.choose()

CCSM.607 <- stack(bio8CCSM.607,bio9CCSM.607,bio13CCSM.607,bio15CCSM.607,bio18CCSM.607) 
CCSM.607 <- crop(CCSM.607,extent(5.8,10.6,45.5,47.9))

CCSM.607MaxEnt <- maxent(CCSM.607, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.60')

Pred.CCSM.607.Full <- predict(CCSM.607MaxEnt, CCSM.607)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.607.Full)

#CCSM.85
bio8CCSM.857 <- file.choose()
bio9CCSM.857 <- file.choose()
bio13CCSM.857 <- file.choose()
bio15CCSM.857 <- file.choose()
bio18CCSM.857 <- file.choose()

CCSM.857 <- stack(bio8CCSM.857,bio9CCSM.857,bio13CCSM.857,bio15CCSM.857,bio18CCSM.857) 
CCSM.857 <- crop(CCSM.857,extent(5.8,10.6,45.5,47.9))

CCSM.857MaxEnt <- maxent(CCSM.857, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.80')

Pred.CCSM.857.Full <- predict(CCSM.857MaxEnt, CCSM.857)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.857.Full)

#CSM1.26
bio8CSM1.267 <- file.choose()
bio9CSM1.267 <- file.choose()
bio13CSM1.267 <- file.choose()
bio15CSM1.267 <- file.choose()
bio18CSM1.267 <- file.choose()

CSM1.267 <- stack(bio8CSM1.267,bio9CSM1.267,bio13CSM1.267,bio15CSM1.267,bio18CSM1.267) 
CSM1.267 <- crop(CSM1.267,extent(5.8,10.6,45.5,47.9))

CSM1.267MaxEnt <- maxent(CSM1.267, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.26')

Pred.CSM1.267.Full <- predict(CSM1.267MaxEnt, CSM1.267)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.267.Full)

#CSM1.45
bio8CSM1.457 <- file.choose()
bio9CSM1.457 <- file.choose()
bio13CSM1.457 <- file.choose()
bio15CSM1.457 <- file.choose()
bio18CSM1.457 <- file.choose()

CSM1.457 <- stack(bio8CSM1.457,bio9CSM1.457,bio13CSM1.457,bio15CSM1.457,bio18CSM1.457) 
CSM1.457 <- crop(CSM1.457,extent(5.8,10.6,45.5,47.9))

CSM1.457MaxEnt <- maxent(CSM1.457, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.45')

Pred.CSM1.457.Full <- predict(CSM1.457MaxEnt, CSM1.457)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.457.Full)

#CSM1.60
bio8CSM1.607 <- file.choose()
bio9CSM1.607 <- file.choose()
bio13CSM1.607 <- file.choose()
bio15CSM1.607 <- file.choose()
bio18CSM1.607 <- file.choose()

CSM1.607 <- stack(bio8CSM1.607,bio9CSM1.607,bio13CSM1.607,bio15CSM1.607,bio18CSM1.607) 
CSM1.607 <- crop(CSM1.607,extent(5.8,10.6,45.5,47.9))

CSM1.607MaxEnt <- maxent(CSM1.607, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.60')

Pred.CSM1.607.Full <- predict(CSM1.607MaxEnt, CSM1.607)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.607.Full)

#CSM1.85
bio8CSM1.857 <- file.choose()
bio9CSM1.857 <- file.choose()
bio13CSM1.857 <- file.choose()
bio15CSM1.857 <- file.choose()
bio18CSM1.857 <- file.choose()

CSM1.857 <- stack(bio8CSM1.857,bio9CSM1.857,bio13CSM1.857,bio15CSM1.857,bio18CSM1.857) 
CSM1.857 <- crop(CSM1.607,extent(5.8,10.6,45.5,47.9))

CSM1.857MaxEnt <- maxent(CSM1.857, SW, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.85')

Pred.CSM1.857.Full <- predict(CSM1.857MaxEnt, CSM1.857)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.857.Full)

```



#### Split by Genomic region

North present
```
# Genomic Cropped
#Genomic Group- North

#Present
North.Pres <- crop(SW, extent(6.1,10.1,46.7,47.5)) #Presence
Present.Bio.N <- crop(Subset.Bio, extent(6.1,10.1,46.7,47.5))

Present.North <- maxent(Present.Bio.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/Present')

Pred.Present.North <- predict(Present.North, Present.Bio.N)#What goes here??? Why doesn't it work!?
plot(Pred.Present.North)
```

North Future 2050
```
#CCSM1.26
CCSM.26.N <- crop(CCSM.26, extent(6.1,10.1,46.7,47.5)) #Bioclim
CCSM.26MaxEnt.N <- maxent(CCSM.26.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.26')

Pred.CCSM.26.North <- predict(CCSM.26MaxEnt.N, CCSM.26.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.26.North)

#CCSM1.45
CCSM.45.N <- crop(CCSM.45,extent(6.1,10.1,46.7,47.5))
CCSM.45MaxEnt.N <- maxent(CCSM.45.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.45')

Pred.CCSM.45.North <- predict(CCSM.45MaxEnt.N, CCSM.45.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.45.North)

#CCSM1.60
CCSM.60.N <- crop(CCSM.60,extent(6.1,10.1,46.7,47.5))
CCSM.60MaxEnt.N <- maxent(CCSM.60.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.60')

Pred.CCSM.60.North <- predict(CCSM.60MaxEnt.N, CCSM.60.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.45.North)

#CCSM1.85
CCSM.85.N <- crop(CCSM.85,extent(6.1,10.1,46.7,47.5))
CCSM.85MaxEnt.N <- maxent(CCSM.85.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.80')

Pred.CCSM.85.North <- predict(CCSM.85MaxEnt.N, CCSM.85.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.85.North)

#CSM1.26
CSM1.26.N <- crop(CSM1.26,extent(6.1,10.1,46.7,47.5))
CSM1.26MaxEnt.N <- maxent(CSM1.26.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.26')

Pred.CSM1.26.North <- predict(CSM1.26MaxEnt.N, CSM1.26.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.26.North)

#CSM1.45
CSM1.45.N <- crop(CSM1.45,extent(6.1,10.1,46.7,47.5))
CSM1.45MaxEnt.N <- maxent(CSM1.45.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.45')

Pred.CSM1.45.North <- predict(CSM1.45MaxEnt.N, CSM1.45.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.45.North)

#CSM1.60
CSM1.60.N <- crop(CSM1.60,extent(6.1,10.1,46.7,47.5))
CSM1.60MaxEnt.N <- maxent(CSM1.60.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.60')

Pred.CSM1.60.North <- predict(CSM1.60MaxEnt.N, CSM1.60.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.60.North)

#CSM1.85
CSM1.85.N <- crop(CSM1.85,extent(6.1,10.1,46.7,47.5))
CSM1.85MaxEnt.N <- maxent(CSM1.85.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.85')

Pred.CSM1.85.North <- predict(CSM1.85MaxEnt.N, CSM1.85.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.85.North)
########################################################################################
```


North Future 2070
```
#CCSM.267
CCSM.267.N <- crop(CCSM.267,extent(6.1,10.1,46.7,47.5))
CCSM.267MaxEnt.N <- maxent(CCSM.267.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.26')

Pred.CCSM.267.North <- predict(CCSM.267MaxEnt.N, CCSM.267.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.267.North)

#CCSM.457
CCSM.457.N <- crop(CCSM.457,extent(6.1,10.1,46.7,47.5))
CCSM.457MaxEnt.N <- maxent(CCSM.457.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.45')

Pred.CCSM.457.North <- predict(CCSM.457MaxEnt.N, CCSM.457.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.457.North)

#CCSM.607
CCSM.607.N <- crop(CCSM.607,extent(6.1,10.1,46.7,47.5))
CCSM.607MaxEnt.N <- maxent(CCSM.607.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.60')

Pred.CCSM.607.North <- predict(CCSM.607MaxEnt.N, CCSM.607.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.607.North)

#CCSM.857
CCSM.857.N <- crop(CCSM.857,extent(6.1,10.1,46.7,47.5))
CCSM.857MaxEnt.N <- maxent(CCSM.857.N, North.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt North Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.80')

Pred.CCSM.857.North <- predict(CCSM.857MaxEnt.N, CCSM.857.N)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.857.North)

#CSM1.267
CSM1.267.N <- crop(CSM1.267,extent(6.1,10.1,46.7,47.5))
CSM1.267MaxEnt.N <- maxent(CSM1.267.N, North.Pres, a=NULL, factors=NULL, nbg=10000)

Pred.CSM1.267.North <- predict(CSM1.267MaxEnt.N, CSM1.267.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.267.North)

#CSM1.457
CSM1.457.N <- crop(CSM1.457,extent(6.1,10.1,46.7,47.5))
CSM1.457MaxEnt.N <- maxent(CSM1.457.N, North.Pres, a=NULL, factors=NULL, nbg=10000)

Pred.CSM1.457.North <- predict(CSM1.457MaxEnt.N, CSM1.457.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.457.North)

#CSM1.607
CSM1.607.N <- crop(CSM1.607,extent(6.1,10.1,46.7,47.5))
CSM1.607MaxEnt.N <- maxent(CSM1.607.N, North.Pres, a=NULL, factors=NULL, nbg=10000)

Pred.CSM1.607.North <- predict(CSM1.607MaxEnt.N, CSM1.607.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.607.North)

#CSM1.857
CSM1.857.N <- crop(CSM1.857,extent(6.1,10.1,46.7,47.5))
CSM1.857MaxEnt.N <- maxent(CSM1.857.N, North.Pres, a=NULL, factors=NULL, nbg=10000)

Pred.CSM1.857.North <- predict(CSM1.857MaxEnt.N, CSM1.857.N)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.857.North)
```


South West present
```
SW.Pres <- crop(SW, extent(6.1,8,46.7,47.5)) #Presence
Present.Bio.SW <- crop(Subset.Bio, extent(6.1,8,46.7,47.5))
SW.pres <- maxent(Present.Bio.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/Present')

Pred.CCSM.SW <- predict(SW.pres, Present.Bio.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.SW)
```

South West Future 2050
```
#CCSM1.26
CCSM.26.SW <- crop(CCSM.26, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.26MaxEnt.SW <- maxent(CCSM.26.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.26')

Pred.CCSM.26.SW <- predict(CCSM.26MaxEnt.SW, CCSM.26.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.26.SW)

#CCSM1.45
CCSM.45.SW <- crop(CCSM.45, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.45MaxEnt.SW <- maxent(CCSM.45.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.45')

Pred.CCSM.45.SW <- predict(CCSM.45MaxEnt.SW, CCSM.45.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.45.SW)

#CCSM1.60
CCSM.60.SW <- crop(CCSM.60, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.60MaxEnt.SW <- maxent(CCSM.60.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.60')

Pred.CCSM.60.SW <- predict(CCSM.60MaxEnt.SW, CCSM.60.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.60.SW)

#CCSM1.85
CCSM.85.SW <- crop(CCSM.85, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.85MaxEnt.SW <- maxent(CCSM.85.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.80')

Pred.CCSM.85.SW <- predict(CCSM.85MaxEnt.SW, CCSM.85.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.85.SW)

#CSM.26
CSM1.26.SW <- crop(CSM1.26, extent(6.1,8,46.7,47.5))
CSM1.26MaxEnt.SW <- maxent(CSM1.26.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.26')

Pred.CSM1.26.SW <- predict(CSM1.26MaxEnt.SW, CSM1.26.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.26.SW)

#CSM.45
CSM1.45.SW <- crop(CSM1.45, extent(6.1,8,46.7,47.5))
CSM1.45MaxEnt.SW <- maxent(CSM1.45.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.45')

Pred.CSM1.45.SW <- predict(CSM1.45MaxEnt.SW, CSM1.45.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.45.SW)

#CSM.60
CSM1.60.SW <- crop(CSM1.60, extent(6.1,8,46.7,47.5))
CSM1.60MaxEnt.SW <- maxent(CSM1.60.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.60')

Pred.CSM1.60.SW <- predict(CSM1.60MaxEnt.SW, CSM1.60.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.60.SW)

#CSM.85
CSM1.85.SW <- crop(CSM1.85, extent(6.1,8,46.7,47.5))
CSM1.85MaxEnt.SW <- maxent(CSM1.85.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.85')

Pred.CSM1.85.SW <- predict(CSM1.85MaxEnt.SW, CSM1.85.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.85.SW)
```


South West Future 2070
```
#CSM1.26
CSM1.267.SW <- crop(CSM1.267,extent(6.1,8,46.7,47.5))
CSM1.267MaxEnt.SW <- maxent(CSM1.267.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.26')

Pred.CSM1.267.SW <- predict(CSM1.267MaxEnt.SW, CSM1.267.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.267.SW)

#CSM1.45
CSM1.457.SW <- crop(CSM1.457,extent(6.1,8,46.7,47.5))
CSM1.457MaxEnt.SW <- maxent(CSM1.457.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.45')

Pred.CSM1.457.SW <- predict(CSM1.457MaxEnt.SW, CSM1.457.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.457.SW)

#CSM1.60
CSM1.607.SW <- crop(CSM1.607,extent(6.1,8,46.7,47.5))
CSM1.607MaxEnt.SW <- maxent(CSM1.607.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.60')

Pred.CSM1.607.SW <- predict(CSM1.607MaxEnt.SW, CSM1.607.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.607.SW)

#CSM1.85
CSM1.857.SW <- crop(CSM1.857,extent(6.1,8,46.7,47.5))
CSM1.857MaxEnt.SW <- maxent(CSM1.857.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.85')

Pred.CSM1.857.SW <- predict(CSM1.857MaxEnt.SW, CSM1.857.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.857.SW)

#CCSM.26
CCSM.267.SW <- crop(CCSM.267, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.267MaxEnt.SW <- maxent(CCSM.267.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.26')

Pred.CCSM.267.SW <- predict(CCSM.267MaxEnt.SW, CCSM.267.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.267.SW)

#CCSM.45
CCSM.457.SW <- crop(CCSM.457, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.457MaxEnt.SW <- maxent(CCSM.457.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.45')

Pred.CCSM.457.SW <- predict(CCSM.457MaxEnt.SW, CCSM.457.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.457.SW)

#CCSM.60
CCSM.607.SW <- crop(CCSM.607, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.607MaxEnt.SW <- maxent(CCSM.607.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.60')

Pred.CCSM.607.SW <- predict(CCSM.607MaxEnt.SW, CCSM.607.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.607.SW)

#CCSM.85
CCSM.857.SW <- crop(CCSM.857, extent(6.1,8,46.7,47.5)) #Bioclim
CCSM.857MaxEnt.SW <- maxent(CCSM.857.SW, SW.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South West Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.80')

Pred.CCSM.857.SW <- predict(CCSM.857MaxEnt.SW, CCSM.857.SW)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.857.SW)
```

South East present

** Is there a mistake in the code here?? crop SW i.s.o SE
```
#South West
SE.Pres <- crop(SW, extent(8.0,10.1,46,47.5)) #Presence
Present.Bio.SE <- crop(Subset.Bio, extent(8.0,10.1,46,47.5))
SE.pres <- maxent(Present.Bio.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/Present')

Pred.East <- predict(SE.pres, Present.Bio.SE)#What goes here??? Why doesn't it work!?
plot(Pred.East)
```

South East Future 2050
```
#CCSM.26
CCSM.26.SE <- crop(CCSM.26, extent(8.0,10.1,46,47.5)) #Bioclim
CCSM.26MaxEnt.SE <- maxent(CCSM.26.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.26')

Pred.CCSM.26.SE <- predict(CCSM.26MaxEnt.SE, CCSM.26.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.26.SE)

#CCSM.45
CCSM.45.SE <- crop(CCSM.45, extent(8.0,10.1,46,47.5)) #Bioclim
CCSM.45MaxEnt.SE <- maxent(CCSM.45.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.45')

Pred.CCSM.45.SE <- predict(CCSM.45MaxEnt.SE, CCSM.45.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.45.SE)

#CCSM.60
CCSM.60.SE <- crop(CCSM.60, extent(8.0,10.1,46,47.5)) #Bioclim
CCSM.60MaxEnt.SE <- maxent(CCSM.60.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.60')

Pred.CCSM.60.SE <- predict(CCSM.60MaxEnt.SE, CCSM.60.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.60.SE)

#CCSM.85
CCSM.85.SE <- crop(CCSM.85, extent(8.0,10.1,46,47.5)) #Bioclim
CCSM.85MaxEnt.SE <- maxent(CCSM.85.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CCSM2.80')

Pred.CCSM.85.SE <- predict(CCSM.85MaxEnt.SE, CCSM.85.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.85.SE)

#CSM1.26
CSM1.26.SE <- crop(CSM1.26,extent(8.0,10.1,46,47.5))
CSM1.26MaxEnt.SE <- maxent(CSM1.26.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.26')

Pred.CSM1.26.SE <- predict(CSM1.26MaxEnt.SE, CSM1.26.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.26.SE)

#CSM1.45#### Issue?
CSM1.45.SE <- crop(CSM1.45,extent(8.0,10.1,46,47.5))
CSM1.45MaxEnt.SE <- maxent(CSM1.26.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.45')

Pred.CSM1.45.SE <- predict(CSM1.45MaxEnt.SE, CSM1.45.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.45.SE)

#CSM1.60
CSM1.60.SE <- crop(CSM1.60,extent(8.0,10.1,46,47.5))
CSM1.60MaxEnt.SE <- maxent(CSM1.60.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.60')

Pred.CSM1.60.SE <- predict(CSM1.60MaxEnt.SE, CSM1.60.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.60.SE)

#CSM1.85
CSM1.85.SE <- crop(CSM1.85,extent(8.0,10.1,46,47.5))
CSM1.85MaxEnt.SE <- maxent(CSM1.85.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2050 Switerland/CSM1.85')

Pred.CSM1.85.SE <- predict(CSM1.85MaxEnt.SE, CSM1.85.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.85.SE)
```

South East Future 2070
```
#CCSM.267
CCSM.267.SE <- crop(CCSM.267,extent(8.0,10.1,46,47.5))
CCSM.267MaxEnt.SE <- maxent(CCSM.267.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.26')

Pred.CCSM.267.SE <- predict(CCSM.267MaxEnt.SE, CCSM.267.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.267.SE)

#CCSM.457
CCSM.457.SE <- crop(CCSM.457,extent(8.0,10.1,46,47.5))
CCSM.457MaxEnt.SE <- maxent(CCSM.457.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.45')

Pred.CCSM.457.SE <- predict(CCSM.457MaxEnt.SE, CCSM.457.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.457.SE)

#CCSM.607
CCSM.607.SE <- crop(CCSM.607,extent(8.0,10.1,46,47.5))
CCSM.607MaxEnt.SE <- maxent(CCSM.607.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.60')

Pred.CCSM.607.SE <- predict(CCSM.607MaxEnt.SE, CCSM.607.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.607.SE)

#CCSM.857
CCSM.857.SE <- crop(CCSM.857,extent(8.0,10.1,46,47.5))
CCSM.857MaxEnt.SE <- maxent(CCSM.857.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CCSM2.80')

Pred.CCSM.857.SE <- predict(CCSM.857MaxEnt.SE, CCSM.857.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CCSM.857.SE)

#CSM1.26
CSM1.267.SE <- crop(CSM1.267,extent(8.0,10.1,46,47.5))
CSM1.267MaxEnt.SE <- maxent(CSM1.267.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.26')

Pred.CSM1.267.SE <- predict(CSM1.267MaxEnt.SE, CSM1.267.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.267.SE)

#CSM1.45
CSM1.457.SE <- crop(CSM1.457,extent(8.0,10.1,46,47.5))
CSM1.457MaxEnt.SE <- maxent(CSM1.457.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.45')

Pred.CSM1.457.SE <- predict(CSM1.457MaxEnt.SE, CSM1.457.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.457.SE)

#CSM1.60
CSM1.607.SE <- crop(CSM1.607,extent(8.0,10.1,46,47.5))
CSM1.607MaxEnt.SE <- maxent(CSM1.607.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.60')

Pred.CSM1.607.SE <- predict(CSM1.607MaxEnt.SE, CSM1.607.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.607.SE)

#CSM1.85
CSM1.857.SE <- crop(CSM1.857,extent(8.0,10.1,46,47.5))
CSM1.857MaxEnt.SE <- maxent(CSM1.857.SE, SE.Pres, a=NULL, factors=NULL, nbg=10000, path='C:/Users/Student/Desktop/Third Year Project/MaxEnt input+output/MaxEnt South East Subset/MaxEnt Cropped Bios Full Future 2070 Switerland/CSM1.85')

Pred.CSM1.857.SE <- predict(CSM1.857MaxEnt.SE, CSM1.857.SE)#What goes here??? Why doesn't it work!?
plot(Pred.CSM1.857.SE)

CSM1.45MaxEnt.SE
```

### Combine Future results

Stack Predictors

North
```
North.Present <- Pred.Present.North
North.2050 <- stack(Pred.CCSM.26.North,Pred.CCSM.45.North,Pred.CCSM.60.North,Pred.CCSM.85.North,Pred.CSM1.26.North,Pred.CSM1.45.North,Pred.CSM1.60.North,Pred.CSM1.85.North)
plot(North.2050)
North.2070 <- stack(Pred.CCSM.267.North,Pred.CCSM.457.North,Pred.CCSM.607.North,Pred.CCSM.857.North,Pred.CSM1.267.North,Pred.CSM1.457.North,Pred.CSM1.607.North,Pred.CSM1.857.North)
```

South East
```
SE.Present <- Pred.East 
#SE.2050 <- stack(Pred.CCSM.26.SE,Pred.CCSM.45.SE,Pred.CCSM.60.SE,Pred.CCSM.85.SE,Pred.CSM1.26.SE,Pred.CSM1.45.SE,Pred.CSM1.60.SE,Pred.CSM1.85.SE)
SE.2070 <- stack(Pred.CCSM.267.SE,Pred.CCSM.457.SE,Pred.CCSM.607.SE,Pred.CCSM.857.SE,Pred.CSM1.267.SE,Pred.CSM1.457.SE,Pred.CSM1.607.SE,Pred.CSM1.857.SE)
```


South West
```
SW.Present <- Pred.CCSM.SW
SW.2050 <- stack(Pred.CCSM.26.SW,Pred.CCSM.45.SW,Pred.CCSM.60.SW,Pred.CCSM.85.SW,Pred.CSM1.26.SW,Pred.CSM1.45.SW,Pred.CSM1.60.SW,Pred.CSM1.85.SW)
SW.2070 <- stack(Pred.CCSM.267.SW,Pred.CCSM.457.SW,Pred.CCSM.607.SW,Pred.CCSM.857.SW,Pred.CSM1.267.SW,Pred.CSM1.457.SW,Pred.CSM1.607.SW,Pred.CSM1.857.SW)
```

Full Switzerland
```
#Full
Full.Present <- Pred.Pres.Full
Full.2050 <- stack(Pred.CCSM.26.Full,Pred.CCSM.45.Full,Pred.CCSM.60.Full,Pred.CCSM.85.Full,Pred.CSM1.26.Full,Pred.CSM1.45.Full,Pred.CSM1.60.Full,Pred.CSM1.85.Full)
Full.2070 <- stack(Pred.CCSM.267.Full,Pred.CCSM.457.Full,Pred.CCSM.607.Full,Pred.CCSM.857.Full,Pred.CSM1.267.Full,Pred.CSM1.457.Full,Pred.CSM1.607.Full,Pred.CSM1.857.Full)
```


## Niche similarity

```
#Schoener’s D test Code- present only 
Full.SW <- modOverlap(SE.Present, SW.Present, na.rm = TRUE) 
{
  p1 <- North.Present/sum(North.Present)
  p2 <- SE.Present/sum(SE.Present)
  SchoenerD <- 1 - 0.5 * sum(abs(p1 - p2))
  HellingerDist <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
  WarrenI <- 1 - ((HellingerDist^2)/2)
  list(SchoenerD = SchoenerD, WarrenI = WarrenI, HellingerDist = HellingerDist)
}

Full.SW <- niche.overlap.pair(Full.Present, SW.Present, method = c("schoener"))

install.packages('fuzzySim')
```




North
```
North.Present <- Pred.Present.North
North.2050 <- stack(Pred.CCSM.26.North,Pred.CCSM.45.North,Pred.CCSM.60.North,Pred.CCSM.85.North,Pred.CSM1.26.North,Pred.CSM1.45.North,Pred.CSM1.60.North,Pred.CSM1.85.North)
plot(North.2050)
North.2070 <- stack(Pred.CCSM.267.North,Pred.CCSM.457.North,Pred.CCSM.607.North,Pred.CCSM.857.North,Pred.CSM1.267.North,Pred.CSM1.457.North,Pred.CSM1.607.North,Pred.CSM1.857.North)
```


#South East
```
SE.Present <- Pred.East 
#SE.2050 <- stack(Pred.CCSM.26.SE,Pred.CCSM.45.SE,Pred.CCSM.60.SE,Pred.CCSM.85.SE,Pred.CSM1.26.SE,Pred.CSM1.45.SE,Pred.CSM1.60.SE,Pred.CSM1.85.SE)
SE.2070 <- stack(Pred.CCSM.267.SE,Pred.CCSM.457.SE,Pred.CCSM.607.SE,Pred.CCSM.857.SE,Pred.CSM1.267.SE,Pred.CSM1.457.SE,Pred.CSM1.607.SE,Pred.CSM1.857.SE)
```

#South West
```
SW.Present <- Pred.CCSM.SW
SW.2050 <- stack(Pred.CCSM.26.SW,Pred.CCSM.45.SW,Pred.CCSM.60.SW,Pred.CCSM.85.SW,Pred.CSM1.26.SW,Pred.CSM1.45.SW,Pred.CSM1.60.SW,Pred.CSM1.85.SW)
SW.2070 <- stack(Pred.CCSM.267.SW,Pred.CCSM.457.SW,Pred.CCSM.607.SW,Pred.CCSM.857.SW,Pred.CSM1.267.SW,Pred.CSM1.457.SW,Pred.CSM1.607.SW,Pred.CSM1.857.SW)
```


#Full
```
Full.Present <- Pred.Pres.Full
Full.2050 <- stack(Pred.CCSM.26.Full,Pred.CCSM.45.Full,Pred.CCSM.60.Full,Pred.CCSM.85.Full,Pred.CSM1.26.Full,Pred.CSM1.45.Full,Pred.CSM1.60.Full,Pred.CSM1.85.Full)
Full.2070 <- stack(Pred.CCSM.267.Full,Pred.CCSM.457.Full,Pred.CCSM.607.Full,Pred.CCSM.857.Full,Pred.CSM1.267.Full,Pred.CSM1.457.Full,Pred.CSM1.607.Full,Pred.CSM1.857.Full)
```

#Schoener’s D test Code- present only 
```
<- modOverlap(SE.Present, SW.Present, na.rm = TRUE) 

class       : RasterLayer 
dimensions  : 36, 50, 1800  (nrow, ncol, ncell)
resolution  : 0.04166667, 0.04166667  (x, y)
extent      : 8, 10.08333, 46, 47.5  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : in memory
names       : layer 
values      : 1, 1  (min, max)

##########################################################################
> p1 <- N.Present/sum(N.Present)
Error: object 'N.Present' not found
> p1 <- North.Present/sum(North.Present)
Warning message:
In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
> p2 <- SE.Present/sum(SE.Present)
Warning message:
In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
> SchoenerD <- 1 - 0.5 * sum(abs(p1 - p2))
Warning messages:
1: In p1 - p2 :
  Raster objects have different extents. Result for their intersection is returned
2: In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
> HellingerDist <- sqrt(sum((sqrt(p1) - sqrt(p2))^2))
Warning messages:
1: In sqrt(p1) - sqrt(p2) :
  Raster objects have different extents. Result for their intersection is returned
2: In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
> WarrenI <- 1 - ((HellingerDist^2)/2)
> list(SchoenerD = SchoenerD, WarrenI = WarrenI, HellingerDist = HellingerDist)
$SchoenerD
class       : RasterLayer 
dimensions  : 19, 50, 950  (nrow, ncol, ncell)
resolution  : 0.04166667, 0.04166667  (x, y)
extent      : 8, 10.08333, 46.70833, 47.5  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : in memory
names       : layer 
values      : 1, 1  (min, max)


$WarrenI
class       : RasterLayer 
dimensions  : 19, 50, 950  (nrow, ncol, ncell)
resolution  : 0.04166667, 0.04166667  (x, y)
extent      : 8, 10.08333, 46.70833, 47.5  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : in memory
names       : layer 
values      : 1, 1  (min, max)


$HellingerDist
class       : RasterLayer 
dimensions  : 19, 50, 950  (nrow, ncol, ncell)
resolution  : 0.04166667, 0.04166667  (x, y)
extent      : 8, 10.08333, 46.70833, 47.5  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : in memory
names       : layer 
values      : 0, 0  (min, max)

> Full.SW <- niche.overlap.pair(Full.Present, SW.Present, method = c("schoener"))
Warning messages:
1: In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
2: In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
3: In pij - pkj :
  Raster objects have different extents. Result for their intersection is returned
4: In sum(<S4 object of class "RasterLayer">, na.rm = FALSE) :
  Nothing to summarize if you provide a single RasterLayer; see cellStats
> Full.SW
class       : RasterLayer 
dimensions  : 19, 46, 874  (nrow, ncol, ncell)
resolution  : 0.04166667, 0.04166667  (x, y)
extent      : 6.083333, 8, 46.70833, 47.5  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : in memory
names       : layer 
values      : 1, 1  (min, max)

```
