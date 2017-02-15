#Gradient Forest Analysis

Here I aim to: 

1)determine the environmental and spatial determinants of variation in allele frequency across the Swedish latitudinal gradient. 

2)map the current and future patterns of genomic variation in relation to climate

I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013


##Data

I'm using 3 different datasets: 

1. Fst outliers as determined by BayEnv2 XtX, and RDA (RD1 & 2)

2. Loci associated with environment as determined by LFMM & BayEnv2 

3. "Neutral" loci (all - Data1+2)


##Presenting Results

As in the Fitzpatric paper.

1. Table1 = summary of the datasets

2. Fig1 = R2 heatmap of reletive importance of different predictor variables

3. Graph of the average allele turnover for each of the three datasets (on one Figure)

4. Graph of the SNPs for BIO5 for ENV & Fst datasets

5. PCA of env variables with sample sites on (PCA centered but not scaled)

6. Geographic map of PC space: colour according to PC space

7. Map difference between datasets (Neutral vs Fst, Neutral vs Env)

8. Pop-level vulnerability to climate change: Map mean genetic offset between current and future climates. 


I had a lot of problems creating the geographic maps. I got scripts from Karina (from Victoria Sork's lab) which created "greengrid" raster files with data extracted from 2-10k points. However, the resolution for this was always really low. It turns out I can create the maps using the scripts from Fitzpatrick & Keller 2015 by just using the raster layers directly (i.e. all the points rather than a subset of random points). 


###Table 1

The R2 can be obtained from the importance of each "Species" in the gf model:  see getAnywhere(performance.plot)
```
perf.NEUTRAL <- importance(gf.NEUTRAL.model.SEtemp, type="Species")
perf.Fst <- importance(gf.Fst.model.SEtemp, type="Species")
perf.ENV <- importance(gf.ENV.model.SEtemp, type="Species")

##mean and range
summary(perf.NEUTRAL)
summary(perf.Fst)
summary(perf.ENV)

##count the number of loci above R2 of x (here 0.5)
length(perf.NEUT.df[which(perf.NEUT.df$perf.NEUTRAL>0.5),])
length(perf.Fst.df[which(perf.Fst.df>0.5),])
length(perf.ENV.df[which(perf.ENV.df>0.5),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.NEUTRAL.model.SEtemp
gf.Fst.model.SEtemp
gf.ENVcandidates.model.SEtemp

#How many loci were originally included?


##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(2,2))
hist(perf.NEUTRAL)
hist(perf.ENV)
hist(perf.Fst)
```




###Fig1: R2 heatmap
```
#R

#Retrieve results from the individual outputs for the three GF models: 

R.sq.alldatasets <- (rowMeans(gf.ENVcandidates.model.SEtemp$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "ENV"
R.sq.Fst <- (rowMeans(gf.Fst.model.SEtemp$imp.rsq, na.rm=T))
R.sq.Fst <- as.data.frame(R.sq.Fst)
R.sq.Fst
R.sq.alldatasets$Fst <- R.sq.Fst$R.sq.Fst

R.sq.Neutral <- (rowMeans(gf.NEUTRAL.model.SEtemp$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)
R.sq.alldatasets$Neutral <- R.sq.Neutral$R.sq.Neutral

R.sq.alldatasets <- as.matrix(R.sq.alldatasets)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.alldatasets)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
```

![alt_txt][Fig1]
[Fig1]:https://cloud.githubusercontent.com/assets/12142475/22826395/85140a42-ef92-11e6-8dea-5c885369b22b.png



###Fig2: Average allelic turnover across all three datasets. 

It took me a while to find the code for the gf plots so that I could figure out what was being plotted, and so combine the plots. 
```
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "bio5") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "bio5")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "bio5")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

ymax=1.2 #set the ymax so that it's the same for the the plots to be overlaid

##bio5
pdf("bio5.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="bio5: Max temp during the warmest month", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##bio15
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "bio15") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "bio15")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "bio15")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("bio15.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="bio15: Precipitation Seasonality", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##bio13
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "bio13") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "bio13")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "bio13")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("bio13.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="bio13: Precipitation in wettest month", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##bio18
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "bio18") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "bio18")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "bio18")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("bio18.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="bio18: Precipitation during warmest quarter", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##bio2
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "bio2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "bio2")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "bio2")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))


pdf("bio2.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="bio2: Mean diurnal range", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##MEM1
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "MEM1") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "MEM1")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "MEM1")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("MEM1.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM1", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM2
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "MEM2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "MEM2")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "MEM2")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("MEM2.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="MEM2", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM3
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "MEM3") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "MEM3")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "MEM3")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("MEM3.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM3", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist.km
CU.Fst <- cumimp(gf.Fst.model.SEtemp, "dist.km") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.NEUTRAL.model.SEtemp, "dist.km")
CU.ENV <- cumimp(gf.ENVcandidates.model.SEtemp, "dist.km")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))
isub.ENV <- seq(1, length(CU.ENV$x), len = pmin(500, length(CU.ENV$x)))

pdf("dist.CumImp.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="Distance (km)", ylim=c(0,0.12), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=3, lwd=1.5, xaxt='n', yaxt='n')
par(new=T)
plot(CU.ENV$x[isub.ENV], CU.ENV$y[isub.ENV], type = "s", ylab = "", xlab="", ylim=c(0,0.12), lty=2, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##Create the combined plot in AdobeIllustrator
```

![alt_txt][Fig2]
[Fig2]:https://cloud.githubusercontent.com/assets/12142475/22886423/13fc2d1c-f1fe-11e6-9f5b-2fbad3e0ca54.png


#4 Individual SNPs
Plot of individual SNPs for a specific dataset and environmental variable
```
plot(gf.Fst.model.SEtemp, plot.type="C", imp.vars="bio5", show.species=T, ylim=c(0,0.5))

```
![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png

#5. PCA for each map

Plot the first two PCs for each dataset. Include samples on the PCA, and the first two axes of differentiation. 

```



```




#6. Geographic plots

I'm using the Fitzpatric & Keller 2015 code for best results. 

1. gf model output (generated above)

gf.ENVcandidates.model.SEtemp

gf.Fst.model.SEtemp

gf.NEUTRAL.model.SEtemp

2. env_Trns.SE.complete

The extracted bioclim variables (same as run in the gf model) and all the cell IDs (i.e. locations in the raster stack). 
The data needs to be complete; NA needs to be removed. 

3. predicted maps

create a raster file which is the product of the bioclim data and the predicted allele frequencies

4. Define functions

Define the functions for the mapping (as in Fitzpatrick&Keller). 

5. rast & coordinates

A raster base layer to be mapped on. Here I create a srtm elevation map in grey scale for the base map. 

6. Mapping



###2. env_Trns.SE.complete

extract all the information from the raster layers. Remove unwanted layers. Remove missing data. 
```
library(rgdal)
#install.packages("raster")
library(raster)
climate <- getData('worldclim', var='bio', res=2.5)
climate  ##make sure it's a RasterStack
names(climate)  ##lists bio1-19
plot(climate$bio1) ##This should be of the whole world
climate2 <- crop(climate, extent(9,23,52,70)) ##crop to map extent
climate2
plot(climate2)

env_trns <- extract(climate2,1:ncell(climate2), df=T)   ##extract all the bioclim data and the cell IDs
env_trns.SE <- subset(env_trns, select=c("ID", "bio5", "bio15", "bio18", "bio2", "bio13"))

env_trns.SE.complete <- env_trns.SE[complete.cases(env_trns.SE),]  ##remove missing data
```

###3. Predictor maps
```
#library(gradientForest)
# transform env using gf models, see ?predict.gradientForest

pred.NEUTRAL <- predict(gf.NEUTRAL.model.SEtemp, env_trns.SE[,-1]) ##remove ID column
pred.Fst <- predict(gf.Fst.model.SEtemp, env_trns.SE[,-1])
pred.ENV <- predict(gf.ENV.model.SEtemp, env_trns.SE[,-1])

pred.NEUTRAL.complete <- pred.NEUTRAL[complete.cases(pred.NEUTRAL),]  ##remove all missing data
pred.Fst.complete <- pred.Fst[complete.cases(pred.Fst),]
pred.ENV.complete <- pred.ENV[complete.cases(pred.ENV),]
```

###4. Define functions

```
# Mapping spatial genetic variation --------------------------------------------
###### functions to support mapping #####
# builds RGB raster from transformed environment
# pred* = dataframe of transformed variables from gf or gdm model
# rast = a raster mask to which RGB values are to be mapped  ##This is just a single layer!
# cellNums = cell IDs to which RGB values should be assigned  ##all this can be extracted from raster stack
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
    
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}
```

###5. Rast & sampling coordinates

I just used one of the raster layers from climate2. 

```
raster <- climate2$bio2
```

coordinates
```
setwd("/Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/GradientForest")
SE.coords <- read.table("SE.coords", header=T)
coordinates(SE.coords) <- c("Long", "Lat")
names(SE.coords)
plot(SE.coords)
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  ##make sure the projection is the same as for the raster files
proj4string(SE.coords) <- crs.geo

pch.SE <- c(0,0,0,1,1,1,2,2,2,5,5,5,6,6,6) ##make pch.SE symbols different for each genetic cluster as in DAPC
```



I tried the following first to get a grey-scale base map. But it doesn't work because I didn't get the same number of columns. The resolution is also way too high, so it takes really long to load. The first solution works fine. 

Create a raster mask to map onto. Here I'm creating an srtm elevation map. Since the area is so large, I need to download several tiles
and mosaic them together. 
Then I want the map to be in grey scale. 

Download all the tiles from here: http://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm
```
#extract .zip files in Downloads. Navigate to each folder and import the list of raster files. 

setwd("/Users/alexjvr/Downloads/P32")  ##go to the directory
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.P32 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.P32) <- NULL
srtm.SE.P32$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos <- do.call(mosaic, srtm.SE.P32)   ##mosaic these together

srtm.SEall <- mosaic(srtm123, mos, mos2, mos3, mos4, mos5, mos6, fun=mean)  ##mosaic all the individual rasters together. 
plot(srtm.SEall)  ##check visually
srtm.SEall  ##check extent
srtm.SEall <- mosaic(srtm.SEall, mos7, mos8, mos9, mos10, mos11, mos12, mos13, fun=mean) ##I had to add more tiles

srtm.SEall.crop <- crop(srtm.SEall, extent(9,23,52,70)) ##crop to the same extent as the bioclim data


```

###6. Mapping

```
# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written


# map continuous variation - NEUTRAL SNPs
NEUTRAL.RGBmap <- pcaToRaster(pred.NEUTRAL.complete, rast, env_trns.SE.complete$ID)
plotRGB(NEUTRAL.RGBmap)
plot(SE.coords, pch=20, cex=1, add=T)
writeRaster(NEUTRAL.RGBmap, "/.../NEUTRAL.RGBmap_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - Fst SNPs
Fst.RGBmap <- pcaToRaster(pred.Fst.complete, rast, env_trns.SE.complete$ID)
plotRGB(Fst.RGBmap)
plot(SE.coords, pch=20, cex=1, add=T)
writeRaster(Fst.RGBmap, "/.../Fst.RGBmap_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - ENV SNPs
ENV.RGBmap <- pcaToRaster(pred.ENV.complete, rast, env_trns.SE.complete$ID)
plotRGB(ENV.RGBmap)
plot(SE.coords, pch=20, cex=1, add=T)
writeRaster(ENV.RGBmap, "/.../ENV.RGBmap_map.tif", format="GTiff", overwrite=TRUE)
```

####Differences between maps
```

# Difference between maps (NEUTRAL and Fst) 
diffNEUTRAL.Fst <- RGBdiffMap(pred.NEUTRAL.complete, pred.Fst.complete, rast, env_trns.SE.complete$ID)
plot(diffNEUTRAL.Fst[[2]])
plot(SE.coords, pch=20, cex=1, add=T)
writeRaster(diffNEUTRAL.Fst[[2]], "/.../diffNEUTRAL.Fst.tif", format="GTiff", overwrite=TRUE)


# Difference between maps (NEUTRAL and ENV) 
diffNEUTRAL.ENV <- RGBdiffMap(pred.NEUTRAL.complete, pred.ENV.complete, rast, env_trns.SE.complete$ID)
plot(diffNEUTRAL.ENV[[2]])
writeRaster(diffNEUTRAL.ENV[[2]], "/.../diffNEUTRAL.ENV.tif", format="GTiff", overwrite=TRUE)


# Difference between maps (Fst and ENV) 
diffFst.ENV <- RGBdiffMap(pred.Fst.complete, pred.ENV.complete, rast, env_trns.SE.complete$ID)
plot(diffFst.ENV[[2]])
writeRaster(diffFst.ENV[[2]], "/.../diffFst.ENV.tif", format="GTiff", overwrite=TRUE)

```



###Code for plot.gradientForest

I couldn't find the code for the gradient forest plots, so I couldn't figure out what was actually being plotted! 

In R getAnywhere() will find the code for any function. 


getAnywhere(plot.gradientForest)
```
function (x, plot.type = c("Overall.Importance", "Split.Density", 
    "Cumulative.Importance", "Performance")[1], par.args = NULL, 
    plot.args = NULL, ...) 
{
    if (!inherits(x, "gradientForest")) 
        stop(paste("'x' must be a gradientForest object"))
    plot.options <- c("Overall.Importance", "Split.Density", 
        "Cumulative.Importance", "Performance")
    if (is.na(plot.option <- pmatch(plot.type, plot.options))) 
        stop(paste("Unmatched plot.type \"", plot.type, "\". Expecting one of \"Overall.Importance\", \"Split.Density\", \"Cumulative.Importance\" or \"Performance\"", 
            sep = ""))
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    amend.args <- function(default.args, new.args) {
        for (arg in intersect(names(default.args), names(new.args))) default.args[[arg]] <- new.args[[arg]]
        extra <- new.args[is.na(match(names(new.args), names(default.args)))]
        c(default.args, extra)
    }
    if (plot.options[plot.option] == "Overall.Importance") {
        plot.args.def <- amend.args(list(cex.axis = 0.7, cex.names = 0.7, 
            las = 1, horiz = TRUE), plot.args)
        plot.args.def <- amend.args(plot.args.def, list(...))
        par.args.def <- amend.args(list(mfrow = c(1, 2), mar = c(4, 
            6, 2, 1)), par.args)
        par(par.args.def)
        do.call("overall.importance.plot", c(list(obj = quote(x)), 
            plot.args.def))
    }
    if (plot.options[plot.option] == "Split.Density") {
        plot.args.def <- amend.args(list(leg.posn = "topright", 
            bin = F, nbin = 101, leg.panel = 1, barwidth = 1, 
            cex.legend = 0.8, line.ylab = 1.5), plot.args)
        plot.args.def <- amend.args(plot.args.def, list(...))
        par.args.def <- amend.args(list(mar = c(4.5, 1.5, 0.5, 
            4.5), omi = c(0.1, 0.25, 0.1, 0.1)), par.args)
        par(par.args.def)
        do.call("Split.density.plot.method2", c(list(obj = quote(x)), 
            plot.args.def))
    }
    if (plot.options[plot.option] == "Cumulative.Importance") {
        plot.args.def <- amend.args(list(leg.posn = "topleft", 
            legend = TRUE, common.scale = F, line.ylab = 1, cex.legend = 0.75, 
            show.species = TRUE, show.overall = TRUE, leg.nspecies = 10), 
            plot.args)
        plot.args.def <- amend.args(plot.args.def, list(...))
        par.args.def <- amend.args(list(mar = c(0, 2.1, 1.1, 
            0), omi = c(0.75, 0.75, 0.1, 0.1)), par.args)
        par(par.args.def)
        do.call("species.cumulative.plot", c(list(obj = quote(x)), 
            plot.args.def))
    }
    if (plot.options[plot.option] == "Performance") {
        plot.args.def <- amend.args(list(horizontal = FALSE, 
            show.names = FALSE, las = 2, cex.axis = 0.7, cex.labels = 0.7, 
            line = 2), plot.args)
        plot.args.def <- amend.args(plot.args.def, list(...))
        par.args.def <- amend.args(list(mfrow = c(1, 1)), par.args)
        par(par.args.def)
        do.call("performance.plot", c(list(obj = quote(x)), plot.args.def))
    }
    invisible()
}


```




getAnywhere(species.cumulative.plot)

```
species.cumulative.plot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", 
    leg.nspecies = 10, legend = TRUE, mfrow = rev(n2mfrow(length(imp.vars) * 
        (show.species + show.overall))), show.species = TRUE, 
    show.overall = TRUE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 
        0.75, 0.1, 0.1), common.scale = F, line.ylab = 1, cex.legend = 0.75, 
    ...) 
{
    if (is.null(imp.vars)) 
        imp.vars <- imp.var.names <- names(importance(obj))[1:2]
    par(mfrow = mfrow)
    cols <- rainbow(length(levels(obj$res.u$spec)))
    names(cols) <- levels(obj$res.u$spec)
    xaxt <- if (show.overall) 
        "n"
    else "s"
    if (show.species) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX, "Species")
            xlim <- range(sapply(CU, "[[", "x"))
            ylim <- c(0,1.2)   ###I'm editing this so that I can plot on one plot
            plot(xlim, ylim, type = "n", xlab = if (show.overall) 
                ""
            else imp.vars.names[imp.vars == varX], ylab = "", 
                xaxt = xaxt, ...)
            for (species in names(CU)) {
                isub <- seq(1, length(CU[[species]]$x), len = pmin(500, 
                  length(CU[[species]]$x)))
                lines(CU[[species]]$x[isub], CU[[species]]$y[isub], 
                  type = "s", col = cols[species])
            }
            no.species <- length(names(cols))
            imp.sp <- sapply(CU, function(cu) max(cu$y))
            best <- order(-imp.sp)[1:min(leg.nspecies, length(imp.sp))]
            if (legend) 
                legend(x = leg.posn, legend = names(cols)[best], 
                  pch = rep(1, no.species)[best], col = cols[best], 
                  bty = "n", cex = cex.legend, pt.lwd = 2)
        }
    }
    if (show.overall) {
        for (varX in imp.vars) {
            CU <- cumimp(obj, varX)
            ymax <- 1.2 ##max(CU$y)
            if (varX == imp.vars[1]) 
                ymax1 <- 1.2  ##change this from ymax 
            isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
            plot(CU$x[isub], CU$y[isub], type = "s", ylab = "", 
                xlab = imp.vars.names[imp.vars == varX], ylim = c(0, 
                  if (common.scale) ymax1 else ymax), ...)
        }
    }
    mtext("Cumulative importance", side = 2, line = line.ylab, 
        outer = TRUE)
}


```
