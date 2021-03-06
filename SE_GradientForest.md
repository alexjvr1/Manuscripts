# Gradient Forest Analysis

Here I aim to: 

1)determine the environmental and spatial determinants of variation in allele frequency across the Swedish latitudinal gradient. 

2)map the current and future patterns of genomic variation in relation to climate

I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013


## Data

I'm using 3 different datasets: 

1. Fst outliers as determined by BayEnv2 XtX, and RDA (RD1 & 2)

2. Loci associated with environment as determined by LFMM & BayEnv2 

3. "Neutral" loci (all - Data1+2)


## Presenting Results

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


### Table 1

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




### Fig1: R2 heatmap
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



### Fig2: Average allelic turnover across all three datasets. 

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


# 4 Individual SNPs
Plot of individual SNPs for a specific dataset and environmental variable
```
plot(gf.Fst.model.SEtemp, plot.type="C", imp.vars="bio5", show.species=T, ylim=c(0,0.5))

```
![alt_txt][Fig3]

[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png

# 5. PCA for each map

Plot the first two PCs for each dataset. Include samples on the PCA, and the first two axes of differentiation. 

Neutral data
```
pdf(file="Neutral.PCA.pdf")
Trns.grid.NEUTRAL.mask <- cbind(env.trns.mask.SE.complete$ID, 
predict(gf.NEUTRAL.model.SEtemp, env.trns.mask.SE.complete[,vec]))  ##define the Trns.grid
PCs.NEUTRAL <- prcomp(Trns.grid.NEUTRAL.mask[,vec])   ##calculate PCs
a1 <- PCs.NEUTRAL$x[,1]   ##assign colour
a2 <- PCs.NEUTRAL$x[,2]
a3 <- PCs.NEUTRAL$x[,3]

nvs <- dim(PCs.NEUTRAL$rotation)[1]
lv <- length(vec)
vind <- rownames(PCs.NEUTRAL$rotation) %in% vec
scal <- 20
xrng <- range(PCs.NEUTRAL$x[,1], PCs.NEUTRAL$rotation[,1]/scal)*1.1
yrng <- range(PCs.NEUTRAL$x[,2], PCs.NEUTRAL$rotation[,2]/scal)*1.1

plot((PCs.NEUTRAL$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=4, col=rgb(r,g,b,max=255), asp=1)

arrows(rep(0,lv), rep(0,lv), PCs.NEUTRAL$rotation[vec,1]/scal, 
PCs.NEUTRAL$rotation[vec,2]/scal, length=0.1)   ##add arrows for all env variables

text(PCs.NEUTRAL$rotation[vec,1]/scal + jit *sign(PCs.NEUTRAL$rotation[vec,
1]), PCs.NEUTRAL$rotation[vec,2]/scal+jit*sign(PCs.NEUTRAL$rotation[vec,2]), 
labels=vec)   ##label the arrows

##for the sampling locations:

NEUTRAL.Trns.site <- predict(gf.NEUTRAL.model.SEtemp)
PCsites <- predict(PCs.NEUTRAL, NEUTRAL.Trns.site[, vec])

points(PCsites[,1:2], pch=pch.SE, cex=4)

dev.off()
```

Fst.outliers PCA
```
pdf(file="Fstoutliers.PCA.pdf")
Trns.grid.Fst.mask <- cbind(env.trns.mask.SE.complete$ID, 
predict(gf.Fst.model.SEtemp, env.trns.mask.SE.complete[,vec]))  ##define the Trns.grid
PCs.Fst <- prcomp(Trns.grid.Fst.mask[,vec])   ##calculate PCs
a1 <- PCs.Fst$x[,1]   ##assign colour
a2 <- PCs.Fst$x[,2]
a3 <- PCs.Fst$x[,3]

nvs <- dim(PCs.Fst$rotation)[1]
lv <- length(vec)
vind <- rownames(PCs.Fst$rotation) %in% vec
scal <- 20
xrng <- range(PCs.Fst$x[,1], PCs.Fst$rotation[,1]/scal)*1.1
yrng <- range(PCs.Fst$x[,2], PCs.Fst$rotation[,2]/scal)*1.1

plot((PCs.Fst$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=4, col=rgb(r,g,b,max=255), asp=1)

arrows(rep(0,lv), rep(0,lv), PCs.Fst$rotation[vec,1]/scal, 
PCs.Fst$rotation[vec,2]/scal, length=0.1)   ##add arrows for all env variables

text(PCs.Fst$rotation[vec,1]/scal + jit *sign(PCs.Fst$rotation[vec,
1]), PCs.Fst$rotation[vec,2]/scal+jit*sign(PCs.Fst$rotation[vec,2]), 
labels=vec)   ##label the arrows

##for the sampling locations:

Fst.Trns.site <- predict(gf.Fst.model.SEtemp)
PCsites <- predict(PCs.Fst, Fst.Trns.site[, vec])

points(PCsites[,1:2], pch=pch.SE, cex=4)

dev.off()
```

ENV PCA
```
pdf(file="ENV.PCA.pdf")
Trns.grid.ENV.mask <- cbind(env.trns.mask.SE.complete$ID, 
predict(gf.ENVcandidates.model.SEtemp, env.trns.mask.SE.complete[,vec]))  ##define the Trns.grid
PCs.ENV <- prcomp(Trns.grid.ENV.mask[,vec])   ##calculate PCs
a1 <- PCs.ENV$x[,1]   ##assign colour
a2 <- PCs.ENV$x[,2]
a3 <- PCs.ENV$x[,3]

nvs <- dim(PCs.ENV$rotation)[1]
lv <- length(vec)
vind <- rownames(PCs.ENV$rotation) %in% vec
scal <- 20
xrng <- range(PCs.ENV$x[,1], PCs.ENV$rotation[,1]/scal)*1.1
yrng <- range(PCs.ENV$x[,2], PCs.ENV$rotation[,2]/scal)*1.1

plot((PCs.ENV$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=4, col=rgb(r,g,b,max=255), asp=1)

arrows(rep(0,lv), rep(0,lv), PCs.ENV$rotation[vec,1]/scal, 
PCs.ENV$rotation[vec,2]/scal, length=0.1)   ##add arrows for all env variables

text(PCs.ENV$rotation[vec,1]/scal + jit *sign(PCs.ENV$rotation[vec,
1]), PCs.ENV$rotation[vec,2]/scal+jit*sign(PCs.ENV$rotation[vec,2]), 
labels=vec)   ##label the arrows

##for the sampling locations:

ENV.Trns.site <- predict(gf.ENVcandidates.model.SEtemp)
PCsites <- predict(PCs.ENV, ENV.Trns.site[, vec])

#> PCsites
#               PC1          PC2           PC3           PC4           PC5
# [1,]  0.060915792 -0.011263055  0.0074592278 -0.0048529781  1.009606e-03
# [2,]  0.064153852 -0.003130803  0.0056775950 -0.0004598817  1.197192e-03
# [3,]  0.058275119 -0.008496364  0.0048235124 -0.0048495990  1.418894e-03
# [4,]  0.014993741  0.007954388 -0.0097689320 -0.0123684209  7.103178e-04
# [5,]  0.008043463  0.003129808 -0.0104306636 -0.0089215196  2.963193e-04
# [6,]  0.014888719  0.007225827 -0.0095453781 -0.0110190268 -1.169056e-05
# [7,] -0.004472003  0.009989169 -0.0045379793  0.0100325083 -1.676996e-03
# [8,] -0.004264044  0.010920251 -0.0048521015  0.0106849330 -1.979413e-03
# [9,] -0.003696640  0.011779877 -0.0047879752  0.0109177834 -1.755980e-03
#[10,] -0.007437736  0.018814787 -0.0005391422  0.0140316586 -1.868844e-03
#[11,] -0.011178884  0.015202306 -0.0008776371  0.0140705200 -1.579113e-03
#[12,] -0.007528816  0.018291778 -0.0005751241  0.0137010842 -1.945908e-03
#[13,] -0.020064543  0.005070986  0.0049694497 -0.0084781652  1.454736e-03
#[14,] -0.015354907  0.011452598  0.0052569741  0.0016329271  1.612022e-03
#[15,] -0.003796294  0.009527998  0.0237646832 -0.0077324617 -1.655614e-03

points(PCsites[,1:2], pch=pch.SE, cex=4)

dev.off()
```






# 6. Geographic plots

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



### 2. env_Trns.SE.complete

extract all the information from the raster layers. Remove unwanted layers. Remove missing data. 

First I extract all the bioclim data for the extent of the bioclim rasterbrick. 

Second I show a way to create a polygon and then extract data only for this region. 

Lastly I show how to mask a rast with a polygon, and then extract raster data from this. This is what I will use for the analysis. 
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


###From polygon

##Get coordinates of the polygon from here: http://www.birdtheme.org/useful/v3tool.html
##paste them into a file 

SE.polygon.coords <- read.table("SE.polygon.coords", header=F, sep=",") ##read into R
SE.polygon2 <- SE.polygon.coords[,1:2]  ##select the first two columns (long & Lat)
summary(SE.polygon2)  
head(SE.polygon2)
colnames(SE.polygon2) <- c("Long", "Lat")  
head(SE.polygon2)
SE.polygon3 <- Polygon(SE.polygon2) 
plot(SE.polygon3)
SE.polygon3 <- SpatialPolygons(SE.polygon3)

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  #define the projection (make sure this is the same as the raster to be used)

SE.polygon5 <- SpatialPolygons(list(Polygons(list(SE.polygon3), ID="a")), 
proj4string=crs.geo)
plot(SE.polygon5, axes=T)  ##make sure its the right shape and that all the samples are inside the polygon
plotRGB(NEUTRAL.RGBmap)
plot(SE.polygon5, add=T)
plot(SE.coords, pch=20, cex=1, add=T)   

env_trns.SE.polygon <- extract(climate2, SE.polygon5, df=T)  ##extract all variables inside the polygon
env_trns.SE.polygon2 <- subset(env_trns.SE.polygon, select=c("ID", "bio5", "bio15", "bio18", "bio2", "bio13"))  ##select only the important variables
env_trns.SE.polygon2.complete <- env_trns.SE.polygon2[complete.cases(env_trns.SE.polygon2),] #remove all missing data. 
```

DataMask: the method I use in the end: 
```
mask.test <- mask(climate2, SE.polygon5)  ##use the polygon created above to mask the raster stack created above (i.e. already subset).
env.trns.mask <- extract(mask.test, 1:ncell(mask.test), df=T) ##extract the data to a dataframe
env.trns.mask.SE <- subset(env.trns.mask, select=c("ID", "bio5", "bio15", "bio18", "bio2", "bio13"))  ##make sure the correct variables are included
env.trns.mask.SE.complte <- env.trns.mask.SE[complete.cases(env.trns.mask.SE),]  ##only complete cases

###see below how this is used to create new predictor maps. 
```


### 3. Predictor maps
```
#library(gradientForest)
# transform env using gf models, see ?predict.gradientForest

pred.NEUTRAL <- predict(gf.NEUTRAL.model.SEtemp, env_trns.SE[,-1]) ##remove ID column
pred.Fst <- predict(gf.Fst.model.SEtemp, env_trns.SE[,-1])
pred.ENV <- predict(gf.ENVcandidates.model.SEtemp, env_trns.SE[,-1])

pred.NEUTRAL.complete <- pred.NEUTRAL[complete.cases(pred.NEUTRAL),]  ##remove all missing data
pred.Fst.complete <- pred.Fst[complete.cases(pred.Fst),]
pred.ENV.complete <- pred.ENV[complete.cases(pred.ENV),]
```

Predictor maps with raster mask
```
#library(gradientForest)
# transform env using gf models, see ?predict.gradientForest

pred.NEUTRAL.mask <- predict(gf.NEUTRAL.model.SEtemp, env.trns.mask.SE[,-1]) ##remove ID column
pred.Fst.mask <- predict(gf.Fst.model.SEtemp, env.trns.mask.SE[,-1])
pred.ENV.mask <- predict(gf.ENVcandidates.model.SEtemp, env.trns.mask.SE[,-1])

pred.NEUTRAL.mask.complete <- pred.NEUTRAL.mask[complete.cases(pred.NEUTRAL.mask),]  ##remove all missing data
pred.Fst.mask.complete <- pred.Fst.mask[complete.cases(pred.Fst.mask),]
pred.ENV.mask.complete <- pred.ENV.mask[complete.cases(pred.ENV.mask),]
```


### 4. Define functions

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

### 5. Rast & sampling coordinates

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
#I'm using R33, R34, Q32, Q33, Q34, P32, P33, P34, O32, O33, O34, N33

#make mos from R33
setwd("/Users/alexjvr/Downloads/R33")  ##go to the directory
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.R33 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.R33) <- NULL
srtm.SE.R33$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos <- do.call(mosaic, srtm.SE.R33)   ##mosaic these together

##mos1 from R34
setwd("/Users/alexjvr/Downloads/R34")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.R34 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.R34) <- NULL
srtm.SE.R34$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos1 <- do.call(mosaic, srtm.SE.R34)   ##mosaic these together

##mos2 from Q33
setwd("/Users/alexjvr/Downloads/Q33")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.Q33 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.Q33) <- NULL
srtm.SE.Q33$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos2 <- do.call(mosaic, srtm.SE.Q33)   ##mosaic these together

##mos3 from Q34
setwd("/Users/alexjvr/Downloads/Q34")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.Q34 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.Q34) <- NULL
srtm.SE.Q34$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos3 <- do.call(mosaic, srtm.SE.Q34)   ##mosaic these together

##mos4 from P33
setwd("/Users/alexjvr/Downloads/P33")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.P33 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.P33) <- NULL
srtm.SE.P33$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos4 <- do.call(mosaic, srtm.SE.P33)   ##mosaic these together

##mos5 from P34
setwd("/Users/alexjvr/Downloads/P34")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.P34 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.P34) <- NULL
srtm.SE.P34$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos5 <- do.call(mosaic, srtm.SE.P34)   ##mosaic these together


##mos6 from O33
setwd("/Users/alexjvr/Downloads/O33")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.O33 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.O33) <- NULL
srtm.SE.O33$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos6 <- do.call(mosaic, srtm.SE.O33)   ##mosaic these together

##mos7 from O34
setwd("/Users/alexjvr/Downloads/O34")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.O34 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.O34) <- NULL
srtm.SE.O34$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos7 <- do.call(mosaic, srtm.SE.O34)   ##mosaic these together


##mos8 from N33
setwd("/Users/alexjvr/Downloads/N33")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.N33 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.N33) <- NULL
srtm.SE.N33$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos8 <- do.call(mosaic, srtm.SE.N33)   ##mosaic these together

##mos9 from O32
setwd("/Users/alexjvr/Downloads/O32")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.O32 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.O32) <- NULL
srtm.SE.O32$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos9 <- do.call(mosaic, srtm.SE.O32)   ##mosaic these together

##mos10 from P32
setwd("/Users/alexjvr/Downloads/P32")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.P32 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.P32) <- NULL
srtm.SE.P32$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos10 <- do.call(mosaic, srtm.SE.P32)   ##mosaic these together

##mos11 from Q32
setwd("/Users/alexjvr/Downloads/Q32")  ##go to the directory
SE2 <- NULL
SE2 <- list.files(pattern="hgt")  ##list all the files of a certain type
srtm.SE.Q32 <- sapply(SE2, function(filename){
r=raster(filename)
})   ##function to import all these as a raster list into R

names(srtm.SE.Q32) <- NULL
srtm.SE.Q32$fun <- mean  ##function with which they will be joined together. For only two layers: raster.mos <- mosaic(ras1, ras2, fun=mean)
mos11 <- do.call(mosaic, srtm.SE.Q32)   ##mosaic these together

##########
###Mosaic all the rasters together
srtm.SEall <- mosaic(mos, mos1, mos2, mos3, mos4, mos5, mos6, mos7, mos8, mos9, mos10, mos11, fun=mean)  ##mosaic all the individual rasters together. 
plot(srtm.SEall)  ##check visually
srtm.SEall  ##check extent
##srtm.SEall <- mosaic(srtm.SEall, mos7, mos8, mos9, mos10, mos11, mos12, mos13, fun=mean) ##I had to add more tiles

srtm.SEall.crop <- crop(srtm.SEall, extent(9,23,52,70)) ##crop to the same extent as the bioclim data


```

### 6. Mapping

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

MAPS using the raster mask

```
# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written

setwd("/Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/GradientForest")

# map continuous variation - NEUTRAL SNPs
NEUTRAL.RGBmap.mask <- pcaToRaster(pred.NEUTRAL.mask.complete, rast, env.trns.mask.SE.complete$ID)
pdf(file="NEUTRAL.RGBmap.mask.pdf")
plotRGB(NEUTRAL.RGBmap.mask)
plot(SE.coords, pch=pch.SE, cex=1, add=T)
dev.off()
writeRaster(NEUTRAL.RGBmap.mask, "NEUTRAL.RGBmap.mask_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - Fst SNPs
Fst.RGBmap.mask <- pcaToRaster(pred.Fst.mask.complete, rast, env.trns.mask.SE.complete$ID)
pdf(file="FST.RGBmap.mask.pdf")
plotRGB(Fst.RGBmap.mask)
plot(SE.coords, pch=pch.SE, cex=1, add=T)
dev.off()
writeRaster(Fst.RGBmap.mask, "Fst.RGBmap.mask_map.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - ENV SNPs
ENV.RGBmap.mask <- pcaToRaster(pred.ENV.mask.complete, rast, env.trns.mask.SE.complete$ID)
pdf(file="ENV.RGBmap.mask.pdf")
plotRGB(ENV.RGBmap.mask)
plot(SE.coords, pch=pch.SE, cex=1, add=T)
dev.off()
writeRaster(ENV.RGBmap.mask, "ENV.RGBmap.mask_map.tif", format="GTiff", overwrite=TRUE)
```

Left to right: Neutral, Fst, ENV

![alt_txt][GF.maps.PCA]

[GF.maps.PCA]:https://cloud.githubusercontent.com/assets/12142475/23069483/eb4e0c82-f527-11e6-98c7-3da7fb9eeecc.png




#### Differences between maps


Neutral vs Fst
```
pdf(file="diffNEUTRAL.FST.map.pdf")

diff.NEUTRAL.FST.mask <- RGBdiffMap(pred.NEUTRAL.mask.complete, pred.Fst.mask.complete, mask.test$bio2, env.trns.mask.SE.complete$ID)   ##I have to make sure that the rast.mask is used - i.e. only some cells. And same with the env.trns file. 

##Now I have to normalise the raster from 0-1.

r <- diff.NEUTRAL.FST.mask[[2]]  ##get the raster layer on its own
r.min =cellStats(r, "min")
r.max = cellStats(r, "max")  ##determine the min and max values
r.scale <- ((r-r.min)/(r.max-r.min))  ##rescale (0-1)

##Define a colour palette that will best show the difference
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


##plot
plot(climate2$bio2, col="grey30", legend=F) ##Then I plot a base map in grey (make sure its the same scale as the raster file
plot(r.scale, col=mypalette, add=T) 
plot(SE.coords, pch=pch.SE, cex=1, add=T)

dev.off()
```


Neutral vs ENV
```
pdf(file="diffNEUTRAL.ENV.map.pdf")
diff.NEUTRAL.ENV.mask <- RGBdiffMap(pred.NEUTRAL.mask.complete, pred.ENV.mask.complete, mask.test$bio2, env.trns.mask.SE.complete$ID)   ##I have to make sure that the rast.mask is used - i.e. only some cells. And same with the env.trns file. 

##Now I have to normalise the raster from 0-1.

r <- diff.NEUTRAL.ENV.mask[[2]]  ##get the raster layer on its own
r.min =cellStats(r, "min")
r.max = cellStats(r, "max")  ##determine the min and max values
r.scale <- ((r-r.min)/(r.max-r.min))  ##rescale (0-1)

##Define a colour palette that will best show the difference
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


##plot
plot(climate2$bio2, col="grey30", legend=F) ##Then I plot a base map in grey (make sure its the same scale as the raster file
plot(r.scale, col=mypalette, add=T) 
plot(SE.coords, pch=pch.SE, cex=1, add=T)

dev.off()
```

ENV vs Fst

```
pdf(file="diffENV.Fst.map.pdf")
diff.ENV.Fst.mask <- RGBdiffMap(pred.ENV.mask.complete, pred.Fst.mask.complete, mask.test$bio2, env.trns.mask.SE.complete$ID)
##I have to make sure that the rast.mask is used - i.e. only some cells. And same with the env.trns file. 

##Now I have to normalise the raster from 0-1.

r <- diff.ENV.Fst.mask[[2]]
r.min =cellStats(r, "min")
r.max = cellStats(r, "max")  ##determine the min and max values
r.scale <- ((r-r.min)/(r.max-r.min))  ##rescale (0-1)

##Define a colour palette that will best show the difference

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)

#plot

plot(climate2$bio2, col="grey30", legend=F) ##Then I plot a base map in grey (make sure its the same scale as the raster file
plot(r.scale, col=mypalette, add=T) 
plot(SE.coords, pch=pch.SE, cex=1, add=T)

dev.off()
```



Old script
```
# Difference between maps (NEUTRAL and Fst) 
diffNEUTRAL.Fst <- RGBdiffMap(pred.NEUTRAL.complete, pred.Fst.complete, rast, env.trns.SE.complete$ID)
plot(diffNEUTRAL.Fst[[2]])
plot(SE.coords, pch=20, cex=1, add=T)
writeRaster(diffNEUTRAL.Fst[[2]], "/.../diffNEUTRAL.Fst.tif", format="GTiff", overwrite=TRUE)


# Difference between maps (NEUTRAL and ENV) 
diffNEUTRAL.ENV <- RGBdiffMap(pred.NEUTRAL.complete, pred.ENV.complete, rast, env.trns.SE.complete$ID)
plot(diffNEUTRAL.ENV[[2]])
writeRaster(diffNEUTRAL.ENV[[2]], "/.../diffNEUTRAL.ENV.tif", format="GTiff", overwrite=TRUE)


# Difference between maps (Fst and ENV) 
diffFst.ENV <- RGBdiffMap(pred.Fst.complete, pred.ENV.complete, rast, env.trns.SE.complete$ID)
plot(diffFst.ENV[[2]])
writeRaster(diffFst.ENV[[2]], "/.../diffFst.ENV.tif", format="GTiff", overwrite=TRUE)

```




### Code for plot.gradientForest

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


##Figure: heatmap of change in allele frequency

Heatmap of the 3 datasets showing the change in allele frequency across the gradient. 
I've generated the files containing the MAF of each dataset before. 


```
getwd()
#/Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/GradientForest

library(gplots)

SE.MAF3.Fst.matrix <- data.matrix(SE.MAF3.Fst.matrix)
#SE.MAF3.Fst.matrix <- SE.MAF3.Fst.matrix[,order(SE.MAF3.Fst.matrix[1,])]  #order by row 1 to exaggerate the effects

SE.MAF3.Env.matrix <- data.matrix(SE.MAF3.Env)
SE.MAF3.Neutral.matrix <- data.matrix(SE.MAF3.Neutral)
SE.MAF3.Neutral.matrix <- SE.MAF3.Neutral.matrix[,order(SE.MAF3.Neutral.matrix[1,])]

pdf(file="heatmap.ENV.MAF.pdf")
heatmap.ENV <- heatmap.2(SE.MAF3.Env.matrix, col=shadesOfGrey, trace="none", scale="column", labCol=F, Rowv=F)
dev.off()

pdf(file="heatmap.Ref.MAF.pdf")
heatmap.Ref <- heatmap.2(SE.MAF3.Neutral.matrix, col=shadesOfGrey, trace="none", scale="column", labCol=F, Rowv=F)
dev.off()

pdf(file="heatmap.Fst.MAF.pdf")
heatmap.Fst <- heatmap.2(SE.MAF3.Fst.matrix2, col=shadesOfGrey, trace="none", scale="column", labCol=F, Rowv=F)
dev.off()
```

Write to pdf and import into Adobe illustrator to edit the images. 




############
#############
Code from one of the notebooks

Trns.grid.NEUTRAL.mask <- cbind(env.trns.mask.SE.complete$ID, 
predict(gf.NEUTRAL.model.SEtemp, env.trns.mask.SE.complete[,vec]))  ##define the Trns.grid
PCs.NEUTRAL <- prcomp(Trns.grid.NEUTRAL.mask[,vec])   ##calculate PCs
a1 <- PCs.NEUTRAL$x[,1]   ##assign colour
a2 <- PCs.NEUTRAL$x[,2]
a3 <- PCs.NEUTRAL$x[,3]

nvs <- dim(PCs.NEUTRAL$rotation)[1]
lv <- length(vec)
vind <- rownames(PCs.NEUTRAL$rotation) %in% vec
scal <- 20
xrng <- range(PCs.NEUTRAL$x[,1], PCs.NEUTRAL$rotation[,1]/scal)*1.1
yrng <- range(PCs.NEUTRAL$x[,2], PCs.NEUTRAL$rotation[,2]/scal)*1.1

plot((PCs.NEUTRAL$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=4, col=rgb(r,g,b,max=255), asp=1)

arrows(rep(0,lv), rep(0,lv), PCs.NEUTRAL$rotation[vec,1]/scal, 
PCs.NEUTRAL$rotation[vec,2]/scal, length=0.1)   ##add arrows for all env variables

text(PCs.NEUTRAL$rotation[vec,1]/scal + jit *sign(PCs.NEUTRAL$rotation[vec,
1]), PCs.NEUTRAL$rotation[vec,2]/scal+jit*sign(PCs.NEUTRAL$rotation[vec,2]), 
labels=vec)   ##label the arrows

##for the sampling locations:

NEUTRAL.Trns.site <- predict(gf.NEUTRAL.model.SEtemp)
PCsites <- predict(PCs.NEUTRAL, NEUTRAL.Trns.site[, vec])

> PCsites
               PC1          PC2           PC3           PC4           PC5
 [1,]  0.060915792 -0.011263055  0.0074592278 -0.0048529781  1.009606e-03
 [2,]  0.064153852 -0.003130803  0.0056775950 -0.0004598817  1.197192e-03
 [3,]  0.058275119 -0.008496364  0.0048235124 -0.0048495990  1.418894e-03
 [4,]  0.014993741  0.007954388 -0.0097689320 -0.0123684209  7.103178e-04
 [5,]  0.008043463  0.003129808 -0.0104306636 -0.0089215196  2.963193e-04
 [6,]  0.014888719  0.007225827 -0.0095453781 -0.0110190268 -1.169056e-05
 [7,] -0.004472003  0.009989169 -0.0045379793  0.0100325083 -1.676996e-03
 [8,] -0.004264044  0.010920251 -0.0048521015  0.0106849330 -1.979413e-03
 [9,] -0.003696640  0.011779877 -0.0047879752  0.0109177834 -1.755980e-03
[10,] -0.007437736  0.018814787 -0.0005391422  0.0140316586 -1.868844e-03
[11,] -0.011178884  0.015202306 -0.0008776371  0.0140705200 -1.579113e-03
[12,] -0.007528816  0.018291778 -0.0005751241  0.0137010842 -1.945908e-03
[13,] -0.020064543  0.005070986  0.0049694497 -0.0084781652  1.454736e-03
[14,] -0.015354907  0.011452598  0.0052569741  0.0016329271  1.612022e-03
[15,] -0.003796294  0.009527998  0.0237646832 -0.0077324617 -1.655614e-03

points(PCsites[,1:2], pch=pch.SE)







diff.NEUTRAL.FST.mask <- RGBdiffMap(pred.NEUTRAL.mask.complete, pred.Fst.mask.complete, mask.test$bio2, env.trns.mask.SE.complete$ID)   ##I have to make sure that the rast.mask is used - i.e. only some cells. And same with the env.trns file. 


##Now I have to normalise the raster from 0-1.

r <- diff.NEUTRAL.FST.mask[[2]]  ##get the raster layer on its own
r.min =cellStats(r, "min")
r.max = cellStats(r, "max")  ##determine the min and max values
r.scale <- ((r-r.min)/(r.max-r.min))  ##rescale (0-1)

##Define a colour palette that will best show the difference
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


##plot
plot(climate2$bio2, col="grey30, legend=F) ##Then I plot a base map in grey (make sure its the same scale as the raster file
plot(r.scale, col=mypalette, add=T) 
plot(SE.coords, pch=pch.SE, cex=1, add=T)





