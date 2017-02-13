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



###Table 1
```


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


Plot of individual SNPs for a specific dataset and environmental variable
```
plot(gf.Fst.model.SEtemp, plot.type="C", imp.vars="bio5", show.species=T, ylim=c(0,0.5))

```
![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png




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
