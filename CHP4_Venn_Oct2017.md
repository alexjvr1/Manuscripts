# Venn Diagrams: CHP4 

Calculating the overlap in loci identified for all the different datasets. 


# MAC filtered data

pwd
/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/VENN.Chp4



I'm showing: 

Venn Diagram within analysis for 
  
  
  1.1. LFMM
  
  1.2. bayEnv2


overlap between LFMM, Bayenv2, PCadapt, XtX, and Bayescan within each region (i.e. 6 plots)

  2.1. temp
  
  2.2. season
  

3. Overlap between all candidate loci identified by 2+ methods for CHN, CHS, CZ


4. Overlap between all candidate loci identified by 2+ methods for CHS, CHS.TI, CHS.VS



### 1.1 Venn Diagram for loci within LFMM

These have already been drawn. 

See: https://github.com/alexjvr1/Manuscripts/blob/aaea5c1b483068f975923b87b7b73b184f72a180/Chp4_ComparisonSEvsCH.md



Venn.CHN.CHP4.n2.LFMMonly_20171009.pdf

Venn.CHS.TI.n2.CHP4.LFMMonly_20171009.pdf

Venn.CHS.VS.n2.CHP4.LFMMonly_20171009.pdf

Venn.CHS.n2.CHP4.LFMMonly_20171009.pdf

Venn.CZ.n2.CHP4.LFMMonly_20171009.pdf

Venn.SE.CHP4.n2.LFMMonly_20171009.pdf



### 1.2 Venn Diagram for loci within BayEnv2


##### SE

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

SE.temp.bayenv.candidates.names
SE.season.days.bayenv.candidates.names

d1 <- length(SE.temp.bayenv.candidates.names)
d2 <- length(SE.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(SE.temp.bayenv.candidates.names, SE.season.days.bayenv.candidates.names)))


pdf("SE.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("SE.mean.temp", "SE.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```




##### CHN

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHN.temp.bayenv.candidates.names
CHN.season.days.bayenv.candidates.names

d1 <- length(CHN.temp.bayenv.candidates.names)
d2 <- length(CHN.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(CHN.temp.bayenv.candidates.names, CHN.season.days.bayenv.candidates.names)))


library(VennDiagram)
pdf("CHN.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("CHN.mean.temp", "CHN.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```



##### CHS

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.temp.bayenv.candidates.names
CHS.season.days.bayenv.candidates.names

d1 <- length(CHS.temp.bayenv.candidates.names)
d2 <- length(CHS.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(CHS.temp.bayenv.candidates.names, CHS.season.days.bayenv.candidates.names)))


pdf("CHS.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("CHS.mean.temp", "CHS.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```



##### CZ

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CZ.temp.bayenv.candidates.names
CZ.season.days.bayenv.candidates.names

d1 <- length(CZ.temp.bayenv.candidates.names)
d2 <- length(CZ.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(CZ.temp.bayenv.candidates.names, CZ.season.days.bayenv.candidates.names)))


pdf("CZ.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("CZ.mean.temp", "CZ.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```



##### CHS.VS

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.VS.temp.bayenv.candidates.names
CHS.VS.season.days.bayenv.candidates.names

d1 <- length(CHS.VS.temp.bayenv.candidates.names)
d2 <- length(CHS.VS.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(CHS.VS.temp.bayenv.candidates.names, CHS.VS.season.days.bayenv.candidates.names)))


pdf("CHS.VS.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("CHS.VS.mean.temp", "CHS.VS.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```


##### CHS.TI

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.TI.temp.bayenv.candidates.names
CHS.TI.season.days.bayenv.candidates.names

d1 <- length(CHS.TI.temp.bayenv.candidates.names)
d2 <- length(CHS.TI.season.days.bayenv.candidates.names)

d12 <- length(Reduce(intersect, list(CHS.TI.temp.bayenv.candidates.names, CHS.TI.season.days.bayenv.candidates.names)))


pdf("CHS.TI.bayenv2.20171009.CHP4.VENN.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("CHS.TI.mean.temp", "CHS.TI.season"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()


```


### 2.1 overlap between LFMM, Bayenv2, PCadapt, XtX, and bayescan within each region



/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/VENN.Chp4/


##### CHN

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("CHN.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHN.pcadapt.outliers")
temp.bayenv.outliers <- read.table("CHN.CHP4.temp.outlier.names")
XtX.outliers <- read.table("CHN.XtX.100outliers")
bayescan.outliers <- read.table("CHN.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHN.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHN.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHN.temp.duplicated.outliers <- CHN.temp.alloutliers[duplicated(CHN.temp.alloutliers),]  ##select only loci occurring more than once (here 94 of 816)

write.table(CHN.temp.duplicated.outliers, "CHN.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("CHN.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHN.pcadapt.outliers")
season.bayenv.outliers <- read.table("CHN.CHP4.season.outlier.names")
XtX.outliers <- read.table("CHN.XtX.100outliers")
bayescan.outliers <- read.table("CHN.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHN.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHN.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHN.season.duplicated.outliers <- CHN.season.alloutliers[duplicated(CHN.season.alloutliers),]  ##select only loci occurring more than once (here 124 of 971)

write.table(CHN.season.duplicated.outliers, "CHN.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```


![alt_txt][CHN.Venn.allmethods]

[CHN.Venn.allmethods]:https://user-images.githubusercontent.com/12142475/31176663-8c5b0aca-a90b-11e7-9f7d-25d1fbdac1d7.png


##### CHS

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("CHS.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.pcadapt.outliers")
temp.bayenv.outliers <- read.table("CHS.CHP4.temp.outlier.names")
XtX.outliers <- read.table("CHS.XtX.100outliers")
bayescan.outliers <- read.table("CHS.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.temp.duplicated.outliers <- CHS.temp.alloutliers[duplicated(CHS.temp.alloutliers),]  ##select only loci occurring more than once (here 153 of 1223)

length(CHS.temp.alloutliers$loci)
length(CHS.temp.duplicated.outliers)

write.table(CHS.temp.duplicated.outliers, "CHS.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("CHS.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.pcadapt.outliers")
season.bayenv.outliers <- read.table("CHS.CHP4.season.outlier.names")
XtX.outliers <- read.table("CHS.XtX.100outliers")
bayescan.outliers <- read.table("CHS.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.season.duplicated.outliers <- CHS.season.alloutliers[duplicated(CHS.season.alloutliers),]  ##select only loci occurring more than once (here 140 of 1200)

length(CHS.season.alloutliers$loci)
length(CHS.season.duplicated.outliers)

write.table(CHS.season.duplicated.outliers, "CHS.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




##### CZ

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("CZ.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("CZ.pcadapt.outliers")
temp.bayenv.outliers <- read.table("CZ.CHP4.temp.outlier.names")
XtX.outliers <- read.table("CZ.XtX.100outliers")
bayescan.outliers <- read.table("CZ.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CZ.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CZ.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CZ.temp.duplicated.outliers <- CZ.temp.alloutliers[duplicated(CZ.temp.alloutliers),]  ##select only loci occurring more than once (here 153 of 1223)

length(CZ.temp.alloutliers$loci)
length(CZ.temp.duplicated.outliers)

write.table(CZ.temp.duplicated.outliers, "CZ.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("CZ.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("CZ.pcadapt.outliers")
season.bayenv.outliers <- read.table("CZ.CHP4.season.outlier.names")
XtX.outliers <- read.table("CZ.XtX.100outliers")
bayescan.outliers <- read.table("CZ.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CZ.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CZ.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CZ.season.duplicated.outliers <- CZ.season.alloutliers[duplicated(CZ.season.alloutliers),]  ##select only loci occurring more than once (here 130 of 1267)

length(CZ.season.alloutliers$loci)
length(CZ.season.duplicated.outliers)

write.table(CZ.season.duplicated.outliers, "CZ.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```


![alt_txt][CZ.Venn.allmethods]

[CZ.Venn.allmethods]:https://user-images.githubusercontent.com/12142475/31231271-4a451eea-a9de-11e7-8a4a-7c17db951078.png



##### CHS.VS

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("CHS.VS.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.VS.pcadapt.outliers")
temp.bayenv.outliers <- read.table("CHS.VS.CHP4.temp.outlier.names")
XtX.outliers <- read.table("CHS.VS.XtX.100outliers")
bayescan.outliers <- read.table("CHS.VS.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.VS.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.VS.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.VS.temp.duplicated.outliers <- CHS.VS.temp.alloutliers[duplicated(CHS.VS.temp.alloutliers),]  ##select only loci occurring more than once (here 187 of 1342)

length(CHS.VS.temp.alloutliers$loci)
length(CHS.VS.temp.duplicated.outliers)

write.table(CHS.VS.temp.duplicated.outliers, "CHS.VS.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("CHS.VS.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.VS.pcadapt.outliers")
season.bayenv.outliers <- read.table("CHS.VS.CHP4.season.outlier.names")
XtX.outliers <- read.table("CHS.VS.XtX.100outliers")
bayescan.outliers <- read.table("CHS.VS.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.VS.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.VS.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.VS.season.duplicated.outliers <- CHS.VS.season.alloutliers[duplicated(CHS.VS.season.alloutliers),]  ##select only loci occurring more than once (here 130 of 1267)

length(CHS.VS.season.alloutliers$loci)
length(CHS.VS.season.duplicated.outliers)

write.table(CHS.VS.season.duplicated.outliers, "CHS.VS.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```


![alt_txt][CHS.VS.Venn.allmethods]

[CHS.VS.Venn.allmethods]:https://user-images.githubusercontent.com/12142475/31231537-0d92898c-a9df-11e7-8890-bec1fef7e6fd.png






##### CHS.TI

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("CHS.TI.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.TI.pcadapt.outliers")
temp.bayenv.outliers <- read.table("CHS.TI.CHP4.temp.outlier.names")
XtX.outliers <- read.table("CHS.TI.XtX.100outliers")
bayescan.outliers <- read.table("CHS.TI.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.TI.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.TI.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.TI.temp.duplicated.outliers <- CHS.TI.temp.alloutliers[duplicated(CHS.TI.temp.alloutliers),]  ##select only loci occurring more than once (here 141 of 1059)

length(CHS.TI.temp.alloutliers$loci)
length(CHS.TI.temp.duplicated.outliers)

write.table(CHS.TI.temp.duplicated.outliers, "CHS.TI.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("CHS.TI.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("CHS.TI.pcadapt.outliers")
season.bayenv.outliers <- read.table("CHS.TI.CHP4.season.outlier.names")
XtX.outliers <- read.table("CHS.TI.XtX.100outliers")
bayescan.outliers <- read.table("CHS.TI.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("CHS.TI.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
CHS.TI.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


CHS.TI.season.duplicated.outliers <- CHS.TI.season.alloutliers[duplicated(CHS.TI.season.alloutliers),]  ##select only loci occurring more than once (here 165 of 1185)

length(CHS.TI.season.alloutliers$loci)
length(CHS.TI.season.duplicated.outliers)

write.table(CHS.TI.season.duplicated.outliers, "CHS.TI.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```

![alt_txt][CHS.TI.allmethods]

[CHS.TI.allmethods]:https://user-images.githubusercontent.com/12142475/31176809-0bdc9db8-a90c-11e7-8779-9d784057db2a.png


##### SE

######## 1. temp

prep data
```
temp.lfmm.outliers <- read.table("SE.CHP4.temp.LFMM.alloutliers")
pcadapt.outliers <- read.table("SE.pcadapt.outliers")
temp.bayenv.outliers <- read.table("SE.CHP4.temp.outlier.names")
XtX.outliers <- read.table("SE.XtX.100outliers")
bayescan.outliers <- read.table("SE.bayescan.outliers.FDR0.01")


temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

temp.lfmm.outliers.new <- gsub("X\\.", "", temp.lfmm.outliers$V1)
temp.lfmm.outliers <- NULL
temp.lfmm.outliers <- as.data.frame(temp.lfmm.outliers.new)

temp.bayenv.outliers.new <- gsub("X\\.", "", temp.bayenv.outliers$V1)
temp.bayenv.outliers <- NULL
temp.bayenv.outliers <- as.data.frame(temp.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(temp.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(temp.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

temp.bayenv.outliers.names <- temp.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
temp.lfmm.outliers.names <- temp.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(temp.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(temp.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,temp.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, temp.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
pcadapt.outliers.names, temp.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, temp.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(temp.lfmm.outliers.names, temp.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(temp.bayenv.outliers.names, XtX.outliers.names, temp.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("SE.Venn.temp.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("temp.bayenv", "XtX", "temp.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
SE.temp.alloutliers <- rbind(bayescan.outliers, temp.lfmm.outliers, temp.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


SE.temp.duplicated.outliers <- SE.temp.alloutliers[duplicated(SE.temp.alloutliers),]  ##select only loci occurring more than once (here 88 of 430)

length(SE.temp.alloutliers$loci)
length(SE.temp.duplicated.outliers)

write.table(SE.temp.duplicated.outliers, "SE.temp.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```




######## 2. season

prep data
```
season.lfmm.outliers <- read.table("SE.CHP4.season.LFMM.alloutliers")
pcadapt.outliers <- read.table("SE.pcadapt.outliers")
season.bayenv.outliers <- read.table("SE.CHP4.season.outlier.names")
XtX.outliers <- read.table("SE.XtX.100outliers")
bayescan.outliers <- read.table("SE.bayescan.outliers.FDR0.01")


season.lfmm.outliers <- as.data.frame(season.lfmm.outliers)
pcadapt.outliers <- as.data.frame(pcadapt.outliers)
XtX.outliers <- as.data.frame(XtX.outliers)
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers)
bayescan.outliers <- as.data.frame(bayescan.outliers)

season.lfmm.outliers.new <- gsub("X\\.", "", season.lfmm.outliers$V1)
season.lfmm.outliers <- NULL
season.lfmm.outliers <- as.data.frame(season.lfmm.outliers.new)

season.bayenv.outliers.new <- gsub("X\\.", "", season.bayenv.outliers$V1)
season.bayenv.outliers <- NULL
season.bayenv.outliers <- as.data.frame(season.bayenv.outliers.new)

XtX.outliers.new <- gsub("X\\.", "", XtX.outliers$V1)
XtX.outliers <- NULL
XtX.outliers <- as.data.frame(XtX.outliers.new)

pcadapt.new <- gsub("X\\.", "", pcadapt.outliers$V1)
pcadapt.new <- as.data.frame(pcadapt.new)
pcadapt.new <- gsub(":", "\\.", pcadapt.new$pcadapt.new)
pcadapt.outliers <- as.data.frame(pcadapt.new)


colnames(season.lfmm.outliers) <- "loci"
colnames(XtX.outliers) <- "loci"
colnames(pcadapt.outliers) <- "loci"
colnames(season.bayenv.outliers) <- "loci"
colnames(bayescan.outliers) <- "loci"

```


and draw Venn
```
library(VennDiagram)

season.bayenv.outliers.names <- season.bayenv.outliers$loci
XtX.outliers.names <- XtX.outliers$loci
season.lfmm.outliers.names <- season.lfmm.outliers$loci
bayescan.outliers.names <- bayescan.outliers$loci
pcadapt.outliers.names <- pcadapt.outliers$loci

d1 <- length(season.bayenv.outliers.names)
d2 <- length(XtX.outliers.names)
d3 <- length(season.lfmm.outliers.names)
d4 <- length(bayescan.outliers.names)
d5 <- length(pcadapt.outliers.names)


d12 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names)))
d13 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names)))
d14 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names)))
d15 <- length(Reduce(intersect, list(season.bayenv.outliers.names, pcadapt.outliers.names)))
d23 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names)))
d24 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names)))
d25 <- length(Reduce(intersect, list(XtX.outliers.names, pcadapt.outliers.names)))
d34 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names)))
d35 <- length(Reduce(intersect, list(season.lfmm.outliers.names, pcadapt.outliers.names)))
d45 <- length(Reduce(intersect, list(bayescan.outliers.names, pcadapt.outliers.names)))

d123 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,season.lfmm.outliers.names)))
d124 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,bayescan.outliers.names)))
d125 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names,pcadapt.outliers.names)))
d234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d134 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,bayescan.outliers.names)))
d135 <- length(Reduce(intersect, list(season.bayenv.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d145 <- length(Reduce(intersect, list(season.bayenv.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,pcadapt.outliers.names)))
d245 <- length(Reduce(intersect, list(XtX.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))
d345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, bayescan.outliers.names,pcadapt.outliers.names)))


d1234 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, season.bayenv.outliers.names)))
d1235 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
pcadapt.outliers.names, season.bayenv.outliers.names)))

d2345 <- length(Reduce(intersect, list(XtX.outliers.names, season.lfmm.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1245 <- length(Reduce(intersect, list(XtX.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d1345 <- length(Reduce(intersect, list(season.lfmm.outliers.names, season.bayenv.outliers.names,
bayescan.outliers.names, pcadapt.outliers.names)))
d12345 <- length(Reduce(intersect, list(season.bayenv.outliers.names, XtX.outliers.names, season.lfmm.outliers.names, 
bayescan.outliers.names, pcadapt.outliers.names)))

pdf("SE.Venn.season.outliers_20171009.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("season.bayenv", "XtX", "season.lfmm", "bayescan", "pcadapt"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


Write all duplicated loci to file

```
SE.season.alloutliers <- rbind(bayescan.outliers, season.lfmm.outliers, season.bayenv.outliers, pcadapt.outliers, XtX.outliers)  ##Join all data.frames by “loci” column. This only works of colnames are the same (at least one column name)


SE.season.duplicated.outliers <- SE.season.alloutliers[duplicated(SE.season.alloutliers),]  ##select only loci occurring more than once (here 83 of 408)

length(SE.season.alloutliers$loci)
length(SE.season.duplicated.outliers)

write.table(SE.season.duplicated.outliers, "SE.season.duplicated.outliers.20171009", col.names=F, row.names=F, quote=F)

```


### 3. Overlap between duplicated candidate loci identified for CHN, CHS, CZ

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/VENN.Chp4


#### 1. temp

```
library(VennDiagram)


CHN.temp.duplicated.outliers <- read.table("CHN/CHN.temp.duplicated.outliers.20171009", header=F)
colnames(CHN.temp.duplicated.outliers) <- "loci"
CHN.temp.duplicated.outliers <- as.character(CHN.temp.duplicated.outliers$loci)

CHS.temp.duplicated.outliers <- read.table("CHS/CHS.temp.duplicated.outliers.20171009", header=F)
colnames(CHS.temp.duplicated.outliers) <- "loci"
CHS.temp.duplicated.outliers <- as.character(CHS.temp.duplicated.outliers$loci)

CZ.temp.duplicated.outliers <- read.table("CZ/CZ.temp.duplicated.outliers.20171009", header=F)
colnames(CZ.temp.duplicated.outliers) <- "loci"
CZ.temp.duplicated.outliers <- as.character(CZ.temp.duplicated.outliers$loci)


d1 <- length(CHN.temp.duplicated.outliers)
d2 <- length(CHS.temp.duplicated.outliers)
d3 <- length(CZ.temp.duplicated.outliers)

d12 <- length(Reduce(intersect, list(CHN.temp.duplicated.outliers, CHS.temp.duplicated.outliers)))
d13 <- length(Reduce(intersect, list(CHN.temp.duplicated.outliers, CZ.temp.duplicated.outliers)))
d23 <- length(Reduce(intersect, list(CHS.temp.duplicated.outliers, CZ.temp.duplicated.outliers)))

d123 <- length(Reduce(intersect, list(CHN.temp.duplicated.outliers, CHS.temp.duplicated.outliers,CZ.temp.duplicated.outliers)))


pdf("Venn.CHN.CHS.CZ.temp.duplicated.outliers.20171009.pdf")
draw.triple.venn(area1=d1, area2=d2, area3=d3, n12=d12, n13=d13, n23=d23, n123=d123, category=c("CHN", "CHS", "CZ"), lty="blank", fill=c("yellow", "orange", "skyblue1"))
dev.off()

```

#### 2. season

```
library(VennDiagram)


CHN.season.duplicated.outliers <- read.table("CHN/CHN.season.duplicated.outliers.20171009", header=F)
colnames(CHN.season.duplicated.outliers) <- "loci"
CHN.season.duplicated.outliers <- as.character(CHN.season.duplicated.outliers$loci)

CHS.season.duplicated.outliers <- read.table("CHS/CHS.season.duplicated.outliers.20171009", header=F)
colnames(CHS.season.duplicated.outliers) <- "loci"
CHS.season.duplicated.outliers <- as.character(CHS.season.duplicated.outliers$loci)

CZ.season.duplicated.outliers <- read.table("CZ/CZ.season.duplicated.outliers.20171009", header=F)
colnames(CZ.season.duplicated.outliers) <- "loci"
CZ.season.duplicated.outliers <- as.character(CZ.season.duplicated.outliers$loci)


d1 <- length(CHN.season.duplicated.outliers)
d2 <- length(CHS.season.duplicated.outliers)
d3 <- length(CZ.season.duplicated.outliers)

d12 <- length(Reduce(intersect, list(CHN.season.duplicated.outliers, CHS.season.duplicated.outliers)))
d13 <- length(Reduce(intersect, list(CHN.season.duplicated.outliers, CZ.season.duplicated.outliers)))
d23 <- length(Reduce(intersect, list(CHS.season.duplicated.outliers, CZ.season.duplicated.outliers)))

d123 <- length(Reduce(intersect, list(CHN.season.duplicated.outliers, CHS.season.duplicated.outliers,CZ.season.duplicated.outliers)))


pdf("Venn.CHN.CHS.CZ.season.duplicated.outliers.20171009.pdf")
draw.triple.venn(area1=d1, area2=d2, area3=d3, n12=d12, n13=d13, n23=d23, n123=d123, category=c("CHN", "CHS", "CZ"), lty="blank", fill=c("yellow", "orange", "skyblue1"))
dev.off()

```





### 4. Overlap between all candidate loci identified for CHS, CHS.TI, CHS.VS

#### 1. temp

```
library(VennDiagram)


CHS.VS.temp.duplicated.outliers <- read.table("CHS.VS/CHS.VS.temp.duplicated.outliers.20171009", header=F)
colnames(CHS.VS.temp.duplicated.outliers) <- "loci"
CHS.VS.temp.duplicated.outliers <- as.character(CHS.VS.temp.duplicated.outliers$loci)

CHS.temp.duplicated.outliers <- read.table("CHS/CHS.temp.duplicated.outliers.20171009", header=F)
colnames(CHS.temp.duplicated.outliers) <- "loci"
CHS.temp.duplicated.outliers <- as.character(CHS.temp.duplicated.outliers$loci)

CHS.TI.temp.duplicated.outliers <- read.table("CHS.TI/CHS.TI.temp.duplicated.outliers.20171009", header=F)
colnames(CHS.TI.temp.duplicated.outliers) <- "loci"
CHS.TI.temp.duplicated.outliers <- as.character(CHS.TI.temp.duplicated.outliers$loci)


d1 <- length(CHS.VS.temp.duplicated.outliers)
d2 <- length(CHS.temp.duplicated.outliers)
d3 <- length(CHS.TI.temp.duplicated.outliers)

d12 <- length(Reduce(intersect, list(CHS.VS.temp.duplicated.outliers, CHS.temp.duplicated.outliers)))
d13 <- length(Reduce(intersect, list(CHS.VS.temp.duplicated.outliers, CHS.TI.temp.duplicated.outliers)))
d23 <- length(Reduce(intersect, list(CHS.temp.duplicated.outliers, CHS.TI.temp.duplicated.outliers)))

d123 <- length(Reduce(intersect, list(CHS.VS.temp.duplicated.outliers, CHS.temp.duplicated.outliers,CHS.TI.temp.duplicated.outliers)))


pdf("Venn.CHS.VS.TI.temp.duplicated.outliers.20171009.pdf")
draw.triple.venn(area1=d1, area2=d2, area3=d3, n12=d12, n13=d13, n23=d23, n123=d123, category=c("CHS.VS", "CHS", "CHS.TI"), lty="blank", fill=c("yellow", "orange", "skyblue1"))
dev.off()

```

#### 2. season

```
library(VennDiagram)


CHS.VS.season.duplicated.outliers <- read.table("CHS.VS/CHS.VS.season.duplicated.outliers.20171009", header=F)
colnames(CHS.VS.season.duplicated.outliers) <- "loci"
CHS.VS.season.duplicated.outliers <- as.character(CHS.VS.season.duplicated.outliers$loci)

CHS.season.duplicated.outliers <- read.table("CHS/CHS.season.duplicated.outliers.20171009", header=F)
colnames(CHS.season.duplicated.outliers) <- "loci"
CHS.season.duplicated.outliers <- as.character(CHS.season.duplicated.outliers$loci)

CHS.TI.season.duplicated.outliers <- read.table("CHS.TI/CHS.TI.season.duplicated.outliers.20171009", header=F)
colnames(CHS.TI.season.duplicated.outliers) <- "loci"
CHS.TI.season.duplicated.outliers <- as.character(CHS.TI.season.duplicated.outliers$loci)


d1 <- length(CHS.VS.season.duplicated.outliers)
d2 <- length(CHS.season.duplicated.outliers)
d3 <- length(CHS.TI.season.duplicated.outliers)

d12 <- length(Reduce(intersect, list(CHS.VS.season.duplicated.outliers, CHS.season.duplicated.outliers)))
d13 <- length(Reduce(intersect, list(CHS.VS.season.duplicated.outliers, CHS.TI.season.duplicated.outliers)))
d23 <- length(Reduce(intersect, list(CHS.season.duplicated.outliers, CHS.TI.season.duplicated.outliers)))

d123 <- length(Reduce(intersect, list(CHS.VS.season.duplicated.outliers, CHS.season.duplicated.outliers,CHS.TI.season.duplicated.outliers)))


pdf("Venn.CHS.VS.TI.season.duplicated.outliers.20171009.pdf")
draw.triple.venn(area1=d1, area2=d2, area3=d3, n12=d12, n13=d13, n23=d23, n123=d123, category=c("CHS.VS", "CHS", "CHS.TI"), lty="blank", fill=c("yellow", "orange", "skyblue1"))
dev.off()

```
