# Generalised dissimilarity model

Aim: Compare the most important environmental drivers of adaptive divergence between gradients. 

This model allows ranking of the most important env. variables. 

I'll use only the candidate adaptive loci. 


/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GDM


The vignette and pdf can both be found here: 

https://cran.r-project.org/web/packages/gdm/index.html



Input files: 

1. Response matrix: Fst matrix based on all candidate loci for each transect. The first column should be all the sites names. i.e. not square

2. Predictor matrix: site by predictor table. 



### CHN

I'm using a dissimilarity matrix, which is what is recommended for Fst data. (see gdm.pdf)

Response matrix: I'm using a site by variable matrix with the 5 important Env Variables and coordinate information

#### CHN
```
install.packages("Gdm01", repos="http://R-Forge.R-project.org")
library(Gdm01)
library(gdm)
library(hierfstat)
library(adegenet)
library(reshape)

CHN.229.25 <- read.structure("CHN.AdaptiveLoci.str")
CHN229.pop <- read.table("CHN229.pop", header=T)
CHN.pop.factor <- as.factor(CHN229.pop$pop)
CHN.229.25@pop <- CHN.pop.factor
CHN.229.25@pop

CHN229.25.fst <- pairwise.fst(CHN.229.25, pop=NULL, res.type=c("dist", "matrix"))  ##this creates the lower half of a pairwise distance matrix
m.CHN229 <- as.matrix(CHN229.25.fst)  ##make this a full matrix
m.CHN229.scaled <- apply(m.CHN229[,1:ncol(m.CHN229)], MARGIN=2, FUN=function(X) (X-min(X))/diff(range(X))) ##scale the data 
CHN.pop.names <- unique(CHN.pop.factor)
rownames(m.CHN229.scaled) <- CHN.pop.names  #select only the unique pop names
write.table(m.CHN229.scaled, "CHN229.25.Fst.table", row.names=T, col.names=F, quote=F)


#####Read in the EnvironData
CHN.EnvData <- read.table("CHN.EnvData", header=T)
CHN.EnvData.importantVariables <- CHN.EnvData[,c("site", "long", "lat", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm")]  ##subset the data to use only the important variables


CHN.fst <- read.table("CHN229.25.Fst.table")
names(CHN.fst)[1] <- "site"  ##make sure site column is named the same in both the Env and fst files


##### GDM
###Create pairwise gdm input format


###Run GDM
CHN.gdm <- gdm(CHN.gdmData, geo=T)

plot(CHN.gdm, plot.layout=c(3,3))

```

