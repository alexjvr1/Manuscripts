# Gradient Forest Analysis

Aim: determine whether mean temp during development or growth season are the most important environmental drivers of adaptive divergence. 

I've identified candidate loci associated with mean temp 60 days and number of days above 6deg C. I have two sets of candidate loci for 
each of the transects: 

1. temp

any locus that was found in at least two of the following tests: 

bayenv2 for temp, lfmm for temp, pcadapt, XtX 100 top loci, and bayescan

2. season

any locus that was found in at least two of the following tests: 

bayenv2 for season, lfmm for season, pcadapt, XtX 100 top loci, and bayescan

3. Neutral loci

A set of 1000 loci from ALL - Candidate loci. 
It's important to subsample the data, since the datasets run for a long time if they have more than 1000 loci. 


Transects: 

CHN, CHS, CHS.VS, CHS.TI, CZ, and SE


# Gradient Forest CHall

### Datasets:

2 Datasets for each regions: Neutral, Potential Adaptive Loci (identified by at least 2 programs)

1. CHN

2. CZ

3. CHS

4. CHS.VS

5. CHS.TI




I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013



## Input files

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/GradientForest

MEM and distance have been calculated for Swiss pops in Chapter 2, so can be copied directly from the .csv input files here: 

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/GradientForest/Oct2017/*Input.csv

### CHN

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHN.TempLoci <- read.table("CHN.temp.duplicated.outliers.20171009")
head(CHN.TempLoci)
CHN.TempLoci <- gsub("\\.", ":", CHN.TempLoci$V1)
head(CHN.TempLoci)
write.table(CHN.TempLoci, "CHN.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.TempLoci.names --recode --recode-INFO-all --out CHN.229.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.TempLoci.recode.vcf --plink --out CHN.229.TempLoci.plink

plink --file CHN.229.TempLoci.plink --noweb --recode --recodeA --out CHN.229.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster

```
##in R

CHN.nosex <- read.table("CHN.229.TempLoci.plink.nosex", header=F)
head(CHN.nosex)
CHN.pop <- gsub("_\\d+", "", CHN.nosex$V1)
CHN.pop
CHN.nosex$V3 <- CHN.pop
head(CHN.nosex)
write.table(CHN.nosex, "CHN.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHN.SeasonLoci <- read.table("CHN.Season.duplicated.outliers.20171009")
head(CHN.SeasonLoci)
CHN.SeasonLoci <- gsub("\\.", ":", CHN.SeasonLoci$V1)
head(CHN.SeasonLoci)
write.table(CHN.SeasonLoci, "CHN.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --snps CHN.SeasonLoci.names --recode --recode-INFO-all --out CHN.229.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.SeasonLoci.recode.vcf --plink --out CHN.229.SeasonLoci.plink

plink --file CHN.229.SeasonLoci.plink --noweb --recode --recodeA --out CHN.229.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHN.229.SeasonLoci.plink --within CHN.PlinkCluster --freq --noweb --out CHN.229.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.SeasonLoci.MAF <- read.table("CHN.229.SeasonLoci.frq.strat", header=T)
head(CHN.SeasonLoci.MAF)

CHN.SeasonLoci.MAF <- CHN.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.SeasonLoci.MAF2 <- melt(CHN.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.SeasonLoci.MAF2)
head(CHN.SeasonLoci.MAF2)


CHN.SeasonLoci.MAF3 <- dcast(CHN.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHN.SeasonLoci.MAF3)
colnames(CHN.SeasonLoci.MAF3) <- paste("X", colnames(CHN.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.SeasonLoci.MAF3, file="CHN.229.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHN.229.5265.recode.vcf --exclude CHN.SeasonLoci.names --recode --recode-INFO-all --out CHN.229.Neutral.temporary
vcftools --vcf CHN.229.Neutral.temporary.recode.vcf --exclude CHN.TempLoci.names --recode --recode-INFO-all --out CHN.229.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHN.229.Neutral.recode.vcf --plink --out CHN.229.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHN.Neutral.loci.names.map <- read.table("CHN.229.Neutral.plink.map", header=F)
CHN.Neutral.loci.names <-  CHN.Neutral.loci.names.map$V2
CHN.Neutral.loci.names <- as.data.frame(CHN.Neutral.loci.names)
CHN.1000.loci.names <- CHN.Neutral.loci.names[sample(nrow(CHN.Neutral.loci.names),1000),]
CHN.1000.loci.names <- as.data.frame(CHN.1000.loci.names)
summary(CHN.1000.loci.names)
write.table(CHN.1000.loci.names, "CHN.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHN.229.Neutral.recode.vcf --snps CHN.1000.loci.names --recode --recode-INFO-all --out CHN.1000NeutralLoci

vcftools --vcf CHN.1000NeutralLoci.recode.vcf --plink --out CHN.1000NeutralLoci.plink
plink --file CHN.1000NeutralLoci.plink --noweb --recode --recodeA --out CHN.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHN.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHN.1000NeutralLoci.plink --within CHN.PlinkCluster --noweb --freq --out CHN.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHN.Neutral.MAF <- read.table("CHN.1000NeutralLoci.frq.strat", header=T)
head(CHN.Neutral.MAF)

CHN.Neutral.MAF <- CHN.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHN.Neutral.MAF2 <- melt(CHN.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHN.Neutral.MAF2)
head(CHN.Neutral.MAF2)


CHN.Neutral.MAF3 <- dcast(CHN.Neutral.MAF2, formula= CLST ~ SNP)
head(CHN.Neutral.MAF3)
colnames(CHN.Neutral.MAF3) <- paste("X", colnames(CHN.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHN.Neutral.MAF3, file="CHN.229.Neutral.MAF.csv")
```




### CZ

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CZ.TempLoci <- read.table("CZ.temp.duplicated.outliers.20171009")
head(CZ.TempLoci)
CZ.TempLoci <- gsub("\\.", ":", CZ.TempLoci$V1)
head(CZ.TempLoci)
write.table(CZ.TempLoci, "CZ.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CZ.404.5265.recode.vcf --snps CZ.TempLoci.names --recode --recode-INFO-all --out CZ.404.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.TempLoci.recode.vcf --plink --out CZ.404.TempLoci.plink

plink --file CZ.404.TempLoci.plink --noweb --recode --recodeA --out CZ.404.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster

```
##in R

CZ.nosex <- read.table("CZ.404.TempLoci.plink.nosex", header=F)
head(CZ.nosex)
CZ.pop <- gsub("_\\d+", "", CZ.nosex$V1)
CZ.pop
CZ.nosex$V3 <- CZ.pop
head(CZ.nosex)
write.table(CZ.nosex, "CZ.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CZ.SeasonLoci <- read.table("CZ.Season.duplicated.outliers.20171009")
head(CZ.SeasonLoci)
CZ.SeasonLoci <- gsub("\\.", ":", CZ.SeasonLoci$V1)
head(CZ.SeasonLoci)
write.table(CZ.SeasonLoci, "CZ.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --snps CZ.SeasonLoci.names --recode --recode-INFO-all --out CZ.404.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.SeasonLoci.recode.vcf --plink --out CZ.404.SeasonLoci.plink

plink --file CZ.404.SeasonLoci.plink --noweb --recode --recodeA --out CZ.404.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CZ.404.SeasonLoci.plink --within CZ.PlinkCluster --freq --noweb --out CZ.404.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.SeasonLoci.MAF <- read.table("CZ.404.SeasonLoci.frq.strat", header=T)
head(CZ.SeasonLoci.MAF)

CZ.SeasonLoci.MAF <- CZ.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.SeasonLoci.MAF2 <- melt(CZ.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.SeasonLoci.MAF2)
head(CZ.SeasonLoci.MAF2)


CZ.SeasonLoci.MAF3 <- dcast(CZ.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CZ.SeasonLoci.MAF3)
colnames(CZ.SeasonLoci.MAF3) <- paste("X", colnames(CZ.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.SeasonLoci.MAF3, file="CZ.404.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CZ.404.7288.recode.vcf --exclude CZ.SeasonLoci.names --recode --recode-INFO-all --out CZ.404.Neutral.temporary
vcftools --vcf CZ.404.Neutral.temporary.recode.vcf --exclude CZ.TempLoci.names --recode --recode-INFO-all --out CZ.404.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CZ.404.Neutral.recode.vcf --plink --out CZ.404.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CZ.Neutral.loci.names.map <- read.table("CZ.404.Neutral.plink.map", header=F)
CZ.Neutral.loci.names <-  CZ.Neutral.loci.names.map$V2
CZ.Neutral.loci.names <- as.data.frame(CZ.Neutral.loci.names)
CZ.1000.loci.names <- CZ.Neutral.loci.names[sample(nrow(CZ.Neutral.loci.names),1000),]
CZ.1000.loci.names <- as.data.frame(CZ.1000.loci.names)
summary(CZ.1000.loci.names)
write.table(CZ.1000.loci.names, "CZ.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CZ.404.Neutral.recode.vcf --snps CZ.1000.loci.names --recode --recode-INFO-all --out CZ.1000NeutralLoci

vcftools --vcf CZ.1000NeutralLoci.recode.vcf --plink --out CZ.1000NeutralLoci.plink
plink --file CZ.1000NeutralLoci.plink --noweb --recode --recodeA --out CZ.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CZ.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CZ.1000NeutralLoci.plink --within CZ.PlinkCluster --noweb --freq --out CZ.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CZ.Neutral.MAF <- read.table("CZ.1000NeutralLoci.frq.strat", header=T)
head(CZ.Neutral.MAF)

CZ.Neutral.MAF <- CZ.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CZ.Neutral.MAF2 <- melt(CZ.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CZ.Neutral.MAF2)
head(CZ.Neutral.MAF2)


CZ.Neutral.MAF3 <- dcast(CZ.Neutral.MAF2, formula= CLST ~ SNP)
head(CZ.Neutral.MAF3)
colnames(CZ.Neutral.MAF3) <- paste("X", colnames(CZ.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CZ.Neutral.MAF3, file="CZ.404.Neutral.MAF.csv")
```


### CHS

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TempLoci <- read.table("CHS.temp.duplicated.outliers.20171009")
head(CHS.TempLoci)
CHS.TempLoci <- gsub("\\.", ":", CHS.TempLoci$V1)
head(CHS.TempLoci)
write.table(CHS.TempLoci, "CHS.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --snps CHS.TempLoci.names --recode --recode-INFO-all --out CHS.275.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.TempLoci.recode.vcf --plink --out CHS.275.TempLoci.plink

plink --file CHS.275.TempLoci.plink --noweb --recode --recodeA --out CHS.275.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.PlinkCluster

```
##in R

CHS.nosex <- read.table("CHS.275.TempLoci.plink.nosex", header=F)
head(CHS.nosex)
CHS.pop <- gsub("_\\d+", "", CHS.nosex$V1)
CHS.pop
CHS.nosex$V3 <- CHS.pop
head(CHS.nosex)
write.table(CHS.nosex, "CHS.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.SeasonLoci <- read.table("CHS.Season.duplicated.outliers.20171009")
head(CHS.SeasonLoci)
CHS.SeasonLoci <- gsub("\\.", ":", CHS.SeasonLoci$V1)
head(CHS.SeasonLoci)
write.table(CHS.SeasonLoci, "CHS.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --snps CHS.SeasonLoci.names --recode --recode-INFO-all --out CHS.275.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.SeasonLoci.recode.vcf --plink --out CHS.275.SeasonLoci.plink

plink --file CHS.275.SeasonLoci.plink --noweb --recode --recodeA --out CHS.275.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.275.SeasonLoci.plink --within CHS.PlinkCluster --freq --noweb --out CHS.275.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.SeasonLoci.MAF <- read.table("CHS.275.SeasonLoci.frq.strat", header=T)
head(CHS.SeasonLoci.MAF)

CHS.SeasonLoci.MAF <- CHS.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.SeasonLoci.MAF2 <- melt(CHS.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.SeasonLoci.MAF2)
head(CHS.SeasonLoci.MAF2)


CHS.SeasonLoci.MAF3 <- dcast(CHS.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.SeasonLoci.MAF3)
colnames(CHS.SeasonLoci.MAF3) <- paste("X", colnames(CHS.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.SeasonLoci.MAF3, file="CHS.275.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.275.6339.recode.vcf --exclude CHS.SeasonLoci.names --recode --recode-INFO-all --out CHS.275.Neutral.temporary
vcftools --vcf CHS.275.Neutral.temporary.recode.vcf --exclude CHS.TempLoci.names --recode --recode-INFO-all --out CHS.275.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.275.Neutral.recode.vcf --plink --out CHS.275.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.Neutral.loci.names.map <- read.table("CHS.275.Neutral.plink.map", header=F)
CHS.Neutral.loci.names <-  CHS.Neutral.loci.names.map$V2
CHS.Neutral.loci.names <- as.data.frame(CHS.Neutral.loci.names)
CHS.1000.loci.names <- CHS.Neutral.loci.names[sample(nrow(CHS.Neutral.loci.names),1000),]
CHS.1000.loci.names <- as.data.frame(CHS.1000.loci.names)
summary(CHS.1000.loci.names)
write.table(CHS.1000.loci.names, "CHS.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.275.Neutral.recode.vcf --snps CHS.1000.loci.names --recode --recode-INFO-all --out CHS.1000NeutralLoci

vcftools --vcf CHS.1000NeutralLoci.recode.vcf --plink --out CHS.1000NeutralLoci.plink
plink --file CHS.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.1000NeutralLoci.plink --within CHS.PlinkCluster --noweb --freq --out CHS.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.Neutral.MAF <- read.table("CHS.1000NeutralLoci.frq.strat", header=T)
head(CHS.Neutral.MAF)

CHS.Neutral.MAF <- CHS.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.Neutral.MAF2 <- melt(CHS.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.Neutral.MAF2)
head(CHS.Neutral.MAF2)


CHS.Neutral.MAF3 <- dcast(CHS.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.Neutral.MAF3)
colnames(CHS.Neutral.MAF3) <- paste("X", colnames(CHS.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.Neutral.MAF3, file="CHS.275.Neutral.MAF.csv")
```




### CHS.TI

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TI.TempLoci <- read.table("CHS.TI.temp.duplicated.outliers.20171009")
head(CHS.TI.TempLoci)
CHS.TI.TempLoci <- gsub("\\.", ":", CHS.TI.TempLoci$V1)
head(CHS.TI.TempLoci)
write.table(CHS.TI.TempLoci, "CHS.TI.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.TempLoci.names --recode --recode-INFO-all --out CHS.TI.140.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.TempLoci.recode.vcf --plink --out CHS.TI.140.TempLoci.plink

plink --file CHS.TI.140.TempLoci.plink --noweb --recode --recodeA --out CHS.TI.140.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster

```
##in R

CHS.TI.nosex <- read.table("CHS.TI.140.TempLoci.plink.nosex", header=F)
head(CHS.TI.nosex)
CHS.TI.pop <- gsub("_\\d+", "", CHS.TI.nosex$V1)
CHS.TI.pop
CHS.TI.nosex$V3 <- CHS.TI.pop
head(CHS.TI.nosex)
write.table(CHS.TI.nosex, "CHS.TI.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.TI.SeasonLoci <- read.table("CHS.TI.Season.duplicated.outliers.20171009")
head(CHS.TI.SeasonLoci)
CHS.TI.SeasonLoci <- gsub("\\.", ":", CHS.TI.SeasonLoci$V1)
head(CHS.TI.SeasonLoci)
write.table(CHS.TI.SeasonLoci, "CHS.TI.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --snps CHS.TI.SeasonLoci.names --recode --recode-INFO-all --out CHS.TI.140.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.SeasonLoci.recode.vcf --plink --out CHS.TI.140.SeasonLoci.plink

plink --file CHS.TI.140.SeasonLoci.plink --noweb --recode --recodeA --out CHS.TI.140.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.TI.140.SeasonLoci.plink --within CHS.TI.PlinkCluster --freq --noweb --out CHS.TI.140.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.SeasonLoci.MAF <- read.table("CHS.TI.140.SeasonLoci.frq.strat", header=T)
head(CHS.TI.SeasonLoci.MAF)

CHS.TI.SeasonLoci.MAF <- CHS.TI.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.SeasonLoci.MAF2 <- melt(CHS.TI.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.SeasonLoci.MAF2)
head(CHS.TI.SeasonLoci.MAF2)


CHS.TI.SeasonLoci.MAF3 <- dcast(CHS.TI.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.TI.SeasonLoci.MAF3)
colnames(CHS.TI.SeasonLoci.MAF3) <- paste("X", colnames(CHS.TI.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.SeasonLoci.MAF3, file="CHS.TI.140.SeasonLoci.MAF.csv")
```


### SE

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.TempLoci <- read.table("SE.temp.duplicated.outliers.20171009")
head(SE.TempLoci)
SE.TempLoci <- gsub("\\.", ":", SE.TempLoci$V1)
head(SE.TempLoci)
write.table(SE.TempLoci, "SE.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.TempLoci.names --recode --recode-INFO-all --out SE.132.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.TempLoci.recode.vcf --plink --out SE.132.TempLoci.plink

plink --file SE.132.TempLoci.plink --noweb --recode --recodeA --out SE.132.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster

```
##in R

SE.nosex <- read.table("SE.132.TempLoci.plink.nosex", header=F)
head(SE.nosex)
SE.pop <- gsub("_\\d+", "", SE.nosex$V1)
SE.pop
SE.nosex$V3 <- SE.pop
head(SE.nosex)
write.table(SE.nosex, "SE.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

SE.SeasonLoci <- read.table("SE.Season.duplicated.outliers.20171009")
head(SE.SeasonLoci)
SE.SeasonLoci <- gsub("\\.", ":", SE.SeasonLoci$V1)
head(SE.SeasonLoci)
write.table(SE.SeasonLoci, "SE.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --snps SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.SeasonLoci.recode.vcf --plink --out SE.132.SeasonLoci.plink

plink --file SE.132.SeasonLoci.plink --noweb --recode --recodeA --out SE.132.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file SE.132.SeasonLoci.plink --within SE.PlinkCluster --freq --noweb --out SE.132.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.SeasonLoci.MAF <- read.table("SE.132.SeasonLoci.frq.strat", header=T)
head(SE.SeasonLoci.MAF)

SE.SeasonLoci.MAF <- SE.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.SeasonLoci.MAF2 <- melt(SE.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.SeasonLoci.MAF2)
head(SE.SeasonLoci.MAF2)


SE.SeasonLoci.MAF3 <- dcast(SE.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(SE.SeasonLoci.MAF3)
colnames(SE.SeasonLoci.MAF3) <- paste("X", colnames(SE.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.SeasonLoci.MAF3, file="SE.132.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf SE132.2027.recode.vcf --exclude SE.SeasonLoci.names --recode --recode-INFO-all --out SE.132.Neutral.temporary
vcftools --vcf SE.132.Neutral.temporary.recode.vcf --exclude SE.TempLoci.names --recode --recode-INFO-all --out SE.132.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf SE.132.Neutral.recode.vcf --plink --out SE.132.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
SE.Neutral.loci.names.map <- read.table("SE.132.Neutral.plink.map", header=F)
SE.Neutral.loci.names <-  SE.Neutral.loci.names.map$V2
SE.Neutral.loci.names <- as.data.frame(SE.Neutral.loci.names)
SE.1000.loci.names <- SE.Neutral.loci.names[sample(nrow(SE.Neutral.loci.names),1000),]
SE.1000.loci.names <- as.data.frame(SE.1000.loci.names)
summary(SE.1000.loci.names)
write.table(SE.1000.loci.names, "SE.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf SE.132.Neutral.recode.vcf --snps SE.1000.loci.names --recode --recode-INFO-all --out SE.1000NeutralLoci

vcftools --vcf SE.1000NeutralLoci.recode.vcf --plink --out SE.1000NeutralLoci.plink
plink --file SE.1000NeutralLoci.plink --noweb --recode --recodeA --out SE.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > SE.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file SE.1000NeutralLoci.plink --within SE.PlinkCluster --noweb --freq --out SE.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



SE.Neutral.MAF <- read.table("SE.1000NeutralLoci.frq.strat", header=T)
head(SE.Neutral.MAF)

SE.Neutral.MAF <- SE.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

SE.Neutral.MAF2 <- melt(SE.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.Neutral.MAF2)
head(SE.Neutral.MAF2)


SE.Neutral.MAF3 <- dcast(SE.Neutral.MAF2, formula= CLST ~ SNP)
head(SE.Neutral.MAF3)
colnames(SE.Neutral.MAF3) <- paste("X", colnames(SE.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(SE.Neutral.MAF3, file="SE.132.Neutral.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.TI.140.5692.recode.vcf --exclude CHS.TI.SeasonLoci.names --recode --recode-INFO-all --out CHS.TI.140.Neutral.temporary
vcftools --vcf CHS.TI.140.Neutral.temporary.recode.vcf --exclude CHS.TI.TempLoci.names --recode --recode-INFO-all --out CHS.TI.140.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.TI.140.Neutral.recode.vcf --plink --out CHS.TI.140.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.TI.Neutral.loci.names.map <- read.table("CHS.TI.140.Neutral.plink.map", header=F)
CHS.TI.Neutral.loci.names <-  CHS.TI.Neutral.loci.names.map$V2
CHS.TI.Neutral.loci.names <- as.data.frame(CHS.TI.Neutral.loci.names)
CHS.TI.1000.loci.names <- CHS.TI.Neutral.loci.names[sample(nrow(CHS.TI.Neutral.loci.names),1000),]
CHS.TI.1000.loci.names <- as.data.frame(CHS.TI.1000.loci.names)
summary(CHS.TI.1000.loci.names)
write.table(CHS.TI.1000.loci.names, "CHS.TI.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.TI.140.Neutral.recode.vcf --snps CHS.TI.1000.loci.names --recode --recode-INFO-all --out CHS.TI.1000NeutralLoci

vcftools --vcf CHS.TI.1000NeutralLoci.recode.vcf --plink --out CHS.TI.1000NeutralLoci.plink
plink --file CHS.TI.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.TI.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.TI.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.TI.1000NeutralLoci.plink --within CHS.TI.PlinkCluster --noweb --freq --out CHS.TI.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.TI.Neutral.MAF <- read.table("CHS.TI.1000NeutralLoci.frq.strat", header=T)
head(CHS.TI.Neutral.MAF)

CHS.TI.Neutral.MAF <- CHS.TI.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.TI.Neutral.MAF2 <- melt(CHS.TI.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.TI.Neutral.MAF2)
head(CHS.TI.Neutral.MAF2)


CHS.TI.Neutral.MAF3 <- dcast(CHS.TI.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.TI.Neutral.MAF3)
colnames(CHS.TI.Neutral.MAF3) <- paste("X", colnames(CHS.TI.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.TI.Neutral.MAF3, file="CHS.TI.140.Neutral.MAF.csv")
```



### CHS.VS

##### 1. TempLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.VS.TempLoci <- read.table("CHS.VS.temp.duplicated.outliers.20171009")
head(CHS.VS.TempLoci)
CHS.VS.TempLoci <- gsub("\\.", ":", CHS.VS.TempLoci$V1)
head(CHS.VS.TempLoci)
write.table(CHS.VS.TempLoci, "CHS.VS.TempLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Temp loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.TempLoci.names --recode --recode-INFO-all --out CHS.VS.135.TempLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.TempLoci.recode.vcf --plink --out CHS.VS.135.TempLoci.plink

plink --file CHS.VS.135.TempLoci.plink --noweb --recode --recodeA --out CHS.VS.135.TempLoci.plink

```

Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster

```
##in R

CHS.VS.nosex <- read.table("CHS.VS.135.TempLoci.plink.nosex", header=F)
head(CHS.VS.nosex)
CHS.VS.pop <- gsub("_\\d+", "", CHS.VS.nosex$V1)
CHS.VS.pop
CHS.VS.nosex$V3 <- CHS.VS.pop
head(CHS.VS.nosex)
write.table(CHS.VS.nosex, "CHS.VS.PlinkCluster", quote=F, row.names=F, col.names=F)

```

##### 2. SeasonLoci
```

###Change Duplicated Adaptive Loci names to vcf format

##R

CHS.VS.SeasonLoci <- read.table("CHS.VS.Season.duplicated.outliers.20171009")
head(CHS.VS.SeasonLoci)
CHS.VS.SeasonLoci <- gsub("\\.", ":", CHS.VS.SeasonLoci$V1)
head(CHS.VS.SeasonLoci)
write.table(CHS.VS.SeasonLoci, "CHS.VS.SeasonLoci.names", quote=F, row.names=F, col.names=F)


#Filter the Season loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS.SeasonLoci.names --recode --recode-INFO-all --out CHS.VS.135.SeasonLoci

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.SeasonLoci.recode.vcf --plink --out CHS.VS.135.SeasonLoci.plink

plink --file CHS.VS.135.SeasonLoci.plink --noweb --recode --recodeA --out CHS.VS.135.SeasonLoci.plink

```


And calculate MAF with Plink

```

plink --file CHS.VS.135.SeasonLoci.plink --within CHS.VS.PlinkCluster --freq --noweb --out CHS.VS.135.SeasonLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.SeasonLoci.MAF <- read.table("CHS.VS.135.SeasonLoci.frq.strat", header=T)
head(CHS.VS.SeasonLoci.MAF)

CHS.VS.SeasonLoci.MAF <- CHS.VS.SeasonLoci.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.SeasonLoci.MAF2 <- melt(CHS.VS.SeasonLoci.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.SeasonLoci.MAF2)
head(CHS.VS.SeasonLoci.MAF2)


CHS.VS.SeasonLoci.MAF3 <- dcast(CHS.VS.SeasonLoci.MAF2, formula= CLST ~ SNP)
head(CHS.VS.SeasonLoci.MAF3)
colnames(CHS.VS.SeasonLoci.MAF3) <- paste("X", colnames(CHS.VS.SeasonLoci.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.SeasonLoci.MAF3, file="CHS.VS.135.SeasonLoci.MAF.csv")
```




##### 3. Neutral loci


```

#Filter the Neutral loci from the vcf file

vcftools --vcf CHS.VS.135.5835.recode.vcf --exclude CHS.VS.SeasonLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral.temporary
vcftools --vcf CHS.VS.135.Neutral.temporary.recode.vcf --exclude CHS.VS.TempLoci.names --recode --recode-INFO-all --out CHS.VS.135.Neutral

#Create a plink file from the vcf file. 

vcftools --vcf CHS.VS.135.Neutral.recode.vcf --plink --out CHS.VS.135.Neutral.plink

```

Import the .map file into R to get a list of the SNP names and to subset this to 1000 random names
```
CHS.VS.Neutral.loci.names.map <- read.table("CHS.VS.135.Neutral.plink.map", header=F)
CHS.VS.Neutral.loci.names <-  CHS.VS.Neutral.loci.names.map$V2
CHS.VS.Neutral.loci.names <- as.data.frame(CHS.VS.Neutral.loci.names)
CHS.VS.1000.loci.names <- CHS.VS.Neutral.loci.names[sample(nrow(CHS.VS.Neutral.loci.names),1000),]
CHS.VS.1000.loci.names <- as.data.frame(CHS.VS.1000.loci.names)
summary(CHS.VS.1000.loci.names)
write.table(CHS.VS.1000.loci.names, "CHS.VS.1000.loci.names", quote=F, row.names=F, col.names=F)
```

Subset the vcf file to get 1000 loci and convert to plink
```
vcftools --vcf CHS.VS.135.Neutral.recode.vcf --snps CHS.VS.1000.loci.names --recode --recode-INFO-all --out CHS.VS.1000NeutralLoci

vcftools --vcf CHS.VS.1000NeutralLoci.recode.vcf --plink --out CHS.VS.1000NeutralLoci.plink
plink --file CHS.VS.1000NeutralLoci.plink --noweb --recode --recodeA --out CHS.VS.1000NeutralLoci.plink
```


Find the sample names in the *nosex file, and add pop names (i.e. 3 columns) to create a file for specifying clusters > CHS.VS.PlinkCluster (done before)

And calculate MAF with Plink

```

plink --file CHS.VS.1000NeutralLoci.plink --within CHS.VS.PlinkCluster --noweb --freq --out CHS.VS.1000NeutralLoci

```

Import into R to reformat the output - by population and loci as columns
```
######Reformat PLINK output
###For Gradient Forest
###MAF for each locus -> melt and reformat rows as pops, and columns as loci. 



CHS.VS.Neutral.MAF <- read.table("CHS.VS.1000NeutralLoci.frq.strat", header=T)
head(CHS.VS.Neutral.MAF)

CHS.VS.Neutral.MAF <- CHS.VS.Neutral.MAF[,c(3,2,6)]

library("ggplot2")
library("reshape2")

CHS.VS.Neutral.MAF2 <- melt(CHS.VS.Neutral.MAF, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(CHS.VS.Neutral.MAF2)
head(CHS.VS.Neutral.MAF2)


CHS.VS.Neutral.MAF3 <- dcast(CHS.VS.Neutral.MAF2, formula= CLST ~ SNP)
head(CHS.VS.Neutral.MAF3)
colnames(CHS.VS.Neutral.MAF3) <- paste("X", colnames(CHS.VS.Neutral.MAF3), sep=".")  ##Change colnames, so that excel doesn't change the SNP names
write.csv(CHS.VS.Neutral.MAF3, file="CHS.VS.135.Neutral.MAF.csv")
```







#### ENV data

Variables used: 

Same as for the EAA: 

1. sol.rad.60d

2. temp.laying.date

3. pcpt.60d

4. shadow.days

5. day10cm


Create a .csv input file with

1. pop names

2. distance

Use the last population (trpa) as distance 0. So all distances are measured from here.

##### Input Files

CHN.GF.AdaptiveLoci.Input.csv

CHN.GF.NeutralLoci.Input.csv

CZ.GF.AdaptiveLoci.Input.csv

CZ.GF.NeutralLoci.Input.csv

CHS.TI.GF.AdaptiveLoci.Input.csv

CHS.TI.GF.NeutralLoci.Input.csv

CHS.VS.GF.AdaptiveLoci.Input.csv

CHS.VS.GF.NeutralLoci.Input.csv


### Run GF

1. Adaptive Loci
```
library(gradientForest)

###CHN

gf.CHN.Adaptive <- read.csv("CHN.GF.AdaptiveLoci.Input.csv", header=T)
envGF.CHN.Adaptive <- gf.CHN.Adaptive[,-1]
colnames(envGF.CHN.Adaptive)

CHN.AdaptiveSNPS <- CHN.Adaptive.MAF3[,grep("X.", colnames(CHN.Adaptive.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHN.Adaptive)/2)
maxLevel

gf.CHN.Adaptive.SNPs <- gradientForest(cbind(envGF.CHN.Adaptive, CHN.AdaptiveSNPS), predictor.vars=colnames(envGF.CHN.Adaptive), response.vars=colnames(CHN.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##12 warnings (loci variable in less than 5 populations)


###CZ

gf.CZ.Adaptive <- read.csv("CZ.GF.AdaptiveLoci.Input.csv", header=T)
envGF.CZ.Adaptive <- gf.CZ.Adaptive[,-1]
colnames(envGF.CZ.Adaptive)

CZ.AdaptiveSNPS <- CZ.Adaptive.MAF3[,grep("X.", colnames(CZ.Adaptive.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CZ.Adaptive)/2)
maxLevel

gf.CZ.Adaptive.SNPs <- gradientForest(cbind(envGF.CZ.Adaptive, CZ.AdaptiveSNPS), predictor.vars=colnames(envGF.CZ.Adaptive), response.vars=colnames(CZ.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

gf.CHS.TI.Adaptive <- read.csv("CHS.TI.GF.AdaptiveLoci.Input.csv", header=T)
envGF.CHS.TI.Adaptive <- gf.CHS.TI.Adaptive[,-1]
colnames(envGF.CHS.TI.Adaptive)

CHS.TI.AdaptiveSNPS <- CHS.TI.Adaptive.MAF3[,grep("X.", colnames(CHS.TI.Adaptive.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.Adaptive)/2)
maxLevel

gf.CHS.TI.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.TI.Adaptive, CHS.TI.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.TI.Adaptive), response.vars=colnames(CHS.TI.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

##26 warnings - loci variable in less than 5 pops

###CHS.VS

gf.CHS.VS.Adaptive <- read.csv("CHS.VS.GF.AdaptiveLoci.Input.csv", header=T)
envGF.CHS.VS.Adaptive <- gf.CHS.VS.Adaptive[,-1]
colnames(envGF.CHS.VS.Adaptive)

CHS.VS.AdaptiveSNPS <- CHS.VS.Adaptive.MAF3[,grep("X.", colnames(CHS.VS.Adaptive.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Adaptive)/2)
maxLevel

gf.CHS.VS.Adaptive.SNPs <- gradientForest(cbind(envGF.CHS.VS.Adaptive, CHS.VS.AdaptiveSNPS), predictor.vars=colnames(envGF.CHS.VS.Adaptive), response.vars=colnames(CHS.VS.AdaptiveSNPS), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

###more than 50 warnings about non-polymorphic loci

```

2. Neutral Loci

```
###CHN

gf.CHN.Neutral <- read.csv("CHN.GF.Neutral.MAF.csv", header=T)
envGF.CHN.Neutral <- gf.CHN.Neutral[,-1]
colnames(envGF.CHN.Neutral)

CHN.Neutral.SNPs <- CHN.Neutral.MAF3[,grep("X.", colnames(CHN.Neutral.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHN.Neutral)/2)
maxLevel

gf.CHN.Neutral.SNPs <- gradientForest(cbind(envGF.CHN.Neutral, CHN.Neutral.SNPs), predictor.vars=colnames(envGF.CHN.Neutral), response.vars=colnames(CHN.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)

###CZ

gf.CZ.Neutral <- read.csv("CZ.GF.Neutral.MAF.csv", header=T)
envGF.CZ.Neutral <- gf.CZ.Neutral[,-1]
colnames(envGF.CZ.Neutral)

CZ.Neutral.SNPs <- CZ.Neutral.MAF3[,grep("X.", colnames(CZ.Neutral.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CZ.Neutral)/2)
maxLevel

gf.CZ.Neutral.SNPs <- gradientForest(cbind(envGF.CZ.Neutral, CZ.Neutral.SNPs), predictor.vars=colnames(envGF.CZ.Neutral), response.vars=colnames(CZ.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.TI

gf.CHS.TI.Neutral <- read.csv("CHS.TI.GF.NeutralLoci.Input.csv", header=T)
envGF.CHS.TI.Neutral <- gf.CHS.TI.Neutral[,-1]
colnames(envGF.CHS.TI.Neutral)

CHS.TI.Neutral.SNPs <- CHS.TI.Neutral.MAF3[,grep("X.", colnames(CHS.TI.Neutral.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.TI.Neutral)/2)
maxLevel

gf.CHS.TI.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.TI.Neutral, CHS.TI.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.TI.Neutral), response.vars=colnames(CHS.TI.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)


###CHS.VS

gf.CHS.VS.Neutral <- read.csv("CHS.VS.GF.NeutralLoci.Input.csv", header=T)
envGF.CHS.VS.Neutral <- gf.CHS.VS.Neutral[,-1]
colnames(envGF.CHS.VS.Neutral)

CHS.VS.Neutral.SNPs <- CHS.VS.Neutral.MAF3[,grep("X.", colnames(CHS.VS.Neutral.MAF3))]
maxLevel <- log2(0.368*nrow(envGF.CHS.VS.Neutral)/2)
maxLevel

gf.CHS.VS.Neutral.SNPs <- gradientForest(cbind(envGF.CHS.VS.Neutral, CHS.VS.Neutral.SNPs), predictor.vars=colnames(envGF.CHS.VS.Neutral), response.vars=colnames(CHS.VS.Neutral.SNPs), ntree=2000, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
```



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
###CHS.VS

perf.NEUTRAL <- importance(gf.CHS.VS.Neutral.SNPs, type="Species")
perf.Fst <- importance(gf.CHS.VS.Adaptive.SNPs, type="Species")

##mean and range
summary(perf.NEUTRAL)
summary(perf.Fst)

perf.NEUT.df <- as.data.frame(perf.NEUTRAL)
perf.Fst.df <- as.data.frame(perf.Fst)

##count the number of loci above R2 of x (here 0.25)
length(perf.NEUT.df[which(perf.NEUT.df$perf.NEUTRAL>0.25),])
length(perf.Fst.df[which(perf.Fst.df>0.25),])

#variable in more than 5 pops: this will be the number of loci run in the final model. Can be seen with: 
gf.CHS.VS.Neutral.SNPs
gf.CHS.VS.Adaptive.SNPs

#How many loci were originally included?

##And we can plot the frequency of R2 for each dataset: 
par(mfrow=c(1,2))
hist(perf.NEUT.df$perf.NEUTRAL)
hist(perf.Fst.df$perf.Fst)
```




### Fig1: R2 heatmap
```
#R

###CHN

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHN.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHN.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHN.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHN") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()


###CZ
#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CZ.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CZ.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CZ.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CZ") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()



###CHS.TI

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHS.TI.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHS.TI.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.TI.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.TI") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()


###CHS.VS

#Retrieve results from the individual outputs for both GF models: 

R.sq.alldatasets <- (rowMeans(gf.CHS.VS.Adaptive.SNPs$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "Adaptive"

R.sq.Neutral <- (rowMeans(gf.CHS.VS.Neutral.SNPs$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)

##we're only interested in the first 8 variables
R.sq.EnvVariables <- R.sq.alldatasets[1:8,]  #get the first 8 rows from the Adaptive dataset
R.sq.EnvVariables <- as.data.frame(R.sq.EnvVariables)
colnames(R.sq.EnvVariables) <- "Adaptive"
row.names(R.sq.EnvVariables) <- c("dist", "sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm", "MEM1", "MEM2")
R.sq.EnvVariables$Neutral <- R.sq.Neutral[1:8,] ##add the first 8 rows from the Neutral dataset. Make sure the order is the same in both. 

R.sq.EnvVariables <- as.matrix(R.sq.EnvVariables)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.EnvVariables)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

R.sq.all.melt

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

pdf("CHS.VS.R2plot.pdf")
ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2 CHS.VS") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
dev.off()
```

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/31220345-07976010-a9b8-11e7-83de-bf481fa4cc4e.png



### Fig2: Average allelic turnover across all three datasets. 

It took me a while to find the code for the gf plots so that I could figure out what was being plotted, and so combine the plots. 
```
leg.txt <- c("Neutral", "Adaptive")  ##define the text that will be added to the legend
```






##### CHS.TI Cumulative plot
```
##sol.rad.60d

CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "sol.rad.60d") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "sol.rad.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

ymax=0.0015 #set the ymax so that it's the same for the the plots to be overlaid

pdf("sol.rad.60d.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="sol.rad.60d: mean solar radiation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##temp.laying.date
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "temp.laying.date") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "temp.laying.date")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("temp.laying.date.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="temp.laying.date", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##pcpt.60d
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "pcpt.60d") 
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "pcpt.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("pcpt.60d.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="pcpt.60d: Precipitation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##shadow.days
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "shadow.days") ##find the cumulative importance for each gf.model output. 
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "shadow.days")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("shadow.days.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="shadow.days", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##day10cm
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "day10cm") ##find the cumulative importance for each gf.model output. Correcting label (see above)
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "day10cm")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("day10cm.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="day10cm: calendar day with 10+cm snow", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##MEM1
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "MEM1") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "MEM1")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("MEM1.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM1", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM2
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "MEM2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "MEM2")


isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("MEM2.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="MEM2", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist
CU.Fst <- cumimp(gf.CHS.TI.Adaptive.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.TI.Neutral.SNPs, "dist")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("dist.CumImp.CHS.TI.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="Distance (km)", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##Create the combined plot in AdobeIllustrator
```



##### CHS.VS Cumulative plot
```
##sol.rad.60d

CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "sol.rad.60d") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "sol.rad.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

ymax=0.0015 #set the ymax so that it's the same for the the plots to be overlaid

pdf("sol.rad.60d.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="sol.rad.60d: mean solar radiation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##temp.laying.date
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "temp.laying.date") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "temp.laying.date")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("temp.laying.date.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="temp.laying.date", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##pcpt.60d
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "day10cm") ##I mislabeled these in the original excel sheet (now corrected). Only for CHS.VS. The rest of the transects are correct
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "day10cm")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("pcpt.60d.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="pcpt.60d: Precipitation 60days after laying", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()



##shadow.days
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "pcpt.60d") ##find the cumulative importance for each gf.model output. Correcting the excel mislabelling here: see above
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "pcpt.60d")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("shadow.days.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="shadow.days", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##day10cm
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "shadow.days") ##find the cumulative importance for each gf.model output. Correcting label (see above)
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "shadow.days")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("day10cm.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="day10cm: calendar day with 10+cm snow", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##MEM1
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "MEM1") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "MEM1")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))

pdf("MEM1.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="MEM1", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##MEM2
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "MEM2") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "MEM2")


isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("MEM2.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="Cumulative Importance", xlab="MEM2", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()


##dist
CU.Fst <- cumimp(gf.CHS.VS.Adaptive.SNPs, "dist") ##find the cumulative importance for each gf.model output
CU.Neutral <- cumimp(gf.CHS.VS.Neutral.SNPs, "dist")

isub.Fst <- seq(1, length(CU.Fst$x), len = pmin(500, length(CU.Fst$x)))
isub.Neutral <- seq(1, length(CU.Neutral$x), len = pmin(500, length(CU.Neutral$x)))


pdf("dist.CumImp.CHS.VS.pdf")
plot(CU.Neutral$x[isub.Neutral], CU.Neutral$y[isub.Neutral], type = "s", ylab ="", xlab="Distance (km)", ylim=c(0,ymax), lty=1, lwd=1.5)
par(new=T)  ##allows you to overlay the plots
plot(CU.Fst$x[isub.Fst], CU.Fst$y[isub.Fst], type = "s", ylab = "", xlab="", ylim=c(0,ymax), lty=3, lwd=1.5, xaxt='n', yaxt='n')
#legend("topleft", leg.txt, col="black", lty=c(1,3,2), lwd=1.5, bty = "n")  ##only in the first box. bty removes border
dev.off()

##Create the combined plot in AdobeIllustrator
```

![alt_txt][Fig2]

[Fig2]:https://cloud.githubusercontent.com/assets/12142475/22886423/13fc2d1c-f1fe-11e6-9f5b-2fbad3e0ca54.png





# 4 Individual SNPs

###First I need to redefine the species.cumulative.plot function that comes with gradient forest:
I've changed the y-limit so everythings on the same axis. And I've definede the imp.var.names. 
I've also change show.overall to FALSE so that I can see the individual SNPs

getAnywhere(species.cumulative.plot)

```
species.cumulative.plot <- function (obj, imp.vars = NULL, imp.vars.names = imp.vars, leg.posn = "topleft", 
    leg.nspecies = 10, legend = TRUE, mfrow = rev(n2mfrow(length(imp.vars) * 
        (show.species + show.overall))), show.species = TRUE, 
    show.overall = FALSE, mar = c(0, 2.1, 1.1, 0), omi = c(0.75, 
        0.75, 0.1, 0.1), common.scale = F, line.ylab = 1, cex.legend = 0.75, 
    ...) 
{
    if (is.null(imp.vars)) 
        imp.vars <- imp.var.names <- c("sol.rad.60d", "temp.laying.date", "pcpt.60d", "shadow.days", "day10cm")
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
            ylim <- c(0,0.02)   ###I'm editing this so that I can plot on one plot
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



Plot of individual SNPs for a specific dataset and environmental variable
```
### CHN
pdf("CHN.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.Adaptive.SNPs)
dev.off()

pdf("CHN.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHN.Neutral.SNPs)
dev.off()

### CZ
pdf("CZ.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.Adaptive.SNPs)
dev.off()

pdf("CZ.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CZ.Neutral.SNPs)
dev.off()


### CHS
pdf("CHS.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.Adaptive.SNPs)
dev.off()

pdf("CHS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.Neutral.SNPs)
dev.off()


### CHS.VS
pdf("CHS.VS.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.Adaptive.SNPs)
dev.off()

pdf("CHS.VS.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.VS.Neutral.SNPs)
dev.off()

### CHS.TI
pdf("CHS.TI.individual.AdaptiveSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.Adaptive.SNPs)
dev.off()

pdf("CHS.TI.individual.NeutralSNPs.GF.pdf")
species.cumulative.plot(gf.CHS.TI.Neutral.SNPs)
dev.off()

```
![alt_txt][Fig3]

[Fig3]:https://cloud.githubusercontent.com/assets/12142475/22886425/19e50622-f1fe-11e6-9b8e-8d91624e3b87.png




