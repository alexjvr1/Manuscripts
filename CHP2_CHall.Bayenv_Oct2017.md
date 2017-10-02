# Bayenv2: CHall data


I want to identify loci associated with environment for the 6 CH datasets. stba is removed from all datasets (due to missing env data),
hence 8 less individuals in CHall, CHS, and CHS.TI. And I've also filtered all the datasets for a Minor allele count of ~5%. 

1. CHall.932.7744.recode.vcf

2. CHN.229.5265.recode.vcf

3. CHS.275.6339.recode.vcf

4. CHS.VS.135.5835.recode.vcf

5. CHS.TI.140.5692.recode.vcf

6. CZ.404.7288.recode.vcf


####Env Variables

Identified in https://github.com/alexjvr1/Manuscripts/blob/master/5.CHP2_CH_LandscapeGenomics.md

I've chosen a different 5 environmental parameters for which to run BayEnv2 and LFMM. Based on the site-specific parameters that were calculated by Josh, and the reduction based on corrolation of 0.8.

In the order that they are in the lfmm.env file

1. shadow.days

2. solar.rad.60d

3. pcpt.60d

4. day10cm

5. temp.laying.date

#Input files:

https://bitbucket.org/tguenther/bayenv2_public/src/8e4039f64d61?at=default

1. Input file with subset of 1000 SNPs for co-variance matrix (pop structure)

2. Input file with all SNPs for the association analysis

3. ENV input file with normalised environmental parameters.

NB: population order in the input files should all be the same.


## 1. Covariance matrix estimation

Use a random set of 1000SNPs for matrix estimation. 

I estimated Fst distributions per locus for each of the six datasets. To make sure I'm excluding loci under selection for the 
estimation of the neutral covariance matrix, I will remove loci in the upper Fst tail from each dataset:

1. Identify neutral loci (ie. exclude upper Fst tail) and select 1000 at random

2. Select 1000 loci from vcf file

3. Convert vcf to bayenv2 format

4. Calculate covariance matrix1 + 2


#### CHall

First calculate the per locus Fst for CHall. Draw the distribution and decide on a cut-off for the Fst to exclude any outliers from the dataset. 
```
library(adegenet)
library(hierfstat)
CHall <- read.structure("CHall.932.7744.str")
CHall.popnames <- read.table("CHall.popnames", header=F)
CHall.pops <- gsub("_\\d+", "", CHall.popnames$V1)
CHall.pops
CHall.pops <- as.data.frame(CHall.pops)
CHall.popnames$V2 <- CHall.pops$CHall.pops
colnames(CHall.popnames) <- c("indiv", "site")
CHall@pop <- CHall.popnames$site
CHall@pop
hier.CHall <- genind2hierfstat(CHall)
stats.hier.CHall <- basic.stats(hier.CHall)
perloc.stats.CHall <- stats.hier.CHall$perloc

hist(perloc.stats.CHall$Fst, xlim=c(-0.2, 1.0), breaks=120)
```


1. Fst<0.6
```
CHall.FstNeutral <- subset(perloc.stats.CHall, Fst<0.6)
summary(CHall.FstNeutral)
CHall932.FstNeutral.names <- row.names(CHall.FstNeutral)  ##get the locus names
CHall932.FstNeutral.names <- gsub("X", "", CHall932.FstNeutral.names)
CHall932.FstNeutral.names <- gsub("\\.", ":", CHall932.FstNeutral.names)
CHall932.FstNeutral.names 
CHall932.FstNeutral.names <- as.data.frame(CHall932.FstNeutral.names )
CHall932.1000loci.names <- CHall932.FstNeutral.names[sample(nrow(CHall932.FstNeutral.names), 1000), ]
write.table(CHall932.1000loci.names, "CHall.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHall.932.7744.recode.vcf --snps CHall.1000loci.names --recode --recode-INFO-all --out CHall932.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHall932.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHall/CHall932.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHall/CHall.pops

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHall/CHall932.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHall/CHall932.1000.loci



#### CHN

```
library(adegenet)
library(hierfstat)
CHN <- read.structure("CHN.229.5265.str")

CHN.popnames <- read.table("CHN.229.popnames") ##popnames from bcftools query -l *vcf 
CHN.pops <- gsub("_\\d+", "", CHN.popnames$V1) ##get popnames on their own
CHN.pops <- as.data.frame(CHN.pops)
CHN.popnames$V2 <- CHN.pops$CHN.pops
colnames(CHN.popnames) <- c("indiv", "site")

CHN@pop <- CHN.popnames$site  ##popnames to genind object
hier.CHN <- genind2hierfstat(CHN)  #change to hierfstat
stats.hier.CHN <- basic.stats(hier.CHN) #calculate basic stats
perloc.stats.CHN <- stats.hier.CHN$perloc  #per locus fst is found here

hist(perloc.stats.CHN$Fst, xlim=c(-0.2, 1.0), breaks=120)  ##decide on cut-off based on the histogram
```


1. Fst<0.4
```
CHN.FstNeutral <- subset(perloc.stats.CHN, Fst<0.4)
summary(CHN.FstNeutral)
CHN229.FstNeutral.names <- row.names(CHN.FstNeutral)  ##get the locus names
CHN229.FstNeutral.names <- gsub("X", "", CHN229.FstNeutral.names)
CHN229.FstNeutral.names <- gsub("\\.", ":", CHN229.FstNeutral.names)
CHN229.FstNeutral.names 
CHN229.FstNeutral.names <- as.data.frame(CHN229.FstNeutral.names )
CHN229.1000loci.names <- CHN229.FstNeutral.names[sample(nrow(CHN229.FstNeutral.names), 1000), ]
write.table(CHN229.1000loci.names, "CHN229.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHN.229.5265.recode.vcf --snps CHN229.1000loci.names --recode --recode-INFO-all --out CHN229.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHN/CHN229.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHN/CHN229.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHN/CHN.pops



#### CHS

```
library(adegenet)
library(hierfstat)
CHS <- read.structure("CHS.275.6339.str")
CHS.popnames <- read.table("CHS.popnames", header=F)
CHS.pops <- gsub("_\\d+", "", CHS.popnames$V1)
CHS.pops <- as.data.frame(CHS.pops)
CHS.popnames$V2 <- CHS.pops$CHS.pops
colnames(CHS.popnames) <- c("indiv", "site")
head(CHS.popnames)
CHS@pop <- CHS.popnames$site
CHS@pop
hier.CHS <- genind2hierfstat(CHS)
stats.hier.CHS <- basic.stats(hier.CHS)
perloc.stats.CHS <- stats.hier.CHS$perloc
hist(perloc.stats.CHS$Fst, xlim=c(-0.2, 1.0), breaks=120)
```



1. Fst<0.5
```
CHS.FstNeutral <- subset(perloc.stats.CHS, Fst<0.5)
summary(CHS.FstNeutral)
CHS275.FstNeutral.names <- row.names(CHS.FstNeutral)  ##get the locus names
CHS275.FstNeutral.names <- gsub("X", "", CHS275.FstNeutral.names)
CHS275.FstNeutral.names <- gsub("\\.", ":", CHS275.FstNeutral.names)
CHS275.FstNeutral.names 
CHS275.FstNeutral.names <- as.data.frame(CHS275.FstNeutral.names )
CHS275.1000loci.names <- CHS275.FstNeutral.names[sample(nrow(CHS275.FstNeutral.names), 1000), ]
write.table(CHS275.1000loci.names, "CHS275.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.275.6339.recode.vcf --snps CHS275.1000loci.names --recode --recode-INFO-all --out CHS275.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHS/CHS275.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHS/CHS275.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHS/CHS.pops


#### CHS.TI

```
library(hierfstat)
library(adegenet)
CHS.TI <- read.structure("CHS.TI.140.5692.str")
CHS.TI.popnames <- read.table("CHS.TI.popnames", header=F)
CHS.TI.pops <- gsub("_\\d+", "", CHS.TI.popnames$V1)
CHS.TI.pops <- as.data.frame(CHS.TI.pops)
CHS.TI.popnames$V2 <- CHS.TI.pops$CHS.TI.pops
colnames(CHS.TI.popnames) <- c("indiv", "site")
head(CHS.TI.popnames)
CHS.TI@pop <- CHS.TI.popnames$site
CHS.TI@pop
hier.CHS.TI <- genind2hierfstat(CHS.TI)
stats.hier.CHS.TI <- basic.stats(hier.CHS.TI)
perloc.stats.CHS.TI <- stats.hier.CHS.TI$perloc

hist(perloc.stats.CHS.TI$Fst, xlim=c(-0.2, 1.0), breaks=120)
```

1. Fst<0.5
```
CHS.TI.FstNeutral <- subset(perloc.stats.CHS.TI, Fst<0.5)
summary(CHS.FstNeutral)
CHS.TI140.FstNeutral.names <- row.names(CHS.TI.FstNeutral)  ##get the locus names
CHS.TI140.FstNeutral.names <- gsub("X", "", CHS.TI140.FstNeutral.names)
CHS.TI140.FstNeutral.names <- gsub("\\.", ":", CHS.TI140.FstNeutral.names)
CHS.TI140.FstNeutral.names 
CHS.TI140.FstNeutral.names <- as.data.frame(CHS.TI140.FstNeutral.names )
CHS.TI140.1000loci.names <- CHS.TI140.FstNeutral.names[sample(nrow(CHS.TI140.FstNeutral.names), 1000), ]
write.table(CHS.TI140.1000loci.names, "CHS.TI140.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.TI.148.9608.recode.vcf --snps CHS.TI148.1000loci.names --recode --recode-INFO-all --out CHS.TI148.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHS.TI

#### CHS.VS

```
library(adegenet)
library(hierfstat)
CHS.VS <- read.structure("CHS.VS.135.5835.str")
CHS.VS.popnames <- read.table("CHS.VS.popnames", header=F)
CHS.VS.pops <- gsub("_\\d+", "", CHS.VS.popnames$V1)
CHS.VS.pops <- as.data.frame(CHS.VS.pops)
CHS.VS.popnames$V2 <- CHS.VS.pops$CHS.VS.pops
colnames(CHS.VS.popnames) <- c("indiv","site")
head(CHS.VS.popnames)
CHS.VS@pop <- CHS.VS.popnames$site
hier.CHS.VS <- genind2hierfstat(CHS.VS)
stats.hier.CHS.VS <- basic.stats(hier.CHS.VS)
perloc.stats.CHS.VS <- stats.hier.CHS.VS$perloc

hist(perloc.stats.CHS.VS$Fst, xlim=c(-0.2, 1.0), breaks=120)
```

1. Fst<0.5
```
CHS.VS.FstNeutral <- subset(perloc.stats.CHS.VS, Fst<0.5)
summary(CHS.VS.FstNeutral)
CHS.VS135.FstNeutral.names <- row.names(CHS.VS.FstNeutral)  ##get the locus names
CHS.VS135.FstNeutral.names <- gsub("X", "", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names <- gsub("\\.", ":", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names 
CHS.VS135.FstNeutral.names <- as.data.frame(CHS.VS135.FstNeutral.names )
CHS.VS135.1000loci.names <- CHS.VS135.FstNeutral.names[sample(nrow(CHS.VS135.FstNeutral.names), 1000), ]
write.table(CHS.VS135.1000loci.names, "CHS.VS135.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.VS.135.5835.recode.vcf --snps CHS.VS135.1000loci.names --recode --recode-INFO-all --out CHS.VS135.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CHS.VS




#### CZ

```
library(adegenet)
library(hierfstat)
CZ <- read.structure("CZ.404.7288.str")
CZ.popnames <- read.table("CZ.popnames", header=F)
CZ.pops <- gsub("_\\d+", "", CZ.popnames$V1)
CZ.pops <- as.data.frame(CZ.pops)
CZ.popnames$V2 <- CZ.pops$CZ.pops
colnames(CZ.popnames) <- c("indiv", "site")
CZ@pop <- CZ.popnames$site
CZ@pop
hier.CZ <- genind2hierfstat(CZ)
stats.hier.CZ <- basic.stats(hier.CZ)
perloc.stats.CZ <- stats.hier.CZ$perloc

hist(perloc.stats.CZ$Fst, xlim=c(-0.2, 1.0), breaks=120)
```


1. Fst<0.4
```
CZ404.FstNeutral <- subset(perloc.stats.CZ, Fst<0.4)
summary(CZ404.FstNeutral)
CZ404.FstNeutral.names <- row.names(CZ404.FstNeutral)  ##get the locus names
CZ404.FstNeutral.names <- gsub("X", "", CZ404.FstNeutral.names)
CZ404.FstNeutral.names <- gsub("\\.", ":", CZ404.FstNeutral.names)
CZ404.FstNeutral.names 
CZ404.FstNeutral.names <- as.data.frame(CZ404.FstNeutral.names )
CZ404.1000loci.names <- CZ404.FstNeutral.names[sample(nrow(CZ404.FstNeutral.names), 1000), ]
write.table(CZ404.1000loci.names, "CZ404.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CZ.404.9608.recode.vcf --snps CZ404.1000loci.names --recode --recode-INFO-all --out CZ404.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_inputfiles/MACFiltered_Oct2017/CZ


####4. Calculate covariance matrix

Copy all input files to the bayenv2 directory on the fgcz server: 

/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/MACfilter/MatrixEstimation

Calculate 2 matrices for each dataset. This takes a long time to run. Started at 10:50 on Monday 2 Oct2017
```
#check how far the run has progressed: 
cat CHS.VS.matrix1.out |grep "ITER"

./bayenv2 -i CHall932.1000Neutral.bayenv.txt -p 81 -k 100000 -r $RANDOM > CHall.matrix1.out
./bayenv2 -i CHall932.1000Neutral.bayenv.txt -p 81 -k 100000 -r $RANDOM > CHall.matrix2.out

./bayenv2 -i CHN.229.1000Neutral.bayenv.txt -p 19 -k 100000 -r $RANDOM > CHN229.matrix1.out
./bayenv2 -i CHN.229.1000Neutral.bayenv.txt -p 19 -k 100000 -r $RANDOM > CHN229.matrix2.out

./bayenv2 -i CZ404.1000Neutral.bayenv.txt -p 38 -k 100000 -r $RANDOM > CZ404.matrix1.out
./bayenv2 -i CZ404.1000Neutral.bayenv.txt -p 38 -k 100000 -r $RANDOM > CZ404.matrix2.out

./bayenv2 -i CHS275.1000Neutral.bayenv.txt -p 24 -k 100000 -r $RANDOM > CHS275.matrix1.out
./bayenv2 -i CHS275.1000Neutral.bayenv.txt -p 24 -k 100000 -r $RANDOM > CHS275.matrix2.out

./bayenv2 -i CHS.TI140.1000Neutral.bayenv.txt -p 14 -k 100000 -r $RANDOM > CHS.TI140.matrix1.out
./bayenv2 -i CHS.TI140.1000Neutral.bayenv.txt -p 14 -k 100000 -r $RANDOM > CHS.TI140.matrix2.out

./bayenv2 -i CHS.VS135.1000Neutral.bayenv.txt -p 10 -k 100000 -r $RANDOM > CHS.VS.matrix1.out
./bayenv2 -i CHS.VS135.1000Neutral.bayenv.txt -p 10 -k 100000 -r $RANDOM > CHS.VS.matrix2.out
```


##### Estimate convergence within and between matrices

###### CHall

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CHall
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CHall.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CHall.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CHall.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CHall.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CHall.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CHall.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CHall.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CHall.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.989983 
0.9915622 
0.9929452 


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.9874641 
0.9939042 
0.9927469


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

0.996168 

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```

![alt_txt][CHall.tree]

[CHall.tree]:https://cloud.githubusercontent.com/assets/12142475/23836525/0d7a437a-077a-11e7-8d01-e4582d5a05a9.png

###### CHN

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CHN
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CHN.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CHN.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CHN.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CHN.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CHN.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CHN.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CHN.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CHN.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.817711  
0.9464372  
0.9704053  


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.7877207 
0.9695955 
0.9758337


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

0.9661045 

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```
![alt_txt][CHN.tree]

[CHN.tree]:https://cloud.githubusercontent.com/assets/12142475/23836421/20e1f586-0778-11e7-9577-7e48b87b5ccb.png



###### CHS

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CHS
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CHS.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CHS.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CHS.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CHS.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CHS.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CHS.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CHS.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CHS.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.940169 
0.9899528 
0.991355 


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.950828  
0.9875196  
0.9884865 


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

0.9874994 

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```
![alt_txt][CHS.tree]

[CHS.tree]:https://cloud.githubusercontent.com/assets/12142475/23836579/d39e3d5e-077a-11e7-84be-862d29565441.png



###### CHS.VS

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CHS.VS
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CHS.VS.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CHS.VS.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CHS.VS.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CHS.VS.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CHS.VS.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CHS.VS.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CHS.VS.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CHS.VS.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.5097888 
0.9651997 
0.968304 


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.6663041  
0.9018106 
0.9105807 


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

0.9462934 

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```

![alt_txt][CHS.VS.tree]

[CHS.VS.tree]:https://cloud.githubusercontent.com/assets/12142475/23826814/525b433e-06a5-11e7-84ed-948b4f9f15da.png



###### CHS.TI

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CHS.TI
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CHS.TI.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CHS.TI.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CHS.TI.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CHS.TI.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CHS.TI.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CHS.TI.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CHS.TI.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CHS.TI.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.6511991 
0.9105099 
0.9568654 


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.6663041 
0.9018106 
0.9105807 


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

1.0 (I double checked that I was using the correct matrices). 

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```

######CZ

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac/CZ
```
#linux
#copy matrices 500, 50k, 75k, and 100k from matrix1 and matrix2

#read into R

matrix1.500 <- read.table("CZ.500")
matrix1.500 <- as.matrix(matrix1.500)

matrix1.500.2 <- read.table("CZ.500.2")
matrix1.500.2 <- as.matrix(matrix1.500.2)

matrix1.50k <- read.table("CZ.50k")
matrix1.50k <- as.matrix(matrix1.50k)

matrix1.50k.2 <- read.table("CZ.50k.2")
matrix1.50k.2 <- as.matrix(matrix1.50k.2)

matrix1.100k <- read.table("CZ.100k")
matrix1.100k <- as.matrix(matrix1.100k)

matrix1.100k.2 <- read.table("CZ.100k.2")
matrix1.100k.2 <- as.matrix(matrix1.100k.2)

matrix1.75k <- read.table("CZ.75k")
matrix1.75k <- as.matrix(matrix1.75k)

matrix1.75k.2 <- read.table("CZ.75k.2")
matrix1.75k.2 <- as.matrix(matrix1.75k.2)

#par(mfrow=c(1,2))
#image(matrix1.10)
#image(matrix2.10)
```




```
#Compare distance trees
##calculate a distance
diss500 <- 1-cor(matrix1.500)
diss50k <- 1-cor(matrix1.50k)
diss75k <- 1-cor(matrix1.75k)
diss100k <- 1-cor(matrix1.100k)

diss500.2 <- 1-cor(matrix1.500.2)
diss50k.2 <- 1-cor(matrix1.50k.2)
diss75k.2 <- 1-cor(matrix1.75k.2)
diss100k.2 <- 1-cor(matrix1.100k.2)


##change the df type
diss500 <- as.dist(diss500)
diss50k <- as.dist(diss50k)
diss75k <- as.dist(diss75k)
diss100k <- as.dist(diss100k)

diss500.2 <- as.dist(diss500.2)
diss50k.2 <- as.dist(diss50k.2)
diss75k.2 <- as.dist(diss75k.2)
diss100k.2 <- as.dist(diss100k.2)


##plot and inspect
par(mfrow=c(2,4))

plot(hclust(diss500), main="matrix500", xlab="")
plot(hclust(diss50k), main="matrix50k", xlab="")
plot(hclust(diss75k), main="matrix75k", xlab="")
plot(hclust(diss100k), main="matrix100k", xlab="")

plot(hclust(diss500.2), main="matrix500.2", xlab="")
plot(hclust(diss50k.2), main="matrix50k.2", xlab="")
plot(hclust(diss75k.2), main="matrix75k.2", xlab="")
plot(hclust(diss100k.2), main="matrix100k.2", xlab="")




#Calculate correlation between matrices within a run

#Run1

cor.test(diss500, diss100k)
cor.test(diss50k, diss100k)
cor.test(diss75k, diss100k)

#Pearson's correlation. Results

0.9633814  
0.9742568 
0.9611082 


#Run2

cor.test(diss500.2, diss100k.2)
cor.test(diss50k.2, diss100k.2)
cor.test(diss75k.2, diss100k.2)

0.942456 
0.9802239  
0.9747231 


#Between runs
cor.test(diss100k, diss100k.2)

#Pearson's correlation result: 

0.9793455

##plot
par(mfrow=c(1,2))
image(matrix1.10)
image(matrix2.10)

```

![alt_txt][CZ.tree]

[CZ.tree]:https://cloud.githubusercontent.com/assets/12142475/23836376/5f027792-0777-11e7-96ed-b82309e1c4c3.png


# BayEnv2 Run

Once all the input files have been created, copy them to the fgcz47 server and fun BayEnv2: 

### 1. Normalised .env file 

These have been created before and are already in the BayEnv2 folder on my computer.

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2

### 2. MATRIX file: copy the last matrix from one of the covariance matrices. 

Make sure to copy the file to excel first to check the spacing between all numbers. Also to remove stba from the CHS, CHall, and CHS.VS 
datasets. Copy and nano into a .env file on the server. 

### 3. BayEnv2 input file 

Created using pgdspider. These were all created before. And are in the same folder as above. For code, see: https://github.com/alexjvr1/Manuscripts/blob/master/5.CHP2_CH_LandscapeGenomics.md


Start the various runs: 
```
##1. split the bayenv2 large input into single loci
##2. loop over the input with bayenv2

#CHall
split -d -a 10 -l 2 CHall.932.9608.BayEnv.txt snp_batch
for f in $(ls snp_batch*); do ./bayenv2 -i $f -e CHall.5n.ENV -m CHall.March.MATRIX -k 100000 -p 81 -n 5 -r $RANDOM -t -c -X; done


#CHN
split -d -a 10 -l 2 CHN.229.8951.bayENV2.txt snp_batch
for f in $(ls snp_batch*); do ./bayenv2 -i $f -e ENVIRONFILE.CHN.p19.env -m CHN.March.MATRIX -k 100000 -p 19 -n 5 -r $RANDOM -t -c -X; done

#CHS
split -d -a 10 -l 2 CHS.275.9608.bayenv2.txt snp_batch
for f in $(ls snp_batch*); do ./bayenv2 -i $f -e ENVIRONFILE.CHS.n5.env -m CHS.March.MATRIX -p 24 -k 100000 -n 5 -r $RANDOM -t -c -X; done

#CHS.VS
split -d -a 10 -l 2 CHS.VS.135.8779.bayEnv2.txt snp_batch
for f in $(ls snp_batch*); do ./bayenv2 -i $f -e standardised.5env.env -m CHS.VS.MATRIX -k 100000 -p 10 -n 5 -r $RANDOM -t -c -X; done

#CHS.TI


#CZ
split -d -a 10 -l 2 CZ.404.9528.bayENV2.txt snp_batch
for f in $(ls snp_batch*); do ./bayenv2 -i $f -e ENVIRON.CZ.5n.env -m CZ.March.MATRIX -k 100000 -p 38 -n 6 -r $RANDOM -t -c -X; done
```

Run analysis three times for all datasets. 

Something went wrong with the output here. My output files for the three runs is not always the same. I can figure out which snps are missing with the following command: 

```
##this prints the first line of each snp output file 
awk '{print $1}' CHS.Run1/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run1.names   
awk '{print $1}' CHS.Run2/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run2.names
awk '{print $1}' CHS.Run3/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run3.names

##this compares the files and prints the differences
diff CHS.Run1.names CHS.Run2.names 
diff CHS.Run1.names CHS.Run3.names 

###to print the last word in a line use NF
##awk 'NF>{print $NF}' CHS.Run1/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run1.names   
```

After rerunning several of the analyes here: 

/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHP2

I have to rerun individual loci that seem to skip during the analysis. 

CHS.Run3 needs snp_batch0000004291

CHS.VS.run3 needs snp_batch0000006240

CZ.run2 needs snp_batch0000005442

CZ.run3 needs snp_batch0000001203

Once this was done, I could concatenate the results as before 


I need to calculate convergence of BF, p, and ranking 1000 loci between the 3 runs for each of the 5 environmental variables. (i.e. 45 correlations)


In R:
```
###CHall
BF.CHall.Run1 <- read.table("CHall.Run1/bf_environ.CHall.5n.ENV", header=F)
colnames(BF.CHall.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHall.Run2 <- read.table("CHall.Run2/bf_environ.CHall.5n.ENV", header=F)
colnames(BF.CHall.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHall.Run3 <- read.table("CHall.Run3/bf_environ.CHall.5n.ENV", header=F)
colnames(BF.CHall.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")


##CHN
BF.CHN.Run1 <- read.table("CHN.Run1/bf_environ.ENVIRONFILE.CHN.p19.env", header=F)
colnames(BF.CHN.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHN.Run2 <- read.table("CHN.Run2/bf_environ.ENVIRONFILE.CHN.p19.env", header=F)
colnames(BF.CHN.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHN.Run3 <- read.table("CHN.Run3/bf_environ.ENVIRONFILE.CHN.p19.env", header=F)
colnames(BF.CHN.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")


###CHS
BF.CHS.Run1 <- read.table("CHS.Run1/bf_environ.ENVIRONFILE.CHS.n5.env", header=F)
colnames(BF.CHS.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.Run2 <- read.table("CHS.Run2/bf_environ.ENVIRONFILE.CHS.n5.env", header=F)
colnames(BF.CHS.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.Run3 <- read.table("CHS.Run3/bf_environ.ENVIRONFILE.CHS.n5.env", header=F)
colnames(BF.CHS.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")


###CZ
BF.CZ.Run1 <- read.table("CZ.Run1/bf_environ.ENVIRON.CZ.5n.env", header=F)
colnames(BF.CZ.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CZ.Run2 <- read.table("CZ.Run2/bf_environ.ENVIRON.CZ.5n.env", header=F)
colnames(BF.CZ.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CZ.Run3 <- read.table("CZ.Run3/bf_environ.ENVIRON.CZ.5n.env", header=F)
colnames(BF.CZ.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

###CHS.VS
BF.CHS.VS.Run1 <- read.table("CHS.VS.Run1/bf_environ.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)
colnames(BF.CHS.VS.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.VS.Run2 <- read.table("CHS.VS.Run2/bf_environ.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)
colnames(BF.CHS.VS.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.VS.Run3 <- read.table("CHS.VS.Run3/bf_environ.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)
colnames(BF.CHS.VS.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

###CHS.TI
BF.CHS.TI.Run1 <- read.table("CHS.TI.Run1/bf_environ.standardised.5env.env", header=F)
colnames(BF.CHS.TI.Run1) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.TI.Run2 <- read.table("CHS.TI.Run2/bf_environ.standardised.5env.env", header=F)
colnames(BF.CHS.TI.Run2) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")

BF.CHS.TI.Run3 <- read.table("CHS.TI.Run3/bf_environ.standardised.5env.env", header=F)
colnames(BF.CHS.TI.Run3) <- c("snp", "rad.BF", "rad.rho", "rad.r", "shadow.days.BF", "shadow.days.rho", "shadow.days.r", "temp.BF", "temp.rho", "temp.r", "pcpt.BF", "pcpt.rho", "pcpt.r", "day.10cm.BF", "day.10cm.rho", "day.10cm.r")
```




Compare ranking SNPs 5% at a time

####CHall

1. rad
```
CHall.rad.run1.top100 <- BF.CHall.Run1[order(-BF.CHall.Run1$rad.rho),]  ##descending order by rho
CHall.rad.run2.top100 <- BF.CHall.Run2[order(-BF.CHall.Run2$rad.rho),]  ##descending order by rho
CHall.rad.run3.top100 <- BF.CHall.Run3[order(-BF.CHall.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHall.rad
###############
```
CHall.rad.run1.top100.set1 <- head(CHall.rad.run1.top100, 100)
CHall.rad.run2.top100.set1 <- head(CHall.rad.run2.top100, 100)
CHall.rad.run3.top100.set1 <- head(CHall.rad.run3.top100, 100)
CHall.rad.12 <- Reduce(intersect, list(CHall.rad.run2.top100.set1$snp, CHall.rad.run1.top100.set1$snp))
CHall.rad.13 <- Reduce(intersect, list(CHall.rad.run3.top100.set1$snp, CHall.rad.run1.top100.set1$snp))
CHall.rad.23 <- Reduce(intersect, list(CHall.rad.run2.top100.set1$snp, CHall.rad.run3.top100.set1$snp))
CHall.rad.123 <- Reduce(intersect, list(CHall.rad.run2.top100.set1$snp, CHall.rad.run3.top100.set1$snp, CHall.rad.run1.top100.set1$snp))

summary(CHall.rad.123)
#70
```

##overlap in 6-10% of outlier loci for CHall.rad
###############
```
CHall.rad.run1.top100.set2 <- CHall.rad.run1.top100[101:200,]
CHall.rad.run2.top100.set2 <- CHall.rad.run2.top100[101:200,]
CHall.rad.run3.top100.set2 <- CHall.rad.run3.top100[101:200,]
CHall.rad.12 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run1.top100.set2$snp))
CHall.rad.13 <- Reduce(intersect, list(CHall.rad.run3.top100.set2$snp, CHall.rad.run1.top100.set2$snp))
CHall.rad.23 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run3.top100.set2$snp))
CHall.rad.123 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run3.top100.set2$snp, CHall.rad.run1.top100.set2$snp))

summary(CHall.rad.123)
#26
```

##overlap in 11-15% of outlier loci for CHall.rad
###############
```
CHall.rad.run1.top100.set2 <- CHall.rad.run1.top100[201:300,]
CHall.rad.run2.top100.set2 <- CHall.rad.run2.top100[201:300,]
CHall.rad.run3.top100.set2 <- CHall.rad.run3.top100[201:300,]
CHall.rad.12 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run1.top100.set2$snp))
CHall.rad.13 <- Reduce(intersect, list(CHall.rad.run3.top100.set2$snp, CHall.rad.run1.top100.set2$snp))
CHall.rad.23 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run3.top100.set2$snp))
CHall.rad.123 <- Reduce(intersect, list(CHall.rad.run2.top100.set2$snp, CHall.rad.run3.top100.set2$snp, CHall.rad.run1.top100.set2$snp))

summary(CHall.rad.123)
#10
```

2. shadow.days

```
CHall.shadow.days.run1.top100 <- BF.CHall.Run1[order(-BF.CHall.Run1$shadow.days.rho),]  ##descending order by rho
CHall.shadow.days.run2.top100 <- BF.CHall.Run2[order(-BF.CHall.Run2$shadow.days.rho),]  ##descending order by rho
CHall.shadow.days.run3.top100 <- BF.CHall.Run3[order(-BF.CHall.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHall.shadow.days
###############
```
CHall.shadow.days.run1.top100.set1 <- head(CHall.shadow.days.run1.top100, 100)
CHall.shadow.days.run2.top100.set1 <- head(CHall.shadow.days.run2.top100, 100)
CHall.shadow.days.run3.top100.set1 <- head(CHall.shadow.days.run3.top100, 100)
CHall.shadow.days.12 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set1$snp, CHall.shadow.days.run1.top100.set1$snp))
CHall.shadow.days.13 <- Reduce(intersect, list(CHall.shadow.days.run3.top100.set1$snp, CHall.shadow.days.run1.top100.set1$snp))
CHall.shadow.days.23 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set1$snp, CHall.shadow.days.run3.top100.set1$snp))
CHall.shadow.days.123 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set1$snp, CHall.shadow.days.run3.top100.set1$snp, CHall.shadow.days.run1.top100.set1$snp))

summary(CHall.shadow.days.123)
#69
```

##overlap in 6-10% of outlier loci for CHall.shadow.days
###############
```
CHall.shadow.days.run1.top100.set2 <- CHall.shadow.days.run1.top100[101:200,]
CHall.shadow.days.run2.top100.set2 <- CHall.shadow.days.run2.top100[101:200,]
CHall.shadow.days.run3.top100.set2 <- CHall.shadow.days.run3.top100[101:200,]
CHall.shadow.days.12 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))
CHall.shadow.days.13 <- Reduce(intersect, list(CHall.shadow.days.run3.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))
CHall.shadow.days.23 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run3.top100.set2$snp))
CHall.shadow.days.123 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run3.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))

summary(CHall.shadow.days.123)
#29
```

##overlap in 11-15% of outlier loci for CHall.shadow.days
###############
```
CHall.shadow.days.run1.top100.set2 <- CHall.shadow.days.run1.top100[201:300,]
CHall.shadow.days.run2.top100.set2 <- CHall.shadow.days.run2.top100[201:300,]
CHall.shadow.days.run3.top100.set2 <- CHall.shadow.days.run3.top100[201:300,]
CHall.shadow.days.12 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))
CHall.shadow.days.13 <- Reduce(intersect, list(CHall.shadow.days.run3.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))
CHall.shadow.days.23 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run3.top100.set2$snp))
CHall.shadow.days.123 <- Reduce(intersect, list(CHall.shadow.days.run2.top100.set2$snp, CHall.shadow.days.run3.top100.set2$snp, CHall.shadow.days.run1.top100.set2$snp))

summary(CHall.shadow.days.123)
#15
```

3. temp

```
CHall.temp.run1.top100 <- BF.CHall.Run1[order(-BF.CHall.Run1$temp.rho),]  ##descending order by rho
CHall.temp.run2.top100 <- BF.CHall.Run2[order(-BF.CHall.Run2$temp.rho),]  ##descending order by rho
CHall.temp.run3.top100 <- BF.CHall.Run3[order(-BF.CHall.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHall.temp
###############
```
CHall.temp.run1.top100.set1 <- head(CHall.temp.run1.top100, 100)
CHall.temp.run2.top100.set1 <- head(CHall.temp.run2.top100, 100)
CHall.temp.run3.top100.set1 <- head(CHall.temp.run3.top100, 100)
CHall.temp.12 <- Reduce(intersect, list(CHall.temp.run2.top100.set1$snp, CHall.temp.run1.top100.set1$snp))
CHall.temp.13 <- Reduce(intersect, list(CHall.temp.run3.top100.set1$snp, CHall.temp.run1.top100.set1$snp))
CHall.temp.23 <- Reduce(intersect, list(CHall.temp.run2.top100.set1$snp, CHall.temp.run3.top100.set1$snp))
CHall.temp.123 <- Reduce(intersect, list(CHall.temp.run2.top100.set1$snp, CHall.temp.run3.top100.set1$snp, CHall.temp.run1.top100.set1$snp))

summary(CHall.temp.123)
#77
```

##overlap in 6-10% of outlier loci for CHall.temp
###############
```
CHall.temp.run1.top100.set2 <- CHall.temp.run1.top100[101:200,]
CHall.temp.run2.top100.set2 <- CHall.temp.run2.top100[101:200,]
CHall.temp.run3.top100.set2 <- CHall.temp.run3.top100[101:200,]
CHall.temp.12 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run1.top100.set2$snp))
CHall.temp.13 <- Reduce(intersect, list(CHall.temp.run3.top100.set2$snp, CHall.temp.run1.top100.set2$snp))
CHall.temp.23 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run3.top100.set2$snp))
CHall.temp.123 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run3.top100.set2$snp, CHall.temp.run1.top100.set2$snp))
summary(CHall.temp.123)
#34
```

##overlap in 11-15% of outlier loci for CHall.temp
###############
```
CHall.temp.run1.top100.set2 <- CHall.temp.run1.top100[201:300,]
CHall.temp.run2.top100.set2 <- CHall.temp.run2.top100[201:300,]
CHall.temp.run3.top100.set2 <- CHall.temp.run3.top100[201:300,]
CHall.temp.12 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run1.top100.set2$snp))
CHall.temp.13 <- Reduce(intersect, list(CHall.temp.run3.top100.set2$snp, CHall.temp.run1.top100.set2$snp))
CHall.temp.23 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run3.top100.set2$snp))
CHall.temp.123 <- Reduce(intersect, list(CHall.temp.run2.top100.set2$snp, CHall.temp.run3.top100.set2$snp, CHall.temp.run1.top100.set2$snp))

summary(CHall.temp.123)
#12
```

4. pcpt

```
CHall.pcpt.run1.top100 <- BF.CHall.Run1[order(-BF.CHall.Run1$pcpt.rho),]  ##descending order by rho
CHall.pcpt.run2.top100 <- BF.CHall.Run2[order(-BF.CHall.Run2$pcpt.rho),]  ##descending order by rho
CHall.pcpt.run3.top100 <- BF.CHall.Run3[order(-BF.CHall.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHall.pcpt
###############
```
CHall.pcpt.run1.top100.set1 <- head(CHall.pcpt.run1.top100, 100)
CHall.pcpt.run2.top100.set1 <- head(CHall.pcpt.run2.top100, 100)
CHall.pcpt.run3.top100.set1 <- head(CHall.pcpt.run3.top100, 100)
CHall.pcpt.12 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set1$snp, CHall.pcpt.run1.top100.set1$snp))
CHall.pcpt.13 <- Reduce(intersect, list(CHall.pcpt.run3.top100.set1$snp, CHall.pcpt.run1.top100.set1$snp))
CHall.pcpt.23 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set1$snp, CHall.pcpt.run3.top100.set1$snp))
CHall.pcpt.123 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set1$snp, CHall.pcpt.run3.top100.set1$snp, CHall.pcpt.run1.top100.set1$snp))

summary(CHall.pcpt.123)
#68
```

##overlap in 6-10% of outlier loci for CHall.pcpt
###############
```
CHall.pcpt.run1.top100.set2 <- CHall.pcpt.run1.top100[101:200,]
CHall.pcpt.run2.top100.set2 <- CHall.pcpt.run2.top100[101:200,]
CHall.pcpt.run3.top100.set2 <- CHall.pcpt.run3.top100[101:200,]
CHall.pcpt.12 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))
CHall.pcpt.13 <- Reduce(intersect, list(CHall.pcpt.run3.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))
CHall.pcpt.23 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run3.top100.set2$snp))
CHall.pcpt.123 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run3.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))

summary(CHall.pcpt.123)
#25
```

##overlap in 11-15% of outlier loci for CHall.pcpt
###############
```
CHall.pcpt.run1.top100.set2 <- CHall.pcpt.run1.top100[201:300,]
CHall.pcpt.run2.top100.set2 <- CHall.pcpt.run2.top100[201:300,]
CHall.pcpt.run3.top100.set2 <- CHall.pcpt.run3.top100[201:300,]
CHall.pcpt.12 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))
CHall.pcpt.13 <- Reduce(intersect, list(CHall.pcpt.run3.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))
CHall.pcpt.23 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run3.top100.set2$snp))
CHall.pcpt.123 <- Reduce(intersect, list(CHall.pcpt.run2.top100.set2$snp, CHall.pcpt.run3.top100.set2$snp, CHall.pcpt.run1.top100.set2$snp))

summary(CHall.pcpt.123)
#12
```

5. day.10cm
```
CHall.day.10cm.run1.top100 <- BF.CHall.Run1[order(-BF.CHall.Run1$day.10cm.rho),]  ##descending order by rho
CHall.day.10cm.run2.top100 <- BF.CHall.Run2[order(-BF.CHall.Run2$day.10cm.rho),]  ##descending order by rho
CHall.day.10cm.run3.top100 <- BF.CHall.Run3[order(-BF.CHall.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHall.day.10cm
###############
```
CHall.day.10cm.run1.top100.set1 <- head(CHall.day.10cm.run1.top100, 100)
CHall.day.10cm.run2.top100.set1 <- head(CHall.day.10cm.run2.top100, 100)
CHall.day.10cm.run3.top100.set1 <- head(CHall.day.10cm.run3.top100, 100)
CHall.day.10cm.12 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set1$snp, CHall.day.10cm.run1.top100.set1$snp))
CHall.day.10cm.13 <- Reduce(intersect, list(CHall.day.10cm.run3.top100.set1$snp, CHall.day.10cm.run1.top100.set1$snp))
CHall.day.10cm.23 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set1$snp, CHall.day.10cm.run3.top100.set1$snp))
CHall.day.10cm.123 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set1$snp, CHall.day.10cm.run3.top100.set1$snp, CHall.day.10cm.run1.top100.set1$snp))

summary(CHall.day.10cm.123)
#73
```

##overlap in 6-10% of outlier loci for CHall.day.10cm
###############
```
CHall.day.10cm.run1.top100.set2 <- CHall.day.10cm.run1.top100[101:200,]
CHall.day.10cm.run2.top100.set2 <- CHall.day.10cm.run2.top100[101:200,]
CHall.day.10cm.run3.top100.set2 <- CHall.day.10cm.run3.top100[101:200,]
CHall.day.10cm.12 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))
CHall.day.10cm.13 <- Reduce(intersect, list(CHall.day.10cm.run3.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))
CHall.day.10cm.23 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run3.top100.set2$snp))
CHall.day.10cm.123 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run3.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))
summary(CHall.day.10cm.123)
#31
```

##overlap in 11-15% of outlier loci for CHall.day.10cm
###############
```
CHall.day.10cm.run1.top100.set2 <- CHall.day.10cm.run1.top100[201:300,]
CHall.day.10cm.run2.top100.set2 <- CHall.day.10cm.run2.top100[201:300,]
CHall.day.10cm.run3.top100.set2 <- CHall.day.10cm.run3.top100[201:300,]
CHall.day.10cm.12 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))
CHall.day.10cm.13 <- Reduce(intersect, list(CHall.day.10cm.run3.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))
CHall.day.10cm.23 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run3.top100.set2$snp))
CHall.day.10cm.123 <- Reduce(intersect, list(CHall.day.10cm.run2.top100.set2$snp, CHall.day.10cm.run3.top100.set2$snp, CHall.day.10cm.run1.top100.set2$snp))

summary(CHall.day.10cm.123)
#10
```

####CHN

1. rad
```
CHN.rad.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$rad.rho),]  ##descending order by rho
CHN.rad.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$rad.rho),]  ##descending order by rho
CHN.rad.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.rad
###############
```
CHN.rad.run1.top100.set1 <- head(CHN.rad.run1.top100, 100)
CHN.rad.run2.top100.set1 <- head(CHN.rad.run2.top100, 100)
CHN.rad.run3.top100.set1 <- head(CHN.rad.run3.top100, 100)
CHN.rad.12 <- Reduce(intersect, list(CHN.rad.run2.top100.set1$snp, CHN.rad.run1.top100.set1$snp))
CHN.rad.13 <- Reduce(intersect, list(CHN.rad.run3.top100.set1$snp, CHN.rad.run1.top100.set1$snp))
CHN.rad.23 <- Reduce(intersect, list(CHN.rad.run2.top100.set1$snp, CHN.rad.run3.top100.set1$snp))
CHN.rad.123 <- Reduce(intersect, list(CHN.rad.run2.top100.set1$snp, CHN.rad.run3.top100.set1$snp, CHN.rad.run1.top100.set1$snp))

summary(CHN.rad.123)
#59
```

##overlap in 6-10% of outlier loci for CHN.rad
###############
```
CHN.rad.run1.top100.set2 <- CHN.rad.run1.top100[101:200,]
CHN.rad.run2.top100.set2 <- CHN.rad.run2.top100[101:200,]
CHN.rad.run3.top100.set2 <- CHN.rad.run3.top100[101:200,]
CHN.rad.12 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run1.top100.set2$snp))
CHN.rad.13 <- Reduce(intersect, list(CHN.rad.run3.top100.set2$snp, CHN.rad.run1.top100.set2$snp))
CHN.rad.23 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run3.top100.set2$snp))
CHN.rad.123 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run3.top100.set2$snp, CHN.rad.run1.top100.set2$snp))

summary(CHN.rad.123)
#17
```

##overlap in 11-15% of outlier loci for CHN.rad
###############
```
CHN.rad.run1.top100.set2 <- CHN.rad.run1.top100[201:300,]
CHN.rad.run2.top100.set2 <- CHN.rad.run2.top100[201:300,]
CHN.rad.run3.top100.set2 <- CHN.rad.run3.top100[201:300,]
CHN.rad.12 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run1.top100.set2$snp))
CHN.rad.13 <- Reduce(intersect, list(CHN.rad.run3.top100.set2$snp, CHN.rad.run1.top100.set2$snp))
CHN.rad.23 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run3.top100.set2$snp))
CHN.rad.123 <- Reduce(intersect, list(CHN.rad.run2.top100.set2$snp, CHN.rad.run3.top100.set2$snp, CHN.rad.run1.top100.set2$snp))

summary(CHN.rad.123)
#7
```

2. shadow.days

```
CHN.shadow.days.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$shadow.days.rho),]  ##descending order by rho
CHN.shadow.days.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$shadow.days.rho),]  ##descending order by rho
CHN.shadow.days.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.shadow.days
###############
```
CHN.shadow.days.run1.top100.set1 <- head(CHN.shadow.days.run1.top100, 100)
CHN.shadow.days.run2.top100.set1 <- head(CHN.shadow.days.run2.top100, 100)
CHN.shadow.days.run3.top100.set1 <- head(CHN.shadow.days.run3.top100, 100)
CHN.shadow.days.12 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set1$snp, CHN.shadow.days.run1.top100.set1$snp))
CHN.shadow.days.13 <- Reduce(intersect, list(CHN.shadow.days.run3.top100.set1$snp, CHN.shadow.days.run1.top100.set1$snp))
CHN.shadow.days.23 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set1$snp, CHN.shadow.days.run3.top100.set1$snp))
CHN.shadow.days.123 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set1$snp, CHN.shadow.days.run3.top100.set1$snp, CHN.shadow.days.run1.top100.set1$snp))

summary(CHN.shadow.days.123)
#35
```

##overlap in 6-10% of outlier loci for CHN.shadow.days
###############
```
CHN.shadow.days.run1.top100.set2 <- CHN.shadow.days.run1.top100[101:200,]
CHN.shadow.days.run2.top100.set2 <- CHN.shadow.days.run2.top100[101:200,]
CHN.shadow.days.run3.top100.set2 <- CHN.shadow.days.run3.top100[101:200,]
CHN.shadow.days.12 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))
CHN.shadow.days.13 <- Reduce(intersect, list(CHN.shadow.days.run3.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))
CHN.shadow.days.23 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run3.top100.set2$snp))
CHN.shadow.days.123 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run3.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))

summary(CHN.shadow.days.123)
#3
```

##overlap in 11-15% of outlier loci for CHN.shadow.days
###############
```
CHN.shadow.days.run1.top100.set2 <- CHN.shadow.days.run1.top100[201:300,]
CHN.shadow.days.run2.top100.set2 <- CHN.shadow.days.run2.top100[201:300,]
CHN.shadow.days.run3.top100.set2 <- CHN.shadow.days.run3.top100[201:300,]
CHN.shadow.days.12 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))
CHN.shadow.days.13 <- Reduce(intersect, list(CHN.shadow.days.run3.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))
CHN.shadow.days.23 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run3.top100.set2$snp))
CHN.shadow.days.123 <- Reduce(intersect, list(CHN.shadow.days.run2.top100.set2$snp, CHN.shadow.days.run3.top100.set2$snp, CHN.shadow.days.run1.top100.set2$snp))

summary(CHN.shadow.days.123)
#2
```

3. temp

```
CHN.temp.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$temp.rho),]  ##descending order by rho
CHN.temp.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$temp.rho),]  ##descending order by rho
CHN.temp.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.temp
###############
```
CHN.temp.run1.top100.set1 <- head(CHN.temp.run1.top100, 100)
CHN.temp.run2.top100.set1 <- head(CHN.temp.run2.top100, 100)
CHN.temp.run3.top100.set1 <- head(CHN.temp.run3.top100, 100)
CHN.temp.12 <- Reduce(intersect, list(CHN.temp.run2.top100.set1$snp, CHN.temp.run1.top100.set1$snp))
CHN.temp.13 <- Reduce(intersect, list(CHN.temp.run3.top100.set1$snp, CHN.temp.run1.top100.set1$snp))
CHN.temp.23 <- Reduce(intersect, list(CHN.temp.run2.top100.set1$snp, CHN.temp.run3.top100.set1$snp))
CHN.temp.123 <- Reduce(intersect, list(CHN.temp.run2.top100.set1$snp, CHN.temp.run3.top100.set1$snp, CHN.temp.run1.top100.set1$snp))

summary(CHN.temp.123)
#35
```

##overlap in 6-10% of outlier loci for CHN.temp
###############
```
CHN.temp.run1.top100.set2 <- CHN.temp.run1.top100[101:200,]
CHN.temp.run2.top100.set2 <- CHN.temp.run2.top100[101:200,]
CHN.temp.run3.top100.set2 <- CHN.temp.run3.top100[101:200,]
CHN.temp.12 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run1.top100.set2$snp))
CHN.temp.13 <- Reduce(intersect, list(CHN.temp.run3.top100.set2$snp, CHN.temp.run1.top100.set2$snp))
CHN.temp.23 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run3.top100.set2$snp))
CHN.temp.123 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run3.top100.set2$snp, CHN.temp.run1.top100.set2$snp))
summary(CHN.temp.123)
#4
```

##overlap in 11-15% of outlier loci for CHN.temp
###############
```
CHN.temp.run1.top100.set2 <- CHN.temp.run1.top100[201:300,]
CHN.temp.run2.top100.set2 <- CHN.temp.run2.top100[201:300,]
CHN.temp.run3.top100.set2 <- CHN.temp.run3.top100[201:300,]
CHN.temp.12 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run1.top100.set2$snp))
CHN.temp.13 <- Reduce(intersect, list(CHN.temp.run3.top100.set2$snp, CHN.temp.run1.top100.set2$snp))
CHN.temp.23 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run3.top100.set2$snp))
CHN.temp.123 <- Reduce(intersect, list(CHN.temp.run2.top100.set2$snp, CHN.temp.run3.top100.set2$snp, CHN.temp.run1.top100.set2$snp))

summary(CHN.temp.123)
#4
```

4. pcpt

```
CHN.pcpt.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$pcpt.rho),]  ##descending order by rho
CHN.pcpt.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$pcpt.rho),]  ##descending order by rho
CHN.pcpt.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.pcpt
###############
```
CHN.pcpt.run1.top100.set1 <- head(CHN.pcpt.run1.top100, 100)
CHN.pcpt.run2.top100.set1 <- head(CHN.pcpt.run2.top100, 100)
CHN.pcpt.run3.top100.set1 <- head(CHN.pcpt.run3.top100, 100)
CHN.pcpt.12 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set1$snp, CHN.pcpt.run1.top100.set1$snp))
CHN.pcpt.13 <- Reduce(intersect, list(CHN.pcpt.run3.top100.set1$snp, CHN.pcpt.run1.top100.set1$snp))
CHN.pcpt.23 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set1$snp, CHN.pcpt.run3.top100.set1$snp))
CHN.pcpt.123 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set1$snp, CHN.pcpt.run3.top100.set1$snp, CHN.pcpt.run1.top100.set1$snp))

summary(CHN.pcpt.123)
#48
```

##overlap in 6-10% of outlier loci for CHN.pcpt
###############
```
CHN.pcpt.run1.top100.set2 <- CHN.pcpt.run1.top100[101:200,]
CHN.pcpt.run2.top100.set2 <- CHN.pcpt.run2.top100[101:200,]
CHN.pcpt.run3.top100.set2 <- CHN.pcpt.run3.top100[101:200,]
CHN.pcpt.12 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))
CHN.pcpt.13 <- Reduce(intersect, list(CHN.pcpt.run3.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))
CHN.pcpt.23 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run3.top100.set2$snp))
CHN.pcpt.123 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run3.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))

summary(CHN.pcpt.123)
#9
```

##overlap in 11-15% of outlier loci for CHN.pcpt
###############
```
CHN.pcpt.run1.top100.set2 <- CHN.pcpt.run1.top100[201:300,]
CHN.pcpt.run2.top100.set2 <- CHN.pcpt.run2.top100[201:300,]
CHN.pcpt.run3.top100.set2 <- CHN.pcpt.run3.top100[201:300,]
CHN.pcpt.12 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))
CHN.pcpt.13 <- Reduce(intersect, list(CHN.pcpt.run3.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))
CHN.pcpt.23 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run3.top100.set2$snp))
CHN.pcpt.123 <- Reduce(intersect, list(CHN.pcpt.run2.top100.set2$snp, CHN.pcpt.run3.top100.set2$snp, CHN.pcpt.run1.top100.set2$snp))

summary(CHN.pcpt.123)
#5
```

5. day.10cm
```
CHN.day.10cm.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$day.10cm.rho),]  ##descending order by rho
CHN.day.10cm.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$day.10cm.rho),]  ##descending order by rho
CHN.day.10cm.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.day.10cm
###############
```
CHN.day.10cm.run1.top100.set1 <- head(CHN.day.10cm.run1.top100, 100)
CHN.day.10cm.run2.top100.set1 <- head(CHN.day.10cm.run2.top100, 100)
CHN.day.10cm.run3.top100.set1 <- head(CHN.day.10cm.run3.top100, 100)
CHN.day.10cm.12 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set1$snp, CHN.day.10cm.run1.top100.set1$snp))
CHN.day.10cm.13 <- Reduce(intersect, list(CHN.day.10cm.run3.top100.set1$snp, CHN.day.10cm.run1.top100.set1$snp))
CHN.day.10cm.23 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set1$snp, CHN.day.10cm.run3.top100.set1$snp))
CHN.day.10cm.123 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set1$snp, CHN.day.10cm.run3.top100.set1$snp, CHN.day.10cm.run1.top100.set1$snp))

summary(CHN.day.10cm.123)
#32
```

##overlap in 6-10% of outlier loci for CHN.day.10cm
###############
```
CHN.day.10cm.run1.top100.set2 <- CHN.day.10cm.run1.top100[101:200,]
CHN.day.10cm.run2.top100.set2 <- CHN.day.10cm.run2.top100[101:200,]
CHN.day.10cm.run3.top100.set2 <- CHN.day.10cm.run3.top100[101:200,]
CHN.day.10cm.12 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))
CHN.day.10cm.13 <- Reduce(intersect, list(CHN.day.10cm.run3.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))
CHN.day.10cm.23 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run3.top100.set2$snp))
CHN.day.10cm.123 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run3.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))
summary(CHN.day.10cm.123)
#5
```

##overlap in 11-15% of outlier loci for CHN.day.10cm
###############
```
CHN.day.10cm.run1.top100.set2 <- CHN.day.10cm.run1.top100[201:300,]
CHN.day.10cm.run2.top100.set2 <- CHN.day.10cm.run2.top100[201:300,]
CHN.day.10cm.run3.top100.set2 <- CHN.day.10cm.run3.top100[201:300,]
CHN.day.10cm.12 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))
CHN.day.10cm.13 <- Reduce(intersect, list(CHN.day.10cm.run3.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))
CHN.day.10cm.23 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run3.top100.set2$snp))
CHN.day.10cm.123 <- Reduce(intersect, list(CHN.day.10cm.run2.top100.set2$snp, CHN.day.10cm.run3.top100.set2$snp, CHN.day.10cm.run1.top100.set2$snp))

summary(CHN.day.10cm.123)
#0
```

####CHS

1. rad
```
CHS.rad.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$rad.rho),]  ##descending order by rho
CHS.rad.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$rad.rho),]  ##descending order by rho
CHS.rad.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.rad
###############
```
CHS.rad.run1.top100.set1 <- head(CHS.rad.run1.top100, 100)
CHS.rad.run2.top100.set1 <- head(CHS.rad.run2.top100, 100)
CHS.rad.run3.top100.set1 <- head(CHS.rad.run3.top100, 100)
CHS.rad.12 <- Reduce(intersect, list(CHS.rad.run2.top100.set1$snp, CHS.rad.run1.top100.set1$snp))
CHS.rad.13 <- Reduce(intersect, list(CHS.rad.run3.top100.set1$snp, CHS.rad.run1.top100.set1$snp))
CHS.rad.23 <- Reduce(intersect, list(CHS.rad.run2.top100.set1$snp, CHS.rad.run3.top100.set1$snp))
CHS.rad.123 <- Reduce(intersect, list(CHS.rad.run2.top100.set1$snp, CHS.rad.run3.top100.set1$snp, CHS.rad.run1.top100.set1$snp))

summary(CHS.rad.123)
#82
```

##overlap in 6-10% of outlier loci for CHS.rad
###############
```
CHS.rad.run1.top100.set2 <- CHS.rad.run1.top100[101:200,]
CHS.rad.run2.top100.set2 <- CHS.rad.run2.top100[101:200,]
CHS.rad.run3.top100.set2 <- CHS.rad.run3.top100[101:200,]
CHS.rad.12 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run1.top100.set2$snp))
CHS.rad.13 <- Reduce(intersect, list(CHS.rad.run3.top100.set2$snp, CHS.rad.run1.top100.set2$snp))
CHS.rad.23 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run3.top100.set2$snp))
CHS.rad.123 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run3.top100.set2$snp, CHS.rad.run1.top100.set2$snp))

summary(CHS.rad.123)
#44
```

##overlap in 11-15% of outlier loci for CHS.rad
###############
```
CHS.rad.run1.top100.set2 <- CHS.rad.run1.top100[201:300,]
CHS.rad.run2.top100.set2 <- CHS.rad.run2.top100[201:300,]
CHS.rad.run3.top100.set2 <- CHS.rad.run3.top100[201:300,]
CHS.rad.12 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run1.top100.set2$snp))
CHS.rad.13 <- Reduce(intersect, list(CHS.rad.run3.top100.set2$snp, CHS.rad.run1.top100.set2$snp))
CHS.rad.23 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run3.top100.set2$snp))
CHS.rad.123 <- Reduce(intersect, list(CHS.rad.run2.top100.set2$snp, CHS.rad.run3.top100.set2$snp, CHS.rad.run1.top100.set2$snp))

summary(CHS.rad.123)
#21
```

2. shadow.days

```
CHS.shadow.days.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$shadow.days.rho),]  ##descending order by rho
CHS.shadow.days.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$shadow.days.rho),]  ##descending order by rho
CHS.shadow.days.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.shadow.days
###############
```
CHS.shadow.days.run1.top100.set1 <- head(CHS.shadow.days.run1.top100, 100)
CHS.shadow.days.run2.top100.set1 <- head(CHS.shadow.days.run2.top100, 100)
CHS.shadow.days.run3.top100.set1 <- head(CHS.shadow.days.run3.top100, 100)
CHS.shadow.days.12 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set1$snp, CHS.shadow.days.run1.top100.set1$snp))
CHS.shadow.days.13 <- Reduce(intersect, list(CHS.shadow.days.run3.top100.set1$snp, CHS.shadow.days.run1.top100.set1$snp))
CHS.shadow.days.23 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set1$snp, CHS.shadow.days.run3.top100.set1$snp))
CHS.shadow.days.123 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set1$snp, CHS.shadow.days.run3.top100.set1$snp, CHS.shadow.days.run1.top100.set1$snp))

summary(CHS.shadow.days.123)
#81
```

##overlap in 6-10% of outlier loci for CHS.shadow.days
###############
```
CHS.shadow.days.run1.top100.set2 <- CHS.shadow.days.run1.top100[101:200,]
CHS.shadow.days.run2.top100.set2 <- CHS.shadow.days.run2.top100[101:200,]
CHS.shadow.days.run3.top100.set2 <- CHS.shadow.days.run3.top100[101:200,]
CHS.shadow.days.12 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))
CHS.shadow.days.13 <- Reduce(intersect, list(CHS.shadow.days.run3.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))
CHS.shadow.days.23 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run3.top100.set2$snp))
CHS.shadow.days.123 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run3.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))

summary(CHS.shadow.days.123)
#39
```

##overlap in 11-15% of outlier loci for CHS.shadow.days
###############
```
CHS.shadow.days.run1.top100.set2 <- CHS.shadow.days.run1.top100[201:300,]
CHS.shadow.days.run2.top100.set2 <- CHS.shadow.days.run2.top100[201:300,]
CHS.shadow.days.run3.top100.set2 <- CHS.shadow.days.run3.top100[201:300,]
CHS.shadow.days.12 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))
CHS.shadow.days.13 <- Reduce(intersect, list(CHS.shadow.days.run3.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))
CHS.shadow.days.23 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run3.top100.set2$snp))
CHS.shadow.days.123 <- Reduce(intersect, list(CHS.shadow.days.run2.top100.set2$snp, CHS.shadow.days.run3.top100.set2$snp, CHS.shadow.days.run1.top100.set2$snp))

summary(CHS.shadow.days.123)
#20
```

3. temp

```
CHS.temp.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$temp.rho),]  ##descending order by rho
CHS.temp.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$temp.rho),]  ##descending order by rho
CHS.temp.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.temp
###############
```
CHS.temp.run1.top100.set1 <- head(CHS.temp.run1.top100, 100)
CHS.temp.run2.top100.set1 <- head(CHS.temp.run2.top100, 100)
CHS.temp.run3.top100.set1 <- head(CHS.temp.run3.top100, 100)
CHS.temp.12 <- Reduce(intersect, list(CHS.temp.run2.top100.set1$snp, CHS.temp.run1.top100.set1$snp))
CHS.temp.13 <- Reduce(intersect, list(CHS.temp.run3.top100.set1$snp, CHS.temp.run1.top100.set1$snp))
CHS.temp.23 <- Reduce(intersect, list(CHS.temp.run2.top100.set1$snp, CHS.temp.run3.top100.set1$snp))
CHS.temp.123 <- Reduce(intersect, list(CHS.temp.run2.top100.set1$snp, CHS.temp.run3.top100.set1$snp, CHS.temp.run1.top100.set1$snp))

summary(CHS.temp.123)
#83
```

##overlap in 6-10% of outlier loci for CHS.temp
###############
```
CHS.temp.run1.top100.set2 <- CHS.temp.run1.top100[101:200,]
CHS.temp.run2.top100.set2 <- CHS.temp.run2.top100[101:200,]
CHS.temp.run3.top100.set2 <- CHS.temp.run3.top100[101:200,]
CHS.temp.12 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run1.top100.set2$snp))
CHS.temp.13 <- Reduce(intersect, list(CHS.temp.run3.top100.set2$snp, CHS.temp.run1.top100.set2$snp))
CHS.temp.23 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run3.top100.set2$snp))
CHS.temp.123 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run3.top100.set2$snp, CHS.temp.run1.top100.set2$snp))
summary(CHS.temp.123)
#44
```

##overlap in 11-15% of outlier loci for CHS.temp
###############
```
CHS.temp.run1.top100.set2 <- CHS.temp.run1.top100[201:300,]
CHS.temp.run2.top100.set2 <- CHS.temp.run2.top100[201:300,]
CHS.temp.run3.top100.set2 <- CHS.temp.run3.top100[201:300,]
CHS.temp.12 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run1.top100.set2$snp))
CHS.temp.13 <- Reduce(intersect, list(CHS.temp.run3.top100.set2$snp, CHS.temp.run1.top100.set2$snp))
CHS.temp.23 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run3.top100.set2$snp))
CHS.temp.123 <- Reduce(intersect, list(CHS.temp.run2.top100.set2$snp, CHS.temp.run3.top100.set2$snp, CHS.temp.run1.top100.set2$snp))

summary(CHS.temp.123)
#31
```

4. pcpt

```
CHS.pcpt.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$pcpt.rho),]  ##descending order by rho
CHS.pcpt.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$pcpt.rho),]  ##descending order by rho
CHS.pcpt.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.pcpt
###############
```
CHS.pcpt.run1.top100.set1 <- head(CHS.pcpt.run1.top100, 100)
CHS.pcpt.run2.top100.set1 <- head(CHS.pcpt.run2.top100, 100)
CHS.pcpt.run3.top100.set1 <- head(CHS.pcpt.run3.top100, 100)
CHS.pcpt.12 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set1$snp, CHS.pcpt.run1.top100.set1$snp))
CHS.pcpt.13 <- Reduce(intersect, list(CHS.pcpt.run3.top100.set1$snp, CHS.pcpt.run1.top100.set1$snp))
CHS.pcpt.23 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set1$snp, CHS.pcpt.run3.top100.set1$snp))
CHS.pcpt.123 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set1$snp, CHS.pcpt.run3.top100.set1$snp, CHS.pcpt.run1.top100.set1$snp))

summary(CHS.pcpt.123)
#77
```

##overlap in 6-10% of outlier loci for CHS.pcpt
###############
```
CHS.pcpt.run1.top100.set2 <- CHS.pcpt.run1.top100[101:200,]
CHS.pcpt.run2.top100.set2 <- CHS.pcpt.run2.top100[101:200,]
CHS.pcpt.run3.top100.set2 <- CHS.pcpt.run3.top100[101:200,]
CHS.pcpt.12 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))
CHS.pcpt.13 <- Reduce(intersect, list(CHS.pcpt.run3.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))
CHS.pcpt.23 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run3.top100.set2$snp))
CHS.pcpt.123 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run3.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))

summary(CHS.pcpt.123)
#34
```

##overlap in 11-15% of outlier loci for CHS.pcpt
###############
```
CHS.pcpt.run1.top100.set2 <- CHS.pcpt.run1.top100[201:300,]
CHS.pcpt.run2.top100.set2 <- CHS.pcpt.run2.top100[201:300,]
CHS.pcpt.run3.top100.set2 <- CHS.pcpt.run3.top100[201:300,]
CHS.pcpt.12 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))
CHS.pcpt.13 <- Reduce(intersect, list(CHS.pcpt.run3.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))
CHS.pcpt.23 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run3.top100.set2$snp))
CHS.pcpt.123 <- Reduce(intersect, list(CHS.pcpt.run2.top100.set2$snp, CHS.pcpt.run3.top100.set2$snp, CHS.pcpt.run1.top100.set2$snp))

summary(CHS.pcpt.123)
#12
```

5. day.10cm
```
CHS.day.10cm.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$day.10cm.rho),]  ##descending order by rho
CHS.day.10cm.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$day.10cm.rho),]  ##descending order by rho
CHS.day.10cm.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.day.10cm
###############
```
CHS.day.10cm.run1.top100.set1 <- head(CHS.day.10cm.run1.top100, 100)
CHS.day.10cm.run2.top100.set1 <- head(CHS.day.10cm.run2.top100, 100)
CHS.day.10cm.run3.top100.set1 <- head(CHS.day.10cm.run3.top100, 100)
CHS.day.10cm.12 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set1$snp, CHS.day.10cm.run1.top100.set1$snp))
CHS.day.10cm.13 <- Reduce(intersect, list(CHS.day.10cm.run3.top100.set1$snp, CHS.day.10cm.run1.top100.set1$snp))
CHS.day.10cm.23 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set1$snp, CHS.day.10cm.run3.top100.set1$snp))
CHS.day.10cm.123 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set1$snp, CHS.day.10cm.run3.top100.set1$snp, CHS.day.10cm.run1.top100.set1$snp))

summary(CHS.day.10cm.123)
#78
```

##overlap in 6-10% of outlier loci for CHS.day.10cm
###############
```
CHS.day.10cm.run1.top100.set2 <- CHS.day.10cm.run1.top100[101:200,]
CHS.day.10cm.run2.top100.set2 <- CHS.day.10cm.run2.top100[101:200,]
CHS.day.10cm.run3.top100.set2 <- CHS.day.10cm.run3.top100[101:200,]
CHS.day.10cm.12 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))
CHS.day.10cm.13 <- Reduce(intersect, list(CHS.day.10cm.run3.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))
CHS.day.10cm.23 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run3.top100.set2$snp))
CHS.day.10cm.123 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run3.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))
summary(CHS.day.10cm.123)
#39
```

##overlap in 11-15% of outlier loci for CHS.day.10cm
###############
```
CHS.day.10cm.run1.top100.set2 <- CHS.day.10cm.run1.top100[201:300,]
CHS.day.10cm.run2.top100.set2 <- CHS.day.10cm.run2.top100[201:300,]
CHS.day.10cm.run3.top100.set2 <- CHS.day.10cm.run3.top100[201:300,]
CHS.day.10cm.12 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))
CHS.day.10cm.13 <- Reduce(intersect, list(CHS.day.10cm.run3.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))
CHS.day.10cm.23 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run3.top100.set2$snp))
CHS.day.10cm.123 <- Reduce(intersect, list(CHS.day.10cm.run2.top100.set2$snp, CHS.day.10cm.run3.top100.set2$snp, CHS.day.10cm.run1.top100.set2$snp))

summary(CHS.day.10cm.123)
#22
```

####CZ

1. rad
```
CZ.rad.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$rad.rho),]  ##descending order by rho
CZ.rad.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$rad.rho),]  ##descending order by rho
CZ.rad.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.rad
###############
```
CZ.rad.run1.top100.set1 <- head(CZ.rad.run1.top100, 100)
CZ.rad.run2.top100.set1 <- head(CZ.rad.run2.top100, 100)
CZ.rad.run3.top100.set1 <- head(CZ.rad.run3.top100, 100)
CZ.rad.12 <- Reduce(intersect, list(CZ.rad.run2.top100.set1$snp, CZ.rad.run1.top100.set1$snp))
CZ.rad.13 <- Reduce(intersect, list(CZ.rad.run3.top100.set1$snp, CZ.rad.run1.top100.set1$snp))
CZ.rad.23 <- Reduce(intersect, list(CZ.rad.run2.top100.set1$snp, CZ.rad.run3.top100.set1$snp))
CZ.rad.123 <- Reduce(intersect, list(CZ.rad.run2.top100.set1$snp, CZ.rad.run3.top100.set1$snp, CZ.rad.run1.top100.set1$snp))

summary(CZ.rad.123)
#80
```

##overlap in 6-10% of outlier loci for CZ.rad
###############
```
CZ.rad.run1.top100.set2 <- CZ.rad.run1.top100[101:200,]
CZ.rad.run2.top100.set2 <- CZ.rad.run2.top100[101:200,]
CZ.rad.run3.top100.set2 <- CZ.rad.run3.top100[101:200,]
CZ.rad.12 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run1.top100.set2$snp))
CZ.rad.13 <- Reduce(intersect, list(CZ.rad.run3.top100.set2$snp, CZ.rad.run1.top100.set2$snp))
CZ.rad.23 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run3.top100.set2$snp))
CZ.rad.123 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run3.top100.set2$snp, CZ.rad.run1.top100.set2$snp))

summary(CZ.rad.123)
#48
```

##overlap in 11-15% of outlier loci for CZ.rad
###############
```
CZ.rad.run1.top100.set2 <- CZ.rad.run1.top100[201:300,]
CZ.rad.run2.top100.set2 <- CZ.rad.run2.top100[201:300,]
CZ.rad.run3.top100.set2 <- CZ.rad.run3.top100[201:300,]
CZ.rad.12 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run1.top100.set2$snp))
CZ.rad.13 <- Reduce(intersect, list(CZ.rad.run3.top100.set2$snp, CZ.rad.run1.top100.set2$snp))
CZ.rad.23 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run3.top100.set2$snp))
CZ.rad.123 <- Reduce(intersect, list(CZ.rad.run2.top100.set2$snp, CZ.rad.run3.top100.set2$snp, CZ.rad.run1.top100.set2$snp))

summary(CZ.rad.123)
#26
```

2. shadow.days

```
CZ.shadow.days.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$shadow.days.rho),]  ##descending order by rho
CZ.shadow.days.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$shadow.days.rho),]  ##descending order by rho
CZ.shadow.days.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.shadow.days
###############
```
CZ.shadow.days.run1.top100.set1 <- head(CZ.shadow.days.run1.top100, 100)
CZ.shadow.days.run2.top100.set1 <- head(CZ.shadow.days.run2.top100, 100)
CZ.shadow.days.run3.top100.set1 <- head(CZ.shadow.days.run3.top100, 100)
CZ.shadow.days.12 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set1$snp, CZ.shadow.days.run1.top100.set1$snp))
CZ.shadow.days.13 <- Reduce(intersect, list(CZ.shadow.days.run3.top100.set1$snp, CZ.shadow.days.run1.top100.set1$snp))
CZ.shadow.days.23 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set1$snp, CZ.shadow.days.run3.top100.set1$snp))
CZ.shadow.days.123 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set1$snp, CZ.shadow.days.run3.top100.set1$snp, CZ.shadow.days.run1.top100.set1$snp))

summary(CZ.shadow.days.123)
#81
```

##overlap in 6-10% of outlier loci for CZ.shadow.days
###############
```
CZ.shadow.days.run1.top100.set2 <- CZ.shadow.days.run1.top100[101:200,]
CZ.shadow.days.run2.top100.set2 <- CZ.shadow.days.run2.top100[101:200,]
CZ.shadow.days.run3.top100.set2 <- CZ.shadow.days.run3.top100[101:200,]
CZ.shadow.days.12 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))
CZ.shadow.days.13 <- Reduce(intersect, list(CZ.shadow.days.run3.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))
CZ.shadow.days.23 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run3.top100.set2$snp))
CZ.shadow.days.123 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run3.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))

summary(CZ.shadow.days.123)
#43
```

##overlap in 11-15% of outlier loci for CZ.shadow.days
###############
```
CZ.shadow.days.run1.top100.set2 <- CZ.shadow.days.run1.top100[201:300,]
CZ.shadow.days.run2.top100.set2 <- CZ.shadow.days.run2.top100[201:300,]
CZ.shadow.days.run3.top100.set2 <- CZ.shadow.days.run3.top100[201:300,]
CZ.shadow.days.12 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))
CZ.shadow.days.13 <- Reduce(intersect, list(CZ.shadow.days.run3.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))
CZ.shadow.days.23 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run3.top100.set2$snp))
CZ.shadow.days.123 <- Reduce(intersect, list(CZ.shadow.days.run2.top100.set2$snp, CZ.shadow.days.run3.top100.set2$snp, CZ.shadow.days.run1.top100.set2$snp))

summary(CZ.shadow.days.123)
#30
```

3. temp

```
CZ.temp.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$temp.rho),]  ##descending order by rho
CZ.temp.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$temp.rho),]  ##descending order by rho
CZ.temp.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.temp
###############
```
CZ.temp.run1.top100.set1 <- head(CZ.temp.run1.top100, 100)
CZ.temp.run2.top100.set1 <- head(CZ.temp.run2.top100, 100)
CZ.temp.run3.top100.set1 <- head(CZ.temp.run3.top100, 100)
CZ.temp.12 <- Reduce(intersect, list(CZ.temp.run2.top100.set1$snp, CZ.temp.run1.top100.set1$snp))
CZ.temp.13 <- Reduce(intersect, list(CZ.temp.run3.top100.set1$snp, CZ.temp.run1.top100.set1$snp))
CZ.temp.23 <- Reduce(intersect, list(CZ.temp.run2.top100.set1$snp, CZ.temp.run3.top100.set1$snp))
CZ.temp.123 <- Reduce(intersect, list(CZ.temp.run2.top100.set1$snp, CZ.temp.run3.top100.set1$snp, CZ.temp.run1.top100.set1$snp))

summary(CZ.temp.123)
#77
```

##overlap in 6-10% of outlier loci for CZ.temp
###############
```
CZ.temp.run1.top100.set2 <- CZ.temp.run1.top100[101:200,]
CZ.temp.run2.top100.set2 <- CZ.temp.run2.top100[101:200,]
CZ.temp.run3.top100.set2 <- CZ.temp.run3.top100[101:200,]
CZ.temp.12 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run1.top100.set2$snp))
CZ.temp.13 <- Reduce(intersect, list(CZ.temp.run3.top100.set2$snp, CZ.temp.run1.top100.set2$snp))
CZ.temp.23 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run3.top100.set2$snp))
CZ.temp.123 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run3.top100.set2$snp, CZ.temp.run1.top100.set2$snp))
summary(CZ.temp.123)
#33
```

##overlap in 11-15% of outlier loci for CZ.temp
###############
```
CZ.temp.run1.top100.set2 <- CZ.temp.run1.top100[201:300,]
CZ.temp.run2.top100.set2 <- CZ.temp.run2.top100[201:300,]
CZ.temp.run3.top100.set2 <- CZ.temp.run3.top100[201:300,]
CZ.temp.12 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run1.top100.set2$snp))
CZ.temp.13 <- Reduce(intersect, list(CZ.temp.run3.top100.set2$snp, CZ.temp.run1.top100.set2$snp))
CZ.temp.23 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run3.top100.set2$snp))
CZ.temp.123 <- Reduce(intersect, list(CZ.temp.run2.top100.set2$snp, CZ.temp.run3.top100.set2$snp, CZ.temp.run1.top100.set2$snp))

summary(CZ.temp.123)
#17
```

4. pcpt

```
CZ.pcpt.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$pcpt.rho),]  ##descending order by rho
CZ.pcpt.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$pcpt.rho),]  ##descending order by rho
CZ.pcpt.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.pcpt
###############
```
CZ.pcpt.run1.top100.set1 <- head(CZ.pcpt.run1.top100, 100)
CZ.pcpt.run2.top100.set1 <- head(CZ.pcpt.run2.top100, 100)
CZ.pcpt.run3.top100.set1 <- head(CZ.pcpt.run3.top100, 100)
CZ.pcpt.12 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set1$snp, CZ.pcpt.run1.top100.set1$snp))
CZ.pcpt.13 <- Reduce(intersect, list(CZ.pcpt.run3.top100.set1$snp, CZ.pcpt.run1.top100.set1$snp))
CZ.pcpt.23 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set1$snp, CZ.pcpt.run3.top100.set1$snp))
CZ.pcpt.123 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set1$snp, CZ.pcpt.run3.top100.set1$snp, CZ.pcpt.run1.top100.set1$snp))

summary(CZ.pcpt.123)
#80
```

##overlap in 6-10% of outlier loci for CZ.pcpt
###############
```
CZ.pcpt.run1.top100.set2 <- CZ.pcpt.run1.top100[101:200,]
CZ.pcpt.run2.top100.set2 <- CZ.pcpt.run2.top100[101:200,]
CZ.pcpt.run3.top100.set2 <- CZ.pcpt.run3.top100[101:200,]
CZ.pcpt.12 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))
CZ.pcpt.13 <- Reduce(intersect, list(CZ.pcpt.run3.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))
CZ.pcpt.23 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run3.top100.set2$snp))
CZ.pcpt.123 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run3.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))

summary(CZ.pcpt.123)
#37
```

##overlap in 11-15% of outlier loci for CZ.pcpt
###############
```
CZ.pcpt.run1.top100.set2 <- CZ.pcpt.run1.top100[201:300,]
CZ.pcpt.run2.top100.set2 <- CZ.pcpt.run2.top100[201:300,]
CZ.pcpt.run3.top100.set2 <- CZ.pcpt.run3.top100[201:300,]
CZ.pcpt.12 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))
CZ.pcpt.13 <- Reduce(intersect, list(CZ.pcpt.run3.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))
CZ.pcpt.23 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run3.top100.set2$snp))
CZ.pcpt.123 <- Reduce(intersect, list(CZ.pcpt.run2.top100.set2$snp, CZ.pcpt.run3.top100.set2$snp, CZ.pcpt.run1.top100.set2$snp))

summary(CZ.pcpt.123)
#20
```

5. day.10cm
```
CZ.day.10cm.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$day.10cm.rho),]  ##descending order by rho
CZ.day.10cm.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$day.10cm.rho),]  ##descending order by rho
CZ.day.10cm.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.day.10cm
###############
```
CZ.day.10cm.run1.top100.set1 <- head(CZ.day.10cm.run1.top100, 100)
CZ.day.10cm.run2.top100.set1 <- head(CZ.day.10cm.run2.top100, 100)
CZ.day.10cm.run3.top100.set1 <- head(CZ.day.10cm.run3.top100, 100)
CZ.day.10cm.12 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set1$snp, CZ.day.10cm.run1.top100.set1$snp))
CZ.day.10cm.13 <- Reduce(intersect, list(CZ.day.10cm.run3.top100.set1$snp, CZ.day.10cm.run1.top100.set1$snp))
CZ.day.10cm.23 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set1$snp, CZ.day.10cm.run3.top100.set1$snp))
CZ.day.10cm.123 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set1$snp, CZ.day.10cm.run3.top100.set1$snp, CZ.day.10cm.run1.top100.set1$snp))

summary(CZ.day.10cm.123)
#73
```

##overlap in 6-10% of outlier loci for CZ.day.10cm
###############
```
CZ.day.10cm.run1.top100.set2 <- CZ.day.10cm.run1.top100[101:200,]
CZ.day.10cm.run2.top100.set2 <- CZ.day.10cm.run2.top100[101:200,]
CZ.day.10cm.run3.top100.set2 <- CZ.day.10cm.run3.top100[101:200,]
CZ.day.10cm.12 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))
CZ.day.10cm.13 <- Reduce(intersect, list(CZ.day.10cm.run3.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))
CZ.day.10cm.23 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run3.top100.set2$snp))
CZ.day.10cm.123 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run3.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))
summary(CZ.day.10cm.123)
#31
```

##overlap in 11-15% of outlier loci for CZ.day.10cm
###############
```
CZ.day.10cm.run1.top100.set2 <- CZ.day.10cm.run1.top100[201:300,]
CZ.day.10cm.run2.top100.set2 <- CZ.day.10cm.run2.top100[201:300,]
CZ.day.10cm.run3.top100.set2 <- CZ.day.10cm.run3.top100[201:300,]
CZ.day.10cm.12 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))
CZ.day.10cm.13 <- Reduce(intersect, list(CZ.day.10cm.run3.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))
CZ.day.10cm.23 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run3.top100.set2$snp))
CZ.day.10cm.123 <- Reduce(intersect, list(CZ.day.10cm.run2.top100.set2$snp, CZ.day.10cm.run3.top100.set2$snp, CZ.day.10cm.run1.top100.set2$snp))

summary(CZ.day.10cm.123)
#24
```

####CHS.VS

1. rad
```
CHS.VS.rad.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$rad.rho),]  ##descending order by rho
CHS.VS.rad.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$rad.rho),]  ##descending order by rho
CHS.VS.rad.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.rad
###############
```
CHS.VS.rad.run1.top100.set1 <- head(CHS.VS.rad.run1.top100, 100)
CHS.VS.rad.run2.top100.set1 <- head(CHS.VS.rad.run2.top100, 100)
CHS.VS.rad.run3.top100.set1 <- head(CHS.VS.rad.run3.top100, 100)
CHS.VS.rad.12 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set1$snp, CHS.VS.rad.run1.top100.set1$snp))
CHS.VS.rad.13 <- Reduce(intersect, list(CHS.VS.rad.run3.top100.set1$snp, CHS.VS.rad.run1.top100.set1$snp))
CHS.VS.rad.23 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set1$snp, CHS.VS.rad.run3.top100.set1$snp))
CHS.VS.rad.123 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set1$snp, CHS.VS.rad.run3.top100.set1$snp, CHS.VS.rad.run1.top100.set1$snp))

summary(CHS.VS.rad.123)
#77
```

##overlap in 6-10% of outlier loci for CHS.VS.rad
###############
```
CHS.VS.rad.run1.top100.set2 <- CHS.VS.rad.run1.top100[101:200,]
CHS.VS.rad.run2.top100.set2 <- CHS.VS.rad.run2.top100[101:200,]
CHS.VS.rad.run3.top100.set2 <- CHS.VS.rad.run3.top100[101:200,]
CHS.VS.rad.12 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))
CHS.VS.rad.13 <- Reduce(intersect, list(CHS.VS.rad.run3.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))
CHS.VS.rad.23 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run3.top100.set2$snp))
CHS.VS.rad.123 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run3.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))

summary(CHS.VS.rad.123)
#38
```

##overlap in 11-15% of outlier loci for CHS.VS.rad
###############
```
CHS.VS.rad.run1.top100.set2 <- CHS.VS.rad.run1.top100[201:300,]
CHS.VS.rad.run2.top100.set2 <- CHS.VS.rad.run2.top100[201:300,]
CHS.VS.rad.run3.top100.set2 <- CHS.VS.rad.run3.top100[201:300,]
CHS.VS.rad.12 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))
CHS.VS.rad.13 <- Reduce(intersect, list(CHS.VS.rad.run3.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))
CHS.VS.rad.23 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run3.top100.set2$snp))
CHS.VS.rad.123 <- Reduce(intersect, list(CHS.VS.rad.run2.top100.set2$snp, CHS.VS.rad.run3.top100.set2$snp, CHS.VS.rad.run1.top100.set2$snp))

summary(CHS.VS.rad.123)
#21
```

2. shadow.days

```
CHS.VS.shadow.days.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$shadow.days.rho),]  ##descending order by rho
CHS.VS.shadow.days.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$shadow.days.rho),]  ##descending order by rho
CHS.VS.shadow.days.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.shadow.days
###############
```
CHS.VS.shadow.days.run1.top100.set1 <- head(CHS.VS.shadow.days.run1.top100, 100)
CHS.VS.shadow.days.run2.top100.set1 <- head(CHS.VS.shadow.days.run2.top100, 100)
CHS.VS.shadow.days.run3.top100.set1 <- head(CHS.VS.shadow.days.run3.top100, 100)
CHS.VS.shadow.days.12 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set1$snp, CHS.VS.shadow.days.run1.top100.set1$snp))
CHS.VS.shadow.days.13 <- Reduce(intersect, list(CHS.VS.shadow.days.run3.top100.set1$snp, CHS.VS.shadow.days.run1.top100.set1$snp))
CHS.VS.shadow.days.23 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set1$snp, CHS.VS.shadow.days.run3.top100.set1$snp))
CHS.VS.shadow.days.123 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set1$snp, CHS.VS.shadow.days.run3.top100.set1$snp, CHS.VS.shadow.days.run1.top100.set1$snp))

summary(CHS.VS.shadow.days.123)
#78
```

##overlap in 6-10% of outlier loci for CHS.VS.shadow.days
###############
```
CHS.VS.shadow.days.run1.top100.set2 <- CHS.VS.shadow.days.run1.top100[101:200,]
CHS.VS.shadow.days.run2.top100.set2 <- CHS.VS.shadow.days.run2.top100[101:200,]
CHS.VS.shadow.days.run3.top100.set2 <- CHS.VS.shadow.days.run3.top100[101:200,]
CHS.VS.shadow.days.12 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))
CHS.VS.shadow.days.13 <- Reduce(intersect, list(CHS.VS.shadow.days.run3.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))
CHS.VS.shadow.days.23 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run3.top100.set2$snp))
CHS.VS.shadow.days.123 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run3.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))

summary(CHS.VS.shadow.days.123)
#44
```

##overlap in 11-15% of outlier loci for CHS.VS.shadow.days
###############
```
CHS.VS.shadow.days.run1.top100.set2 <- CHS.VS.shadow.days.run1.top100[201:300,]
CHS.VS.shadow.days.run2.top100.set2 <- CHS.VS.shadow.days.run2.top100[201:300,]
CHS.VS.shadow.days.run3.top100.set2 <- CHS.VS.shadow.days.run3.top100[201:300,]
CHS.VS.shadow.days.12 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))
CHS.VS.shadow.days.13 <- Reduce(intersect, list(CHS.VS.shadow.days.run3.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))
CHS.VS.shadow.days.23 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run3.top100.set2$snp))
CHS.VS.shadow.days.123 <- Reduce(intersect, list(CHS.VS.shadow.days.run2.top100.set2$snp, CHS.VS.shadow.days.run3.top100.set2$snp, CHS.VS.shadow.days.run1.top100.set2$snp))

summary(CHS.VS.shadow.days.123)
#27
```

3. temp

```
CHS.VS.temp.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$temp.rho),]  ##descending order by rho
CHS.VS.temp.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$temp.rho),]  ##descending order by rho
CHS.VS.temp.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.temp
###############
```
CHS.VS.temp.run1.top100.set1 <- head(CHS.VS.temp.run1.top100, 100)
CHS.VS.temp.run2.top100.set1 <- head(CHS.VS.temp.run2.top100, 100)
CHS.VS.temp.run3.top100.set1 <- head(CHS.VS.temp.run3.top100, 100)
CHS.VS.temp.12 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set1$snp, CHS.VS.temp.run1.top100.set1$snp))
CHS.VS.temp.13 <- Reduce(intersect, list(CHS.VS.temp.run3.top100.set1$snp, CHS.VS.temp.run1.top100.set1$snp))
CHS.VS.temp.23 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set1$snp, CHS.VS.temp.run3.top100.set1$snp))
CHS.VS.temp.123 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set1$snp, CHS.VS.temp.run3.top100.set1$snp, CHS.VS.temp.run1.top100.set1$snp))

summary(CHS.VS.temp.123)
#87
```

##overlap in 6-10% of outlier loci for CHS.VS.temp
###############
```
CHS.VS.temp.run1.top100.set2 <- CHS.VS.temp.run1.top100[101:200,]
CHS.VS.temp.run2.top100.set2 <- CHS.VS.temp.run2.top100[101:200,]
CHS.VS.temp.run3.top100.set2 <- CHS.VS.temp.run3.top100[101:200,]
CHS.VS.temp.12 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))
CHS.VS.temp.13 <- Reduce(intersect, list(CHS.VS.temp.run3.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))
CHS.VS.temp.23 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run3.top100.set2$snp))
CHS.VS.temp.123 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run3.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))
summary(CHS.VS.temp.123)
#51
```

##overlap in 11-15% of outlier loci for CHS.VS.temp
###############
```
CHS.VS.temp.run1.top100.set2 <- CHS.VS.temp.run1.top100[201:300,]
CHS.VS.temp.run2.top100.set2 <- CHS.VS.temp.run2.top100[201:300,]
CHS.VS.temp.run3.top100.set2 <- CHS.VS.temp.run3.top100[201:300,]
CHS.VS.temp.12 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))
CHS.VS.temp.13 <- Reduce(intersect, list(CHS.VS.temp.run3.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))
CHS.VS.temp.23 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run3.top100.set2$snp))
CHS.VS.temp.123 <- Reduce(intersect, list(CHS.VS.temp.run2.top100.set2$snp, CHS.VS.temp.run3.top100.set2$snp, CHS.VS.temp.run1.top100.set2$snp))

summary(CHS.VS.temp.123)
#17
```

4. pcpt

```
CHS.VS.pcpt.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$pcpt.rho),]  ##descending order by rho
CHS.VS.pcpt.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$pcpt.rho),]  ##descending order by rho
CHS.VS.pcpt.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.pcpt
###############
```
CHS.VS.pcpt.run1.top100.set1 <- head(CHS.VS.pcpt.run1.top100, 100)
CHS.VS.pcpt.run2.top100.set1 <- head(CHS.VS.pcpt.run2.top100, 100)
CHS.VS.pcpt.run3.top100.set1 <- head(CHS.VS.pcpt.run3.top100, 100)
CHS.VS.pcpt.12 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set1$snp, CHS.VS.pcpt.run1.top100.set1$snp))
CHS.VS.pcpt.13 <- Reduce(intersect, list(CHS.VS.pcpt.run3.top100.set1$snp, CHS.VS.pcpt.run1.top100.set1$snp))
CHS.VS.pcpt.23 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set1$snp, CHS.VS.pcpt.run3.top100.set1$snp))
CHS.VS.pcpt.123 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set1$snp, CHS.VS.pcpt.run3.top100.set1$snp, CHS.VS.pcpt.run1.top100.set1$snp))

summary(CHS.VS.pcpt.123)
#75
```

##overlap in 6-10% of outlier loci for CHS.VS.pcpt
###############
```
CHS.VS.pcpt.run1.top100.set2 <- CHS.VS.pcpt.run1.top100[101:200,]
CHS.VS.pcpt.run2.top100.set2 <- CHS.VS.pcpt.run2.top100[101:200,]
CHS.VS.pcpt.run3.top100.set2 <- CHS.VS.pcpt.run3.top100[101:200,]
CHS.VS.pcpt.12 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))
CHS.VS.pcpt.13 <- Reduce(intersect, list(CHS.VS.pcpt.run3.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))
CHS.VS.pcpt.23 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run3.top100.set2$snp))
CHS.VS.pcpt.123 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run3.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))

summary(CHS.VS.pcpt.123)
#29
```

##overlap in 11-15% of outlier loci for CHS.VS.pcpt
###############
```
CHS.VS.pcpt.run1.top100.set2 <- CHS.VS.pcpt.run1.top100[201:300,]
CHS.VS.pcpt.run2.top100.set2 <- CHS.VS.pcpt.run2.top100[201:300,]
CHS.VS.pcpt.run3.top100.set2 <- CHS.VS.pcpt.run3.top100[201:300,]
CHS.VS.pcpt.12 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))
CHS.VS.pcpt.13 <- Reduce(intersect, list(CHS.VS.pcpt.run3.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))
CHS.VS.pcpt.23 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run3.top100.set2$snp))
CHS.VS.pcpt.123 <- Reduce(intersect, list(CHS.VS.pcpt.run2.top100.set2$snp, CHS.VS.pcpt.run3.top100.set2$snp, CHS.VS.pcpt.run1.top100.set2$snp))

summary(CHS.VS.pcpt.123)
#21
```

5. day.10cm
```
CHS.VS.day.10cm.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$day.10cm.rho),]  ##descending order by rho
CHS.VS.day.10cm.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$day.10cm.rho),]  ##descending order by rho
CHS.VS.day.10cm.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.day.10cm
###############
```
CHS.VS.day.10cm.run1.top100.set1 <- head(CHS.VS.day.10cm.run1.top100, 100)
CHS.VS.day.10cm.run2.top100.set1 <- head(CHS.VS.day.10cm.run2.top100, 100)
CHS.VS.day.10cm.run3.top100.set1 <- head(CHS.VS.day.10cm.run3.top100, 100)
CHS.VS.day.10cm.12 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set1$snp, CHS.VS.day.10cm.run1.top100.set1$snp))
CHS.VS.day.10cm.13 <- Reduce(intersect, list(CHS.VS.day.10cm.run3.top100.set1$snp, CHS.VS.day.10cm.run1.top100.set1$snp))
CHS.VS.day.10cm.23 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set1$snp, CHS.VS.day.10cm.run3.top100.set1$snp))
CHS.VS.day.10cm.123 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set1$snp, CHS.VS.day.10cm.run3.top100.set1$snp, CHS.VS.day.10cm.run1.top100.set1$snp))

summary(CHS.VS.day.10cm.123)
#71
```

##overlap in 6-10% of outlier loci for CHS.VS.day.10cm
###############
```
CHS.VS.day.10cm.run1.top100.set2 <- CHS.VS.day.10cm.run1.top100[101:200,]
CHS.VS.day.10cm.run2.top100.set2 <- CHS.VS.day.10cm.run2.top100[101:200,]
CHS.VS.day.10cm.run3.top100.set2 <- CHS.VS.day.10cm.run3.top100[101:200,]
CHS.VS.day.10cm.12 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))
CHS.VS.day.10cm.13 <- Reduce(intersect, list(CHS.VS.day.10cm.run3.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))
CHS.VS.day.10cm.23 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run3.top100.set2$snp))
CHS.VS.day.10cm.123 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run3.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))
summary(CHS.VS.day.10cm.123)
#32
```

##overlap in 11-15% of outlier loci for CHS.VS.day.10cm
###############
```
CHS.VS.day.10cm.run1.top100.set2 <- CHS.VS.day.10cm.run1.top100[201:300,]
CHS.VS.day.10cm.run2.top100.set2 <- CHS.VS.day.10cm.run2.top100[201:300,]
CHS.VS.day.10cm.run3.top100.set2 <- CHS.VS.day.10cm.run3.top100[201:300,]
CHS.VS.day.10cm.12 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))
CHS.VS.day.10cm.13 <- Reduce(intersect, list(CHS.VS.day.10cm.run3.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))
CHS.VS.day.10cm.23 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run3.top100.set2$snp))
CHS.VS.day.10cm.123 <- Reduce(intersect, list(CHS.VS.day.10cm.run2.top100.set2$snp, CHS.VS.day.10cm.run3.top100.set2$snp, CHS.VS.day.10cm.run1.top100.set2$snp))

summary(CHS.VS.day.10cm.123)
#18
```

####CHS.TI

1. rad
```
CHS.TI.rad.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$rad.rho),]  ##descending order by rho
CHS.TI.rad.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$rad.rho),]  ##descending order by rho
CHS.TI.rad.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$rad.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.rad
###############
```
CHS.TI.rad.run1.top100.set1 <- head(CHS.TI.rad.run1.top100, 100)
CHS.TI.rad.run2.top100.set1 <- head(CHS.TI.rad.run2.top100, 100)
CHS.TI.rad.run3.top100.set1 <- head(CHS.TI.rad.run3.top100, 100)
CHS.TI.rad.12 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set1$snp, CHS.TI.rad.run1.top100.set1$snp))
CHS.TI.rad.13 <- Reduce(intersect, list(CHS.TI.rad.run3.top100.set1$snp, CHS.TI.rad.run1.top100.set1$snp))
CHS.TI.rad.23 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set1$snp, CHS.TI.rad.run3.top100.set1$snp))
CHS.TI.rad.123 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set1$snp, CHS.TI.rad.run3.top100.set1$snp, CHS.TI.rad.run1.top100.set1$snp))

summary(CHS.TI.rad.123)
#68
```

##overlap in 6-10% of outlier loci for CHS.TI.rad
###############
```
CHS.TI.rad.run1.top100.set2 <- CHS.TI.rad.run1.top100[101:200,]
CHS.TI.rad.run2.top100.set2 <- CHS.TI.rad.run2.top100[101:200,]
CHS.TI.rad.run3.top100.set2 <- CHS.TI.rad.run3.top100[101:200,]
CHS.TI.rad.12 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))
CHS.TI.rad.13 <- Reduce(intersect, list(CHS.TI.rad.run3.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))
CHS.TI.rad.23 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run3.top100.set2$snp))
CHS.TI.rad.123 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run3.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))

summary(CHS.TI.rad.123)
#33
```

##overlap in 11-15% of outlier loci for CHS.TI.rad
###############
```
CHS.TI.rad.run1.top100.set2 <- CHS.TI.rad.run1.top100[201:300,]
CHS.TI.rad.run2.top100.set2 <- CHS.TI.rad.run2.top100[201:300,]
CHS.TI.rad.run3.top100.set2 <- CHS.TI.rad.run3.top100[201:300,]
CHS.TI.rad.12 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))
CHS.TI.rad.13 <- Reduce(intersect, list(CHS.TI.rad.run3.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))
CHS.TI.rad.23 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run3.top100.set2$snp))
CHS.TI.rad.123 <- Reduce(intersect, list(CHS.TI.rad.run2.top100.set2$snp, CHS.TI.rad.run3.top100.set2$snp, CHS.TI.rad.run1.top100.set2$snp))

summary(CHS.TI.rad.123)
#15
```

2. shadow.days

```
CHS.TI.shadow.days.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$shadow.days.rho),]  ##descending order by rho
CHS.TI.shadow.days.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$shadow.days.rho),]  ##descending order by rho
CHS.TI.shadow.days.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$shadow.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.shadow.days
###############
```
CHS.TI.shadow.days.run1.top100.set1 <- head(CHS.TI.shadow.days.run1.top100, 100)
CHS.TI.shadow.days.run2.top100.set1 <- head(CHS.TI.shadow.days.run2.top100, 100)
CHS.TI.shadow.days.run3.top100.set1 <- head(CHS.TI.shadow.days.run3.top100, 100)
CHS.TI.shadow.days.12 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set1$snp, CHS.TI.shadow.days.run1.top100.set1$snp))
CHS.TI.shadow.days.13 <- Reduce(intersect, list(CHS.TI.shadow.days.run3.top100.set1$snp, CHS.TI.shadow.days.run1.top100.set1$snp))
CHS.TI.shadow.days.23 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set1$snp, CHS.TI.shadow.days.run3.top100.set1$snp))
CHS.TI.shadow.days.123 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set1$snp, CHS.TI.shadow.days.run3.top100.set1$snp, CHS.TI.shadow.days.run1.top100.set1$snp))

summary(CHS.TI.shadow.days.123)
#72
```

##overlap in 6-10% of outlier loci for CHS.TI.shadow.days
###############
```
CHS.TI.shadow.days.run1.top100.set2 <- CHS.TI.shadow.days.run1.top100[101:200,]
CHS.TI.shadow.days.run2.top100.set2 <- CHS.TI.shadow.days.run2.top100[101:200,]
CHS.TI.shadow.days.run3.top100.set2 <- CHS.TI.shadow.days.run3.top100[101:200,]
CHS.TI.shadow.days.12 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))
CHS.TI.shadow.days.13 <- Reduce(intersect, list(CHS.TI.shadow.days.run3.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))
CHS.TI.shadow.days.23 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run3.top100.set2$snp))
CHS.TI.shadow.days.123 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run3.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))

summary(CHS.TI.shadow.days.123)
#33
```

##overlap in 11-15% of outlier loci for CHS.TI.shadow.days
###############
```
CHS.TI.shadow.days.run1.top100.set2 <- CHS.TI.shadow.days.run1.top100[201:300,]
CHS.TI.shadow.days.run2.top100.set2 <- CHS.TI.shadow.days.run2.top100[201:300,]
CHS.TI.shadow.days.run3.top100.set2 <- CHS.TI.shadow.days.run3.top100[201:300,]
CHS.TI.shadow.days.12 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))
CHS.TI.shadow.days.13 <- Reduce(intersect, list(CHS.TI.shadow.days.run3.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))
CHS.TI.shadow.days.23 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run3.top100.set2$snp))
CHS.TI.shadow.days.123 <- Reduce(intersect, list(CHS.TI.shadow.days.run2.top100.set2$snp, CHS.TI.shadow.days.run3.top100.set2$snp, CHS.TI.shadow.days.run1.top100.set2$snp))

summary(CHS.TI.shadow.days.123)
#20
```

3. temp

```
CHS.TI.temp.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$temp.rho),]  ##descending order by rho
CHS.TI.temp.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$temp.rho),]  ##descending order by rho
CHS.TI.temp.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$temp.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.temp
###############
```
CHS.TI.temp.run1.top100.set1 <- head(CHS.TI.temp.run1.top100, 100)
CHS.TI.temp.run2.top100.set1 <- head(CHS.TI.temp.run2.top100, 100)
CHS.TI.temp.run3.top100.set1 <- head(CHS.TI.temp.run3.top100, 100)
CHS.TI.temp.12 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set1$snp, CHS.TI.temp.run1.top100.set1$snp))
CHS.TI.temp.13 <- Reduce(intersect, list(CHS.TI.temp.run3.top100.set1$snp, CHS.TI.temp.run1.top100.set1$snp))
CHS.TI.temp.23 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set1$snp, CHS.TI.temp.run3.top100.set1$snp))
CHS.TI.temp.123 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set1$snp, CHS.TI.temp.run3.top100.set1$snp, CHS.TI.temp.run1.top100.set1$snp))

summary(CHS.TI.temp.123)
#77
```

##overlap in 6-10% of outlier loci for CHS.TI.temp
###############
```
CHS.TI.temp.run1.top100.set2 <- CHS.TI.temp.run1.top100[101:200,]
CHS.TI.temp.run2.top100.set2 <- CHS.TI.temp.run2.top100[101:200,]
CHS.TI.temp.run3.top100.set2 <- CHS.TI.temp.run3.top100[101:200,]
CHS.TI.temp.12 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))
CHS.TI.temp.13 <- Reduce(intersect, list(CHS.TI.temp.run3.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))
CHS.TI.temp.23 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run3.top100.set2$snp))
CHS.TI.temp.123 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run3.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))
summary(CHS.TI.temp.123)
#28
```

##overlap in 11-15% of outlier loci for CHS.TI.temp
###############
```
CHS.TI.temp.run1.top100.set2 <- CHS.TI.temp.run1.top100[201:300,]
CHS.TI.temp.run2.top100.set2 <- CHS.TI.temp.run2.top100[201:300,]
CHS.TI.temp.run3.top100.set2 <- CHS.TI.temp.run3.top100[201:300,]
CHS.TI.temp.12 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))
CHS.TI.temp.13 <- Reduce(intersect, list(CHS.TI.temp.run3.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))
CHS.TI.temp.23 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run3.top100.set2$snp))
CHS.TI.temp.123 <- Reduce(intersect, list(CHS.TI.temp.run2.top100.set2$snp, CHS.TI.temp.run3.top100.set2$snp, CHS.TI.temp.run1.top100.set2$snp))

summary(CHS.TI.temp.123)
#12
```

4. pcpt

```
CHS.TI.pcpt.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$pcpt.rho),]  ##descending order by rho
CHS.TI.pcpt.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$pcpt.rho),]  ##descending order by rho
CHS.TI.pcpt.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$pcpt.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.pcpt
###############
```
CHS.TI.pcpt.run1.top100.set1 <- head(CHS.TI.pcpt.run1.top100, 100)
CHS.TI.pcpt.run2.top100.set1 <- head(CHS.TI.pcpt.run2.top100, 100)
CHS.TI.pcpt.run3.top100.set1 <- head(CHS.TI.pcpt.run3.top100, 100)
CHS.TI.pcpt.12 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set1$snp, CHS.TI.pcpt.run1.top100.set1$snp))
CHS.TI.pcpt.13 <- Reduce(intersect, list(CHS.TI.pcpt.run3.top100.set1$snp, CHS.TI.pcpt.run1.top100.set1$snp))
CHS.TI.pcpt.23 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set1$snp, CHS.TI.pcpt.run3.top100.set1$snp))
CHS.TI.pcpt.123 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set1$snp, CHS.TI.pcpt.run3.top100.set1$snp, CHS.TI.pcpt.run1.top100.set1$snp))

summary(CHS.TI.pcpt.123)
#67
```

##overlap in 6-10% of outlier loci for CHS.TI.pcpt
###############
```
CHS.TI.pcpt.run1.top100.set2 <- CHS.TI.pcpt.run1.top100[101:200,]
CHS.TI.pcpt.run2.top100.set2 <- CHS.TI.pcpt.run2.top100[101:200,]
CHS.TI.pcpt.run3.top100.set2 <- CHS.TI.pcpt.run3.top100[101:200,]
CHS.TI.pcpt.12 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))
CHS.TI.pcpt.13 <- Reduce(intersect, list(CHS.TI.pcpt.run3.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))
CHS.TI.pcpt.23 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run3.top100.set2$snp))
CHS.TI.pcpt.123 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run3.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))

summary(CHS.TI.pcpt.123)
#26
```

##overlap in 11-15% of outlier loci for CHS.TI.pcpt
###############
```
CHS.TI.pcpt.run1.top100.set2 <- CHS.TI.pcpt.run1.top100[201:300,]
CHS.TI.pcpt.run2.top100.set2 <- CHS.TI.pcpt.run2.top100[201:300,]
CHS.TI.pcpt.run3.top100.set2 <- CHS.TI.pcpt.run3.top100[201:300,]
CHS.TI.pcpt.12 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))
CHS.TI.pcpt.13 <- Reduce(intersect, list(CHS.TI.pcpt.run3.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))
CHS.TI.pcpt.23 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run3.top100.set2$snp))
CHS.TI.pcpt.123 <- Reduce(intersect, list(CHS.TI.pcpt.run2.top100.set2$snp, CHS.TI.pcpt.run3.top100.set2$snp, CHS.TI.pcpt.run1.top100.set2$snp))

summary(CHS.TI.pcpt.123)
#14
```

5. day.10cm
```
CHS.TI.day.10cm.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$day.10cm.rho),]  ##descending order by rho
CHS.TI.day.10cm.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$day.10cm.rho),]  ##descending order by rho
CHS.TI.day.10cm.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$day.10cm.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.day.10cm
###############
```
CHS.TI.day.10cm.run1.top100.set1 <- head(CHS.TI.day.10cm.run1.top100, 100)
CHS.TI.day.10cm.run2.top100.set1 <- head(CHS.TI.day.10cm.run2.top100, 100)
CHS.TI.day.10cm.run3.top100.set1 <- head(CHS.TI.day.10cm.run3.top100, 100)
CHS.TI.day.10cm.12 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set1$snp, CHS.TI.day.10cm.run1.top100.set1$snp))
CHS.TI.day.10cm.13 <- Reduce(intersect, list(CHS.TI.day.10cm.run3.top100.set1$snp, CHS.TI.day.10cm.run1.top100.set1$snp))
CHS.TI.day.10cm.23 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set1$snp, CHS.TI.day.10cm.run3.top100.set1$snp))
CHS.TI.day.10cm.123 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set1$snp, CHS.TI.day.10cm.run3.top100.set1$snp, CHS.TI.day.10cm.run1.top100.set1$snp))

summary(CHS.TI.day.10cm.123)
#70
```

##overlap in 6-10% of outlier loci for CHS.TI.day.10cm
###############
```
CHS.TI.day.10cm.run1.top100.set2 <- CHS.TI.day.10cm.run1.top100[101:200,]
CHS.TI.day.10cm.run2.top100.set2 <- CHS.TI.day.10cm.run2.top100[101:200,]
CHS.TI.day.10cm.run3.top100.set2 <- CHS.TI.day.10cm.run3.top100[101:200,]
CHS.TI.day.10cm.12 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))
CHS.TI.day.10cm.13 <- Reduce(intersect, list(CHS.TI.day.10cm.run3.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))
CHS.TI.day.10cm.23 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run3.top100.set2$snp))
CHS.TI.day.10cm.123 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run3.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))
summary(CHS.TI.day.10cm.123)
#32
```

##overlap in 11-15% of outlier loci for CHS.TI.day.10cm
###############
```
CHS.TI.day.10cm.run1.top100.set2 <- CHS.TI.day.10cm.run1.top100[201:300,]
CHS.TI.day.10cm.run2.top100.set2 <- CHS.TI.day.10cm.run2.top100[201:300,]
CHS.TI.day.10cm.run3.top100.set2 <- CHS.TI.day.10cm.run3.top100[201:300,]
CHS.TI.day.10cm.12 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))
CHS.TI.day.10cm.13 <- Reduce(intersect, list(CHS.TI.day.10cm.run3.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))
CHS.TI.day.10cm.23 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run3.top100.set2$snp))
CHS.TI.day.10cm.123 <- Reduce(intersect, list(CHS.TI.day.10cm.run2.top100.set2$snp, CHS.TI.day.10cm.run3.top100.set2$snp, CHS.TI.day.10cm.run1.top100.set2$snp))

summary(CHS.TI.day.10cm.123)
#18
```

#####Concatenating the data and finding candidate loci: BayENV2
################################


According to Blair et al. 2014, there is a lot of inter-run variation in BayEnv2 results. They suggest that the best approach (if the covariance matrix converges) is to combine results from several different runs.

Im following the methods in Božičević 2016 (Mol Ecol) to get the most robust results:

" To test the convergence of Bayenv2, we used several independent Markov chain Monte Carlo runs with a maximum chain length of 10 000 iterations. We observed convergence after about 5000 iterations (Fig. S1, Supporting information). However, these chains might converge to different solutions. To be most stringent, we used the median results from 10 independent runs (Blair et al. 2014). We then tested for correlations between each single SNP and six environmental variables: geographical latitude, height above mean sea level and four temperature measures (average daily minimum of the coldest and warmest month, and average daily minimum and maximum throughout the year) (Table 2). The results for environmental variables are given as Bayes factors (BFs). A higher BF gives higher support to the model where the environmental variable has a significant effect on allele frequency distribution over an alternative model with no effect (Coop et al. 2010). Similar as above, we finally report the median BF of ten independent runs of each SNP that has been described to improve the proportion of false positives (Blair et al. 2014; Lotterhos & Whitlock 2014). With BF values, we also used a resampling approach, analogously to the one applied on FST. We randomly sampled sets of SNPs of the same size from the genomic background and then assessed the null distribution for BF of CCRT-, RSS- and SR-associated SNPs. "

I've run 3 independend BayEnv2 association analyses. now Im combining all the output tables by the median of all the values.


Calculate a final matrix as the median of all the matrices

```
library(data.table)

##CHall
CHall.BF.all <- rbindlist(list(BF.CHall.Run1, BF.CHall.Run2, BF.CHall.Run3))[,lapply(.SD,median), list(snp)]


##CHN
CHN.BF.all <- rbindlist(list(BF.CHN.Run1, BF.CHN.Run2, BF.CHN.Run3))[,lapply(.SD,median), list(snp)]

##CHS
CHS.BF.all <- rbindlist(list(BF.CHS.Run1, BF.CHS.Run2, BF.CHS.Run3))[,lapply(.SD,median), list(snp)]

##CZ
CZ.BF.all <- rbindlist(list(BF.CZ.Run1, BF.CZ.Run2, BF.CZ.Run3))[,lapply(.SD,median), list(snp)]

##CHS.VS
CHS.VS.BF.all <- rbindlist(list(BF.CHS.VS.Run1, BF.CHS.VS.Run2, BF.CHS.VS.Run3))[,lapply(.SD,median), list(snp)]

##CHS.TI
CHS.TI.BF.all <- rbindlist(list(BF.CHS.TI.Run1, BF.CHS.TI.Run2, BF.CHS.TI.Run3))[,lapply(.SD,median), list(snp)]
```


Identify the outliers associated with the different env variables. " Strength of evidence for significant associations was based on the value of the log10 Bayes factor (log10BF), with the following log10BF cut-offs: 0.5–1= substantial evidence; 1–2 = strong evidence; >2 =decisive (Kass & Raftery 1995). The linear model underlying the Bayes factor might not be correct or outliers within our data might misguide the model (bayenv2.0 manual, https://bitbucket.org/tguenther/bayenv2_public/src). To deal with this, bayenv2 also calculates the nonparametric Spearman's rank correlation coefficient, ρ. SNPs with a log10BF >0.5 as well as an absolute value of ρ > 0.3 (where ρ ranges from −1 to 1) were therefore considered as robust candidates demonstrating signatures of selection. " Christmas et al. 2016


###CHall
```
#calculate the log10 of BF for each environmental variable

CHall.BF.all$rad.log10BF <- log10(CHall.BF.all$rad.BF)
CHall.BF.all$shadow.days.log10BF <- log10(CHall.BF.all$shadow.days.BF)
CHall.BF.all$temp.log10BF <- log10(CHall.BF.all$temp.BF)
CHall.BF.all$pcpt.log10BF <- log10(CHall.BF.all$pcpt.BF)
CHall.BF.all$day.10cm.log10BF <- log10(CHall.BF.all$day.10cm.BF)

##Sort and plot

CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$rad.log10BF),]

pdf("CHall.BF.plot.pdf")
par(mfrow=c(3,2))
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$rad.log10BF),]
plot(CHall.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$shadow.days.log10BF),]
plot(CHall.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$temp.log10BF),]
plot(CHall.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$pcpt.log10BF),]
plot(CHall.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$day.10cm.log10BF),]
plot(CHall.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHall.BFvsp.pdf")
par(mfrow=c(3,2))
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$rad.log10BF),]
plot(CHall.BF.all.sort$rad.log10BF~(abs(CHall.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$shadow.days.log10BF),]
plot(CHall.BF.all.sort$shadow.days.log10BF~(abs(CHall.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$temp.log10BF),]
plot(CHall.BF.all.sort$temp.log10BF~(abs(CHall.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$pcpt.log10BF),]
plot(CHall.BF.all.sort$pcpt.log10BF~(abs(CHall.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHall.BF.all.sort <- CHall.BF.all[order(CHall.BF.all$day.10cm.log10BF),]
plot(CHall.BF.all.sort$day.10cm.log10BF~(abs(CHall.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHall.rad.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$rad.log10BF>0.5 & (abs(CHall.BF.all$rad.rho))>0.3),]
CHall.shadow.days.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$shadow.days.log10BF>0.5 & (abs(CHall.BF.all$shadow.days.rho))>0.3),]
CHall.temp.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$temp.log10BF>0.5 & (abs(CHall.BF.all$temp.rho))>0.3),]
CHall.pcpt.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$pcpt.log10BF>0.5 & (abs(CHall.BF.all$pcpt.rho))>0.3),]
CHall.day.10cm.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$day.10cm.log10BF>0.5 & (abs(CHall.BF.all$day.10cm.rho))>0.3),]
```





####CHN
```
#calculate the log10 of BF for each environmental variable

CHN.BF.all$rad.log10BF <- log10(CHN.BF.all$rad.BF)
CHN.BF.all$shadow.days.log10BF <- log10(CHN.BF.all$shadow.days.BF)
CHN.BF.all$temp.log10BF <- log10(CHN.BF.all$temp.BF)
CHN.BF.all$pcpt.log10BF <- log10(CHN.BF.all$pcpt.BF)
CHN.BF.all$day.10cm.log10BF <- log10(CHN.BF.all$day.10cm.BF)

##Sort and plot

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$rad.log10BF),]

pdf("CHN.BF.plot.pdf")
par(mfrow=c(3,2))
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$rad.log10BF),]
plot(CHN.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$shadow.days.log10BF),]
plot(CHN.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$temp.log10BF),]
plot(CHN.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$pcpt.log10BF),]
plot(CHN.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$day.10cm.log10BF),]
plot(CHN.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```


##plot BF vs absolute p value
```
pdf("CHN.BFvsp.pdf")
par(mfrow=c(3,2))
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$rad.log10BF),]
plot(CHN.BF.all.sort$rad.log10BF~(abs(CHN.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$shadow.days.log10BF),]
plot(CHN.BF.all.sort$shadow.days.log10BF~(abs(CHN.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$temp.log10BF),]
plot(CHN.BF.all.sort$temp.log10BF~(abs(CHN.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$pcpt.log10BF),]
plot(CHN.BF.all.sort$pcpt.log10BF~(abs(CHN.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$day.10cm.log10BF),]
plot(CHN.BF.all.sort$day.10cm.log10BF~(abs(CHN.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHN.rad.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$rad.log10BF>0.5 & (abs(CHN.BF.all$rad.rho))>0.3),]
CHN.shadow.days.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$shadow.days.log10BF>0.5 & (abs(CHN.BF.all$shadow.days.rho))>0.3),]
CHN.temp.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$temp.log10BF>0.5 & (abs(CHN.BF.all$temp.rho))>0.3),]
CHN.pcpt.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$pcpt.log10BF>0.5 & (abs(CHN.BF.all$pcpt.rho))>0.3),]
CHN.day.10cm.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$day.10cm.log10BF>0.5 & (abs(CHN.BF.all$day.10cm.rho))>0.3),]
```





####CHS

```
#calculate the log10 of BF for each environmental variable

CHS.BF.all$rad.log10BF <- log10(CHS.BF.all$rad.BF)
CHS.BF.all$shadow.days.log10BF <- log10(CHS.BF.all$shadow.days.BF)
CHS.BF.all$temp.log10BF <- log10(CHS.BF.all$temp.BF)
CHS.BF.all$pcpt.log10BF <- log10(CHS.BF.all$pcpt.BF)
CHS.BF.all$day.10cm.log10BF <- log10(CHS.BF.all$day.10cm.BF)

##Sort and plot

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$rad.log10BF),]

pdf("CHS.BF.plot.pdf")
par(mfrow=c(3,2))
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$rad.log10BF),]
plot(CHS.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$shadow.days.log10BF),]
plot(CHS.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$temp.log10BF),]
plot(CHS.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$pcpt.log10BF),]
plot(CHS.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$day.10cm.log10BF),]
plot(CHS.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHS.BFvsp.pdf")
par(mfrow=c(3,2))
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$rad.log10BF),]
plot(CHS.BF.all.sort$rad.log10BF~(abs(CHS.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$shadow.days.log10BF),]
plot(CHS.BF.all.sort$shadow.days.log10BF~(abs(CHS.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$temp.log10BF),]
plot(CHS.BF.all.sort$temp.log10BF~(abs(CHS.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$pcpt.log10BF),]
plot(CHS.BF.all.sort$pcpt.log10BF~(abs(CHS.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$day.10cm.log10BF),]
plot(CHS.BF.all.sort$day.10cm.log10BF~(abs(CHS.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHS.rad.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$rad.log10BF>0.5 & (abs(CHS.BF.all$rad.rho))>0.3),]
CHS.shadow.days.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.BF.all$shadow.days.rho))>0.3),]
CHS.temp.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$temp.log10BF>0.5 & (abs(CHS.BF.all$temp.rho))>0.3),]
CHS.pcpt.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$pcpt.log10BF>0.5 & (abs(CHS.BF.all$pcpt.rho))>0.3),]
CHS.day.10cm.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.BF.all$day.10cm.rho))>0.3),]
```





####CZ

```
#calculate the log10 of BF for each environmental variable

CZ.BF.all$rad.log10BF <- log10(CZ.BF.all$rad.BF)
CZ.BF.all$shadow.days.log10BF <- log10(CZ.BF.all$shadow.days.BF)
CZ.BF.all$temp.log10BF <- log10(CZ.BF.all$temp.BF)
CZ.BF.all$pcpt.log10BF <- log10(CZ.BF.all$pcpt.BF)
CZ.BF.all$day.10cm.log10BF <- log10(CZ.BF.all$day.10cm.BF)

##Sort and plot

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$rad.log10BF),]

pdf("CZ.BF.plot.pdf")
par(mfrow=c(3,2))
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$rad.log10BF),]
plot(CZ.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$shadow.days.log10BF),]
plot(CZ.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$temp.log10BF),]
plot(CZ.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$pcpt.log10BF),]
plot(CZ.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$day.10cm.log10BF),]
plot(CZ.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CZ.BFvsp.pdf")
par(mfrow=c(3,2))
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$rad.log10BF),]
plot(CZ.BF.all.sort$rad.log10BF~(abs(CZ.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$shadow.days.log10BF),]
plot(CZ.BF.all.sort$shadow.days.log10BF~(abs(CZ.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$temp.log10BF),]
plot(CZ.BF.all.sort$temp.log10BF~(abs(CZ.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$pcpt.log10BF),]
plot(CZ.BF.all.sort$pcpt.log10BF~(abs(CZ.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$day.10cm.log10BF),]
plot(CZ.BF.all.sort$day.10cm.log10BF~(abs(CZ.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CZ.rad.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$rad.log10BF>0.5 & (abs(CZ.BF.all$rad.rho))>0.3),]
CZ.shadow.days.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$shadow.days.log10BF>0.5 & (abs(CZ.BF.all$shadow.days.rho))>0.3),]
CZ.temp.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$temp.log10BF>0.5 & (abs(CZ.BF.all$temp.rho))>0.3),]
CZ.pcpt.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$pcpt.log10BF>0.5 & (abs(CZ.BF.all$pcpt.rho))>0.3),]
CZ.day.10cm.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$day.10cm.log10BF>0.5 & (abs(CZ.BF.all$day.10cm.rho))>0.3),]
```









####CHS.VS
```
#calculate the log10 of BF for each environmental variable

CHS.VS.BF.all$rad.log10BF <- log10(CHS.VS.BF.all$rad.BF)
CHS.VS.BF.all$shadow.days.log10BF <- log10(CHS.VS.BF.all$shadow.days.BF)
CHS.VS.BF.all$temp.log10BF <- log10(CHS.VS.BF.all$temp.BF)
CHS.VS.BF.all$pcpt.log10BF <- log10(CHS.VS.BF.all$pcpt.BF)
CHS.VS.BF.all$day.10cm.log10BF <- log10(CHS.VS.BF.all$day.10cm.BF)

##Sort and plot

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$rad.log10BF),]

pdf("CHS.VS.BF.plot.pdf")
par(mfrow=c(3,2))
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$rad.log10BF),]
plot(CHS.VS.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$shadow.days.log10BF),]
plot(CHS.VS.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$temp.log10BF),]
plot(CHS.VS.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$pcpt.log10BF),]
plot(CHS.VS.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$day.10cm.log10BF),]
plot(CHS.VS.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHS.VS.BFvsp.pdf")
par(mfrow=c(3,2))
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$rad.log10BF),]
plot(CHS.VS.BF.all.sort$rad.log10BF~(abs(CHS.VS.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$shadow.days.log10BF),]
plot(CHS.VS.BF.all.sort$shadow.days.log10BF~(abs(CHS.VS.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$temp.log10BF),]
plot(CHS.VS.BF.all.sort$temp.log10BF~(abs(CHS.VS.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$pcpt.log10BF),]
plot(CHS.VS.BF.all.sort$pcpt.log10BF~(abs(CHS.VS.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$day.10cm.log10BF),]
plot(CHS.VS.BF.all.sort$day.10cm.log10BF~(abs(CHS.VS.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHS.VS.rad.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$rad.log10BF>0.5 & (abs(CHS.VS.BF.all$rad.rho))>0.3),]
CHS.VS.shadow.days.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.VS.BF.all$shadow.days.rho))>0.3),]
CHS.VS.temp.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$temp.log10BF>0.5 & (abs(CHS.VS.BF.all$temp.rho))>0.3),]
CHS.VS.pcpt.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$pcpt.log10BF>0.5 & (abs(CHS.VS.BF.all$pcpt.rho))>0.3),]
CHS.VS.day.10cm.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.VS.BF.all$day.10cm.rho))>0.3),]
```




####CHS.TI
```
#calculate the log10 of BF for each environmental variable

CHS.TI.BF.all$rad.log10BF <- log10(CHS.TI.BF.all$rad.BF)
CHS.TI.BF.all$shadow.days.log10BF <- log10(CHS.TI.BF.all$shadow.days.BF)
CHS.TI.BF.all$temp.log10BF <- log10(CHS.TI.BF.all$temp.BF)
CHS.TI.BF.all$pcpt.log10BF <- log10(CHS.TI.BF.all$pcpt.BF)
CHS.TI.BF.all$day.10cm.log10BF <- log10(CHS.TI.BF.all$day.10cm.BF)

##Sort and plot

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$rad.log10BF),]

pdf("CHS.TI.BF.plot.pdf")
par(mfrow=c(3,2))
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$rad.log10BF),]
plot(CHS.TI.BF.all.sort$rad.log10BF, main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$shadow.days.log10BF),]
plot(CHS.TI.BF.all.sort$shadow.days.log10BF, main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$temp.log10BF),]
plot(CHS.TI.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)


CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$pcpt.log10BF),]
plot(CHS.TI.BF.all.sort$pcpt.log10BF, main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$day.10cm.log10BF),]
plot(CHS.TI.BF.all.sort$day.10cm.log10BF, main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHS.TI.BFvsp.pdf")
par(mfrow=c(3,2))
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$rad.log10BF),]
plot(CHS.TI.BF.all.sort$rad.log10BF~(abs(CHS.TI.BF.all.sort$rad.rho)), main="rad")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$shadow.days.log10BF),]
plot(CHS.TI.BF.all.sort$shadow.days.log10BF~(abs(CHS.TI.BF.all.sort$shadow.days.rho)), main="shadow.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$temp.log10BF),]
plot(CHS.TI.BF.all.sort$temp.log10BF~(abs(CHS.TI.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$pcpt.log10BF),]
plot(CHS.TI.BF.all.sort$pcpt.log10BF~(abs(CHS.TI.BF.all.sort$pcpt.rho)), main="pcpt")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$day.10cm.log10BF),]
plot(CHS.TI.BF.all.sort$day.10cm.log10BF~(abs(CHS.TI.BF.all.sort$day.10cm.rho)), main="day.10cm")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHS.TI.rad.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$rad.log10BF>0.5 & (abs(CHS.TI.BF.all$rad.rho))>0.3),]
CHS.TI.shadow.days.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.TI.BF.all$shadow.days.rho))>0.3),]
CHS.TI.temp.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$temp.log10BF>0.5 & (abs(CHS.TI.BF.all$temp.rho))>0.3),]
CHS.TI.pcpt.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$pcpt.log10BF>0.5 & (abs(CHS.TI.BF.all$pcpt.rho))>0.3),]
CHS.TI.day.10cm.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.TI.BF.all$day.10cm.rho))>0.3),]
```




###Determine the names of all the loci. 

My input was obtained using pgdspider 2.1.0.3 (vcf -> bayenv2)

However, this removes all monomorphic loci from the input file. Thus, to match this back to the vcf file, I need to remove all monomorphic loci from this file. 

```
vcftools --vcf CHall.932.9608.recode.vcf --maf 0.0001 --recode --recode-INFO-all --out CHall.932.

vcftools --vcf CHN.229.9608.recode.vcf --maf 0.000001 --recode --recode-INFO-all --out CHN.229.8951

vcftools --vcf CHS.275.9608.recode.vcf --maf 0.000001 --recode --recode-INFO-all --out CHS.275.9295

vcftools --vcf CZ.404.9608.recode.vcf --maf 0.000001 --recode --recode-INFO-all --out CZ.404.9528

vcftools --vcf CHS.VS.135.9608.recode.vcf --maf 0.000001 --recode --recode-INFO-all --out CHS.VS.135.8779

vcftools --vcf CHS.TI.140.9608.recode.vcf --maf 0.000001 --recode --recode-INFO-all --out CHS.TI.140.8343
```
Check that this is the right number of loci by comparing this to the number of snp_batch files found when the .txt input file is split (see above). Remember that these are numbered from 0. 

Next, convert to plink. 

Once this is done, the names can be read into R and indexed. 

```
#linux

vcftools --vcf CHall.932.9608.recode.vcf --plink --out CHall.932.9608.plink

vcftools --vcf CHN.229.8951.recode.vcf --plink --out CHN.229.8951.plink

vcftools --vcf CHS.275.9295.recode.vcf --plink --out CHS.275.9295.plink

vcftools --vcf CZ.404.9528.recode.vcf --plink --out CZ.404.9528.plink

vcftools --vcf CHS.VS.135.8779.recode.vcf --plink --out CHS.VS.135.8779.plink

vcftools --vcf CHS.TI.140.8343.recode.vcf --plink --out CHS.TI.140.8343.plink

```



```
##R

###CHall
CHall.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHall.932.9608.plink.map", header=F) #import the locus names into R
CHall.locus.names$ID <- seq.int(nrow(CHall.locus.names)) ##index the CHall.locus.names file so that all the loci are numbered in order of appearance

CHall.BF.all$ID <- seq.int(nrow(CHall.BF.all))  ##do the same with the CHall.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHall.BF.all file

##Find all the candidate loci
CHall.rad.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$rad.log10BF>0.5 & (abs(CHall.BF.all$rad.rho))>0.3),]
CHall.shadow.days.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$shadow.days.log10BF>0.5 & (abs(CHall.BF.all$shadow.days.rho))>0.3),]
CHall.temp.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$temp.log10BF>0.5 & (abs(CHall.BF.all$temp.rho))>0.3),]
CHall.pcpt.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$pcpt.log10BF>0.5 & (abs(CHall.BF.all$pcpt.rho))>0.3),]
CHall.day.10cm.bayenv.candidates <- CHall.BF.all[which(CHall.BF.all$day.10cm.log10BF>0.5 & (abs(CHall.BF.all$day.10cm.rho))>0.3),]

CHall.rad.bayenv.candidates$ID <- as.character(CHall.rad.bayenv.candidates$ID)
CHall.rad.bayenv.candidates.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHall.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHall.rad.bayenv.candidates.names <- paste("X", CHall.rad.bayenv.candidates.names$SNP, sep=".")

CHall.shadow.days.bayenv.candidates$ID <- as.character(CHall.shadow.days.bayenv.candidates$ID)
CHall.shadow.days.bayenv.candidates.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHall.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHall.shadow.days.bayenv.candidates.names <- paste("X", CHall.shadow.days.bayenv.candidates.names$SNP, sep=".")

CHall.temp.bayenv.candidates$ID <- as.character(CHall.temp.bayenv.candidates$ID)
CHall.temp.bayenv.candidates.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHall.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHall.temp.bayenv.candidates.names <- paste("X", CHall.temp.bayenv.candidates.names$SNP, sep=".")

CHall.pcpt.bayenv.candidates$ID <- as.character(CHall.pcpt.bayenv.candidates$ID)
CHall.pcpt.bayenv.candidates.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHall.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHall.pcpt.bayenv.candidates.names <- paste("X", CHall.pcpt.bayenv.candidates.names$SNP, sep=".")

CHall.day.10cm.bayenv.candidates$ID <- as.character(CHall.day.10cm.bayenv.candidates$ID)
CHall.day.10cm.bayenv.candidates.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHall.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHall.day.10cm.bayenv.candidates.names <- paste("X", CHall.day.10cm.bayenv.candidates.names$SNP, sep=".")

```
find names and write to file

```
CHall.rad.names <- as.data.frame(CHall.rad.bayenv.candidates.names)
colnames(CHall.rad.names) <- "names"

CHall.shadow.days.names <- as.data.frame(CHall.shadow.days.bayenv.candidates.names)
colnames(CHall.shadow.days.names) <- "names"

CHall.temp.names <- as.data.frame(CHall.temp.bayenv.candidates.names)
colnames(CHall.temp.names) <- "names"

CHall.pcpt.names <- as.data.frame(CHall.pcpt.bayenv.candidates.names)
colnames(CHall.pcpt.names) <- "names"

CHall.day.10cm.names <- as.data.frame(CHall.day.10cm.bayenv.candidates.names)
colnames(CHall.day.10cm.names) <- "names"


CHall.alloutliers.names <- rbind(CHall.rad.names, CHall.shadow.days.names, CHall.temp.names, CHall.pcpt.names, CHall.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHall.alloutliers.names <- lapply(CHall.alloutliers.names, unique)  #select only the unique rows. 

CHall.alloutliers.names <- sub(":", ".", CHall.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHall.alloutliers.names <- as.data.frame(CHall.alloutliers.names)
write.table(CHall.alloutliers.names$CHall.alloutliers.names, "CHall.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CHall.BayEnv.outliers <- read.table("CHall.BayEnv.alloutliers")
colnames(CHall.BayEnv.outliers) <- ("loci")
CHall.BayEnv.outliers <- as.character(CHall.BayEnv.outliers$loci)

```

Write XtX outliers top 100 to file

```
CHall.XtX.run1 <- read.table("CHall.Run1/XtX_out.CHall.5n.ENV", header=F)
CHall.XtX.run2 <- read.table("CHall.Run2/XtX_out.CHall.5n.ENV", header=F)
CHall.XtX.run3 <- read.table("CHall.Run3/XtX_out.CHall.5n.ENV", header=F)


CHall.XtX.all <- rbindlist(list(CHall.XtX.run1, CHall.XtX.run2, CHall.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CHall.XtX.all) <- c("snp", "XtX")
##index the loci 
CHall.XtX.all$ID <- seq.int(nrow(CHall.XtX.all))   ##index the loci in the original output
CHall.XtX.all.sort <- CHall.XtX.all[order(-CHall.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CHall.XtX.top100 <- CHall.XtX.all.sort[1:100,]  #select the first 100 loci
CHall.XtX.top100.names <- CHall.locus.names[CHall.locus.names$ID %in% CHall.XtX.top100$ID,]  ##get their names
summary(CHall.XtX.top100.names)
colnames(CHall.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CHall.XtX.top100.names <- paste("X", CHall.XtX.top100.names$SNP, sep=".")  #rename loci
CHall.XtX.top100.names <- as.data.frame(CHall.XtX.top100.names)
CHall.XtX.top100.names
colnames(CHall.XtX.top100.names) <- "names"
CHall.XtX.top100.names <- sub(":", ".", CHall.XtX.top100.names$names)  ##rename to match Fst table
CHall.XtX.top100.names <- as.data.frame(CHall.XtX.top100.names)
write.table(CHall.XtX.top100.names$CHall.XtX.top100.names, "CHall.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```


###### CHN
```
CHN.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHN.229.8951.plink.map", header=F) #import the locus names into R
CHN.locus.names$ID <- seq.int(nrow(CHN.locus.names)) ##index the CHN.locus.names file so that all the loci are numbered in order of appearance

CHN.BF.all$ID <- seq.int(nrow(CHN.BF.all))  ##do the same with the CHN.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHN.BF.all file

##Find all the candidate loci
CHN.rad.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$rad.log10BF>0.5 & (abs(CHN.BF.all$rad.rho))>0.3),]
CHN.shadow.days.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$shadow.days.log10BF>0.5 & (abs(CHN.BF.all$shadow.days.rho))>0.3),]
CHN.temp.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$temp.log10BF>0.5 & (abs(CHN.BF.all$temp.rho))>0.3),]
CHN.pcpt.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$pcpt.log10BF>0.5 & (abs(CHN.BF.all$pcpt.rho))>0.3),]
CHN.day.10cm.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$day.10cm.log10BF>0.5 & (abs(CHN.BF.all$day.10cm.rho))>0.3),]

CHN.rad.bayenv.candidates$ID <- as.character(CHN.rad.bayenv.candidates$ID)
CHN.rad.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.rad.bayenv.candidates.names <- paste("X", CHN.rad.bayenv.candidates.names$SNP, sep=".")

CHN.shadow.days.bayenv.candidates$ID <- as.character(CHN.shadow.days.bayenv.candidates$ID)
CHN.shadow.days.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.shadow.days.bayenv.candidates.names <- paste("X", CHN.shadow.days.bayenv.candidates.names$SNP, sep=".")

CHN.temp.bayenv.candidates$ID <- as.character(CHN.temp.bayenv.candidates$ID)
CHN.temp.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.temp.bayenv.candidates.names <- paste("X", CHN.temp.bayenv.candidates.names$SNP, sep=".")

CHN.pcpt.bayenv.candidates$ID <- as.character(CHN.pcpt.bayenv.candidates$ID)
CHN.pcpt.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.pcpt.bayenv.candidates.names <- paste("X", CHN.pcpt.bayenv.candidates.names$SNP, sep=".")

CHN.day.10cm.bayenv.candidates$ID <- as.character(CHN.day.10cm.bayenv.candidates$ID)
CHN.day.10cm.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.day.10cm.bayenv.candidates.names <- paste("X", CHN.day.10cm.bayenv.candidates.names$SNP, sep=".")

```

Find names and write to file
```
CHN.rad.names <- as.data.frame(CHN.rad.bayenv.candidates.names)
colnames(CHN.rad.names) <- "names"

CHN.shadow.days.names <- as.data.frame(CHN.shadow.days.bayenv.candidates.names)
colnames(CHN.shadow.days.names) <- "names"

CHN.temp.names <- as.data.frame(CHN.temp.bayenv.candidates.names)
colnames(CHN.temp.names) <- "names"

CHN.pcpt.names <- as.data.frame(CHN.pcpt.bayenv.candidates.names)
colnames(CHN.pcpt.names) <- "names"

CHN.day.10cm.names <- as.data.frame(CHN.day.10cm.bayenv.candidates.names)
colnames(CHN.day.10cm.names) <- "names"


CHN.alloutliers.names <- rbind(CHN.rad.names, CHN.shadow.days.names, CHN.temp.names, CHN.pcpt.names, CHN.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHN.alloutliers.names <- lapply(CHN.alloutliers.names, unique)  #select only the unique rows. 

CHN.alloutliers.names <- sub(":", ".", CHN.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHN.alloutliers.names <- as.data.frame(CHN.alloutliers.names)
write.table(CHN.alloutliers.names$CHN.alloutliers.names, "CHN.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CHN.BayEnv.outliers <- read.table("CHN.BayEnv.alloutliers")
colnames(CHN.BayEnv.outliers) <- ("loci")
CHN.BayEnv.outliers <- as.character(CHN.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHN.XtX.run1 <- read.table("CHN.Run1/XtX_out.ENVIRONFILE.CHN.p19.env", header=F)
CHN.XtX.run2 <- read.table("CHN.Run2/XtX_out.ENVIRONFILE.CHN.p19.env", header=F)
CHN.XtX.run3 <- read.table("CHN.Run3/XtX_out.ENVIRONFILE.CHN.p19.env", header=F)


CHN.XtX.all <- rbindlist(list(CHN.XtX.run1, CHN.XtX.run2, CHN.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CHN.XtX.all) <- c("snp", "XtX")
##index the loci 
CHN.XtX.all$ID <- seq.int(nrow(CHN.XtX.all))   ##index the loci in the original output
CHN.XtX.all.sort <- CHN.XtX.all[order(-CHN.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CHN.XtX.top100 <- CHN.XtX.all.sort[1:100,]  #select the first 100 loci
CHN.XtX.top100.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.XtX.top100$ID,]  ##get their names
summary(CHN.XtX.top100.names)
colnames(CHN.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CHN.XtX.top100.names <- paste("X", CHN.XtX.top100.names$SNP, sep=".")  #rename loci
CHN.XtX.top100.names <- as.data.frame(CHN.XtX.top100.names)
CHN.XtX.top100.names
colnames(CHN.XtX.top100.names) <- "names"
CHN.XtX.top100.names <- sub(":", ".", CHN.XtX.top100.names$names)  ##rename to match Fst table
CHN.XtX.top100.names <- as.data.frame(CHN.XtX.top100.names)
write.table(CHN.XtX.top100.names$CHN.XtX.top100.names, "CHN.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```


########  CHS
```
CHS.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHS.275.9295.plink.map", header=F) #import the locus names into R
CHS.locus.names$ID <- seq.int(nrow(CHS.locus.names)) ##index the CHS.locus.names file so that all the loci are numbered in order of appearance

CHS.BF.all$ID <- seq.int(nrow(CHS.BF.all))  ##do the same with the CHS.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.BF.all file

##Find all the candidate loci
CHS.rad.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$rad.log10BF>0.5 & (abs(CHS.BF.all$rad.rho))>0.3),]
CHS.shadow.days.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.BF.all$shadow.days.rho))>0.3),]
CHS.temp.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$temp.log10BF>0.5 & (abs(CHS.BF.all$temp.rho))>0.3),]
CHS.pcpt.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$pcpt.log10BF>0.5 & (abs(CHS.BF.all$pcpt.rho))>0.3),]
CHS.day.10cm.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.BF.all$day.10cm.rho))>0.3),]

CHS.rad.bayenv.candidates$ID <- as.character(CHS.rad.bayenv.candidates$ID)
CHS.rad.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.rad.bayenv.candidates.names <- paste("X", CHS.rad.bayenv.candidates.names$SNP, sep=".")

CHS.shadow.days.bayenv.candidates$ID <- as.character(CHS.shadow.days.bayenv.candidates$ID)
CHS.shadow.days.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.shadow.days.bayenv.candidates.names <- paste("X", CHS.shadow.days.bayenv.candidates.names$SNP, sep=".")

CHS.temp.bayenv.candidates$ID <- as.character(CHS.temp.bayenv.candidates$ID)
CHS.temp.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.temp.bayenv.candidates.names <- paste("X", CHS.temp.bayenv.candidates.names$SNP, sep=".")

CHS.pcpt.bayenv.candidates$ID <- as.character(CHS.pcpt.bayenv.candidates$ID)
CHS.pcpt.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.pcpt.bayenv.candidates.names <- paste("X", CHS.pcpt.bayenv.candidates.names$SNP, sep=".")

CHS.day.10cm.bayenv.candidates$ID <- as.character(CHS.day.10cm.bayenv.candidates$ID)
CHS.day.10cm.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.day.10cm.bayenv.candidates.names <- paste("X", CHS.day.10cm.bayenv.candidates.names$SNP, sep=".")

```

find names and write to file
```
CHS.rad.names <- as.data.frame(CHS.rad.bayenv.candidates.names)
colnames(CHS.rad.names) <- "names"

CHS.shadow.days.names <- as.data.frame(CHS.shadow.days.bayenv.candidates.names)
colnames(CHS.shadow.days.names) <- "names"

CHS.temp.names <- as.data.frame(CHS.temp.bayenv.candidates.names)
colnames(CHS.temp.names) <- "names"

CHS.pcpt.names <- as.data.frame(CHS.pcpt.bayenv.candidates.names)
colnames(CHS.pcpt.names) <- "names"

CHS.day.10cm.names <- as.data.frame(CHS.day.10cm.bayenv.candidates.names)
colnames(CHS.day.10cm.names) <- "names"


CHS.alloutliers.names <- rbind(CHS.rad.names, CHS.shadow.days.names, CHS.temp.names, CHS.pcpt.names, CHS.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.alloutliers.names <- lapply(CHS.alloutliers.names, unique)  #select only the unique rows. 

CHS.alloutliers.names <- sub(":", ".", CHS.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.alloutliers.names <- as.data.frame(CHS.alloutliers.names)
write.table(CHS.alloutliers.names$CHS.alloutliers.names, "CHS.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CHS.BayEnv.outliers <- read.table("CHS.BayEnv.alloutliers")
colnames(CHS.BayEnv.outliers) <- ("loci")
CHS.BayEnv.outliers <- as.character(CHS.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHS.XtX.run1 <- read.table("CHS.Run1/XtX_out.ENVIRONFILE.CHS.n5.env", header=F)
CHS.XtX.run2 <- read.table("CHS.Run2/XtX_out.ENVIRONFILE.CHS.n5.env", header=F)
CHS.XtX.run3 <- read.table("CHS.Run3/XtX_out.ENVIRONFILE.CHS.n5.env", header=F)


CHS.XtX.all <- rbindlist(list(CHS.XtX.run1, CHS.XtX.run2, CHS.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CHS.XtX.all) <- c("snp", "XtX")
##index the loci 
CHS.XtX.all$ID <- seq.int(nrow(CHS.XtX.all))   ##index the loci in the original output
CHS.XtX.all.sort <- CHS.XtX.all[order(-CHS.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CHS.XtX.top100 <- CHS.XtX.all.sort[1:100,]  #select the first 100 loci
CHS.XtX.top100.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.XtX.top100$ID,]  ##get their names
summary(CHS.XtX.top100.names)
colnames(CHS.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CHS.XtX.top100.names <- paste("X", CHS.XtX.top100.names$SNP, sep=".")  #rename loci
CHS.XtX.top100.names <- as.data.frame(CHS.XtX.top100.names)
CHS.XtX.top100.names
colnames(CHS.XtX.top100.names) <- "names"
CHS.XtX.top100.names <- sub(":", ".", CHS.XtX.top100.names$names)  ##rename to match Fst table
CHS.XtX.top100.names <- as.data.frame(CHS.XtX.top100.names)
write.table(CHS.XtX.top100.names$CHS.XtX.top100.names, "CHS.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```


####### CZ
```
CZ.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CZ.404.9528.plink.map", header=F) #import the locus names into R
CZ.locus.names$ID <- seq.int(nrow(CZ.locus.names)) ##index the CZ.locus.names file so that all the loci are numbered in order of appearance

CZ.BF.all$ID <- seq.int(nrow(CZ.BF.all))  ##do the same with the CZ.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CZ.BF.all file

##Find all the candidate loci
CZ.rad.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$rad.log10BF>0.5 & (abs(CZ.BF.all$rad.rho))>0.3),]
CZ.shadow.days.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$shadow.days.log10BF>0.5 & (abs(CZ.BF.all$shadow.days.rho))>0.3),]
CZ.temp.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$temp.log10BF>0.5 & (abs(CZ.BF.all$temp.rho))>0.3),]
CZ.pcpt.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$pcpt.log10BF>0.5 & (abs(CZ.BF.all$pcpt.rho))>0.3),]
CZ.day.10cm.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$day.10cm.log10BF>0.5 & (abs(CZ.BF.all$day.10cm.rho))>0.3),]

CZ.rad.bayenv.candidates$ID <- as.character(CZ.rad.bayenv.candidates$ID)
CZ.rad.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.rad.bayenv.candidates.names <- paste("X", CZ.rad.bayenv.candidates.names$SNP, sep=".")

CZ.shadow.days.bayenv.candidates$ID <- as.character(CZ.shadow.days.bayenv.candidates$ID)
CZ.shadow.days.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.shadow.days.bayenv.candidates.names <- paste("X", CZ.shadow.days.bayenv.candidates.names$SNP, sep=".")

CZ.temp.bayenv.candidates$ID <- as.character(CZ.temp.bayenv.candidates$ID)
CZ.temp.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.temp.bayenv.candidates.names <- paste("X", CZ.temp.bayenv.candidates.names$SNP, sep=".")

CZ.pcpt.bayenv.candidates$ID <- as.character(CZ.pcpt.bayenv.candidates$ID)
CZ.pcpt.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.pcpt.bayenv.candidates.names <- paste("X", CZ.pcpt.bayenv.candidates.names$SNP, sep=".")

CZ.day.10cm.bayenv.candidates$ID <- as.character(CZ.day.10cm.bayenv.candidates$ID)
CZ.day.10cm.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.day.10cm.bayenv.candidates.names <- paste("X", CZ.day.10cm.bayenv.candidates.names$SNP, sep=".")

```
find names and write to file
```
CZ.rad.names <- as.data.frame(CZ.rad.bayenv.candidates.names)
colnames(CZ.rad.names) <- "names"

CZ.shadow.days.names <- as.data.frame(CZ.shadow.days.bayenv.candidates.names)
colnames(CZ.shadow.days.names) <- "names"

CZ.temp.names <- as.data.frame(CZ.temp.bayenv.candidates.names)
colnames(CZ.temp.names) <- "names"

CZ.pcpt.names <- as.data.frame(CZ.pcpt.bayenv.candidates.names)
colnames(CZ.pcpt.names) <- "names"

CZ.day.10cm.names <- as.data.frame(CZ.day.10cm.bayenv.candidates.names)
colnames(CZ.day.10cm.names) <- "names"


CZ.alloutliers.names <- rbind(CZ.rad.names, CZ.shadow.days.names, CZ.temp.names, CZ.pcpt.names, CZ.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CZ.alloutliers.names <- lapply(CZ.alloutliers.names, unique)  #select only the unique rows. 

CZ.alloutliers.names <- sub(":", ".", CZ.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CZ.alloutliers.names <- as.data.frame(CZ.alloutliers.names)
write.table(CZ.alloutliers.names$CZ.alloutliers.names, "CZ.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CZ.BayEnv.outliers <- read.table("CZ.BayEnv.alloutliers")
colnames(CZ.BayEnv.outliers) <- ("loci")
CZ.BayEnv.outliers <- as.character(CZ.BayEnv.outliers$loci)
```


Write top 100 XtX loci to file

```
CZ.XtX.run1 <- read.table("CZ.Run1/XtX_out.ENVIRON.CZ.5n.env", header=F)
CZ.XtX.run2 <- read.table("CZ.Run2/XtX_out.ENVIRON.CZ.5n.env", header=F)
CZ.XtX.run3 <- read.table("CZ.Run3/XtX_out.ENVIRON.CZ.5n.env", header=F)


CZ.XtX.all <- rbindlist(list(CZ.XtX.run1, CZ.XtX.run2, CZ.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CZ.XtX.all) <- c("snp", "XtX")
##index the loci 
CZ.XtX.all$ID <- seq.int(nrow(CZ.XtX.all))   ##index the loci in the original output
CZ.XtX.all.sort <- CZ.XtX.all[order(-CZ.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CZ.XtX.top100 <- CZ.XtX.all.sort[1:100,]  #select the first 100 loci
CZ.XtX.top100.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.XtX.top100$ID,]  ##get their names
summary(CZ.XtX.top100.names)
colnames(CZ.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CZ.XtX.top100.names <- paste("X", CZ.XtX.top100.names$SNP, sep=".")  #rename loci
CZ.XtX.top100.names <- as.data.frame(CZ.XtX.top100.names)
CZ.XtX.top100.names
colnames(CZ.XtX.top100.names) <- "names"
CZ.XtX.top100.names <- sub(":", ".", CZ.XtX.top100.names$names)  ##rename to match Fst table
CZ.XtX.top100.names <- as.data.frame(CZ.XtX.top100.names)
write.table(CZ.XtX.top100.names$CZ.XtX.top100.names, "CZ.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```




####### CHS.VS
```
CHS.VS.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHS.VS.135.8779.plink.map", header=F) #import the locus names into R
CHS.VS.locus.names$ID <- seq.int(nrow(CHS.VS.locus.names)) ##index the CHS.VS.locus.names file so that all the loci are numbered in order of appearance

CHS.VS.BF.all$ID <- seq.int(nrow(CHS.VS.BF.all))  ##do the same with the CHS.VS.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.VS.BF.all file

##Find all the candidate loci
CHS.VS.rad.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$rad.log10BF>0.5 & (abs(CHS.VS.BF.all$rad.rho))>0.3),]
CHS.VS.shadow.days.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.VS.BF.all$shadow.days.rho))>0.3),]
CHS.VS.temp.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$temp.log10BF>0.5 & (abs(CHS.VS.BF.all$temp.rho))>0.3),]
CHS.VS.pcpt.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$pcpt.log10BF>0.5 & (abs(CHS.VS.BF.all$pcpt.rho))>0.3),]
CHS.VS.day.10cm.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.VS.BF.all$day.10cm.rho))>0.3),]

CHS.VS.rad.bayenv.candidates$ID <- as.character(CHS.VS.rad.bayenv.candidates$ID)
CHS.VS.rad.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.rad.bayenv.candidates.names <- paste("X", CHS.VS.rad.bayenv.candidates.names$SNP, sep=".")

CHS.VS.shadow.days.bayenv.candidates$ID <- as.character(CHS.VS.shadow.days.bayenv.candidates$ID)
CHS.VS.shadow.days.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.shadow.days.bayenv.candidates.names <- paste("X", CHS.VS.shadow.days.bayenv.candidates.names$SNP, sep=".")

CHS.VS.temp.bayenv.candidates$ID <- as.character(CHS.VS.temp.bayenv.candidates$ID)
CHS.VS.temp.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.temp.bayenv.candidates.names <- paste("X", CHS.VS.temp.bayenv.candidates.names$SNP, sep=".")

CHS.VS.pcpt.bayenv.candidates$ID <- as.character(CHS.VS.pcpt.bayenv.candidates$ID)
CHS.VS.pcpt.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.pcpt.bayenv.candidates.names <- paste("X", CHS.VS.pcpt.bayenv.candidates.names$SNP, sep=".")

CHS.VS.day.10cm.bayenv.candidates$ID <- as.character(CHS.VS.day.10cm.bayenv.candidates$ID)
CHS.VS.day.10cm.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.day.10cm.bayenv.candidates.names <- paste("X", CHS.VS.day.10cm.bayenv.candidates.names$SNP, sep=".")

```
find names and write to file

```
CHS.VS.rad.names <- as.data.frame(CHS.VS.rad.bayenv.candidates.names)
colnames(CHS.VS.rad.names) <- "names"

CHS.VS.shadow.days.names <- as.data.frame(CHS.VS.shadow.days.bayenv.candidates.names)
colnames(CHS.VS.shadow.days.names) <- "names"

CHS.VS.temp.names <- as.data.frame(CHS.VS.temp.bayenv.candidates.names)
colnames(CHS.VS.temp.names) <- "names"

CHS.VS.pcpt.names <- as.data.frame(CHS.VS.pcpt.bayenv.candidates.names)
colnames(CHS.VS.pcpt.names) <- "names"

CHS.VS.day.10cm.names <- as.data.frame(CHS.VS.day.10cm.bayenv.candidates.names)
colnames(CHS.VS.day.10cm.names) <- "names"


CHS.VS.alloutliers.names <- rbind(CHS.VS.rad.names, CHS.VS.shadow.days.names, CHS.VS.temp.names, CHS.VS.pcpt.names, CHS.VS.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.VS.alloutliers.names <- lapply(CHS.VS.alloutliers.names, unique)  #select only the unique rows. 

CHS.VS.alloutliers.names <- sub(":", ".", CHS.VS.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.VS.alloutliers.names <- as.data.frame(CHS.VS.alloutliers.names)
write.table(CHS.VS.alloutliers.names$CHS.VS.alloutliers.names, "CHS.VS.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CHS.VS.BayEnv.outliers <- read.table("CHS.VS.BayEnv.alloutliers")
colnames(CHS.VS.BayEnv.outliers) <- ("loci")
CHS.VS.BayEnv.outliers <- as.character(CHS.VS.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHS.VS.XtX.run1 <- read.table("CHS.VS.Run1/XtX_out.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)
CHS.VS.XtX.run2 <- read.table("CHS.VS.Run2/XtX_out.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)
CHS.VS.XtX.run3 <- read.table("CHS.VS.Run3/XtX_out.CHS.VS.p10.n5.ENVIRONFILE.env", header=F)


CHS.VS.XtX.all <- rbindlist(list(CHS.VS.XtX.run1, CHS.VS.XtX.run2, CHS.VS.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CHS.VS.XtX.all) <- c("snp", "XtX")
##index the loci 
CHS.VS.XtX.all$ID <- seq.int(nrow(CHS.VS.XtX.all))   ##index the loci in the original output
CHS.VS.XtX.all.sort <- CHS.VS.XtX.all[order(-CHS.VS.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CHS.VS.XtX.top100 <- CHS.VS.XtX.all.sort[1:100,]  #select the first 100 loci
CHS.VS.XtX.top100.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.XtX.top100$ID,]  ##get their names
summary(CHS.VS.XtX.top100.names)
colnames(CHS.VS.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CHS.VS.XtX.top100.names <- paste("X", CHS.VS.XtX.top100.names$SNP, sep=".")  #rename loci
CHS.VS.XtX.top100.names <- as.data.frame(CHS.VS.XtX.top100.names)
CHS.VS.XtX.top100.names
colnames(CHS.VS.XtX.top100.names) <- "names"
CHS.VS.XtX.top100.names <- sub(":", ".", CHS.VS.XtX.top100.names$names)  ##rename to match Fst table
CHS.VS.XtX.top100.names <- as.data.frame(CHS.VS.XtX.top100.names)
write.table(CHS.VS.XtX.top100.names$CHS.VS.XtX.top100.names, "CHS.VS.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```



####### CHS.TI
```
CHS.TI.locus.names <- read.table("/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHS.TI.140.8343.plink.map", header=F) #import the locus names into R
CHS.TI.locus.names$ID <- seq.int(nrow(CHS.TI.locus.names)) ##index the CHS.TI.locus.names file so that all the loci are numbered in order of appearance

CHS.TI.BF.all$ID <- seq.int(nrow(CHS.TI.BF.all))  ##do the same with the CHS.TI.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.TI.BF.all file

##Find all the candidate loci
CHS.TI.rad.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$rad.log10BF>0.5 & (abs(CHS.TI.BF.all$rad.rho))>0.3),]
CHS.TI.shadow.days.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$shadow.days.log10BF>0.5 & (abs(CHS.TI.BF.all$shadow.days.rho))>0.3),]
CHS.TI.temp.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$temp.log10BF>0.5 & (abs(CHS.TI.BF.all$temp.rho))>0.3),]
CHS.TI.pcpt.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$pcpt.log10BF>0.5 & (abs(CHS.TI.BF.all$pcpt.rho))>0.3),]
CHS.TI.day.10cm.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$day.10cm.log10BF>0.5 & (abs(CHS.TI.BF.all$day.10cm.rho))>0.3),]

CHS.TI.rad.bayenv.candidates$ID <- as.character(CHS.TI.rad.bayenv.candidates$ID)
CHS.TI.rad.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.rad.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.rad.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.rad.bayenv.candidates.names <- paste("X", CHS.TI.rad.bayenv.candidates.names$SNP, sep=".")

CHS.TI.shadow.days.bayenv.candidates$ID <- as.character(CHS.TI.shadow.days.bayenv.candidates$ID)
CHS.TI.shadow.days.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.shadow.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.shadow.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.shadow.days.bayenv.candidates.names <- paste("X", CHS.TI.shadow.days.bayenv.candidates.names$SNP, sep=".")

CHS.TI.temp.bayenv.candidates$ID <- as.character(CHS.TI.temp.bayenv.candidates$ID)
CHS.TI.temp.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.temp.bayenv.candidates.names <- paste("X", CHS.TI.temp.bayenv.candidates.names$SNP, sep=".")

CHS.TI.pcpt.bayenv.candidates$ID <- as.character(CHS.TI.pcpt.bayenv.candidates$ID)
CHS.TI.pcpt.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.pcpt.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.pcpt.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.pcpt.bayenv.candidates.names <- paste("X", CHS.TI.pcpt.bayenv.candidates.names$SNP, sep=".")

CHS.TI.day.10cm.bayenv.candidates$ID <- as.character(CHS.TI.day.10cm.bayenv.candidates$ID)
CHS.TI.day.10cm.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.day.10cm.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.day.10cm.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.day.10cm.bayenv.candidates.names <- paste("X", CHS.TI.day.10cm.bayenv.candidates.names$SNP, sep=".")

```

Find names and write to file
```
CHS.TI.rad.names <- as.data.frame(CHS.TI.rad.bayenv.candidates.names)
colnames(CHS.TI.rad.names) <- "names"

CHS.TI.shadow.days.names <- as.data.frame(CHS.TI.shadow.days.bayenv.candidates.names)
colnames(CHS.TI.shadow.days.names) <- "names"

CHS.TI.temp.names <- as.data.frame(CHS.TI.temp.bayenv.candidates.names)
colnames(CHS.TI.temp.names) <- "names"

CHS.TI.pcpt.names <- as.data.frame(CHS.TI.pcpt.bayenv.candidates.names)
colnames(CHS.TI.pcpt.names) <- "names"

CHS.TI.day.10cm.names <- as.data.frame(CHS.TI.day.10cm.bayenv.candidates.names)
colnames(CHS.TI.day.10cm.names) <- "names"


CHS.TI.alloutliers.names <- rbind(CHS.TI.rad.names, CHS.TI.shadow.days.names, CHS.TI.temp.names, CHS.TI.pcpt.names, CHS.TI.day.10cm.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.TI.alloutliers.names <- lapply(CHS.TI.alloutliers.names, unique)  #select only the unique rows. 

CHS.TI.alloutliers.names <- sub(":", ".", CHS.TI.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.TI.alloutliers.names <- as.data.frame(CHS.TI.alloutliers.names)
write.table(CHS.TI.alloutliers.names$CHS.TI.alloutliers.names, "CHS.TI.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CHS.TI.BayEnv.outliers <- read.table("CHS.TI.BayEnv.alloutliers")
colnames(CHS.TI.BayEnv.outliers) <- ("loci")
CHS.TI.BayEnv.outliers <- as.character(CHS.TI.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHS.TI.XtX.run1 <- read.table("CHS.TI.Run1/XtX_out.standardised.5env.env", header=F)
CHS.TI.XtX.run2 <- read.table("CHS.TI.Run2/XtX_out.standardised.5env.env", header=F)
CHS.TI.XtX.run3 <- read.table("CHS.TI.Run3/XtX_out.standardised.5env.env", header=F)


CHS.TI.XtX.all <- rbindlist(list(CHS.TI.XtX.run1, CHS.TI.XtX.run2, CHS.TI.XtX.run3))[,lapply(.SD,median), list(V1)]   ##combine by SNP and calculate the median XtX value

colnames(CHS.TI.XtX.all) <- c("snp", "XtX")
##index the loci 
CHS.TI.XtX.all$ID <- seq.int(nrow(CHS.TI.XtX.all))   ##index the loci in the original output
CHS.TI.XtX.all.sort <- CHS.TI.XtX.all[order(-CHS.TI.XtX.all$XtX),]  ##order to get the top 100 loci. The minus means its descending order
CHS.TI.XtX.top100 <- CHS.TI.XtX.all.sort[1:100,]  #select the first 100 loci
CHS.TI.XtX.top100.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.XtX.top100$ID,]  ##get their names
summary(CHS.TI.XtX.top100.names)
colnames(CHS.TI.XtX.top100.names) <- c("V1", "SNP", "V3", "V4", "ID")   #rename the columns to match the code
CHS.TI.XtX.top100.names <- paste("X", CHS.TI.XtX.top100.names$SNP, sep=".")  #rename loci
CHS.TI.XtX.top100.names <- as.data.frame(CHS.TI.XtX.top100.names)
CHS.TI.XtX.top100.names
colnames(CHS.TI.XtX.top100.names) <- "names"
CHS.TI.XtX.top100.names <- sub(":", ".", CHS.TI.XtX.top100.names$names)  ##rename to match Fst table
CHS.TI.XtX.top100.names <- as.data.frame(CHS.TI.XtX.top100.names)
write.table(CHS.TI.XtX.top100.names$CHS.TI.XtX.top100.names, "CHS.TI.XtX.100outliers", col.names=F, row.names=F, quote=F)  ##write the table


##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")
```
