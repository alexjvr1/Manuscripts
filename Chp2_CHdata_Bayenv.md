# Bayenv2: CHall data


I want to identify loci associated with environment for the 6 CH datasets. stba is removed from all datasets (due to missing env data),
hence 8 less individuals in CHall, CHS, and CHS.TI.

1. CHall.932.9608

2. CHN.229.9608

3. CHS.283.9608

4. CHS.VS.135.9608

5. CHS.TI.148.9608

6. CZ.404.9608

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

3. ENV input file with normalised environmental parameters. I will first use only temperature.

NB: population order in the input files should all be the same.


## 1. Covariance matrix estimation

Use a random set of 1000SNPs for matrix estimation. 

I estimated Fst distributions per locus for each of the six datasets. To make sure I'm excluding loci under selection for the estimation of the neutral covariance matrix, I will remove loci in the upper Fst tail from each dataset:

1. Identify neutral loci (ie. exclude upper Fst tail) and select 1000 at random

2. Select 1000 loci from vcf file

3. Convert vcf to bayenv2 format

4. Calculate covariance matrix1 + 2

#### CHall

1. Fst<0.6
```
CHall.FstNeutral <- subset(stats.hier.CHall940.perlocus, Fst<0.6)
summary(CHall.FstNeutral)
CHall940.FstNeutral.names <- row.names(CHall.FstNeutral)  ##get the locus names
CHall940.FstNeutral.names <- gsub("X", "", CHall940.FstNeutral.names)
CHall940.FstNeutral.names <- gsub("\\.", ":", CHall940.FstNeutral.names)
CHall940.FstNeutral.names 
CHall940.FstNeutral.names <- as.data.frame(CHall940.FstNeutral.names )
CHall940.1000loci.names <- CHall940.FstNeutral.names[sample(nrow(CHall940.FstNeutral.names), 1000), ]
write.table(CHall940.1000loci.names, "CHall.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHall.940.9608.newnames.vcf --snps CHall.1000loci.names  --recode --recode-INFO-all --out CHall940.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHall940.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/CHall.spidpopfile

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000



#### CHN

1. Fst<0.4
```
CHN.FstNeutral <- subset(stats.hier.CHN229.perlocus, Fst<0.4)
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
vcftools --vcf CHN.229.9608.recode.vcf --snps CHN229.1000loci.names --recode --recode-INFO-all --out CHN229.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHN229.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000


#### CHS

1. Fst<0.5
```
CHS.FstNeutral <- subset(stats.hier.CHS283.perlocus, Fst<0.5)
summary(CHS.FstNeutral)
CHS283.FstNeutral.names <- row.names(CHS.FstNeutral)  ##get the locus names
CHS283.FstNeutral.names <- gsub("X", "", CHS283.FstNeutral.names)
CHS283.FstNeutral.names <- gsub("\\.", ":", CHS283.FstNeutral.names)
CHS283.FstNeutral.names 
CHS283.FstNeutral.names <- as.data.frame(CHS283.FstNeutral.names )
CHS283.1000loci.names <- CHS283.FstNeutral.names[sample(nrow(CHS283.FstNeutral.names), 1000), ]
write.table(CHS283.1000loci.names, "CHS283.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.283.9608.recode.vcf --snps CHS283.1000loci.names --recode --recode-INFO-all --out CHS283.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS283.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000


#### CHS.TI

1. Fst<0.5
```
CHS.TI.FstNeutral <- subset(stats.hier.CHS.TI148.perlocus, Fst<0.5)
summary(CHS.FstNeutral)
CHS.TI148.FstNeutral.names <- row.names(CHS.TI.FstNeutral)  ##get the locus names
CHS.TI148.FstNeutral.names <- gsub("X", "", CHS.TI148.FstNeutral.names)
CHS.TI148.FstNeutral.names <- gsub("\\.", ":", CHS.TI148.FstNeutral.names)
CHS.TI148.FstNeutral.names 
CHS.TI148.FstNeutral.names <- as.data.frame(CHS.TI148.FstNeutral.names )
CHS.TI148.1000loci.names <- CHS.TI148.FstNeutral.names[sample(nrow(CHS.TI148.FstNeutral.names), 1000), ]
write.table(CHS.TI148.1000loci.names, "CHS.TI148.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.TI.148.9608.recode.vcf --snps CHS.TI148.1000loci.names --recode --recode-INFO-all --out CHS.TI148.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS.TI148.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000


#### CHS.VS

1. Fst<0.5
```
CHS.VS.FstNeutral <- subset(stats.hier.CHS.VS135.perlocus, Fst<0.5)
summary(CHS.VS.FstNeutral)
CHS.VS135.FstNeutral.names <- row.names(CHS.TI.FstNeutral)  ##get the locus names
CHS.VS135.FstNeutral.names <- gsub("X", "", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names <- gsub("\\.", ":", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names 
CHS.VS135.FstNeutral.names <- as.data.frame(CHS.VS135.FstNeutral.names )
CHS.VS135.1000loci.names <- CHS.VS135.FstNeutral.names[sample(nrow(CHS.VS135.FstNeutral.names), 1000), ]
write.table(CHS.VS135.1000loci.names, "CHS.VS135.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.VS.135.9608.recode.vcf --snps CHS.VS135.1000loci.names --recode --recode-INFO-all --out CHS.VS135.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS.VS135.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000



#### CZ

1. Fst<0.4
```
CZ404.FstNeutral <- subset(stats.hier.CZ404.perlocus, Fst<0.4)
summary(CZ404.FstNeutral)
CZ404.FstNeutral.names <- row.names(CZ404.FstNeutral)  ##get the locus names
CZ404.FstNeutral.names <- gsub("X", "", CZ404.FstNeutral.names)
CZ404.FstNeutral.names <- gsub("\\.", ":", CZ404.FstNeutral.names)
CZ404.FstNeutral.names 
CZ404.FstNeutral.names <- as.data.frame(CHS.VS135.FstNeutral.names )
CZ404.1000loci.names <- CZ404.FstNeutral.names[sample(nrow(CZ404.FstNeutral.names), 1000), ]
write.table(CZ404.1000loci.names, "CZ404.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CZ.404.9608.recode.vcf --snps CZ404.1000loci.names --recode --recode-INFO-all --out CZ404.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CZ404.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000

####4. Calculate covariance matrix

Copy all input files to the bayenv2 directory: 

/Users/alexjvr/Applications/bayenv2/compiled_on_a_mac

Calculate 2 matrices for each dataset. This takes a long time to run. 
```
#check how far the run has progressed: 
cat CHS.VS.matrix1.out |grep "ITER"

./bayenv2 -i CHall940.1000Neutral.bayenv.txt -p 82 -k 100000 -r 19866 > CHall.matrix1.out
./bayenv2 -i CHall940.1000Neutral.bayenv.txt -p 82 -k 100000 -r 74626 > CHall.matrix2.out

./bayenv2 -i CHN.229.1000Neutral.bayenv.txt -p 19 -k 100000 -r 46843 > CHN229.matrix1.out
./bayenv2 -i CHN.229.1000Neutral.bayenv.txt -p 19 -k 100000 -r 87292 > CHN229.matrix2.out

./bayenv2 -i CZ404.1000Neutral.bayenv.txt -p 38 -k 100000 -r 54839 > CZ404.matrix1.out
./bayenv2 -i CZ404.1000Neutral.bayenv.txt -p 38 -k 100000 -r 27362 > CZ404.matrix2.out

./bayenv2 -i CHS.283.1000Neutral.bayenv.txt -p 25 -k 100000 -r 76290 > CHS283.matrix1.out
./bayenv2 -i CHS.283.1000Neutral.bayenv.txt -p 25 -k 100000 -r 65829 > CHS283.matrix2.out

#./bayenv2 -i CHS.283.1000Neutral.bayenv.txt -p 25 -k 100000 -r 76290 > CHS283.matrix1.out
#./bayenv2 -i CHS.283.1000Neutral.bayenv.txt -p 25 -k 100000 -r 65829 > CHS283.matrix2.out

./bayenv2 -i CHS.VS135.1000Neutral.bayenv.txt -p 10 -k 100000 -r 93863 > CHS.VS.matrix1.out
./bayenv2 -i CHS.VS135.1000Neutral.bayenv.txt -p 10 -k 100000 -r 74952 > CHS.VS.matrix2.out
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

####STILL NEED TO RUN
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
#####STILL NEEDS TO BE DONE
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
#70
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
#26
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
#10
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
#69
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
#29
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
#15
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
#34
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
#68
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
#25
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
#12
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
#73
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
#31
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
#10
```


According to Blair et al. 2014, there is a lot of inter-run variation in BayEnv2 results. They suggest that the best approach (if the covariance matrix converges) is to combine results from several different runs.

Im following the methods in Božičević 2016 (Mol Ecol) to get the most robust results:

" To test the convergence of Bayenv2, we used several independent Markov chain Monte Carlo runs with a maximum chain length of 10 000 iterations. We observed convergence after about 5000 iterations (Fig. S1, Supporting information). However, these chains might converge to different solutions. To be most stringent, we used the median results from 10 independent runs (Blair et al. 2014). We then tested for correlations between each single SNP and six environmental variables: geographical latitude, height above mean sea level and four temperature measures (average daily minimum of the coldest and warmest month, and average daily minimum and maximum throughout the year) (Table 2). The results for environmental variables are given as Bayes factors (BFs). A higher BF gives higher support to the model where the environmental variable has a significant effect on allele frequency distribution over an alternative model with no effect (Coop et al. 2010). Similar as above, we finally report the median BF of ten independent runs of each SNP that has been described to improve the proportion of false positives (Blair et al. 2014; Lotterhos & Whitlock 2014). With BF values, we also used a resampling approach, analogously to the one applied on FST. We randomly sampled sets of SNPs of the same size from the genomic background and then assessed the null distribution for BF of CCRT-, RSS- and SR-associated SNPs. "

Ive run 3 independend BayEnv2 association analyses. now Im combining all the output tables by the median of all the values.


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

####STILL NEEDS TO BE DONE
##CHS.TI
CHS.TI.BF.all <- rbindlist(list(BF.CHS.TI.Run1, BF.CHS.TI.Run2, BF.CHS.TI.Run3))[,lapply(.SD,median), list(snp)]
```


Identify the outliers associated with the different env variables. " Strength of evidence for significant associations was based on the value of the log10 Bayes factor (log10BF), with the following log10BF cut-offs: 0.5–1= substantial evidence; 1–2 = strong evidence; >2 =decisive (Kass & Raftery 1995). The linear model underlying the Bayes factor might not be correct or outliers within our data might misguide the model (bayenv2.0 manual, https://bitbucket.org/tguenther/bayenv2_public/src). To deal with this, bayenv2 also calculates the nonparametric Spearman's rank correlation coefficient, ρ. SNPs with a log10BF >0.5 as well as an absolute value of ρ > 0.3 (where ρ ranges from −1 to 1) were therefore considered as robust candidates demonstrating signatures of selection. " Christmas et al. 2016


```
#calculate the log10 of BF for each environmental variable

BF.all$bio2.log10BF <- log10(BF.all$bio2.BF)
BF.all$bio5.log10BF <- log10(BF.all$bio5.BF)
BF.all$bio13.log10BF <- log10(BF.all$bio13.BF)
BF.all$bio15.log10BF <- log10(BF.all$bio15.BF)
BF.all$bio18.log10BF <- log10(BF.all$bio18.BF)

##Sort and plot

BF.all.sort <- BF.all[order(BF.all$bio2.log10BF),]

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
