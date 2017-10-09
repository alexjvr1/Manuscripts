# CHP4

Adaptation across latitude may be driven by temperature during the growth season, or by the season length. Both of these have been
suggested previously, and in my SE MS I find that adaptation to temperature might be important. 

Here I investigate which of the two variables are the most important, and determine whether adaptation across elevation is also driven
by season length or temperature during the growth season. 

Strategy: 

1. Identify candidate adaptive loci. 
  a)I'm using EAA (LFMM and BayENV) to identify loci associated with season length or temperature during
the breeding season. 
  b)I'm identifying outlier loci using Bayescan, PCAdapt, and XtX from Bayenv2
  
2. Test whether these are important across the gradient compared with neutral loci: 
  a) Gradient forest analysis
  b) GDM
  
3) Identify regions with high genomic turnover -> predict the effects of future climate change. 


#### Datasets

1. CHall.932.7744.recode.vcf

2. CHN.229.5265.recode.vcf

3. CHS.275.6339.recode.vcf

4. CHS.VS.135.5835.recode.vcf

5. CHS.TI.140.5692.recode.vcf

6. CZ.404.7288.recode.vcf

#### EnvVariables

1. mean.temp: mean temperature 60days after spawning

2. season.length: number of days above 6 deg C


## Part1: Identify Candidate adaptive loci

### a) Env. Assoc. Analyses (EAA)

I'm using mean temperature 60 days after laying as temp during development time, and number of days above 6deg C as season length. 

I'm using the same input files for Switzerland as in the Pop Genomics chapter (CHP2)

For SE, I have to filter for MAC as I did for the Swiss populations. 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH
```
vcftools --vcf SE132.FINAL.recode.vcf --mac 13 --recode --recode-INFO-all --out SE132.2027
```
This removes only 54 loci from the previous dataset. 


## EAA

### 1. BayEnv2

#### SE

Input files: 

1. Env file

2. Genotype file

3. Matrix file


##### Env file

I got the new derived environmental data from Josh today (9 Oct 2017) after getting approximate spawning dates from Maria. 

A copy of the per individual data is found here /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/SE.bayenv
 in the same order as the vcf file. 
 
Env data needs to be normalised and transcribed
``` 
##R

library(clusterSim)


file <- read.csv("SE.Derived.EnvData_20171009.csv", header=T) ##make sure the order here is the same as in the vcf file

df <- data.frame(pop=file$pop, mean.temp=file$mean.temp.60.days, season=file$days.6.degrees)

n <- df$pop

df.norm <- data.Normalization(df[,2:ncol(df)], type="n1", normalization="column") ##normalise the data. See --help for data.Normalisation to see other normalisation options

## transpose all standardised data
df.norm.transposed <- as.data.frame(t(df.norm))
colnames(df.norm.transposed) <- n ##just to check the order of the pops in the file
head(df.norm.transposed)
colnames(df.norm.transposed) <- NULL ##and drop the column names
#df.norm.transposed$myfactor <- factor(row.names(df.norm.transposed))

str(df.norm.transposed) # Check the column types
head(df.norm.transposed)

write.table(df.norm.transposed, "SE.ENVIRONFILE.csv", sep=",")##write to file

```


##### Matrix file

We're using the same matrix as before - just check that the pop order is the same in both vcf files. 


##### Genetic data input


converted with pgdspider

##### Run BayEnv

fgcz47:/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/MACfilter/CHP4/

```
split -d -a 10 -l 2 SE132.2027.bayenv.txt snp_batch

for f in $(ls snp_batch*); do ./bayenv2 -i $f -m SE132.MATRIX -e SE.Oct.ENVIRONFILE -p 15 -n 2 -k 100000 -r $RANDOM -t -c -X; done

```



#### CH
Input files for CH have been generated in CHP2. 

Env input is taken from the data that Josh calculated. 

All runs on the GDC server: 


co-variance matrix calculated before: https://github.com/alexjvr1/Manuscripts/blob/master/CHP2_CHall.Bayenv_Oct2017.md

Everything was run on the gdc server. All results are transferred here: 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

First check if all the runs are complete: 

```
## CHN

##this prints the first line of each snp output file 
awk '{print $1}' CHN.Run1/bf_environ.CHN.ENV > names.1      
awk '{print $1}' CHN.Run2/bf_environ.CHN.ENV > names.2   
awk '{print $1}' CHN.Run3/bf_environ.CHN.ENV > names.3   

##this compares the files and prints the differences
diff names.1 names.2
diff names.1 names.3
diff names.2 names.3


## CHS

##this prints the first line of each snp output file 
awk '{print $1}' CHS.Run1/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run1.names   
awk '{print $1}' CHS.Run2/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run2.names
awk '{print $1}' CHS.Run3/bf_environ.ENVIRONFILE.CHS.n5.env > CHS.Run3.names

##this compares the files and prints the differences
diff CHS.Run1.names CHS.Run2.names 
diff CHS.Run1.names CHS.Run3.names 



## CZ

##this prints the first line of each snp output file 
awk '{print $1}' CZ.Run1/bf_environ.CZ.ENV > names.1      
awk '{print $1}' CZ.Run2/bf_environ.CZ.ENV > names.2   
awk '{print $1}' CZ.Run3/bf_environ.CZ.ENV > names.3   

##this compares the files and prints the differences
diff names.1 names.2
diff names.1 names.3
diff names.2 names.3


## CHS.VS

##this prints the first line of each snp output file 
awk '{print $1}' CHS.VS.Run1/bf_environ.CHS.VS.ENV > names.1      
awk '{print $1}' CHS.VS.Run2/bf_environ.CHS.VS.ENV > names.2   
awk '{print $1}' CHS.VS.Run3/bf_environ.CHS.VS.ENV > names.3   

##this compares the files and prints the differences
diff names.1 names.2
diff names.1 names.3
diff names.2 names.3


## CHS.TI

##this prints the first line of each snp output file 
awk '{print $1}' CHS.TI.Run1/bf_environ.CHS.TI.ENV > names.1      
awk '{print $1}' CHS.TI.Run2/bf_environ.CHS.TI.ENV > names.2   
awk '{print $1}' CHS.TI.Run3/bf_environ.CHS.TI.ENV > names.3   

##this compares the files and prints the differences
diff names.1 names.2
diff names.1 names.3
diff names.2 names.3

```


Once this was done, I could concatenate the results as before

I need to calculate convergence of BF, p, and ranking 1000 loci between the 3 runs for each of the 5 environmental variables. (i.e. 45 correlations)


In R
```
##CHN
BF.CHN.Run1 <- read.table("CHN.Run1/bf_environ.CHN.ENV", header=F)
colnames(BF.CHN.Run1) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHN.Run2 <- read.table("CHN.Run2/bf_environ.CHN.ENV", header=F)
colnames(BF.CHN.Run2) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHN.Run3 <- read.table("CHN.Run3/bf_environ.CHN.ENV", header=F)
colnames(BF.CHN.Run3) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")


###CHS
BF.CHS.Run1 <- read.table("CHS.Run1/bf_environ.CHS.ENV", header=F)
colnames(BF.CHS.Run1) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.Run2 <- read.table("CHS.Run2/bf_environ.CHS.ENV", header=F)
colnames(BF.CHS.Run2) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.Run3 <- read.table("CHS.Run3/bf_environ.CHS.ENV", header=F)
colnames(BF.CHS.Run3) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")


###CZ
BF.CZ.Run1 <- read.table("CZ.Run1/bf_environ.CZ.ENV", header=F)
colnames(BF.CZ.Run1) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CZ.Run2 <- read.table("CZ.Run2/bf_environ.CZ.ENV", header=F)
colnames(BF.CZ.Run2) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CZ.Run3 <- read.table("CZ.Run3/bf_environ.CZ.ENV", header=F)
colnames(BF.CZ.Run3) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

###CHS.VS
BF.CHS.VS.Run1 <- read.table("CHS.VS.Run1/bf_environ.CHS.VS.ENV", header=F)
colnames(BF.CHS.VS.Run1) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.VS.Run2 <- read.table("CHS.VS.Run2/bf_environ.CHS.VS.ENV", header=F)
colnames(BF.CHS.VS.Run2) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.VS.Run3 <- read.table("CHS.VS.Run3/bf_environ.CHS.VS.ENV", header=F)
colnames(BF.CHS.VS.Run3) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

###CHS.TI
BF.CHS.TI.Run1 <- read.table("CHS.TI.Run1/bf_environ.CHS.TI.ENV", header=F)
colnames(BF.CHS.TI.Run1) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.TI.Run2 <- read.table("CHS.TI.Run2/bf_environ.CHS.TI.ENV", header=F)
colnames(BF.CHS.TI.Run2) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")

BF.CHS.TI.Run3 <- read.table("CHS.TI.Run2/bf_environ.CHS.TI.ENV", header=F)
colnames(BF.CHS.TI.Run3) <- c("snp", "temp.BF", "temp.rho", "temp.r", "season.days.BF", "season.days.rho", "season.days.r")


```


Compare ranking SNPs 5% at a time



####CHN

1. temp
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
#76
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
#47
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
#28
```

2. season.days

```
CHN.season.days.run1.top100 <- BF.CHN.Run1[order(-BF.CHN.Run1$season.days.rho),]  ##descending order by rho
CHN.season.days.run2.top100 <- BF.CHN.Run2[order(-BF.CHN.Run2$season.days.rho),]  ##descending order by rho
CHN.season.days.run3.top100 <- BF.CHN.Run3[order(-BF.CHN.Run3$season.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHN.season.days
###############
```
CHN.season.days.run1.top100.set1 <- head(CHN.season.days.run1.top100, 100)
CHN.season.days.run2.top100.set1 <- head(CHN.season.days.run2.top100, 100)
CHN.season.days.run3.top100.set1 <- head(CHN.season.days.run3.top100, 100)
CHN.season.days.12 <- Reduce(intersect, list(CHN.season.days.run2.top100.set1$snp, CHN.season.days.run1.top100.set1$snp))
CHN.season.days.13 <- Reduce(intersect, list(CHN.season.days.run3.top100.set1$snp, CHN.season.days.run1.top100.set1$snp))
CHN.season.days.23 <- Reduce(intersect, list(CHN.season.days.run2.top100.set1$snp, CHN.season.days.run3.top100.set1$snp))
CHN.season.days.123 <- Reduce(intersect, list(CHN.season.days.run2.top100.set1$snp, CHN.season.days.run3.top100.set1$snp, CHN.season.days.run1.top100.set1$snp))

summary(CHN.season.days.123)
#81
```

##overlap in 6-10% of outlier loci for CHN.season.days
###############
```
CHN.season.days.run1.top100.set2 <- CHN.season.days.run1.top100[101:200,]
CHN.season.days.run2.top100.set2 <- CHN.season.days.run2.top100[101:200,]
CHN.season.days.run3.top100.set2 <- CHN.season.days.run3.top100[101:200,]
CHN.season.days.12 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))
CHN.season.days.13 <- Reduce(intersect, list(CHN.season.days.run3.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))
CHN.season.days.23 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run3.top100.set2$snp))
CHN.season.days.123 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run3.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))

summary(CHN.season.days.123)
#43
```

##overlap in 11-15% of outlier loci for CHN.season.days
###############
```
CHN.season.days.run1.top100.set2 <- CHN.season.days.run1.top100[201:300,]
CHN.season.days.run2.top100.set2 <- CHN.season.days.run2.top100[201:300,]
CHN.season.days.run3.top100.set2 <- CHN.season.days.run3.top100[201:300,]
CHN.season.days.12 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))
CHN.season.days.13 <- Reduce(intersect, list(CHN.season.days.run3.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))
CHN.season.days.23 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run3.top100.set2$snp))
CHN.season.days.123 <- Reduce(intersect, list(CHN.season.days.run2.top100.set2$snp, CHN.season.days.run3.top100.set2$snp, CHN.season.days.run1.top100.set2$snp))

summary(CHN.season.days.123)
#23
```


####CHS

1. temp
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
#78
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
#41
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
#21
```

2. season.days

```
CHS.season.days.run1.top100 <- BF.CHS.Run1[order(-BF.CHS.Run1$season.days.rho),]  ##descending order by rho
CHS.season.days.run2.top100 <- BF.CHS.Run2[order(-BF.CHS.Run2$season.days.rho),]  ##descending order by rho
CHS.season.days.run3.top100 <- BF.CHS.Run3[order(-BF.CHS.Run3$season.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.season.days
###############
```
CHS.season.days.run1.top100.set1 <- head(CHS.season.days.run1.top100, 100)
CHS.season.days.run2.top100.set1 <- head(CHS.season.days.run2.top100, 100)
CHS.season.days.run3.top100.set1 <- head(CHS.season.days.run3.top100, 100)
CHS.season.days.12 <- Reduce(intersect, list(CHS.season.days.run2.top100.set1$snp, CHS.season.days.run1.top100.set1$snp))
CHS.season.days.13 <- Reduce(intersect, list(CHS.season.days.run3.top100.set1$snp, CHS.season.days.run1.top100.set1$snp))
CHS.season.days.23 <- Reduce(intersect, list(CHS.season.days.run2.top100.set1$snp, CHS.season.days.run3.top100.set1$snp))
CHS.season.days.123 <- Reduce(intersect, list(CHS.season.days.run2.top100.set1$snp, CHS.season.days.run3.top100.set1$snp, CHS.season.days.run1.top100.set1$snp))

summary(CHS.season.days.123)
#78
```

##overlap in 6-10% of outlier loci for CHS.season.days
###############
```
CHS.season.days.run1.top100.set2 <- CHS.season.days.run1.top100[101:200,]
CHS.season.days.run2.top100.set2 <- CHS.season.days.run2.top100[101:200,]
CHS.season.days.run3.top100.set2 <- CHS.season.days.run3.top100[101:200,]
CHS.season.days.12 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))
CHS.season.days.13 <- Reduce(intersect, list(CHS.season.days.run3.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))
CHS.season.days.23 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run3.top100.set2$snp))
CHS.season.days.123 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run3.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))

summary(CHS.season.days.123)
#43
```

##overlap in 11-15% of outlier loci for CHS.season.days
###############
```
CHS.season.days.run1.top100.set2 <- CHS.season.days.run1.top100[201:300,]
CHS.season.days.run2.top100.set2 <- CHS.season.days.run2.top100[201:300,]
CHS.season.days.run3.top100.set2 <- CHS.season.days.run3.top100[201:300,]
CHS.season.days.12 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))
CHS.season.days.13 <- Reduce(intersect, list(CHS.season.days.run3.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))
CHS.season.days.23 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run3.top100.set2$snp))
CHS.season.days.123 <- Reduce(intersect, list(CHS.season.days.run2.top100.set2$snp, CHS.season.days.run3.top100.set2$snp, CHS.season.days.run1.top100.set2$snp))

summary(CHS.season.days.123)
#21
```



####CZ

1. temp
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
#68
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
#23
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
#9
```

2. season.days

```
CZ.season.days.run1.top100 <- BF.CZ.Run1[order(-BF.CZ.Run1$season.days.rho),]  ##descending order by rho
CZ.season.days.run2.top100 <- BF.CZ.Run2[order(-BF.CZ.Run2$season.days.rho),]  ##descending order by rho
CZ.season.days.run3.top100 <- BF.CZ.Run3[order(-BF.CZ.Run3$season.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CZ.season.days
###############
```
CZ.season.days.run1.top100.set1 <- head(CZ.season.days.run1.top100, 100)
CZ.season.days.run2.top100.set1 <- head(CZ.season.days.run2.top100, 100)
CZ.season.days.run3.top100.set1 <- head(CZ.season.days.run3.top100, 100)
CZ.season.days.12 <- Reduce(intersect, list(CZ.season.days.run2.top100.set1$snp, CZ.season.days.run1.top100.set1$snp))
CZ.season.days.13 <- Reduce(intersect, list(CZ.season.days.run3.top100.set1$snp, CZ.season.days.run1.top100.set1$snp))
CZ.season.days.23 <- Reduce(intersect, list(CZ.season.days.run2.top100.set1$snp, CZ.season.days.run3.top100.set1$snp))
CZ.season.days.123 <- Reduce(intersect, list(CZ.season.days.run2.top100.set1$snp, CZ.season.days.run3.top100.set1$snp, CZ.season.days.run1.top100.set1$snp))

summary(CZ.season.days.123)
#69
```

##overlap in 6-10% of outlier loci for CZ.season.days
###############
```
CZ.season.days.run1.top100.set2 <- CZ.season.days.run1.top100[101:200,]
CZ.season.days.run2.top100.set2 <- CZ.season.days.run2.top100[101:200,]
CZ.season.days.run3.top100.set2 <- CZ.season.days.run3.top100[101:200,]
CZ.season.days.12 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))
CZ.season.days.13 <- Reduce(intersect, list(CZ.season.days.run3.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))
CZ.season.days.23 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run3.top100.set2$snp))
CZ.season.days.123 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run3.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))

summary(CZ.season.days.123)
#33
```

##overlap in 11-15% of outlier loci for CZ.season.days
###############
```
CZ.season.days.run1.top100.set2 <- CZ.season.days.run1.top100[201:300,]
CZ.season.days.run2.top100.set2 <- CZ.season.days.run2.top100[201:300,]
CZ.season.days.run3.top100.set2 <- CZ.season.days.run3.top100[201:300,]
CZ.season.days.12 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))
CZ.season.days.13 <- Reduce(intersect, list(CZ.season.days.run3.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))
CZ.season.days.23 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run3.top100.set2$snp))
CZ.season.days.123 <- Reduce(intersect, list(CZ.season.days.run2.top100.set2$snp, CZ.season.days.run3.top100.set2$snp, CZ.season.days.run1.top100.set2$snp))

summary(CZ.season.days.123)
#19
```



####CHS.VS

1. temp
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
#80
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
#46
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
#28
```

2. season.days

```
CHS.VS.season.days.run1.top100 <- BF.CHS.VS.Run1[order(-BF.CHS.VS.Run1$season.days.rho),]  ##descending order by rho
CHS.VS.season.days.run2.top100 <- BF.CHS.VS.Run2[order(-BF.CHS.VS.Run2$season.days.rho),]  ##descending order by rho
CHS.VS.season.days.run3.top100 <- BF.CHS.VS.Run3[order(-BF.CHS.VS.Run3$season.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.VS.season.days
###############
```
CHS.VS.season.days.run1.top100.set1 <- head(CHS.VS.season.days.run1.top100, 100)
CHS.VS.season.days.run2.top100.set1 <- head(CHS.VS.season.days.run2.top100, 100)
CHS.VS.season.days.run3.top100.set1 <- head(CHS.VS.season.days.run3.top100, 100)
CHS.VS.season.days.12 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set1$snp, CHS.VS.season.days.run1.top100.set1$snp))
CHS.VS.season.days.13 <- Reduce(intersect, list(CHS.VS.season.days.run3.top100.set1$snp, CHS.VS.season.days.run1.top100.set1$snp))
CHS.VS.season.days.23 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set1$snp, CHS.VS.season.days.run3.top100.set1$snp))
CHS.VS.season.days.123 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set1$snp, CHS.VS.season.days.run3.top100.set1$snp, CHS.VS.season.days.run1.top100.set1$snp))

summary(CHS.VS.season.days.123)
#81
```

##overlap in 6-10% of outlier loci for CHS.VS.season.days
###############
```
CHS.VS.season.days.run1.top100.set2 <- CHS.VS.season.days.run1.top100[101:200,]
CHS.VS.season.days.run2.top100.set2 <- CHS.VS.season.days.run2.top100[101:200,]
CHS.VS.season.days.run3.top100.set2 <- CHS.VS.season.days.run3.top100[101:200,]
CHS.VS.season.days.12 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))
CHS.VS.season.days.13 <- Reduce(intersect, list(CHS.VS.season.days.run3.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))
CHS.VS.season.days.23 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run3.top100.set2$snp))
CHS.VS.season.days.123 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run3.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))

summary(CHS.VS.season.days.123)
#41
```

##overlap in 11-15% of outlier loci for CHS.VS.season.days
###############
```
CHS.VS.season.days.run1.top100.set2 <- CHS.VS.season.days.run1.top100[201:300,]
CHS.VS.season.days.run2.top100.set2 <- CHS.VS.season.days.run2.top100[201:300,]
CHS.VS.season.days.run3.top100.set2 <- CHS.VS.season.days.run3.top100[201:300,]
CHS.VS.season.days.12 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))
CHS.VS.season.days.13 <- Reduce(intersect, list(CHS.VS.season.days.run3.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))
CHS.VS.season.days.23 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run3.top100.set2$snp))
CHS.VS.season.days.123 <- Reduce(intersect, list(CHS.VS.season.days.run2.top100.set2$snp, CHS.VS.season.days.run3.top100.set2$snp, CHS.VS.season.days.run1.top100.set2$snp))

summary(CHS.VS.season.days.123)
#9
```



####CHS.TI

1. temp
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
#83
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
#56
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
#47
```

2. season.days

```
CHS.TI.season.days.run1.top100 <- BF.CHS.TI.Run1[order(-BF.CHS.TI.Run1$season.days.rho),]  ##descending order by rho
CHS.TI.season.days.run2.top100 <- BF.CHS.TI.Run2[order(-BF.CHS.TI.Run2$season.days.rho),]  ##descending order by rho
CHS.TI.season.days.run3.top100 <- BF.CHS.TI.Run3[order(-BF.CHS.TI.Run3$season.days.rho),]  ##descending order by rho
```

##overlap in top 5% of outliers loci for CHS.TI.season.days
###############
```
CHS.TI.season.days.run1.top100.set1 <- head(CHS.TI.season.days.run1.top100, 100)
CHS.TI.season.days.run2.top100.set1 <- head(CHS.TI.season.days.run2.top100, 100)
CHS.TI.season.days.run3.top100.set1 <- head(CHS.TI.season.days.run3.top100, 100)
CHS.TI.season.days.12 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set1$snp, CHS.TI.season.days.run1.top100.set1$snp))
CHS.TI.season.days.13 <- Reduce(intersect, list(CHS.TI.season.days.run3.top100.set1$snp, CHS.TI.season.days.run1.top100.set1$snp))
CHS.TI.season.days.23 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set1$snp, CHS.TI.season.days.run3.top100.set1$snp))
CHS.TI.season.days.123 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set1$snp, CHS.TI.season.days.run3.top100.set1$snp, CHS.TI.season.days.run1.top100.set1$snp))

summary(CHS.TI.season.days.123)
#84
```

##overlap in 6-10% of outlier loci for CHS.TI.season.days
###############
```
CHS.TI.season.days.run1.top100.set2 <- CHS.TI.season.days.run1.top100[101:200,]
CHS.TI.season.days.run2.top100.set2 <- CHS.TI.season.days.run2.top100[101:200,]
CHS.TI.season.days.run3.top100.set2 <- CHS.TI.season.days.run3.top100[101:200,]
CHS.TI.season.days.12 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))
CHS.TI.season.days.13 <- Reduce(intersect, list(CHS.TI.season.days.run3.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))
CHS.TI.season.days.23 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run3.top100.set2$snp))
CHS.TI.season.days.123 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run3.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))

summary(CHS.TI.season.days.123)
#61
```

##overlap in 11-15% of outlier loci for CHS.TI.season.days
###############
```
CHS.TI.season.days.run1.top100.set2 <- CHS.TI.season.days.run1.top100[201:300,]
CHS.TI.season.days.run2.top100.set2 <- CHS.TI.season.days.run2.top100[201:300,]
CHS.TI.season.days.run3.top100.set2 <- CHS.TI.season.days.run3.top100[201:300,]
CHS.TI.season.days.12 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))
CHS.TI.season.days.13 <- Reduce(intersect, list(CHS.TI.season.days.run3.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))
CHS.TI.season.days.23 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run3.top100.set2$snp))
CHS.TI.season.days.123 <- Reduce(intersect, list(CHS.TI.season.days.run2.top100.set2$snp, CHS.TI.season.days.run3.top100.set2$snp, CHS.TI.season.days.run1.top100.set2$snp))

summary(CHS.TI.season.days.123)
#47
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



####CHN
```
#calculate the log10 of BF for each environmental variable

CHN.BF.all$temp.log10BF <- log10(CHN.BF.all$temp.BF)
CHN.BF.all$season.days.log10BF <- log10(CHN.BF.all$season.days.BF)


##Sort and plot

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$temp.log10BF),]

pdf("CHN.CHP4.BF.plot.pdf")
par(mfrow=c(1,2))
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$temp.log10BF),]
plot(CHN.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$season.days.log10BF),]
plot(CHN.BF.all.sort$season.days.log10BF, main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```


##plot BF vs absolute p value
```
pdf("CHN.CHP4.BFvsp.pdf")
par(mfrow=c(1,2))
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$temp.log10BF),]
plot(CHN.BF.all.sort$temp.log10BF~(abs(CHN.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHN.BF.all.sort <- CHN.BF.all[order(CHN.BF.all$season.days.log10BF),]
plot(CHN.BF.all.sort$season.days.log10BF~(abs(CHN.BF.all.sort$season.days.rho)), main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHN.temp.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$temp.log10BF>0.5 & (abs(CHN.BF.all$temp.rho))>0.3),]
CHN.season.days.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$season.days.log10BF>0.5 & (abs(CHN.BF.all$season.days.rho))>0.3),]
```





####CHS

```
#calculate the log10 of BF for each environmental variable

CHS.BF.all$temp.log10BF <- log10(CHS.BF.all$temp.BF)
CHS.BF.all$season.days.log10BF <- log10(CHS.BF.all$season.days.BF)

##Sort and plot

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$temp.log10BF),]

pdf("CHS.CHP4.BF.plot.pdf")
par(mfrow=c(1,2))
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$temp.log10BF),]
plot(CHS.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$season.days.log10BF),]
plot(CHS.BF.all.sort$season.days.log10BF, main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```


##plot BF vs absolute p value
```
pdf("CHS.CHP4.BFvsp.pdf")
par(mfrow=c(1,2))
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$temp.log10BF),]
plot(CHS.BF.all.sort$temp.log10BF~(abs(CHS.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.BF.all.sort <- CHS.BF.all[order(CHS.BF.all$season.days.log10BF),]
plot(CHS.BF.all.sort$season.days.log10BF~(abs(CHS.BF.all.sort$season.days.rho)), main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHS.temp.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$temp.log10BF>0.5 & (abs(CHS.BF.all$temp.rho))>0.3),]
CHS.season.days.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$season.days.log10BF>0.5 & (abs(CHS.BF.all$season.days.rho))>0.3),]
```





####CZ

```
#calculate the log10 of BF for each environmental variable

CZ.BF.all$temp.log10BF <- log10(CZ.BF.all$temp.BF)
CZ.BF.all$season.days.log10BF <- log10(CZ.BF.all$season.days.BF)

##Sort and plot

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$temp.log10BF),]

pdf("CZ.CHP4.BF.plot.pdf")
par(mfrow=c(1,2))
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$temp.log10BF),]
plot(CZ.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$season.days.log10BF),]
plot(CZ.BF.all.sort$season.days.log10BF, main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CZ.CHP4.BFvsp.pdf")
par(mfrow=c(1,2))
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$temp.log10BF),]
plot(CZ.BF.all.sort$temp.log10BF~(abs(CZ.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CZ.BF.all.sort <- CZ.BF.all[order(CZ.BF.all$season.days.log10BF),]
plot(CZ.BF.all.sort$season.days.log10BF~(abs(CZ.BF.all.sort$season.days.rho)), main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CZ.temp.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$temp.log10BF>0.5 & (abs(CZ.BF.all$temp.rho))>0.3),]
CZ.season.days.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$season.days.log10BF>0.5 & (abs(CZ.BF.all$season.days.rho))>0.3),]
```


####CHS.VS
```
#calculate the log10 of BF for each environmental variable

CHS.VS.BF.all$temp.log10BF <- log10(CHS.VS.BF.all$temp.BF)
CHS.VS.BF.all$season.days.log10BF <- log10(CHS.VS.BF.all$season.days.BF)

##Sort and plot

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$temp.log10BF),]

pdf("CHS.VS.CHP4.BF.plot.pdf")
par(mfrow=c(1,2))
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$temp.log10BF),]
plot(CHS.VS.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$season.days.log10BF),]
plot(CHS.VS.BF.all.sort$season.days.log10BF, main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHS.VS.CHP4.BFvsp.pdf")
par(mfrow=c(1,2))
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$temp.log10BF),]
plot(CHS.VS.BF.all.sort$temp.log10BF~(abs(CHS.VS.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.VS.BF.all.sort <- CHS.VS.BF.all[order(CHS.VS.BF.all$season.days.log10BF),]
plot(CHS.VS.BF.all.sort$season.days.log10BF~(abs(CHS.VS.BF.all.sort$season.days.rho)), main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)

dev.off()
```

Identify the candidates for each env variable
```
CHS.VS.temp.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$temp.log10BF>0.5 & (abs(CHS.VS.BF.all$temp.rho))>0.3),]
CHS.VS.season.days.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$season.days.log10BF>0.5 & (abs(CHS.VS.BF.all$season.days.rho))>0.3),]
```




####CHS.TI
```
#calculate the log10 of BF for each environmental variable

CHS.TI.BF.all$temp.log10BF <- log10(CHS.TI.BF.all$temp.BF)
CHS.TI.BF.all$season.days.log10BF <- log10(CHS.TI.BF.all$season.days.BF)

##Sort and plot

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$temp.log10BF),]

pdf("CHS.TI.CHP4.BF.plot.pdf")
par(mfrow=c(1,2))
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$temp.log10BF),]
plot(CHS.TI.BF.all.sort$temp.log10BF, main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$season.days.log10BF),]
plot(CHS.TI.BF.all.sort$season.days.log10BF, main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)

dev.off()
```

##plot BF vs absolute p value
```
pdf("CHS.TI.CHP4.BFvsp.pdf")
par(mfrow=c(1,2))
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$temp.log10BF),]
plot(CHS.TI.BF.all.sort$temp.log10BF~(abs(CHS.TI.BF.all.sort$temp.rho)), main="temp")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
CHS.TI.BF.all.sort <- CHS.TI.BF.all[order(CHS.TI.BF.all$season.days.log10BF),]
plot(CHS.TI.BF.all.sort$season.days.log10BF~(abs(CHS.TI.BF.all.sort$season.days.rho)), main="season.days")
abline(h=0.5, col=4, lty=2)
abline(h=1.0, col=3, lty=2)
abline(v=0.3, col=3, lty=2)
dev.off()
```

Identify the candidates for each env variable
```
CHS.TI.temp.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$temp.log10BF>0.5 & (abs(CHS.TI.BF.all$temp.rho))>0.3),]
CHS.TI.season.days.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$season.days.log10BF>0.5 & (abs(CHS.TI.BF.all$season.days.rho))>0.3),]
```




###Determine the names of all the loci. 


Check that this is the right number of loci by comparing this to the number of snp_batch files found when the .txt input file is split 
(see above). Remember that these are numbered from 0. 

Next, convert to plink. 

Once this is done, the names can be read into R and indexed. 

```
#linux

vcftools --vcf CHN.229.5265.recode.vcf --plink --out CHN.229.5265.plink

vcftools --vcf CHS.275.6339.recode.vcf --plink --out CHS.275.6339.plink

vcftools --vcf CZ.404.7288.recode.vcf --plink --out CZ.404.7288.plink

vcftools --vcf CHS.VS.135.5835.recode.vcf --plink --out CHS.VS.135.5835.plink

vcftools --vcf CHS.TI.140.5692.recode.vcf --plink --out CHS.TI.140.5692.plink

```


###### CHN
```
CHN.locus.names <- read.table("CHN.229.5265.plink.map", header=F) #import the locus names into R
CHN.locus.names$ID <- seq.int(nrow(CHN.locus.names)) ##index the CHN.locus.names file so that all the loci are numbered in order of appearance

CHN.BF.all$ID <- seq.int(nrow(CHN.BF.all))  ##do the same with the CHN.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHN.BF.all file

##Find all the candidate loci
CHN.temp.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$temp.log10BF>0.5 & (abs(CHN.BF.all$temp.rho))>0.3),]
CHN.season.days.bayenv.candidates <- CHN.BF.all[which(CHN.BF.all$season.days.log10BF>0.5 & (abs(CHN.BF.all$season.days.rho))>0.3),]

CHN.temp.bayenv.candidates$ID <- as.character(CHN.temp.bayenv.candidates$ID)
CHN.temp.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.temp.bayenv.candidates.names <- paste("X", CHN.temp.bayenv.candidates.names$SNP, sep=".")

CHN.season.days.bayenv.candidates$ID <- as.character(CHN.season.days.bayenv.candidates$ID)
CHN.season.days.bayenv.candidates.names <- CHN.locus.names[CHN.locus.names$ID %in% CHN.season.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHN.season.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHN.season.days.bayenv.candidates.names <- paste("X", CHN.season.days.bayenv.candidates.names$SNP, sep=".")

```

Find names and write to file
```
CHN.temp.names <- as.data.frame(CHN.temp.bayenv.candidates.names)
colnames(CHN.temp.names) <- "names"

CHN.temp.outlier.names <- sub(":", ".", CHN.temp.names$names)
CHN.temp.outlier.names <- as.data.frame(CHN.temp.outlier.names)

write.table(CHN.temp.outlier.names$CHN.temp.outlier.names, "CHN.CHP4.temp.outlier.names", col.names=F, row.names=F, quote=F)

CHN.season.days.names <- as.data.frame(CHN.season.days.bayenv.candidates.names)
colnames(CHN.season.days.names) <- "names"

CHN.season.outlier.names <- sub(":", ".", CHN.season.days.names$names)
CHN.season.outlier.names <- as.data.frame(CHN.season.outlier.names)

write.table(CHN.season.outlier.names$CHN.season.outlier.names, "CHN.CHP4.season.outlier.names", col.names=F, row.names=F, quote=F)


CHN.alloutliers.names <- rbind(CHN.temp.names, CHN.season.days.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHN.alloutliers.names <- lapply(CHN.alloutliers.names, unique)  #select only the unique rows. 

CHN.alloutliers.names <- sub(":", ".", CHN.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHN.alloutliers.names <- as.data.frame(CHN.alloutliers.names)
write.table(CHN.alloutliers.names$CHN.alloutliers.names, "CHN.CHP4.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
CHN.BayEnv.outliers <- read.table("CHN.BayEnv.alloutliers")
colnames(CHN.BayEnv.outliers) <- ("loci")
CHN.BayEnv.outliers <- as.character(CHN.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHN.XtX.run1 <- read.table("CHN.Run1/XtX_out.CHN.ENV", header=F)
CHN.XtX.run2 <- read.table("CHN.Run2/XtX_out.CHN.ENV", header=F)
CHN.XtX.run3 <- read.table("CHN.Run3/XtX_out.CHN.ENV", header=F)

library(data.table)

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
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```


###### CHS
```
CHS.locus.names <- read.table("CHS.275.6339.plink.map", header=F) #import the locus names into R
CHS.locus.names$ID <- seq.int(nrow(CHS.locus.names)) ##index the CHS.locus.names file so that all the loci are numbered in order of appearance

CHS.BF.all$ID <- seq.int(nrow(CHS.BF.all))  ##do the same with the CHS.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.BF.all file

##Find all the candidate loci
CHS.temp.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$temp.log10BF>0.5 & (abs(CHS.BF.all$temp.rho))>0.3),]
CHS.season.days.bayenv.candidates <- CHS.BF.all[which(CHS.BF.all$season.days.log10BF>0.5 & (abs(CHS.BF.all$season.days.rho))>0.3),]

CHS.temp.bayenv.candidates$ID <- as.character(CHS.temp.bayenv.candidates$ID)
CHS.temp.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.temp.bayenv.candidates.names <- paste("X", CHS.temp.bayenv.candidates.names$SNP, sep=".")

CHS.season.days.bayenv.candidates$ID <- as.character(CHS.season.days.bayenv.candidates$ID)
CHS.season.days.bayenv.candidates.names <- CHS.locus.names[CHS.locus.names$ID %in% CHS.season.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.season.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.season.days.bayenv.candidates.names <- paste("X", CHS.season.days.bayenv.candidates.names$SNP, sep=".")

```

Find names and write to file
```
CHS.temp.names <- as.data.frame(CHS.temp.bayenv.candidates.names)
colnames(CHS.temp.names) <- "names"

CHS.temp.outlier.names <- sub(":", ".", CHS.temp.names$names)
CHS.temp.outlier.names <- as.data.frame(CHS.temp.outlier.names)

write.table(CHS.temp.outlier.names$CHS.temp.outlier.names, "CHS.CHP4.temp.outlier.names", col.names=F, row.names=F, quote=F)

CHS.season.days.names <- as.data.frame(CHS.season.days.bayenv.candidates.names)
colnames(CHS.season.days.names) <- "names"

CHS.season.outlier.names <- sub(":", ".", CHS.season.days.names$names)
CHS.season.outlier.names <- as.data.frame(CHS.season.outlier.names)

write.table(CHS.season.outlier.names$CHS.season.outlier.names, "CHS.CHP4.season.outlier.names", col.names=F, row.names=F, quote=F)



CHS.alloutliers.names <- rbind(CHS.temp.names, CHS.season.days.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.alloutliers.names <- lapply(CHS.alloutliers.names, unique)  #select only the unique rows. 

CHS.alloutliers.names <- sub(":", ".", CHS.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.alloutliers.names <- as.data.frame(CHS.alloutliers.names)
write.table(CHS.alloutliers.names$CHS.alloutliers.names, "CHS.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
CHS.BayEnv.outliers <- read.table("CHS.BayEnv.alloutliers")
colnames(CHS.BayEnv.outliers) <- ("loci")
CHS.BayEnv.outliers <- as.character(CHS.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
library(data.table)

CHS.XtX.run1 <- read.table("CHS.Run1/XtX_out.CHS.ENV", header=F)
CHS.XtX.run2 <- read.table("CHS.Run2/XtX_out.CHS.ENV", header=F)
CHS.XtX.run3 <- read.table("CHS.Run3/XtX_out.CHS.ENV", header=F)


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
##copy the list over to /Users/alexjvr/2016tempAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```


####### CZ
```
CZ.locus.names <- read.table("CZ.404.7288.plink.map", header=F) #import the locus names into R
CZ.locus.names$ID <- seq.int(nrow(CZ.locus.names)) ##index the CZ.locus.names file so that all the loci are numbered in order of appearance

CZ.BF.all$ID <- seq.int(nrow(CZ.BF.all))  ##do the same with the CZ.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CZ.BF.all file

##Find all the candidate loci
CZ.temp.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$temp.log10BF>0.5 & (abs(CZ.BF.all$temp.rho))>0.3),]
CZ.season.days.bayenv.candidates <- CZ.BF.all[which(CZ.BF.all$season.days.log10BF>0.5 & (abs(CZ.BF.all$season.days.rho))>0.3),]

CZ.temp.bayenv.candidates$ID <- as.character(CZ.temp.bayenv.candidates$ID)
CZ.temp.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.temp.bayenv.candidates.names <- paste("X", CZ.temp.bayenv.candidates.names$SNP, sep=".")

CZ.season.days.bayenv.candidates$ID <- as.character(CZ.season.days.bayenv.candidates$ID)
CZ.season.days.bayenv.candidates.names <- CZ.locus.names[CZ.locus.names$ID %in% CZ.season.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CZ.season.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CZ.season.days.bayenv.candidates.names <- paste("X", CZ.season.days.bayenv.candidates.names$SNP, sep=".")
```

find names and write to file
```
CZ.temp.names <- as.data.frame(CZ.temp.bayenv.candidates.names)
colnames(CZ.temp.names) <- "names"

CZ.temp.outlier.names <- sub(":", ".", CZ.temp.names$names)
CZ.temp.outlier.names <- as.data.frame(CZ.temp.outlier.names)

write.table(CZ.temp.outlier.names$CZ.temp.outlier.names, "CZ.CHP4.temp.outlier.names", col.names=F, row.names=F, quote=F)

CZ.season.days.names <- as.data.frame(CZ.season.days.bayenv.candidates.names)
colnames(CZ.season.days.names) <- "names"

CZ.season.outlier.names <- sub(":", ".", CZ.season.days.names$names)
CZ.season.outlier.names <- as.data.frame(CZ.season.outlier.names)

write.table(CZ.season.outlier.names$CZ.season.outlier.names, "CZ.CHP4.season.outlier.names", col.names=F, row.names=F, quote=F)




CZ.alloutliers.names <- rbind(CZ.temp.names, CZ.season.days.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CZ.alloutliers.names <- lapply(CZ.alloutliers.names, unique)  #select only the unique rows. 

CZ.alloutliers.names <- sub(":", ".", CZ.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CZ.alloutliers.names <- as.data.frame(CZ.alloutliers.names)
write.table(CZ.alloutliers.names$CZ.alloutliers.names, "CZ.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016tempAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
CZ.BayEnv.outliers <- read.table("CZ.BayEnv.alloutliers")
colnames(CZ.BayEnv.outliers) <- ("loci")
CZ.BayEnv.outliers <- as.character(CZ.BayEnv.outliers$loci)
```


Write top 100 XtX loci to file

```
CZ.XtX.run1 <- read.table("CZ.Run1/XtX_out.CZ.ENV", header=F)
CZ.XtX.run2 <- read.table("CZ.Run2/XtX_out.CZ.ENV", header=F)
CZ.XtX.run3 <- read.table("CZ.Run3/XtX_out.CZ.ENV", header=F)


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
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```




####### CHS.VS
```
CHS.VS.locus.names <- read.table("CHS.VS.135.5835.plink.map", header=F) #import the locus names into R
CHS.VS.locus.names$ID <- seq.int(nrow(CHS.VS.locus.names)) ##index the CHS.VS.locus.names file so that all the loci are numbered in order of appearance

CHS.VS.BF.all$ID <- seq.int(nrow(CHS.VS.BF.all))  ##do the same with the CHS.VS.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.VS.BF.all file

##Find all the candidate loci
CHS.VS.temp.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$temp.log10BF>0.5 & (abs(CHS.VS.BF.all$temp.rho))>0.3),]
CHS.VS.season.days.bayenv.candidates <- CHS.VS.BF.all[which(CHS.VS.BF.all$season.days.log10BF>0.5 & (abs(CHS.VS.BF.all$season.days.rho))>0.3),]

CHS.VS.temp.bayenv.candidates$ID <- as.character(CHS.VS.temp.bayenv.candidates$ID)
CHS.VS.temp.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.temp.bayenv.candidates.names <- paste("X", CHS.VS.temp.bayenv.candidates.names$SNP, sep=".")

CHS.VS.season.days.bayenv.candidates$ID <- as.character(CHS.VS.season.days.bayenv.candidates$ID)
CHS.VS.season.days.bayenv.candidates.names <- CHS.VS.locus.names[CHS.VS.locus.names$ID %in% CHS.VS.season.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.VS.season.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.VS.season.days.bayenv.candidates.names <- paste("X", CHS.VS.season.days.bayenv.candidates.names$SNP, sep=".")

```


find names and write to file

```
CHS.VS.temp.names <- as.data.frame(CHS.VS.temp.bayenv.candidates.names)
colnames(CHS.VS.temp.names) <- "names"

CHS.VS.temp.outlier.names <- sub(":", ".", CHS.VS.temp.names$names)
CHS.VS.temp.outlier.names <- as.data.frame(CHS.VS.temp.outlier.names)

write.table(CHS.VS.temp.outlier.names$CHS.VS.temp.outlier.names, "CHS.VS.CHP4.temp.outlier.names", col.names=F, row.names=F, quote=F)

CHS.VS.season.days.names <- as.data.frame(CHS.VS.season.days.bayenv.candidates.names)
colnames(CHS.VS.season.days.names) <- "names"

CHS.VS.season.outlier.names <- sub(":", ".", CHS.VS.season.days.names$names)
CHS.VS.season.outlier.names <- as.data.frame(CHS.VS.season.outlier.names)

write.table(CHS.VS.season.outlier.names$CHS.VS.season.outlier.names, "CHS.VS.CHP4.season.outlier.names", col.names=F, row.names=F, quote=F)




CHS.VS.alloutliers.names <- rbind(CHS.VS.temp.names, CHS.VS.season.days.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.VS.alloutliers.names <- lapply(CHS.VS.alloutliers.names, unique)  #select only the unique rows. 

CHS.VS.alloutliers.names <- sub(":", ".", CHS.VS.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.VS.alloutliers.names <- as.data.frame(CHS.VS.alloutliers.names)
write.table(CHS.VS.alloutliers.names$CHS.VS.alloutliers.names, "CHS.VS.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
CHS.VS.BayEnv.outliers <- read.table("CHS.VS.BayEnv.alloutliers")
colnames(CHS.VS.BayEnv.outliers) <- ("loci")
CHS.VS.BayEnv.outliers <- as.character(CHS.VS.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHS.VS.XtX.run1 <- read.table("CHS.VS.Run1/XtX_out.CHS.VS.ENV", header=F)
CHS.VS.XtX.run2 <- read.table("CHS.VS.Run2/XtX_out.CHS.VS.ENV", header=F)
CHS.VS.XtX.run3 <- read.table("CHS.VS.Run3/XtX_out.CHS.VS.ENV", header=F)


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
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")


```



####### CHS.TI
```
CHS.TI.locus.names <- read.table("CHS.TI.140.5692.plink.map", header=F) #import the locus names into R
CHS.TI.locus.names$ID <- seq.int(nrow(CHS.TI.locus.names)) ##index the CHS.TI.locus.names file so that all the loci are numbered in order of appearance

CHS.TI.BF.all$ID <- seq.int(nrow(CHS.TI.BF.all))  ##do the same with the CHS.TI.BF.all file. Make sure this is the original output from BayEnv, and not a sorted file. 

##find the candidate loci using the indexed .CHS.TI.BF.all file

##Find all the candidate loci
CHS.TI.temp.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$temp.log10BF>0.5 & (abs(CHS.TI.BF.all$temp.rho))>0.3),]
CHS.TI.season.days.bayenv.candidates <- CHS.TI.BF.all[which(CHS.TI.BF.all$season.days.log10BF>0.5 & (abs(CHS.TI.BF.all$season.days.rho))>0.3),]

CHS.TI.temp.bayenv.candidates$ID <- as.character(CHS.TI.temp.bayenv.candidates$ID)
CHS.TI.temp.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.temp.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.temp.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.temp.bayenv.candidates.names <- paste("X", CHS.TI.temp.bayenv.candidates.names$SNP, sep=".")

CHS.TI.season.days.bayenv.candidates$ID <- as.character(CHS.TI.season.days.bayenv.candidates$ID)
CHS.TI.season.days.bayenv.candidates.names <- CHS.TI.locus.names[CHS.TI.locus.names$ID %in% CHS.TI.season.days.bayenv.candidates$ID,]  #Find the actual locus names
colnames(CHS.TI.season.days.bayenv.candidates.names) <- c("V1", "SNP", "V3", "V4", "ID")
CHS.TI.season.days.bayenv.candidates.names <- paste("X", CHS.TI.season.days.bayenv.candidates.names$SNP, sep=".")

```

Find names and write to file
```
CHS.TI.temp.names <- as.data.frame(CHS.TI.temp.bayenv.candidates.names)
colnames(CHS.TI.temp.names) <- "names"

CHS.TI.temp.outlier.names <- sub(":", ".", CHS.TI.temp.names$names)
CHS.TI.temp.outlier.names <- as.data.frame(CHS.TI.temp.outlier.names)

write.table(CHS.TI.temp.outlier.names$CHS.TI.temp.outlier.names, "CHS.TI.CHP4.temp.outlier.names", col.names=F, row.names=F, quote=F)

CHS.TI.season.days.names <- as.data.frame(CHS.TI.season.days.bayenv.candidates.names)
colnames(CHS.TI.season.days.names) <- "names"

CHS.TI.season.outlier.names <- sub(":", ".", CHS.TI.season.days.names$names)
CHS.TI.season.outlier.names <- as.data.frame(CHS.TI.season.outlier.names)

write.table(CHS.TI.season.outlier.names$CHS.TI.season.outlier.names, "CHS.TI.CHP4.season.outlier.names", col.names=F, row.names=F, quote=F)




CHS.TI.alloutliers.names <- rbind(CHS.TI.temp.names, CHS.TI.season.days.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

CHS.TI.alloutliers.names <- lapply(CHS.TI.alloutliers.names, unique)  #select only the unique rows. 

CHS.TI.alloutliers.names <- sub(":", ".", CHS.TI.alloutliers.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

CHS.TI.alloutliers.names <- as.data.frame(CHS.TI.alloutliers.names)
write.table(CHS.TI.alloutliers.names$CHS.TI.alloutliers.names, "CHS.TI.BayEnv.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/BayEnv

#Read into R
CHS.TI.BayEnv.outliers <- read.table("CHS.TI.BayEnv.alloutliers")
colnames(CHS.TI.BayEnv.outliers) <- ("loci")
CHS.TI.BayEnv.outliers <- as.character(CHS.TI.BayEnv.outliers$loci)
```

Write top 100 XtX loci to file

```
CHS.TI.XtX.run1 <- read.table("CHS.TI.Run1/XtX_out.CHS.TI.ENV", header=F)
CHS.TI.XtX.run2 <- read.table("CHS.TI.Run2/XtX_out.CHS.TI.ENV", header=F)
CHS.TI.XtX.run3 <- read.table("CHS.TI.Run3/XtX_out.CHS.TI.ENV", header=F)


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
##copy the list over to /Users/alexjvr/2016tempAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
XtX.outliers <- read.table("XtX.100outliers")
colnames(XtX.outliers) <- ("loci")
```




### 2. LFMM


Prepeare SE data: 

I got the new derived variables from Josh after getting estimated spawning dates from Maria: 

SE.temp.season.Data_no.indivs.20171009.csv   

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/SE.lfmm

```
library(LEA)
env <- read.csv("SE.temp.season.Data_no.indivs.20171009.csv", header=F)
write.env(env, "SE.CHP4.env")
```

Move the .geno and .env files onto the gdc server and run lfmm: 


```
scp * fgcz47:/srv/kenlab/alexjvr_p1795/CHcomplete/lfmm/lfmm.MAFfilter.Oct2017/CHP4/SE

##R
library(LEA)

project=lfmm("SE132.2027.recode.geno", "SE.CHP4.env", K=2, repetitions=5, project="new")

export.lfmmProject("SE132.2027.recode_SE_lfmmProject")
```



Import project into R

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/SE.lfmm

```
##R

project=import.lfmmProject("SE132.2027.recode_SE_lfmmProject.zip")  ##import the lfmm project in the directory



```

###### Postprocessing: 

Based on the lambda calculated from run 1, I have decided to run 3 K for each dataset. Lambda should be close to or just below 1. Mine were all just above 1. 

stba had to be removed from CHS.TI, because of missing data in the env. variables. So this is also removed from CHS. 



First calculate the lambda for each K to determine what the correct K is: 

##### SE

snmf K = 2


Open R and import project. 
```
library(LEA)

project=import.lfmmProject("SE132.2027.recode_SE_lfmmProject.zip")


zs.d1 <- z.scores(project, K=2, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 3.605604

    
lambda=median(zs.d1.median^2)/1.5
lambda
  #  [1] 1.093547

adj.p.values.d1 =pchisq(zs.d1.median^2/1.5, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=2, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1]  3.765671
lambda=median(zs.d2.median^2)/1.7
lambda
   # [1] 1.00773
adj.p.values.d2 =pchisq(zs.d2.median^2/1.7, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]



pdf(file="SE.LFMM.CHP4.hist.pdf")
par(mfrow=c(1,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][SE.hist]

[SE.hist]:https://user-images.githubusercontent.com/12142475/31121188-1e4eff0a-a82f-11e7-8060-e2ae5ec595dc.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("SE132.2017.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```


VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.SE.CHP4.n2.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange"))
dev.off()
```

![alt_txt][SE.Venn]

[SE.Venn]:https://user-images.githubusercontent.com/12142475/31121243-52d11a88-a82f-11e7-81b8-6cce216ca8e9.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "SE.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "SE.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)


colnames(d1.names) <- "names"
colnames(d2.names) <- "names"


d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12.names <- as.data.frame(d12.names)

write.table(d12.names$d12.names, "SE.CHP4.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

SE.LFMM.outliers <- read.table("SE.LFMM.alloutliers", header=F)

colnames(SE.LFMM.outliers) <- ("loci")
SE.LFMM.outliers <- as.character(SE.LFMM.outliers$loci)

```




##### CHN
snmf K = 4


Open R and import project. This CHN project contains K3-5. 5 runs for each K. d1-5. I show only results from K=4 here: 
```
library(LEA)

project=import.lfmmProject("CHN.229.5265.recode_CHN.temp_lfmmProject.zip")
#project=load.lfmmProject("CHN.229.5265.recode_CHN.temp.lfmmProject")

zs.d1 <- z.scores(project, K=4, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.536812

    
lambda=median(zs.d1.median^2)/0.7
lambda
  #  [1] 0.9987883

adj.p.values.d1 =pchisq(zs.d1.median^2/0.7, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=4, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1]  1.68876
lambda=median(zs.d2.median^2)/0.75
lambda
   # [1] 1.024371
adj.p.values.d2 =pchisq(zs.d2.median^2/0.75, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]



pdf(file="CHN.LFMM.CHP4.hist.pdf")
par(mfrow=c(1,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][CHN.hist]

[CHN.hist]:https://user-images.githubusercontent.com/12142475/31121188-1e4eff0a-a82f-11e7-8060-e2ae5ec595dc.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CHN.229.5265.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```


VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.CHN.CHP4.n2.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange"))
dev.off()
```

![alt_txt][CHN.Venn]

[CHN.Venn]:https://user-images.githubusercontent.com/12142475/31121243-52d11a88-a82f-11e7-81b8-6cce216ca8e9.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "CHN.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "CHN.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)



d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12.names <- as.data.frame(d12.names)

write.table(d12.names$d12.names, "CHN.CHP4.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CHN.LFMM.outliers <- read.table("CHN.LFMM.alloutliers", header=F)

colnames(CHN.LFMM.outliers) <- ("loci")
CHN.LFMM.outliers <- as.character(CHN.LFMM.outliers$loci)

```


##### CHS
snmf K = 10

These datasets all have less individuals than the full dataset. I have to find out why. 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/lfmm

Open R and import project. This CHN project contains K3-5. 5 runs for each K. d1-5. I show only results from K=4 here: 
```
CHS.lfmm=import.lfmmProject("CHS.275.6339.recode_CHS.temp_lfmmProject.zip")
#project=load.lfmmProject("CHS.275.6339.recode_CHS.temp.lfmmProject")

project <- CHS.lfmm

zs.d1 <- z.scores(project, K=10, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1]  1.913426
    
lambda=median(zs.d1.median^2)/0.85
lambda
  #  [1]  1.024103

adj.p.values.d1 =pchisq(zs.d1.median^2/0.85, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=10, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.798
lambda=median(zs.d2.median^2)/0.80
lambda
   # [1] 1.022469
adj.p.values.d2 =pchisq(zs.d2.median^2/0.80, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]



pdf(file="CHS.LFMM.CHP4.hist.pdf")
par(mfrow=c(1,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][CHS.hist]

[CHS.hist]:https://user-images.githubusercontent.com/12142475/31120625-00cf6f7a-a82d-11e7-8fce-6b238faf414d.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CHS.275.6339.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.CHS.n5.CHP4.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2,
cross.area=d12, 
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()
```

![alt_txt][CHS.Venn]

[CHS.Venn]:https://user-images.githubusercontent.com/12142475/31120683-408a8dfc-a82d-11e7-8de8-7a78297aceac.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "CHS.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "CHS.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)



colnames(d1.names) <- "names"
colnames(d2.names) <- "names"

d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists
d12.names <- as.data.frame(d12.names)
write.table(d12.names$d12.names, "CHS.LFMM.CHP4.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CHS.LFMM.outlier <- read.table("CHS.LFMM.alloutliers", header=F)
colnames(CHS.LFMM.outliers) <- ("loci")
CHS.LFMM.outliers <- as.character(CHS.LFMM.outliers$loci)
```


##### CHS.VS
snmf K = 4

These datasets all have less individuals than the full dataset. I have to find out why. 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/lfmm

Open R and import project.  
```
project=import.lfmmProject("CHS.VS.135.5835.recode_CHS.VS.temp_lfmmProject.zip")
#project=load.lfmmProject("CHS.VS.135.5835.recode_CHS.VS.temp.lfmmProject")

zs.d1 <- z.scores(project, K=4, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.708365
    
lambda=median(zs.d1.median^2)/0.75
lambda
  #  [1] 1.036263

adj.p.values.d1 =pchisq(zs.d1.median^2/0.75, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=4, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.507618
lambda=median(zs.d2.median^2)/0.65
lambda
   # [1]  1.055185
adj.p.values.d2 =pchisq(zs.d2.median^2/0.65, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]



pdf(file="CHS.VS.LFMM.CHP4.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][CHS.VS.hist]

[CHS.VS.hist]:https://user-images.githubusercontent.com/12142475/31120940-33b1cf9a-a82e-11e7-9aca-a77203ab8fbf.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CHS.VS.135.5835.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.CHS.VS.n2.CHP4.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()
```

![alt_txt][CHS.VS.Venn]

[CHS.VS.Venn]:https://user-images.githubusercontent.com/12142475/31120999-69f379dc-a82e-11e7-855f-97d8d828f30e.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "CHS.VS.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "CHS.VS.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)



colnames(d1.names) <- "names"
colnames(d2.names) <- "names"

d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists
d12.names <- as.data.frame(d12.names)
write.table(d12.names$d12.names, "CHS.VS.LFMM.CHP4.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CHS.VS.LFMM.outliers <- read.table("CHS.VS.LFMM.outliers", header=F)
colnames(CHS.VS.LFMM.outliers) <- ("loci")
CHS.VS.LFMM.outliers <- as.character(CHS.VS.LFMM.outliers$loci)

```


##### CHS.TI
snmf K = 6

These datasets all have less individuals than the full dataset. I have to find out why. 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/lfmm

Open R and import project.  
```
project=import.lfmmProject("CHS.TI.140.5692.recode_CHS.TI.temp_lfmmProject.zip")
#project=load.lfmmProject("CHS.TI.140.5692.recode_CHS.TI.temp.lfmmProject")

zs.d1 <- z.scores(project, K=6, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.246838
    
lambda=median(zs.d1.median^2)/0.55
lambda
  #  [1] 1.031331

adj.p.values.d1 =pchisq(zs.d1.median^2/0.55, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=6, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.243534
lambda=median(zs.d2.median^2)/0.55
lambda
   # [1]  1.028598
adj.p.values.d2 =pchisq(zs.d2.median^2/0.55, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]



pdf(file="CHS.TI.LFMM.CHP4.hist.pdf")
par(mfrow=c(1,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][CHS.TI.hist]

[CHS.TI.hist]:https://user-images.githubusercontent.com/12142475/31121410-00c997a0-a830-11e7-846f-75b8cc943689.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CHS.TI.140.5692.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.CHS.TI.n2.CHP4.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12,
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()
```

![alt_txt][CHS.TI.Venn]

[CHS.TI.Venn]:https://user-images.githubusercontent.com/12142475/31121490-47e83ea2-a830-11e7-8269-b799ae35fbcf.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "CHS.TI.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "CHS.TI.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)



colnames(d1.names) <- "names"
colnames(d2.names) <- "names"

d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12.names <- as.data.frame(d12.names)

write.table(d12.names$d12.names, "CHS.TI.LFMM.CHP4.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CHS.TI.LFMM.outliers <- read.table("CHS.TI.LFMM.alloutliers", header=F)
colnames(CHS.TI.LFMM.outliers) <- ("loci")
CHS.TI.LFMM.outliers <- as.character(CHS.TI.LFMM.outliers$loci)
```


##### CZ
snmf K = 9

These datasets all have less individuals than the full dataset. I have to find out why. 

/Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/lfmm

Open R and import project.  
```
project=import.lfmmProject("CZ.404.7288.recode_CZ.temp_lfmmProject.zip")
#project=load.lfmmProject("CZ.404.7288.recode_CZ.temp.lfmmProject")

zs.d1 <- z.scores(project, K=9, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.269226
    
lambda=median(zs.d1.median^2)/0.55
lambda
  #  [1] 1.049849

adj.p.values.d1 =pchisq(zs.d1.median^2/0.55, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=9, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.469544
lambda=median(zs.d2.median^2)/0.65
lambda
   # [1]  1.028537
adj.p.values.d2 =pchisq(zs.d2.median^2/0.65, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]


pdf(file="CZ.LFMM.CHP4.hist.pdf")
par(mfrow=c(1,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
dev.off()

```


![alt_txt][CZ.hist]

[CZ.hist]:https://user-images.githubusercontent.com/12142475/31121666-170a4a90-a831-11e7-9bd9-660ec8b1bbfd.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CZ.404.7288.plink.map", header=F)
colnames(locus.names) <- c("V1", "SNP", "V3", "V4")
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers, since the LFMM output is a numbered list corresponding to the original genotype input order
candidates.d1.k5 <- as.character(candidates.d1.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d1.k5.names <- locus.names[locus.names$ID %in% candidates.d1.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d1.k5.names <- paste("X", candidates.d1.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel


candidates.d2.k5 <- as.character(candidates.d2.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d2.k5.names <- locus.names[locus.names$ID %in% candidates.d2.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d2.k5.names <- paste("X", candidates.d2.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel

```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))

pdf(file="Venn.CZ.n2.CHP4.LFMMonly_20171009.pdf")
draw.pairwise.venn(area1=d1, area2=d2, cross.area=d12, 
category=c("mean.temp", "season.length"),
lty="blank", 
fill=c("yellow", "orange")
)
dev.off()
```

![alt_txt][CZ.Venn]

[CZ.Venn]:https://user-images.githubusercontent.com/12142475/31121744-612b5f2e-a831-11e7-9af5-bd66a3c3c470.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
temp.names <- d1.names
colnames(temp.names) <- "names"
temp.names <- sub(":", ".", temp.names$names) 
temp.names <- as.data.frame(temp.names)
write.table(temp.names$temp.names, "CZ.CHP4.temp.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

d2.names <- as.data.frame(candidates.d2.k5.names)
season.names <- d2.names
colnames(season.names) <- "names"
season.names <- sub(":", ".", season.names$names) 
season.names <- as.data.frame(season.names)
write.table(season.names$season.names, "CZ.CHP4.season.LFMM.alloutliers", col.names=F, row.names=F, quote=F)



colnames(d1.names) <- "names"
colnames(d2.names) <- "names"

d12.names <- rbind(d1.names, d2.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12.names <- lapply(d12.names, unique)  #select only the unique rows. 

d12.names <- sub(":", ".", d12.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12.names <- as.data.frame(d12.names)

write.table(d12.names$d12.names, "CZ.LFMM.CHP4.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CZ.LFMM.outliers <- read.table("CZ.LFMM.alloutliers", header=F)
colnames(CZ.LFMM.outliers) <- ("loci")
CZ.LFMM.outliers <- as.character(CZ.LFMM.outliers$loci)
```




## Fst outlier analyses

### 1. Bayescan

CH: Results as in CHP2

SE: Ran the new SE dataset on the sork lab computer. It took about 8 hours: 

copied everything here: /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/SE.bayescan

```
source("plot_R.r")
library(coda)

chain <- read.table("SE132.out.sel", header=T)
chain <- mcmc(chain, thin=10)

plot(chain)  ##check for convergence
summary(chain)
autocorr.diag(chain) ## check correlation between the chains. Make sure the chains didn't get stuck
effectiveSize(chain) ##check that this is close to the sample size (here 5000). If there is correlation (chain got stuck) the sample size will be much smaller than the input
geweke.diag(chain, frac1=0.1, frac2=0.5)  ##The diagnostic reports the z-scores for each parameter. For example, with α = 0.05, the critical values of z are – 1.96 and +1.96. We reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96.

heidel.diag(chain, eps=0.1, pvalue=0.05) ##another test whether the chains have reached stationarity. 

Find the outliers and rename them
SE.results <- read.table("SE132.out_fst.txt")

pdf("SE.bayescan.results.pdf")
par(mfrow=c(1,2))
plot_bayescan(SE.results, FDR=0.05, add_text=F)
plot_bayescan(SE.results, FDR=0.01, add_text=F)
dev.off()


SE.results$rownumber <- 1:nrow(SE.results)  ##add an index of the rownumbers, as these correspond to the loci in the input file. 

SE.outliers <- plot_bayescan(SE.results, FDR=0.01) # this is an R script distributed with bayescan for plotting and identifying outliers. In this case I find 101 when FDR=0.05 and 61 at FDR=0.01
SE.outliers.df <- SE.outliers$outliers
SE.outliers.df <- as.data.frame(SE.outliers.df)
colnames(SE.outliers.df) <- "rownames"  ##the column headers need to be the same for dplyr to work

locus.names <- read.table("SE132.2027.plink.map") ##read in the locus names found in the map file generated vcf --plink
locus.names$rownames <- 1:nrow(locus.names) ##index in the same way as the bayescan output

library(dplyr)
##join by rownames
SE.outlier.names <- semi_join(locus.names, SE.outliers.df)
SE.outlier.names <- gsub(":", ".", SE.outlier.names$V2 )

SE.outlier.names <- as.data.frame(SE.outlier.names)
write.table(SE.outlier.names, "SE.bayescan.outliers.FDR0.01", quote=F, row.names=F, col.names=F)

```

### 2. PCAdapt

CH: results as in CHP2. 

SE: /Users/alexjvr/2016RADAnalysis/6_CHP4.SEvsCH/SE.PCAdapt

```
##convert .vcf to plink 
##linux

vcftools --vcf SE132.2027.recode.vcf --plink --out SE132.2027.plink

plink --file SE132.2027.plink --recode --recodeA --out SE132.2027.plink

##R
library(pcadapt)

SE.132 <- read.pcadapt("SE132.2027.plink.ped", type="ped")
Summary:

        - input file      SE132.2027.plink.ped
        - output file     SE132.2027.plink.pcadapt

	- number of individuals detected:	132
	- number of loci detected:		2027

File has been sucessfully converted.

##Check the nr of PCs

x <- pcadapt(SE.132, K=20)

Reading file SE132.2027.plink.pcadapt...
Number of SNPs: 2027
Number of individuals: 132
Number of SNPs with minor allele frequency lower than 0.05 ignored: 3
16994 out of 267564 missing data ignored.

plot(x, option="screeplot")  ##PC for pop structure = on the steep curve

```

![alt_txt][SE.PCA]

[SE.PCA]:https://user-images.githubusercontent.com/12142475/31168411-21289d30-a8ed-11e7-8f19-404a2f8ea14a.png

Based on these results I'm choosing K=4. In the SE chapter I chose K=5, but after re-reading the manual, I realise that I should choose K as the last point before proportion of variance explained flattens out. 

#### PCAdapt outlier loci

Vignette: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

Useful info about interpreting p-value distributions: http://varianceexplained.org/statistics/interpreting-pvalue-histogram/

I did this in R following the vignette quite closely:

I chose K=4 following the scree plot from PCAdapt and DAPC results

```
x.maf0.05 <- pcadapt(SE.132, K=4, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.maf0.1 <- pcadapt(SE.132, K=4, min.maf=0.1) #calculate for maf 0.1

pdf(file="pcadapt.pvalues.pdf")
par(mfrow=c(2,1))
hist(x.maf0.05$pvalues,xlab="p-values maf0.05",main=NULL,breaks=50)
hist(x.maf0.1$pvalues,xlab="p-values maf0.1",main=NULL,breaks=50)
dev.off()
```

p-value distributions: 

![alt_txt][SE.p.val]

[SE.p.val]:https://user-images.githubusercontent.com/12142475/31168584-bebecdda-a8ed-11e7-95c3-a15a2bd20367.png


I used x.maf0.05 since it has a slightly more uniform p-value distribution

use qvalue in R to specify FDR

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.maf0.05$pvalues)$qvalues
outliers.SE <- which(qval<alpha)
outliers.SE
SE.snp_pc <- get.pc(x.maf0.05,outliers.SE) ##see PCs associated with the outliers
```
60 outliers identified.

Rename the loci
```
locus.names <- read.table("SE132.2027.plink.map", header=F)
locus.names
locus.names$ID <- seq.int(nrow(locus.names)) #add an index of the SNP numbers
outliers.pcadapt <- as.character(outliers.SE)
outliers.pcadapt.names <- locus.names[locus.names$ID %in% outliers.pcadapt,]
outliers.pcadapt.names <- paste("X", outliers.pcadapt.names$V2, sep=".")

write.table(outliers.pcadapt.names, "SE.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table 

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R
pcadapt.outliers <- read.table("pcadapt.outliers")
colnames(pcadapt.outliers) <- ("loci")
pcadapt.outliers <- as.character(pcadapt.outliers$loci)
```



## RDA

Here I want to see how much of the variation is explained by the candidate loci

## GDM & GF

The aim here is to

1. determine whether temp or season length are the most important in all transects

2. identify any thresholds and compare between transects

3. predict affect of future climates on adaptation: are there more vulnerable regions? 




