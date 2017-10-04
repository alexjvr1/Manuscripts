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


Input files for CH have been generated in CHP2. 

Env input is taken from the data that Josh calculated. 

All runs on the GDC server: 

```

```


### 2. LFMM




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




