# Bayescan

Since the EAA analyses didn't yield very convincing results, I've decided to identify outlier loci using Bayescan. 

Website: http://cmpg.unibe.ch/software/BayeScan/

This is available on the Sork computers:

soney@sorklab3.eeb.ucla.edu

signem123

Check with Sorel before running anything: sorel@ucla.edu

more programs in ~/data


## Input files

working directory: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/Bayescan/bayescan.Oct2017

Input files are allele frequencies. They can be produced using pgdspider. 

I need allele frequencies per population since this is the level at which my environmental data is sampled. 

I'll convert vcf files to bayescan input files. 

popfiles for pgdspider can be obtained as follows

```
###CHall

##get a list of names from each vcf file
bcftools query -l CHall.932.7744.recode.vcf > CHall.popnames

##read into R
CHall.popnames <- read.table("CHall.popnames", header=F)
CHall.pops <- gsub("_\\d+", "", CHall.popnames$V1)
CHall.pops <- as.data.frame(CHall.pops)
head(CHall.pops)
CHall.popnames$V2 <- CHall.pops$CHall.pops
write.table(CHall.popnames, "CHall.932.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CHall.popnames


###CHN

##get a list of names from each vcf file
bcftools query -l CHN.229.5265.recode.vcf > CHN.popnames

##read into R
CHN.popnames <- read.table("CHN.popnames", header=F)
CHN.pops <- gsub("_\\d+", "", CHN.popnames$V1)
CHN.pops <- as.data.frame(CHN.pops)
head(CHN.pops)
CHN.popnames$V2 <- CHN.pops$CHN.pops
write.table(CHN.popnames, "CHN.229.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CHN.popnames





###CHS

##get a list of names from each vcf file
bcftools query -l CHS.275.6339.recode.vcf > CHS.popnames

##read into R
CHS.popnames <- read.table("CHS.popnames", header=F)
CHS.pops <- gsub("_\\d+", "", CHS.popnames$V1)
CHS.pops <- as.data.frame(CHS.pops)
head(CHS.pops)
CHS.popnames$V2 <- CHS.pops$CHS.pops
write.table(CHS.popnames, "CHS.275.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CHS.popnames




###CHS.VS

##get a list of names from each vcf file
bcftools query -l CHS.VS.229.5265.recode.vcf > CHS.VS.popnames

##read into R
CHS.VS.popnames <- read.table("CHS.VS.popnames", header=F)
CHS.VS.pops <- gsub("_\\d+", "", CHS.VS.popnames$V1)
CHS.VS.pops <- as.data.frame(CHS.VS.pops)
head(CHS.VS.pops)
CHS.VS.popnames$V2 <- CHS.VS.pops$CHS.VS.pops
write.table(CHS.VS.popnames, "CHS.VS.135.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CHS.VS.popnames




###CHS.TI

##get a list of names from each vcf file
bcftools query -l CHS.TI.140.5692.recode.vcf > CHS.TI.popnames

##read into R
CHS.TI.popnames <- read.table("CHS.TI.popnames", header=F)
CHS.TI.pops <- gsub("_\\d+", "", CHS.TI.popnames$V1)
CHS.TI.pops <- as.data.frame(CHS.TI.pops)
head(CHS.TI.pops)
CHS.TI.popnames$V2 <- CHS.TI.pops$CHS.TI.pops
write.table(CHS.TI.popnames, "CHS.TI.140.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CHS.TI.popnames




###CZ
##get a list of names from each vcf file
bcftools query -l CZ.404.7288.recode.vcf > CZ.popnames

##read into R
CZ.popnames <- read.table("CZ.popnames", header=F)
CZ.pops <- gsub("_\\d+", "", CZ.popnames$V1)
CZ.pops <- as.data.frame(CZ.pops)
head(CZ.pops)
CZ.popnames$V2 <- CZ.pops$CZ.pops
write.table(CZ.popnames, "CZ.404.spidpopfile", col.names=F, row.names=F, quote=F)
q()

##linux: remove the popfile from bcftools
rm CZ.popnames

```


## Run Bayescan

All the runs were completed on the UCLA server. Shorter runs took around 8 hours. Longest run was around 16 hours. 

I copied all the output to here: 

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/Bayescan/bayescan.Oct2017/CHP2.bayescan

```
bayescan input.bayescan.txt -o input.out -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 20

```

pr_odds here is 20, i.e. 5%, which is the cut-off I used in all the other methods. 


## Interpretation

1. Check that the runs have converged

2. Identify outliers. 

I'm using the following as guidance: 

https://evomics.org/learning/population-and-speciation-genomics/bayescan-exercise/



```
source("plot_R.r")
library(coda)

chain <- read.table("CHN.out.sel", header=T)
chain <- mcmc(chain, thin=10)

plot(chain)  ##check for convergence
summary(chain)
autocorr.diag(chain) ## check correlation between the chains. Make sure the chains didn't get stuck
effectiveSize(chain) ##check that this is close to the sample size (here 5000). If there is correlation (chain got stuck) the sample size will be much smaller than the input
geweke.diag(chain, frac1=0.1, frac2=0.5)  ##The diagnostic reports the z-scores for each parameter. For example, with α = 0.05, the critical values of z are – 1.96 and +1.96. We reject H0 (equality of means => convergence) if z < -1.96 or z > +1.96.

heidel.diag(chain, eps=0.1, pvalue=0.05) ##another test whether the chains have reached stationarity. 

plot_bayescan("CHN.out_fst.txt", FDR=0.05)
plot_bayescan("CHN.out_fst.txt", FDR=0.01)
```


Find the outliers and rename them
```
CHN.results <- read.table("CHN.out_fst.txt")
CHN.results$rownumber <- 1:nrow(CHN.results)  ##add an index of the rownumbers, as these correspond to the loci in the input file. 

CHN.outliers <- plot_bayescan(CHN.results, FDR=0.05) # this is an R script distributed with bayescan for plotting and identifying outliers. In this case I find 318 when FDR=0.05 and 217 at FDR=0.01
CHN.outliers.df <- CHN.outliers$outliers
CHN.outliers.df <- as.data.frame(CHN.outliers.df)
colnames(CHN.outliers.df) <- "rownames"  ##the column headers need to be the same for dplyr to work

locus.names <- read.table("CHN.229.5265.plink.map") ##read in the locus names found in the map file generated vcf --plink
locus.names$rownames <- 1:nrow(locus.names) ##index in the same way as the bayescan output

library(dplyr)
##join by rownames
CHN.outlier.names <- semi_join(locus.names, CHN.outliers.df)
CHN.outlier.names <- gsub(":", ".", CHN.outlier.names$V2 )

CHN.outlier.names <- as.data.frame(CHN.outlier.names)
write.table(CHN.outlier.names, "CHN.bayenv.outliers.FDR0.05", quote=F, row.names=F, col.names=F)
```

