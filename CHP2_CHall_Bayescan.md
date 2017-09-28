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

