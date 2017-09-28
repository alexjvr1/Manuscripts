# PCAdapt: CHall data

I want to identify all outlier loci for each of the datasets using PCAdapt

1. CHall.932.7744.recode.vcf

2. CHN.229.5265.recode.vcf

3. CHS.275.6339.recode.vcf

4. CHS.VS.135.5835.recode.vcf

5. CHS.TI.140.5692.recode.vcf

6. CZ.404.7288.recode.vcf


input files .ped from plink

copy from input.files: 

```
cp /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/RDA/Sept2017/AnalysisFullDataset/*ped /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/PCadapt/Oct2017/

```


## 1. CHall

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHall <- read.pcadapt("CHall.932.7744.plink.ped", type="ped")
Summary:

        - input file      CHall.932.7744.plink.ped
        - output file     CHall.932.7744.plink.pcadapt

	- number of individuals detected:	932
	- number of loci detected:		7744

File has been sucessfully converted.

##Check the nr of PCs

x.CHall <- pcadapt(CHall, K=20)

Reading file CHall.932.7744.plink.pcadapt...
Number of SNPs: 7744
Number of individuals: 932
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
2148762 out of 7217408 missing data ignored.

plot(x.CHall, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3. The random eigenvalues fall in the straight line, while the ones explaining structure are above this. 

See vignette: https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html on choosing K. 

Plot the PCA using population information

```
pop.CHall <- read.table("CH932.cluster.pop", header=F)
pop.CHall
poplist <- as.character(pop.CHall[,3])
poplist
plot(x.CHall, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=5 following the scree plot from PCAdapt and DAPC results

```
x.CHall.maf0.05 <- pcadapt(CHall, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHall.maf0.1 <- pcadapt(CHall, K=3, min.maf=0.1) #calculate for maf 0.1

pdf(file="pcadapt.pvalues.CHall.pdf")
par(mfrow=c(2,1))
hist(x.CHall.maf0.05$pvalues,xlab="p-values CHall maf0.05",main=NULL,breaks=50)
hist(x.CHall.maf0.1$pvalues,xlab="p-values CHall maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CHall.maf0.05$pvalues)$qvalues
outliers.CHall <- which(qval<alpha)
outliers.CHall
CHall.snp_pc <- get.pc(x.CHall.maf0.05,outliers.CHall) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHall <- read.table("CHall.932.7744.plink.map", header=F)
locus.names.CHall
locus.names.CHall$ID <- seq.int(nrow(locus.names.CHall)) #add an index of the SNP numbers
CHall.outliers.pcadapt <- as.character(outliers.CHall)
CHall.outliers.pcadapt.names <- locus.names.CHall[locus.names.CHall$ID %in% CHall.outliers.pcadapt,]
CHall.outliers.pcadapt.names <- paste("X", CHall.outliers.pcadapt.names$V2, sep=".")

write.table(CHall.outliers.pcadapt.names, "CHall.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```



## 2. CHN

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHN <- read.pcadapt("CHN.229.9608.plink.ped", type="ped")
Summary:

        - input file      CHN.229.5265.plink.ped
        - output file     CHN.229.5265.plink.pcadapt

	- number of individuals detected:	229
	- number of loci detected:		5265

File has been sucessfully converted.

##Check the nr of PCs

x.CHN <- pcadapt(CHN, K=20)

Reading file CHN.229.5265.plink.pcadapt...
Number of SNPs: 5265
Number of individuals: 229
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
365454 out of 1205685 missing data ignored.

plot(x.CHN, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=4

Plot the PCA using population information

```
pop.CHN <- read.table("CHN229.cluster.pop", header=F)
pop.CHN
poplist <- as.character(pop.CHN[,3])
poplist
plot(x.CHN, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=4 following the scree plot from PCAdapt and DAPC results

```
x.CHN.maf0.05 <- pcadapt(CHN, K=4, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHN.maf0.1 <- pcadapt(CHN, K=4, min.maf=0.1) #calculate for maf 0.1
Reading file CHN.229.5265.plink.pcadapt...
Number of SNPs: 5265
Number of individuals: 229
Number of SNPs with minor allele frequency lower than 0.1 ignored: 712
333699 out of 1205685 missing data ignored.

pdf(file="pcadapt.pvalues.CHN.pdf")
par(mfrow=c(2,1))
hist(x.CHN.maf0.05$pvalues,xlab="p-values CHN maf0.05",main=NULL,breaks=50)
hist(x.CHN.maf0.1$pvalues,xlab="p-values CHN maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CHN.maf0.05$pvalues)$qvalues
outliers.CHN <- which(qval<alpha)
outliers.CHN
CHN.snp_pc <- get.pc(x.CHN.maf0.05, outliers.CHN) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHN <- read.table("CHN.229.5265.plink.map", header=F)
locus.names.CHN
locus.names.CHN$ID <- seq.int(nrow(locus.names.CHN)) #add an index of the SNP numbers
CHN.outliers.pcadapt <- as.character(outliers.CHN)
CHN.outliers.pcadapt.names <- locus.names.CHN[locus.names.CHN$ID %in% CHN.outliers.pcadapt,]
CHN.outliers.pcadapt.names <- paste("X", CHN.outliers.pcadapt.names$V2, sep=".")

write.table(CHN.outliers.pcadapt.names, "CHN.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```


## 3.CHS

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHS <- read.pcadapt("CHS.275.6339.plink.ped", type="ped")
Summary:

        - input file      CHS.275.6339.plink.ped
        - output file     CHS.275.6339.plink.pcadapt

	- number of individuals detected:	275
	- number of loci detected:		6339

File has been sucessfully converted.

##Check the nr of PCs

x.CHS <- pcadapt(CHS, K=20)

Reading file CHS.275.6339.plink.pcadapt...
Number of SNPs: 6339
Number of individuals: 275
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
529915 out of 1743225 missing data ignored.

plot(x.CHS, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3

Plot the PCA using population information

```
pop.CHS <- read.table("CHS275.cluster.pop", header=F)
pop.CHS
poplist <- as.character(pop.CHS[,3])
poplist
plot(x.CHS, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=3 following the scree plot from PCAdapt and DAPC results

```
x.CHS.maf0.05 <- pcadapt(CHS, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.maf0.1 <- pcadapt(CHS, K=3, min.maf=0.1) #calculate for maf 0.1

pdf(file="pcadapt.pvalues.CHS.pdf")
par(mfrow=c(2,1))
hist(x.CHS.maf0.05$pvalues,xlab="p-values CHS maf0.05",main=NULL,breaks=50)
hist(x.CHS.maf0.1$pvalues,xlab="p-values CHS maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CHS.maf0.05$pvalues)$qvalues
outliers.CHS <- which(qval<alpha)
outliers.CHS
CHS.snp_pc <- get.pc(x.CHS.maf0.05,outliers.CHS) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS <- read.table("CHS.275.6339.plink.map", header=F)
locus.names.CHS
locus.names.CHS$ID <- seq.int(nrow(locus.names.CHS)) #add an index of the SNP numbers
CHS.outliers.pcadapt <- as.character(outliers.CHS)
CHS.outliers.pcadapt.names <- locus.names.CHS[locus.names.CHS$ID %in% CHS.outliers.pcadapt,]
CHS.outliers.pcadapt.names <- paste("X", CHS.outliers.pcadapt.names$V2, sep=".")

write.table(CHS.outliers.pcadapt.names, "CHS.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```


## 4.CHS.VS

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHS.VS <- read.pcadapt("CHS.VS.135.5835.plink.ped", type="ped")
Summary:

        - input file      CHS.VS.135.5835.plink.ped
        - output file     CHS.VS.135.5835.plink.pcadapt

	- number of individuals detected:	135
	- number of loci detected:		5835

File has been sucessfully converted.

##Check the nr of PCs

x.CHS.VS <- pcadapt(CHS.VS, K=20)
Reading file CHS.VS.135.5835.plink.pcadapt...
Number of SNPs: 5835
Number of individuals: 135
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
243032 out of 787725 missing data ignored.

plot(x.CHS.VS, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=7

Plot the PCA using population information

```
pop.CHS.VS <- read.table("CHS.VS.cluster.pop", header=F)
pop.CHS.VS
poplist <- as.character(pop.CHS.VS[,3])
poplist
plot(x.CHS.VS, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=4 following the scree plot from PCAdapt and DAPC results

```
x.CHS.VS.maf0.05 <- pcadapt(CHS.VS, K=7, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.VS.maf0.1 <- pcadapt(CHS.VS, K=7, min.maf=0.1) #calculate for maf 0.1
Reading file CHS.VS.135.5835.plink.pcadapt...
Number of SNPs: 5835
Number of individuals: 135
Number of SNPs with minor allele frequency lower than 0.1 ignored: 631
227436 out of 787725 missing data ignored.


pdf(file="pcadapt.pvalues.CHS.VS.pdf")
par(mfrow=c(2,1))
hist(x.CHS.VS.maf0.05$pvalues,xlab="p-values CHS.VS maf0.05",main=NULL,breaks=50)
hist(x.CHS.VS.maf0.1$pvalues,xlab="p-values CHS.VS maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CHS.VS.maf0.05$pvalues)$qvalues
outliers.CHS.VS <- which(qval<alpha)
outliers.CHS.VS
CHS.VS.snp_pc <- get.pc(x.CHS.VS.maf0.05,outliers.CHS.VS) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS.VS <- read.table("CHS.VS.135.5835.plink.map", header=F)
locus.names.CHS.VS
locus.names.CHS.VS$ID <- seq.int(nrow(locus.names.CHS.VS)) #add an index of the SNP numbers
CHS.VS.outliers.pcadapt <- as.character(outliers.CHS.VS)
CHS.VS.outliers.pcadapt.names <- locus.names.CHS.VS[locus.names.CHS.VS$ID %in% CHS.VS.outliers.pcadapt,]
CHS.VS.outliers.pcadapt.names <- paste("X", CHS.VS.outliers.pcadapt.names$V2, sep=".")

write.table(CHS.VS.outliers.pcadapt.names, "CHS.VS.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```


## 5.CHS.TI

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHS.TI <- read.pcadapt("CHS.TI.140.5692.plink.ped", type="ped")
Summary:

        - input file      CHS.TI.140.5692.plink.ped
        - output file     CHS.TI.140.5692.plink.pcadapt

	- number of individuals detected:	140
	- number of loci detected:		5692

File has been sucessfully converted.

##Check the nr of PCs

x.CHS.TI <- pcadapt(CHS.TI, K=20)

Reading file CHS.TI.140.5692.plink.pcadapt...
Number of SNPs: 5692
Number of individuals: 140
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
233871 out of 796880 missing data ignored.

plot(x.CHS.TI, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=5

Plot the PCA using population information

```
pop.CHS.TI <- read.table("CHS.TI.148.cluster.pop", header=F)
pop.CHS.TI
poplist <- as.character(pop.CHS.TI[,3])
poplist
plot(x.CHS.TI, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=5 following the scree plot from PCAdapt and DAPC results

```
x.CHS.TI.maf0.05 <- pcadapt(CHS.TI, K=5, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.TI.maf0.1 <- pcadapt(CHS.TI, K=5, min.maf=0.1) #calculate for maf 0.1
Reading file CHS.TI.140.5692.plink.pcadapt...
Number of SNPs: 5692
Number of individuals: 140
Number of SNPs with minor allele frequency lower than 0.1 ignored: 671
214993 out of 796880 missing data ignored.

pdf(file="pcadapt.pvalues.CHS.TI.pdf")
par(mfrow=c(2,1))
hist(x.CHS.TI.maf0.05$pvalues,xlab="p-values CHS.TI maf0.05",main=NULL,breaks=50)
hist(x.CHS.TI.maf0.1$pvalues,xlab="p-values CHS.TI maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CHS.TI.maf0.05$pvalues)$qvalues
outliers.CHS.TI <- which(qval<alpha)
outliers.CHS.TI
CHS.TI.snp_pc <- get.pc(x.CHS.TI.maf0.05,outliers.CHS.TI) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS.TI <- read.table("CHS.TI.140.5692.plink.map", header=F)
locus.names.CHS.TI
locus.names.CHS.TI$ID <- seq.int(nrow(locus.names.CHS.TI)) #add an index of the SNP numbers
CHS.TI.outliers.pcadapt <- as.character(outliers.CHS.TI)
CHS.TI.outliers.pcadapt.names <- locus.names.CHS.TI[locus.names.CHS.TI$ID %in% CHS.TI.outliers.pcadapt,]
CHS.TI.outliers.pcadapt.names <- paste("X", CHS.TI.outliers.pcadapt.names$V2, sep=".")

write.table(CHS.TI.outliers.pcadapt.names, "CHS.TI.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```

## 6.CZ

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CZ <- read.pcadapt("CZ.404.7288.plink.ped", type="ped")
Summary:

        - input file      CZ.404.7288.plink.ped
        - output file     CZ.404.7288.plink.pcadapt

	- number of individuals detected:	404
	- number of loci detected:		7288

File has been sucessfully converted.

##Check the nr of PCs

x.CZ <- pcadapt(CZ, K=20)

Reading file CZ.404.7288.plink.pcadapt...
Number of SNPs: 7288
Number of individuals: 404
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
869474 out of 2944352 missing data ignored.

plot(x.CZ, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3

Plot the PCA using population information

```
pop.CZ <- read.table("CZ404.cluster.pop", header=F)
pop.CZ
poplist <- as.character(pop.CZ[,3])
poplist
plot(x.CZ, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=3 following the scree plot from PCAdapt and DAPC results

```
x.CZ.maf0.05 <- pcadapt(CZ, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CZ.maf0.1 <- pcadapt(CZ, K=3, min.maf=0.1) #calculate for maf 0.1
Reading file CZ.404.7288.plink.pcadapt...
Number of SNPs: 7288
Number of individuals: 404
Number of SNPs with minor allele frequency lower than 0.1 ignored: 1013
774179 out of 2944352 missing data ignored.

pdf(file="pcadapt.pvalues.CZ.pdf")
par(mfrow=c(2,1))
hist(x.CZ.maf0.05$pvalues,xlab="p-values CZ maf0.05",main=NULL,breaks=50)
hist(x.CZ.maf0.1$pvalues,xlab="p-values CZ maf0.1",main=NULL,breaks=50)
dev.off()
```



I used x.maf0.05 since the p-distribution was flat

```
library(qvalue)
alpha <- 0.05  ##FDR
qval <- qvalue(x.CZ.maf0.05$pvalues)$qvalues
outliers.CZ <- which(qval<alpha)
outliers.CZ
CZ.snp_pc <- get.pc(x.CZ.maf0.05,outliers.CZ) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CZ <- read.table("CZ.404.7288.plink.map", header=F)
locus.names.CZ
locus.names.CZ$ID <- seq.int(nrow(locus.names.CZ)) #add an index of the SNP numbers
CZ.outliers.pcadapt <- as.character(outliers.CZ)
CZ.outliers.pcadapt.names <- locus.names.CZ[locus.names.CZ$ID %in% CZ.outliers.pcadapt,]
CZ.outliers.pcadapt.names <- paste("X", CZ.outliers.pcadapt.names$V2, sep=".")

write.table(CZ.outliers.pcadapt.names, "CZ.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```

