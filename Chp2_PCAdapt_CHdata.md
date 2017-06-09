# PCAdapt: CHall data

I want to identify all outlier loci for each of the datasets using PCAdapt

1. CHall.932.9608

2. CHN.229.9608

3. CHS.283.9608

4. CHS.VS.135.9608

5. CHS.TI.148.9608

6. CZ.404.9608


input files .ped from plink

copy from input.files: 

```
cp /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/*ped /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/PCadapt/

```

## 2. CHN

### 1. Calculate the number of Principal Components

```
library(pcadapt)

CHN <- read.pcadapt("CHN.229.9608.plink.ped", type="ped")
Summary:

        - input file      CHN.229.9608.plink.ped
        - output file     CHN.229.9608.plink.pcadapt

	- number of individuals detected:	229
	- number of loci detected:		9608

File has been sucessfully converted.

##Check the nr of PCs

x.CHN <- pcadapt(CHN, K=20)

Reading file CHN.229.9608.plink.pcadapt...
Number of SNPs: 9608
Number of individuals: 229
Number of SNPs with minor allele frequency lower than 0.05 ignored: 3648
439399 out of 2200232 missing data ignored.

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
x.CHN.maf0.05 <- pcadapt(x.CHN, K=5, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x..CHN.maf0.1 <- pcadapt(x.CHN, K=5, min.maf=0.1) #calculate for maf 0.1

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
qval <- qvalue(x.CHN$pvalues)$qvalues
outliers.CHN <- which(qval<alpha)
outliers.CHN
snp_pc <- get.pc(x.CHN.maf0.1,outliers) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHN <- read.table("CHN.229.9608.plink.map", header=F)
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

CHS <- read.pcadapt("CHS.283.9608.plink.ped", type="ped")
Summary:

        - input file      CHS.283.9608.plink.ped
        - output file     CHS.283.9608.plink.pcadapt

	- number of individuals detected:	283
	- number of loci detected:		9608

File has been sucessfully converted.

##Check the nr of PCs

x.CHS <- pcadapt(CHS, K=20)

Reading file CHS.283.9608.plink.pcadapt...
Number of SNPs: 9608
Number of individuals: 283
Number of SNPs with minor allele frequency lower than 0.05 ignored: 2509
640938 out of 2719064 missing data ignored.

plot(x.CHS, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3

Plot the PCA using population information

```
pop.CHS <- read.table("CHS283.cluster.pop", header=F)
pop.CHS
poplist <- as.character(pop.CHS[,3])
poplist
plot(x.CHS, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=3 following the scree plot from PCAdapt and DAPC results

```
x.CHS.maf0.05 <- pcadapt(CHS, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.maf0.1 <- pcadapt(x.CHS, K=3, min.maf=0.1) #calculate for maf 0.1

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
qval <- qvalue(x.CHS$pvalues)$qvalues
outliers.CHS <- which(qval<alpha)
outliers.CHS
snp_pc <- get.pc(x.CHS.maf0.05,outliers.CHS) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS <- read.table("CHS.283.9608.plink.map", header=F)
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

CHS.VS <- read.pcadapt("CHS.VS.135.plink.ped", type="ped")
Summary:

        - input file      CHS.VS.135.plink.ped
        - output file     CHS.VS.135.plink.pcadapt

	- number of individuals detected:	135
	- number of loci detected:		9608

File has been sucessfully converted.

##Check the nr of PCs

x.CHS.VS <- pcadapt(CHS.VS, K=20)
Reading file CHS.VS.135.plink.pcadapt...
Number of SNPs: 9608
Number of individuals: 135
Number of SNPs with minor allele frequency lower than 0.05 ignored: 3015
288950 out of 1297080 missing data ignored.

plot(x.CHS.VS, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3

Plot the PCA using population information

```
pop.CHS.VS <- read.table("CHS.VS.cluster.pop", header=F)
pop.CHS.VS
poplist <- as.character(pop.CHS.VS[,3])
poplist
plot(x.CHS.VS, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=3 following the scree plot from PCAdapt and DAPC results

```
x.CHS.VS.maf0.05 <- pcadapt(CHS.VS, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.VS.maf0.1 <- pcadapt(CHS.VS, K=3, min.maf=0.1) #calculate for maf 0.1

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
qval <- qvalue(x.CHS.VS$pvalues)$qvalues
outliers.CHS.VS <- which(qval<alpha)
outliers.CHS.VS
snp_pc <- get.pc(x.CHS.VS.maf0.05,outliers.CHS.VS) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS.VS <- read.table("CHS.VS.283.9608.plink.map", header=F)
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

CHS.TI <- read.pcadapt("CHS.TI.148.plink.ped", type="ped")
Summary:

        - input file      CHS.TI.148.plink.ped
        - output file     CHS.TI.148.plink.pcadapt

	- number of individuals detected:	148
	- number of loci detected:		9608

File has been sucessfully converted.

##Check the nr of PCs

x.CHS.TI <- pcadapt(CHS.TI, K=20)

Reading file CHS.TI.148.plink.pcadapt...
Number of SNPs: 9608
Number of individuals: 148
Number of SNPs with minor allele frequency lower than 0.05 ignored: 3310
295646 out of 1421984 missing data ignored.

plot(x.CHS.TI, option="screeplot")  ##PC for pop structure = on the steep curve
```

Based on this I choose K=3

Plot the PCA using population information

```
pop.CHS.TI <- read.table("CHS.TI.148.cluster.pop", header=F)
pop.CHS.TI
poplist <- as.character(pop.CHS.TI[,3])
poplist
plot(x.CHS.TI, option="scores", pop=poplist)
```


### 2.Determine outliers

I chose K=3 following the scree plot from PCAdapt and DAPC results

```
x.CHS.TI.maf0.05 <- pcadapt(CHS.TI, K=3, min.maf=0.05)   ##calculate z-statistics and transformed values for chi-squared distribution
x.CHS.TI.maf0.1 <- pcadapt(CHS.TI, K=3, min.maf=0.1) #calculate for maf 0.1

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
qval <- qvalue(x.CHS.TI$pvalues)$qvalues
outliers.CHS.TI <- which(qval<alpha)
outliers.CHS.TI
snp_pc <- get.pc(x.CHS.TI.maf0.05,outliers.CHS.TI) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CHS.TI <- read.table("CHS.TI.148.plink.map", header=F)
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

CZ <- read.pcadapt("CZ.404.9608.plink.ped", type="ped")
Summary:

        - input file      CZ.404.9608.plink.ped
        - output file     CZ.404.9608.plink.pcadapt

	- number of individuals detected:	404
	- number of loci detected:		9608

File has been sucessfully converted.

##Check the nr of PCs

x.CZ <- pcadapt(CZ, K=20)

Reading file CZ.404.9608.plink.pcadapt...
Number of SNPs: 9608
Number of individuals: 404
Number of SNPs with minor allele frequency lower than 0.05 ignored: 1493
988533 out of 3881632 missing data ignored.

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
qval <- qvalue(x.CZ$pvalues)$qvalues
outliers.CZ <- which(qval<alpha)
outliers.CZ
snp_pc <- get.pc(x.CZ.maf0.05,outliers.CZ) ##see PCs associated with the outliers
```


### 3.Rename loci

```
locus.names.CZ <- read.table("CZ.plink.map", header=F)
locus.names.CZ
locus.names.CZ$ID <- seq.int(nrow(locus.names.CZ)) #add an index of the SNP numbers
CZ.outliers.pcadapt <- as.character(outliers.CZ)
CZ.outliers.pcadapt.names <- locus.names.CZ[locus.names.CZ$ID %in% CZ.outliers.pcadapt,]
CZ.outliers.pcadapt.names <- paste("X", CZ.outliers.pcadapt.names$V2, sep=".")

write.table(CZ.outliers.pcadapt.names, "CZ.pcadapt.outliers", col.names=F, row.names=F, quote=F)  ##write the table #write to file
```

