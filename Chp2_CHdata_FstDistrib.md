#CH Data: Fst distributions

1. Neutral Fst distribution 

2. Neutral + RDA

3. Neutral + XtX

4. Neutral + all outliers

5. Neutral + LFMM

6. Neutral + BayEnv2 

7. Neutral + EAA loci

Total: 42 distributions

###Input files

This is to find loci that explain differentiation between populations more than neutral expectation
```
#Convert plink .ped files to structure using pgdspider

/Users/alexjvr/Applications/PGDSpider_2.1.0.3/

java -Xmx1024m -Xms512m -jar PGDSpider2.jar
```

create population files using the *nosex files. 

copy everything to here: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats




###Neutral Frequency distributions of Fst per locus

1. CHall
```
#R
library(adegenet)
library(hierfstat)
CHall940.genind <- read.structure("CHall.940.9608.str")

CHall940.genind@pop
CHall940.pop <- read.table("CHall940.pop", header=T)
CHall940.pop.factor <- as.factor(CHall940.pop$pop)
CHall940.genind@pop <- CHall940.pop.factor
CHall940.genind@pop

hier.CHall940 <- genind2hierfstat(CHall940.genind)
stats.hier.CHall940 <- basic.stats(hier.CHall940)
stats.hier.CHall940.perlocus <- stats.hier.CHall940$perloc

pdf(file="CHall940.NeutralFst.pdf")
hist(stats.hier.CHall940.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

CHall940.neutral Fst

![alt_txt][CH940.Fst]
[CH940.Fst]:https://cloud.githubusercontent.com/assets/12142475/23822706/939bf45e-0652-11e7-93cc-25a8ec1ed89e.png


2. CHN
```
#R
library(adegenet)
library(hierfstat)
CHN229.genind <- read.structure("CHN.229.9608.str")

CHN229.genind@pop
CHN229.pop <- read.table("CHN229.pop", header=T)
CHN229.pop.factor <- as.factor(CHN229.pop$pop)
CHN229.genind@pop <- CHN229.pop.factor
CHN229.genind@pop

hier.CHN229 <- genind2hierfstat(CHN229.genind)
stats.hier.CHN229 <- basic.stats(hier.CHN229)
stats.hier.CHN229.perlocus <- stats.hier.CHN229$perloc

pdf(file="CHN229.NeutralFst.pdf")
hist(stats.hier.CHN229.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

![alt_txt][CHN.Fst]
[CHN.Fst]:https://cloud.githubusercontent.com/assets/12142475/23822918/191f2318-0657-11e7-9baf-1d8c660f0192.png


3. CHS
```
#R
library(adegenet)
library(hierfstat)
CHS283.genind <- read.structure("CHS283.9608.str")

CHS283.genind@pop
CHS283.pop <- read.table("CHS283.pop", header=T)
CHS283.pop.factor <- as.factor(CHS283.pop$pop)
CHS283.genind@pop <- CHS283.pop.factor
CHS283.genind@pop

hier.CHS283 <- genind2hierfstat(CHS283.genind)
stats.hier.CHS283 <- basic.stats(hier.CHS283)
stats.hier.CHS283.perlocus <- stats.hier.CHS283$perloc

pdf(file="CHS283.NeutralFst.pdf")
hist(stats.hier.CHS283.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

![alt_txt][CHS.Fst]
[CHS.Fst]:https://cloud.githubusercontent.com/assets/12142475/23823039/ba436b4e-0659-11e7-859f-e4ede291522d.png



4. CHS.TI
```
#R
library(adegenet)
library(hierfstat)
CHS.TI148.genind <- read.structure("CHS.TI148.9608.str")

CHS.TI148.genind@pop
CHS.TI148.pop <- read.table("CHS.TI148.pop", header=T)
CHS.TI148.pop.factor <- as.factor(CHS.TI148.pop$pop)
CHS.TI148.genind@pop <- CHS.TI148.pop.factor
CHS.TI148.genind@pop

hier.CHS.TI148 <- genind2hierfstat(CHS.TI148.genind)
stats.hier.CHS.TI148 <- basic.stats(hier.CHS.TI148)
stats.hier.CHS.TI148.perlocus <- stats.hier.CHS.TI148$perloc

pdf(file="CHS.TI148.NeutralFst.pdf")
hist(stats.hier.CHS.TI148.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

![alt_txt][CHS.TI.Fst]
[CHS.TI.Fst]:


5. CHS.VS
```
#R
library(adegenet)
library(hierfstat)
CHS.VS135.genind <- read.structure("CHS.VS135.9608.str")

CHS.VS135.genind@pop
CHS.VS135.pop <- read.table("CHS.VS135.pop", header=T)
CHS.VS135.pop.factor <- as.factor(CHS.VS135.pop$pop)
CHS.VS135.genind@pop <- CHS.VS135.pop.factor
CHS.VS135.genind@pop

hier.CHS.VS135 <- genind2hierfstat(CHS.VS135.genind)
stats.hier.CHS.VS135 <- basic.stats(hier.CHS.VS135)
stats.hier.CHS.VS135.perlocus <- stats.hier.CHS.VS135$perloc

pdf(file="CHS.VS135.NeutralFst.pdf")
hist(stats.hier.CHS.VS135.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

![alt_txt][CHS.VS.Fst]
[CHS.VS.Fst]:



6. CZ
```
#R
library(adegenet)
library(hierfstat)
CZ404.genind <- read.structure("CZ404.9608.str")

CZ404.genind@pop
CZ404.pop <- read.table("CZ404.pop", header=T)
CZ404.pop.factor <- as.factor(CZ404.pop$pop)
CZ404.genind@pop <- CZ404.pop.factor
CZ404.genind@pop

hier.CZ404 <- genind2hierfstat(CZ404.genind)
stats.hier.CZ404 <- basic.stats(hier.CZ404)
stats.hier.CZ404.perlocus <- stats.hier.CZ404$perloc

pdf(file="CZ404.NeutralFst.pdf")
hist(stats.hier.CHall940.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
dev.off()
```

![alt_txt][CZ.Fst]
[CZ.Fst]:

