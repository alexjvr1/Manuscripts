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
```
#R
library(adegenet)
library(hierfstat)
SE132.genind <- read.structure("SE132.plink.str")
hier.SE132 <- genind2hierfstat(SE132.genind)
SE132.genind@pop
SE132.pop <- read.table("SE132.pop", header=T)
SE132.pop.factor <- as.factor(SE132.pop$pop)
SE132.genind@pop <- SE132.pop.factor
SE132.genind@pop
hier.SE132 <- genind2hierfstat(SE132.genind)
stats.SE132 <- basic.stats(hier.SE132)
stats.SE132.perlocus <- stats.SE132$perloc

hist(stats.SE132.perlocus$Fst, xlim=c(-0.2, 1.0), breaks=120)
```
