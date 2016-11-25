#SE Full Data Filtering

I'm getting strange results with the 1Million reads/sample subset of the SE dataset. The SFS and LD (r2) are unexpected. 

I am comparing those results to the Full Dataset here. 


Step 1. Filter loci with missingness >20%

    - SFS and LD overall for the full dataset

Step 2. Filter MAF 0.05 

    - SFS and LD overall for the full dataset

Step 3. Max Het 0.5 for the full dataset

    - SFS and LD for full dataset

Step 4. 
  
  - Split into a)populations and b) regions

  - Filter loci with max Het 0.5 and deviating from HWE in > 3 pops / 3 regions

   - SFS and LD

Step 5.1. Remove indivs with >0.5 missing data & remove populations with <5 individuals

  - SFS and LD

Step 5.2. Filter MAF 0.05 overall

Step 6. removed fixed loci per pop/region

  - SFS and LD



##Step 1: Filter for missingness across loci

Remove loci genotyped in <80% individuals. 


```
gdcsrv1.ethz.ch:/gdc_home4/alexjvr/SEFinalSamples/outfiles

vcftools --vcf SEFinalc94d6m4p3.vcf --max-missing 0.8 --recode --recode-INFO-all --out SE.Full_Filtered_20161118/SE.s1.maxmiss0.8

After filtering, kept 15670 out of a possible 1300946 Sites
```

Rename and plot missingness per pop using R
```
bcftools reheader SE.s1.maxmiss0.8.recode.vcf -s newnamesforvcf.s1 -o SE.s1.maxmiss0.8.names.vcf
vcftools --vcf SE.s1.maxmiss0.8.names.vcf --recode --recode-INFO-all --out SE.s1.maxmiss0.8.names

vcftools --vcf SE.s1.maxmiss0.8.names.recode.vcf --missing-indv
 
mawk '!/IN/' out.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


##R
library(ggplot2)
SE.Full.s1 <- read.table("out.imiss", header=T)
pop <- read.table("SE.193.pops", header=F)
SE.Full.s1$pop <- pop$V2
SE.Full.s1$pop.order <- pop$V3

SE.Full.s1.sort <- SE.Full.s1[order(SE.Full.s1$pop.order),]

SE.Full.s1.sort$pop <- factor(SE.Full.s1.sort$pop, levels=SE.Full.s1.sort$pop)   ##sort pop.nr. Numbers from south to North

qplot(pop, F_MISS, data=SE.Full.s1.sort, geom=c("boxplot", "jitter"))
```

![alt_txt][Fig1.s1]
[Fig1.s1]:https://cloud.githubusercontent.com/assets/12142475/20427224/300d40c8-ad83-11e6-88d3-fb0c1ea68ee9.png


![alt_txt][Fig2.s1]
[Fig2.s1]:https://cloud.githubusercontent.com/assets/12142475/20428010/3105762c-ad87-11e6-82c7-f53c1988d402.png


SFS and LD (r2) for the full dataset: 
```
##linux

vcftools --vcf SE.s1.maxmiss0.8.names.recode.vcf --plink --out SE.s1.plink
plink --file SE.s1.plink --recode --recodeA

plink --file SE.s1.plink --freq --out SE.s1
plink --file SE.s1.plink --r2 --out SE.s1

###R

SE.s1.freq <- read.table("SE.s1.frq", header=T)
hist(SE.s1.freq, main="SE s1 (max.miss0.8; 15294; 193) SFS")

SE.s1.r2 <- read.table("SE.s1.ld", header=T)
hist(SE.s1.r2$R2, main="SE s1 (max.miss0.8; 15294; 193) R2")
```

![alt_txt][SFS.s1]
[SFS.s1]:https://cloud.githubusercontent.com/assets/12142475/20428334/c9cbff2e-ad88-11e6-8793-8247b7530f27.png


![alt_txt][R2.s1]
[R2.s1]:https://cloud.githubusercontent.com/assets/12142475/20428410/452de7ae-ad89-11e6-9acf-0f8773e34525.png



##Step 2-4. Filter MAF 0.05, thin, max-missing 0.25 

```
vcftools --vcf SE.s1.maxmiss0.8.names.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SE.s2.0.8.maf0.05
 vcftools --vcf SE.s2.0.8.maf0.05.recode.vcf --plink --out SE.s2.plink
plink --file SE.s2.plink --recode --recodeA
plink --file SE.s2.plink --freq --out SE.s2

vcftools --vcf SE.s2.0.8.maf0.05.recode.vcf --thin 125 --recode --recode-INFO-all --out SE.s3.0.8.maf0.05.thin 
vcftools --vcf SE.s3.0.8.maf0.05.thin.recode.vcf --plink --out SE.s3.plink
plink --file SE.s3.plink --recode --recodeA 
plink --file SE.s3.plink --freq --out SE.s3

vcftools --vcf SE.s3.0.8.maf0.05.thin.recode.vcf --remove lowDP.indiv --recode --recode-INFO-all --out SE.s4.171
vcftools --vcf SE.s4.171.recode.vcf --plink --out SE.s4.plink
plink --file SE.s4.plink --recode --recodeA
plink --file SE.s4.plink --freq --out SE.s4


vcftools --vcf SE.s4.171.recode.vcf --remove lowDP.indiv --recode --recode-INFO-all --out SE.s4.167
vcftools --vcf SE.s4.167.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SE.s5.167.maf0.05

vcftools --vcf SE.s5.167.maf0.05.recode.vcf --plink --out SE.s5.plink
plink --file SE.s5.plink --recode --recodeA
plink --file SE.s5.plink --freq --out SE.s5


##Remove loci identified as out of HWE in 5_Chp5_SE.LandscapeGenomics.md, but these loci have already been removed with the other filters

plink --file SE.s5.plink --exclude HWE.loci.remove.names --recode --recodeA --out SE.s6.plink

```


##Removing Sk_Ho from the analyses

It looks like Sk_Ho is problematic: there is more variation in this population compared with the others. 
All of the individuals in this population are based on concatenated datasets: i.e. data from 2 or more runs that have been added together. I think there might be a problem with how the data were concatenated (perhaps individuals were not correctly combined?) or one of the plates (I think perhaps H22 is problematic - based on the misassignment of CH individuals from this plate) has mislabeled samples. 

These data will have to be re-analysed. But I will do this at a later stage. For now I will remove these samples from the analysis to finish this chapter. 

##Removing Sk.Ho from analyses

![alt_txt][SFS.Sk.pops]
[SFS.Sk.pops]:https://cloud.githubusercontent.com/assets/12142475/20597838/46e279c2-b246-11e6-8ee8-752cc48daa91.png


##remove Sk.Ho from dataset
```
vcftools --vcf SE.s4.171.recode.vcf --remove lowDP.indiv --recode --recode-INFO-all --out SE.s4.148
vcftools --vcf SE.s4.148.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SE.s5.148.maf0.05
vcftools --vcf SE.s5.148.maf0.05.recode.vcf --plink --out SE.s5.148.plink
plink --file SE.s5.148.plink --recode --recodeA
plink --file SE.s5.148.plink --exclude HWE.loci.remove.names --recode --recodeA --out SE.s6.148.plink 
plink --file SE.s6.148.plink --freq --out SE.s6.148

```

Final Dataset:

2118 SNPs (plink; 2134 in vcf)

148 individuals

93.3% genotyping rate

```
plink --file SE.s5.148.plink --freq --out SE.s5.148

##R

SE.s1.freq <- read.table("SE.s1.frq", header=T)
SE.s2.freq <- read.table("SE.s2.frq", header=T)
SE.s3.freq <- read.table("SE.s3.frq", header=T)
SE.s4.freq <- read.table("SE.s4.frq", header=T)
SE.s5.148.freq <- read.table("SE.s5.148.frq", header=T)


my.bin.width <- 0.05

par(mfrow=c(3,2))
hist(SE.s1.freq$MAF, main="SE s1 (max.miss0.8; 15294; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s2.freq$MAF, main="SE s2 (max.miss0.8; maf 0.05; 4519; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s3.freq$MAF, main="SE s3 (max.miss0.8; maf 0.05; thin; 2199; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s4.freq$MAF, main="SE s4 (max.miss0.8; maf 0.05; thin; 2199; 171) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s5.148.freq$MAF, main="SE s5 (max.miss0.8; maf 0.05; thin; 2118; 148) SFS", breaks=seq(0,0.5, by=my.bin.width))
```

![alt_txt][SFS.Final.148]
[SFS.Final.148]:https://cloud.githubusercontent.com/assets/12142475/20598485/93144f2a-b249-11e6-9b6a-0a5edbf95ee4.png


```
##R

SE.s1.freq <- read.table("SE.s1.frq", header=T)
SE.s2.freq <- read.table("SE.s2.frq", header=T)
SE.s3.freq <- read.table("SE.s3.frq", header=T)
SE.s4.freq <- read.table("SE.s4.frq", header=T)
SE.s5.freq <- read.table("SE.s5.frq", header=T)
SE.s6.freq <- read.table("SE.s6.frq", header=T)

my.bin.width <- 0.05

par(mfrow=c(3,2))
hist(SE.s1.freq$MAF, main="SE s1 (max.miss0.8; 15294; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s2.freq$MAF, main="SE s2 (max.miss0.8; maf 0.05; 4519; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s3.freq$MAF, main="SE s3 (max.miss0.8; maf 0.05; thin; 2199; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s4.freq$MAF, main="SE s4 (max.miss0.8; maf 0.05; thin; 2199; 171) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s5.freq$MAF, main="SE s4 (max.miss0.8; maf 0.05; thin; 2167; 167) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s6.freq$MAF, main="SE s4 (max.miss0.8; maf 0.05; thin; 2167; 167) SFS", breaks=seq(0,0.5, by=my.bin.width))


```


![alt_txt][SFS.s1-6]
[SFS.s1-4]:https://cloud.githubusercontent.com/assets/12142475/20595904/36f88ede-b23c-11e6-9f04-9112c431147a.png


Missingness across individuals
```
##s4 with LT2, Sk_Ho, & UT3 still included
##s6 with final dataset (148indivs) 

alexjvr$ vcftools --vcf SE.s4.171.recode.vcf --missing-indv --out SE.s4
 
##R
library(ggplot2)
SE.s4.miss <- read.table("SE.s4.imiss", header=T)
pop.s4 <- read.table("SE.s4.popnames", header=F)
SE.s4.miss$pop <- pop.s4$V1
SE.s4.miss$pop.order <- pop.s4$V2

SE.s4.miss.sort <- SE.s4.miss[order(SE.s4.miss$pop.order),]

SE.s4.miss.sort$pop <- factor(SE.s4.miss.sort$pop, levels=SE.s4.miss.sort$pop)




vcftools --vcf SE.s5.148.maf0.05.recode.vcf --missing-indv --out SE.s6
 
##R
library(ggplot2)
SE.s6.miss <- read.table("SE.s6.imiss", header=T)
pop.s6 <- read.table("SE.s6.popnames", header=F)
SE.s6.miss$pop <- pop.s6$V2
SE.s6.miss$pop.order <- pop.s6$V1

SE.s6.miss.sort <- SE.s6.miss[order(SE.s6.miss$pop.order),]

SE.s6.miss.sort$pop <- factor(SE.s6.miss.sort$pop, levels=SE.s6.miss.sort$pop)

qplot(pop, F_MISS, data=SE.s6.miss.sort, geom=c("boxplot", "jitter"))
```

![alt_txt][miss.s4]
[miss.s4]:https://cloud.githubusercontent.com/assets/12142475/20431948/394aada0-ad9d-11e6-9162-12229591b63e.png

![alt_txt][miss.s6]
[miss.s6]:https://cloud.githubusercontent.com/assets/12142475/20599639/7968fe1c-b24f-11e6-8f66-fbbf2776ee2b.png


##Subset data

into populations and regions


/Users/alexjvr/2016RADAnalysis/5_SE.MS1/SE.FullData.Analyses/sumstats/DataFiltering/SE.Full_Filtered_20161118
```
mkdir subset.data
for i in $(ls popnames.plink.folder/); do plink --file SE.s4.plink --keep popnames.plink.folder/$i --recode --recodeA --out subset.data/$i.plink; done

mkdir subset.data.s6
for i in $(ls popnames.plink.folder/); do plink --file SE.s5.148.plink --keep popnames.plink.folder/$i --recode --recodeA --out subset.data.s6/$i.plink; done
```

and calculate SFS and LD for each 

```
##s4
for i in $(ls popnames.plink.folder/); do plink --file subset.data/$i.plink --freq --out subset.data/$i; done
for i in $(ls popnames.plink.folder/); do plink --file subset.data/$i.plink --r2 --out subset.data/$i; done

##R

DE.3pop.frq <- read.table("subset.data/DE.3pop.frq", header=T)
Sk.3pop.frq <- read.table("subset.data/Sk.3pop.frq", header=T)
Upp.3pop.frq <- read.table("subset.data/Upp.3pop.frq", header=T)
Umea.3pop.frq <- read.table("subset.data/Umea.3pop.frq", header=T)
Lulea.3pop.frq <- read.table("subset.data/Lulea.3pop.frq", header=T)
Kir.2pop.frq <- read.table("subset.data/Kir.2pop.frq", header=T)
FIN.frq <- read.table("subset.data/FIN.frq", header=T)

par(mfrow=c(2, 4))
hist(DE.3pop.frq$MAF, main="DE.3pop (29, 2199, 90%) SFS")
hist(Sk.3pop.frq$MAF, main="Sk.3pop (48, 2199, 90%) SFS")
hist(Upp.3pop.frq$MAF, main="Upp.3pop (28, 2199, 96%) SFS")
hist(Umea.3pop.frq$MAF, main="Umea.3pop (19, 2199, 93%) SFS")
hist(Lulea.3pop.frq$MAF, main="Lulea.3pop (19, 2199, 97%) SFS")
hist(Kir.2pop.frq$MAF, main="Kir.2pop (20, 2199, 95%) SFS")
hist(FIN.frq$MAF, main="FIN (7, 2199, 88%) SFS")

##s6
for i in $(ls popnames.plink.folder/); do plink --file subset.data.s6/$i.plink --freq --out subset.data.s6/$i; done
for i in $(ls popnames.plink.folder/); do plink --file subset.data.s6/$i.plink --r2 --out subset.data.s6/$i; done

##R

DE.3pop.frq <- read.table("subset.data.s6/DE.3pop.frq", header=T)
Sk.2pop.frq <- read.table("subset.data.s6/Sk.3pop.frq", header=T)
Upp.3pop.frq <- read.table("subset.data.s6/Upp.3pop.frq", header=T)
Umea.2pop.frq <- read.table("subset.data.s6/Umea.3pop.frq", header=T)
Lulea.2pop.frq <- read.table("subset.data.s6/Lulea.3pop.frq", header=T)
Kir.2pop.frq <- read.table("subset.data.s6/Kir.2pop.frq", header=T)
FIN.frq <- read.table("subset.data.s6/FIN.frq", header=T)


par(mfrow=c(2, 4))
hist(DE.3pop.frq$MAF, main="DE.3pop (29, 2167, 90%) SFS")
hist(Sk.2pop.frq$MAF, main="Sk.2pop (29, 2167, 93%) SFS")
hist(Upp.3pop.frq$MAF, main="Upp.3pop (28, 2167, 96%) SFS")
hist(Umea.2pop.frq$MAF, main="Umea.2pop (18, 2167, 94%) SFS")
hist(Lulea.2pop.frq$MAF, main="Lulea.2pop (16, 2167, 97%) SFS")
hist(Kir.2pop.frq$MAF, main="Kir.2pop (20, 2167, 94%) SFS")
hist(FIN.frq$MAF, main="FIN (7, 2167, 88%) SFS")




#And with only variable loci: 

DE.3pop.frq.var.new <- subset(DE.3pop.frq, MAF>0.001)
Sk.2pop.frq.var.new <- subset(Sk.2pop.frq, MAF>0.001)
Upp.3pop.frq.var.new <- subset(Upp.3pop.frq, MAF>0.001)
Umea.2pop.frq.var.new <- subset(Umea.2pop.frq, MAF>0.001)
Lulea.2pop.frq.var.new <- subset(Lulea.2pop.frq, MAF>0.001)
Kir.2pop.frq.var.new <- subset(Kir.2pop.frq, MAF>0.001)
FIN.frq.var.new <- subset(FIN.frq, MAF>0.1)  ##I think because of missing data, I have to filter for a MAC of 2 in this dataset

my.bin.width <- 0.01

par(mfrow=c(2, 4))
hist(DE.3pop.frq.var$MAF, main="DE.3pop.var (29, 1673, 22.7% fixed) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Sk.2pop.frq.var$MAF, main="Sk.3pop (29, 1760, 18.8% fixed) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Upp.3pop.frq.var$MAF, main="Upp.3pop (28, 1256, 42.0%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Umea.2pop.frq.var$MAF, main="Umea.3pop (18, 1129, 47.9%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Lulea.2pop.frq.var$MAF, main="Lulea.3pop (16, 1062, 51.0%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Kir.2pop.frq.var$MAF, main="Kir.2pop (20, 1187, 45.2%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(FIN.frq.var$MAF, main="FIN (7, 739, 65.8%) SFS", breaks=seq(0,0.5, by=my.bin.width))


##Loci variable in x populations
SE.var.loci.freq <- do.call(rbind, lapply(ls(pattern=".frq.var.new$"), get))
summary(SE.var.loci.freq)
SE.var.loci.freq.keep <- data.frame(table(SE.var.loci.freq$SNP)) 
#SE.var.loci.freq.keep <- subset(SE.var.loci.freq.keep, Freq>0) 
hist(SE.var.loci.freq.keep$Freq, xlab="Nr pops", ylab="Frequency", main="Frequency of variable loci across regions (s6)", breaks=seq(0.5,7.5, by=1.0))


###Looking only at fixed loci

#And without the fixed loci: (Frequency here is based on 1/(nx2) for MAC >1)

DE.3.pop.frq.fix.new <- subset(DE.3pop.frq, MAF<0.0018)
Sk.2pop.frq.fix.new <- subset(Sk.2pop.frq, MAF<0.0018)
Upp.3pop.frq.fix.new <- subset(Upp.3pop.frq, MAF<0.00179)
Umea.3pop.frq.fix.new <- subset(Umea.2pop.frq, MAF<0.0027)
Lulea.3pop.frq.fix.new <- subset(Lulea.2pop.frq, MAF<0.0031)
Kir.2pop.frq.fix.new <- subset(Kir.2pop.frq, MAF<0.0026)
FIN.frq.fix.new <- subset(FIN.frq, MAF<0.070)



SE.region.frq.fix.table <- do.call(rbind, lapply(ls(pattern="fix.new$"), get))

SE.fix.region.table.keep <- data.frame(table(SE.region.frq.fix.table$SNP))
my.bin.width=1
hist(SE.fix.region.table.keep$Freq, xlab="Number of pops", ylab="Frequency", main="Frequency of fixed loci found in increasing number of populations", breaks=seq(0,7, by=my.bin.width))
```


##SFS.s4.regions
![alt_txt][SFS.s4.regions]
[SFS.s4.regions]:https://cloud.githubusercontent.com/assets/12142475/20433442/df6d3152-ada3-11e6-9a02-4123ad1ff7f4.png


##SFS.s6.regions
![alt_txt][SFS.s6.regions]
[SFS.s6.regions]:https://cloud.githubusercontent.com/assets/12142475/20599759/32c752e6-b250-11e6-92a7-61b9f659ddf5.png


##SFS.s4.variable only
![alt_txt][var.only]
[var.only]:https://cloud.githubusercontent.com/assets/12142475/20435155/4d959fb4-adab-11e6-99d7-cc0fd54519f2.png

##SFS.s6.variable only
![alt_txt][var.only.s6]
[var.only.s6]:https://cloud.githubusercontent.com/assets/12142475/20599880/e3c6013c-b250-11e6-9d94-b3108f2b0269.png



##s4 Frequency fixed loci
![alt_txt][freq.fixed]
[freq.fixed]:https://cloud.githubusercontent.com/assets/12142475/20435161/56ab279a-adab-11e6-8ad5-2a512cb2c689.png


##s4 Frequency variable loci
![alt_txt][freq.var]
[freq.var]:https://cloud.githubusercontent.com/assets/12142475/20495176/452bf2ac-b020-11e6-87bb-ec23ccb29357.png

##s4 Variable sites only
![alt_txt][SFS.var.allregions]
[SFS.var.allregions]:https://cloud.githubusercontent.com/assets/12142475/20478502/ad6c623e-afda-11e6-8607-c630767538ba.png

##s4 Adjusted bin width increments of ~MAC 1 
![alt_txt][SFS.var.allregions.adjbinwidth]
[SFS.var.allregions.adjbinwidth]:https://cloud.githubusercontent.com/assets/12142475/20478501/ad4d9f84-afda-11e6-9e3e-c936431fac6e.png



##s6 Frequency fixed loci
![alt_txt][freq.fixed.s6]
[freq.fixed.s6]:https://cloud.githubusercontent.com/assets/12142475/20600482/25939de2-b254-11e6-9a58-949041539ccd.png


##s6 Frequency variable loci
![alt_txt][freq.var.s6]
[freq.var.s6]:https://cloud.githubusercontent.com/assets/12142475/20600481/25935440-b254-11e6-8840-2c3640e46715.png






##R2 LD

Checking for linkage 
```
##R
write.table(DE.3pop.frq.fix.new$SNP, "DE.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Sk.2pop.frq.fix.new$SNP, "Sk.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Upp.3pop.frq.fix.new$SNP, "Upp.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Umea.3pop.frq.fix.new$SNP, "Umea.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Lulea.3pop.frq.fix.new$SNP, "Lulea.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Kir.2pop.frq.fix.new$SNP, "Kir.2pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(FIN.frq.fix.new$SNP, "FIN.fixedloci", row.names=F, col.names=F, quote=F)


##linux
plink --file subset.data.s6/DE.3pop.plink --exclude DE.3pop.fixedloci --recode --recodeA --out subset.data.s6/DE.3pop.var.plink
plink --file subset.data.s6/DE.3pop.var.plink --r2 --out subset.data.s6/DE.3pop.var
plink --file subset.data.s6/Sk.3pop.plink --exclude Sk.3pop.fixedloci --recode --recodeA --out subset.data.s6/Sk.3pop.var.plink
plink --file subset.data.s6/Sk.3pop.var.plink --r2 --out subset.data.s6/Sk.3pop.var
plink --file subset.data.s6/Upp.3pop.plink --exclude Upp.3pop.fixedloci --recode --recodeA --out subset.data.s6/Upp.3pop.var.plink
plink --file subset.data.s6/Upp.3pop.var.plink --r2 --out subset.data.s6/Upp.3pop.var
plink --file subset.data.s6/Umea.3pop.plink --exclude Umea.3pop.fixedloci --recode --recodeA --out subset.data.s6/Um.3pop.var.plink
plink --file subset.data.s6/Um.3pop.var.plink --r2 --out subset.data.s6/Um.3pop.var
plink --file subset.data.s6/Lulea.3pop.plink --exclude Lulea.3pop.fixedloci --recode --recodeA --out subset.data.s6/Lulea.3pop.var.plink
plink --file subset.data.s6/Lulea.3pop.var.plink --r2 --out subset.data.s6/Lulea.3pop.var
plink --file subset.data.s6/Kir.2pop.plink --exclude Kir.2pop.fixedloci --recode --recodeA --out subset.data.s6/Kir.2pop.var.plink
plink --file subset.data.s6/Kir.2pop.var.plink --r2 --out subset.data.s6/Kir.2pop.var
plink --file subset.data.s6/FIN.plink --exclude FIN.fixedloci --recode --recodeA --out subset.data.s6/FIN.var.plink
plink --file subset.data.s6/FIN.var.plink --r2 --out subset.data.s6/FIN.var




##R
DE.3pop.var.ld <- read.table("subset.data.s6/DE.3pop.var.ld", header=T)
Sk.3pop.var.ld <- read.table("subset.data.s6/Sk.3pop.var.ld", header=T)
Upp.3pop.var.ld <- read.table("subset.data.s6/Upp.3pop.var.ld", header=T)
Lulea.3pop.var.ld <- read.table("subset.data.s6/Lulea.3pop.var.ld", header=T)
Um.3pop.var.ld <- read.table("subset.data.s6/Um.3pop.var.ld", header=T)
Kir.2pop.var.ld <- read.table("subset.data.s6/Kir.2pop.var.ld", header=T)
FIN.var.ld <- read.table("subset.data.s6/FIN.var.ld", header=T)

my.bin.width=0.05

par(mfrow=c(2,4))
hist(DE.3pop.var.ld$R2, main="DE.3pop.var (29, 1623 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Sk.3pop.var.ld$R2, main="Sk.2pop.var (29, 1760 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Upp.3pop.var.ld$R2, main="Upp.3pop (28, 1255 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Um.3pop.var.ld$R2, main="Umea.2pop (18, 1127 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Lulea.3pop.var.ld$R2, main="Lulea.2pop (16, 1065 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Kir.2pop.var.ld$R2,main="Kir.2pop (20, 1190 loci) R2",  breaks=seq(0,1.0, by=my.bin.width))
hist(FIN.var.ld$R2, main="FIN (7, 897 loci) R2", breaks=seq(0,1.0, by=my.bin.width))


##Keep only loci R2>0.8 and plot frequency in R
DE.3pop.ld.remove <- subset(DE.3pop.var.ld, R2>0.799999)
Sk.3pop.ld.remove <- subset(Sk.3pop.var.ld, R2>0.799999)
Upp.3pop.ld.remove <- subset(Upp.3pop.var.ld, R2>0.799999)
Umea.3pop.ld.remove <- subset(Um.3pop.var.ld, R2>0.799999)
Lulea.3pop.ld.remove <- subset(Lulea.3pop.var.ld, R2>0.799999)
Kir.2pop.ld.remove <- subset(Kir.2pop.var.ld, R2>0.799999)
FIN.ld.remove <- subset(FIN.var.ld, R2>0.799999)

SE.ld.remove.file <- do.call(rbind, lapply(ls(pattern="ld.remove$"), get))
test.table.ld <- data.frame(table(SE.ld.remove.file$SNP_A, SE.ld.remove.file$SNP_B))
test.table.ld.subset <- subset(test.table.ld, Freq>0)
hist(test.table.ld.subset$Freq, main="Frequency of loci R2>0.8 across 7 regions")
```

R2 of variable loci only

Starts at 0.2, because: 

Finland still has high R2, because: 

![alt_txt][R2.variableLoci]
[R2.variableLoci]:https://cloud.githubusercontent.com/assets/12142475/20603559/5b809c3a-b262-11e6-9a90-d9c0da234060.png



Frequency of loci R2>0.8
![alt_txt][R2.3pop]
[R2.3pop]:https://cloud.githubusercontent.com/assets/12142475/20603558/5b804d16-b262-11e6-984e-54cd5a71c70a.png



##Population structure

###Fst

Convert plink s5 to structure using pgdspider, and import into R
```
/Users/alexjvr/2016RADAnalysis/5_SE.MS1/SE.FullData.Analyses/pop.str

library(adegenet)
library(hierfstat)
library(reshape)

SEall.171 <- read.structure("SE.s5.148.15pops.str")
SEall.171

pop.SE171 <- read.table("popnames.s5.structure", header=T)  ##make sure the populations are numbered "01.DE.B", etc.
popnames.15.SEall.factor <- as.factor(pop.SE171$pop)
SEall.171@pop <- popnames.15.SEall.factor

hier.SEall <- genind2hierfstat(SEall.171)

SEall.fst <- pairwise.fst(SEall.171, pop=NULL, res.type=c("dist", "matrix"))

m <- as.matrix(SEall.fst)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2)<- c("c1","c2", "distance")

library(gplots)

shadesOfGrey <- colorRampPalette(c("grey100", "grey0"))  ##define the colourpalette. 

Dend <- read.table("heatmap.popcolours", header=T)  ##list of colour names for each population based on R colour palatte. In alphabetical order (as in genind file)
Dend.Colours <- as.character(Dend$colours.pop)

par(oma=c(1,1,2,1))
heatmap.2(as.matrix(SEall.fst), trace="none", RowSideColors=Dend.Colours, ColSideColors=Dend.Colours, col=shadesOfGrey, labRow=F, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst and dendrograms of SEall: 15pops, 7regions, 2167loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input. The pop order is alphabetical in the output. 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

popnames.all <- as.character(Dend$pop)
legend("bottom", popnames.all, xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="o", pch=15, col=Dend.Colours, title="Side Dendrogram:Region")



##Fst between regions
##pop levels are always in alphabetical order

SEall.171.regions <- SEall.171
popnames.7.SEall.factor <- as.factor(pop.SE171$region)
SEall.171.regions@pop <- popnames.7.SEall.factor

SEall.region.Fst <- pairwise.fst(SEall.171.regions, pop=NULL, res.type=c("dist", "matrix"))

#Dend3 <- read.table("heatmap.regioncolours", header=T)  ##list of colour names for each population based on R colour palatte. In #alphabetical order (as in genind file)
#Dend.Colours3 <- as.character(Dend3$colours.region)

Dend3 <- read.table("heatmap.regions.nonumber", header=T)  ##list of colour names for each population based on R colour palatte. In alphabetical order (as in genind file)
Dend.Colours3 <- as.character(Dend3$colours.region)

par(oma=c(1,1,2,1))
heatmap.2(as.matrix(SE.region.Fst), trace="none", RowSideColors=Dend.Colours3, ColSideColors=Dend.Colours3, col=shadesOfGrey, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst and dendrograms of SEall: 15pops, 7regions, 2167loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input.


par(oma=c(1,1,2,1))
heatmap.2(as.matrix(SE.region.Fst), trace="none", RowSideColors=Dend.Colours3, ColSideColors=Dend.Colours3, col=shadesOfGrey, labRow=F, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst and dendrograms of SEall: 15pops, 7regions, 2167loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input.
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottomleft", c("DE", "SK", "Upp", "Um", "Lul", "Kir", "FIN"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="o", pch=15, col=c("gold1","chartreuse2","darkolivegreen1", "darkslategray1", "deepskyblue1", "dodgerblue2", "dodgerblue4"), title="Side Dendrogram:Region")
```


![alt_txt][Fst.17pops]
[Fst.17pops]:https://cloud.githubusercontent.com/assets/12142475/20605840/a3af1652-b26d-11e6-81a9-fffef37154d4.png

![alt_txt][Fst.7regions]
[Fst.7regions]:https://cloud.githubusercontent.com/assets/12142475/20605871/d05a0996-b26d-11e6-80f3-176b2086edff.png




###IBD

Fst/(1-Fst) vs log.dist(km) - according to Rousset et al. 1997, this is correlated with the effective population density (De) and effective dispersal distance (variance)

```
##Fst 15 pops -> Fst/(1-Fst)

library(reshape)
library(fields)


m <- as.matrix(SEall.fst)
m
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("c1", "c2", "distance")
m2
m2$IBD <- m2$distance/(1-m2$distance)


SE.pop.coords <- read.table("SE.coords.15pop", header=T)
SEpop_lon.lat <- cbind(SE.pop.coords$Long, SE.pop.coords$Lat)
distance.matrix.SEpop <- rdist.earth(SEpop_lon.lat, miles=F)  ##great circle dist based on the coordinates
m.dist <- as.matrix(distance.matrix.SEpop)
summary(m.dist)

m2.dist <- melt(m.dist)[melt(upper.tri(m.dist))$value,]
names(m2.dist) <- c("c1", "c2", "distance")
summary(m2.dist)
m2.dist$log.km <- log(m2.dist$distance)


library(MASS)
#dens <- kde2d(m2$IBD,m2.dist$log.km, n=10)
#myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(m2$IBD~m2.dist$log.km, pch=20,cex=.5, xlab="log Geographic distance (km)", ylab="Fst/(1-Fst)")
#image(dens, col=transp(myPal(10),.7), add=TRUE)
abline(fit <- lm(m2$IBD~m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - SEall")


##Remove Germany

DE.only.m2 <- subset(m2, m2$c1>3)
DE.only.m2.dist <- subset(m2.dist, c1>3)
plot(DE.only.m2$IBD~DE.only.m2.dist$log.km, pch=20,cex=.5, xlab="log Geographic distance (km)", ylab="Fst/(1-Fst)")
abline(fit <- lm(DE.only.m2$IBD~DE.only.m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - SE (no DE)")


###Regions only

m.region <- as.matrix(SE.region.Fst)
m.region
m2.region <- melt(m.region)[melt(upper.tri(m.region))$value,]
names(m2.region) <- c("c1", "c2", "distance")
m2.region
m2.region$IBD <- m2.region$distance/(1-m2.region$distance)


SE.region.coords <- read.table("SE.regions.coords", header=T)
SEregion_lon.lat <- cbind(SE.region.coords$Long, SE.region.coords$Lat)
distance.matrix.SEregion <- rdist.earth(SEregion_lon.lat, miles=F)  ##great circle dist based on the coordinates
m.dist.region <- as.matrix(distance.matrix.SEregion)
summary(m.dist.region)

m2.dist.region <- melt(m.dist.region)[melt(upper.tri(m.dist.region))$value,]
names(m2.dist.region) <- c("c1", "c2", "distance")
summary(m2.dist.region)
m2.dist.region$log.km <- log(m2.dist.region$distance)


library(MASS)
#dens <- kde2d(m2.region$IBD,m2.dist.region$log.km, n=300)
#myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(m2.region$IBD~m2.dist.region$log.km, pch=20,cex=.5, xlab="log Geographic distance (km)", ylab="Fst/(1-Fst)")
#image(dens, col=transp(myPal(10),.7), add=TRUE)
abline(fit <- lm(m2.region$IBD~m2.dist.region$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - SE regions")

```

![alt_txt][IBD.SEall]
[IBD.SEall]:https://cloud.githubusercontent.com/assets/12142475/20606328/7690e44a-b270-11e6-97af-4cc2eb089190.png

![alt_txt][IBD.SEonly]
[IBD.SEonly]:https://cloud.githubusercontent.com/assets/12142475/20606327/7690b54c-b270-11e6-806b-9abe07862590.png


![alt_txt][IBD.region]
[IBD.region]:https://cloud.githubusercontent.com/assets/12142475/20490586/cbf78402-b00e-11e6-9866-bc70d20398db.png



###AMOVA

http://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html

```
library(adegenet)
library(poppr)

SE148.strata <- pop.SE171[,1:3]  ##from text file. each column has one hierarchy level specified for all individuals. (indiv, pop, region)
SEall.148.genind@other <- SE148.strata ##renamed the SEall.171 genind to a better name

strata(SEall.148.genind) <- other(SEall.148.genind)

SEall.148.genclone <- as.genclone(SEall.148.genind)

SE148.amova <- poppr.amova(SEall.148.genclone, ~region/pop)

SE148.amova 
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                            Df    Sum Sq    Mean Sq
Between region               6 14249.577 2374.92948
Between pop Within region    8  1267.817  158.47708
Between samples Within pop 133 10411.288   78.28036
Within samples             148 10492.128   70.89276
Total                      295 36420.810  123.46037

$componentsofcovariance
                                            Sigma          %
Variations  Between region              53.634789  40.539819
Variations  Between pop Within region    4.080152   3.083980
Variations  Between samples Within pop   3.693802   2.791958
Variations  Within samples              70.892758  53.584243
Total variations                       132.301501 100.000000

$statphi
                         Phi
Phi-samples-total 0.46415757
Phi-samples-pop   0.04952370
Phi-pop-region    0.05186631
Phi-region-total  0.40539819

SE148.amovatest <- randtest(SE148.amova, nrepet=999)

SE148.amovatest

class: krandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = SE148.amova, nrepet = 999)

Number of tests:   4 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test       Obs    Std.Obs   Alter Pvalue
1  Variations within samples 70.892758 -16.254415    less  0.001
2 Variations between samples  3.693802   2.715691 greater  0.006
3     Variations between pop  4.080152  16.816844 greater  0.001
4  Variations between region 53.634789   4.478215    less  1.000

```

![alt_txt][AMOVA]
[AMOVA]:https://cloud.githubusercontent.com/assets/12142475/20623482/1040a2be-b308-11e6-92ad-570212235f23.png


i.e. significant structure between regions, but not within regions. 

```
##Randomised test
SE148.new <- SEall.148.genclone
set.seed(9001)
strata(SE148.new) <- strata(SEall.148.genclone)[sample(nInd(SEall.148.genclone)), -1]
head(strata(SE148.new))
head(strata(SEall.148.genclone))
SE148.new.amova <- poppr.amova(SE148.new, ~region/pop)

SE148.new.amova   ##now all the variation is within samples and within populations. So no population structure evident. 
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                            Df     Sum Sq   Mean Sq
Between region               6   909.3634 151.56057
Between pop Within region    8  1445.1514 180.64392
Between samples Within pop 133 23574.1669 177.24938
Within samples             148 10492.1281  70.89276
Total                      295 36420.8098 123.46037

$componentsofcovariance
                                             Sigma           %
Variations  Between region              -0.7034904  -0.5694421
Variations  Between pop Within region    0.1727036   0.1397954
Variations  Between samples Within pop  53.1783088  43.0453201
Variations  Within samples              70.8927576  57.3843267
Total variations                       123.5402795 100.0000000

$statphi
                           Phi
Phi-samples-total  0.426156733
Phi-samples-pop    0.428611685
Phi-pop-region     0.001390038
Phi-region-total  -0.005694421


SE148.new.amova.test<- randtest(SE148.new.amova, nrepet=999) 

SE148.new.amova.test

class: krandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = SE148.new.amova, nrepet = 999)

Number of tests:   4 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test        Obs     Std.Obs   Alter Pvalue
1  Variations within samples 70.8927576 -15.8040911    less  0.001
2 Variations between samples 53.1783088  16.1312400 greater  0.001
3     Variations between pop  0.1727036   0.2045010 greater  0.368
4  Variations between region -0.7034904  -0.5337828    less  0.312
```
Randomised samples

![alt_txt][AMOVA.random]
[AMOVA.random]:https://cloud.githubusercontent.com/assets/12142475/20623863/faefc384-b309-11e6-8e13-44fd4d78190e.png



###DAPC

tutorial-dapc: A tutorial for Discriminant Analysis of Principal Components (DAPC) using adegenet 2.0.0

total variance = (variance between groups) + (variance within groups)

or more simply, denoting X the data matrix:

VAR(X) = B(X) +W(X)

Usual approaches such as Principal Component Analysis (PCA) or Principal Coordinates
Analysis (PCoA / MDS) focus on V AR(X). That is, they only describe the global diversity,
possibly overlooking differences between groups. On the contrary, DAPC optimizes B(X)
while minimizing W(X): it seeks synthetic variables, the discriminant functions, which show
differences between groups as best as possible while minimizing variation within clusters.

```
##1. estimate the number of clusters

Using k-means. Which finds the number of clusters with minimises W(X) and maximises B(X). Compare using BIC

Run algorithm on PCA transformed data. I.e. reduce the dataset so that it can run faster. 

grp.SE148 <- find.clusters(SEall.148.genind, max.n.clust=40)

> choose nr of PCs: 200  ##I try to keep all the PCs

> choose k: 5 ##see figure below


names.15 <- c("DE.B", "DE.K", "DE.W", "Sk.SF", "SK.SL", "Upp.Gra", "Upp.K", "Upp.O", "Um.Gr", "Um.Taf", "LT1", "LT3", "Kir.G", "Kir.L", "FIN")
names.15 <- as.character(names.15)
table.value(table(pop(SEall.148.genind), grp.SE148$grp), col.lab=paste("inf", 1:6), row.lab=names.15)

dapc1.SE148 <- dapc(SEall.148.genind, grp.SE148$grp)
scatter(dapc1.SE148)
```

![alt_txt][DAPC.148]
[DAPC.148]:https://cloud.githubusercontent.com/assets/12142475/20624232/0ea393ae-b30c-11e6-8840-3685e3fd023f.png

![alt_txt][DAPC.148.2]
[DAPC.148.2]:https://cloud.githubusercontent.com/assets/12142475/20625975/b2dbf8ee-b316-11e6-8b64-c0bca7d057eb.png

![alt_txt][DAPC.148.3]
[DAPC.148.3]:https://cloud.githubusercontent.com/assets/12142475/20625944/8461363c-b316-11e6-8c10-5e490c976210.png



###PCA


PCAdapt in R:

```
##convert .vcf to plink 
##linux

vcftools --vcf SE.148.recode.vcf --plink --out SE.148.plink

plink --file SE.148.plink --recode --recodeA

##È
SE.148 <- read.pcadapt("SE.148.plink.ped", type="ped")
Summary:

        - input file      SE.148.plink.ped
        - output file     SE.148.plink.pcadapt

	- number of individuals detected:	148
	- number of loci detected:		2118

File has been sucessfully converted.
```

Check the nr of PCs

```
x <- pcadapt(SE148, K=20)

Reading file SE.148.plink.pcadapt...
Number of SNPs: 2118
Number of individuals: 148
Number of SNPs with minor allele frequency lower than 0.05 ignored: 0
20694 out of 313464 missing data ignored.

plot(x, option="screeplot")  ##PC for pop structure = on the steep curve
```

![alt_txt][pcadapt.148.scree]
[pcadapt.148.scree]:https://cloud.githubusercontent.com/assets/12142475/20627427/399c5d0c-b320-11e6-8492-cb21faea7f8c.png


Plot the PCA using population information
```
popnames <- read.table("popnames", header=T)
head(popnames)
   indiv     pop region region2
1 DE_B01 01.DE.B   1.DE      DE
2 DE_B02 01.DE.B   1.DE      DE
3 DE_B03 01.DE.B   1.DE      DE
4 DE_B04 01.DE.B   1.DE      DE
5 DE_B05 01.DE.B   1.DE      DE
6 DE_B06 01.DE.B   1.DE      DE
poplist <- popnames[,2]

plot(x,option="scores",pop=poplist)

poplist <- as.character(popnames[,3]) ##select regions
plot(x,option="scores",pop=poplist)
```

![alt_txt][pcadapt.pop]
[pcadapt.pop]:https://cloud.githubusercontent.com/assets/12142475/20627810/cd7147c0-b322-11e6-95b4-a61f35f7ec7c.png

![alt_txt][pcadapt.regions]
[pcadapt.regions]:https://cloud.githubusercontent.com/assets/12142475/20627811/cd719086-b322-11e6-99ef-531d948b7f0c.png


###fastStructure



###TESS3



##RDA

Input files: 

1. MAF of all loci 

2. Geographic coordinates

3. Climate variables

2199 loci 

```
###1. MAF

#Calculate MAF for the full dataset within region using PLINK

/Users/alexjvr/2016RADAnalysis/5_SE.MS1/SE.FullData.Analyses/sumstats/DataFiltering/SE.Full_Filtered_20161118

###Use the *nosex file to create a file for subsetting the data. Here I've made 2: 1=within region, 2=population

plink --file /Users/alexjvr/2016RADAnalysis/5_SE.MS1/SE.FullData.Analyses/sumstats/DataFiltering/SE.Full_Filtered_20161118/SE.s4.plink --freq --within SE.Full.s4.region.Cluster --out SE.s4


###R
######Reformat PLINK output
###MAF for each locus -> melt and reformat rows as pops, and columns as loci.

SE.MAF <- read.table("SE.s4.frq.strat", header=T)
SE.MAF2 <- SE.MAF[,c(3,2,6)]
> summary(SE.MAF2)
    CLST             SNP             MAF         
 DE   :2199   100096:2 :    7   Min.   :0.00000  
 FIN  :2199   100865:23:    7   1st Qu.:0.00000  
 Kir  :2199   101108:35:    7   Median :0.08889  
 Lulea:2199   101270:3 :    7   Mean   :0.20265  
 Sk   :2199   101367:84:    7   3rd Qu.:0.32500  
 Umea :2199   101697:24:    7   Max.   :1.00000  
 Upp  :2199   (Other)  :15351                    
 
library("ggplot2")
library("reshape2")

SE.MAF3 <- melt(SE.MAF2, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.MAF3)
head(SE.MAF3)


SE.MAF4 <- dcast(SE.MAF3, formula= CLST ~ SNP)
head(SE.MAF4)
write.csv(SE.MAF4, file="SE.MAF.csv")


##Add X infront of all locusnames. 
colnames(SE.MAF4) <- paste("X", colnames(SE.MAF4), sep=".")

```

Run RDA

See this tutorial for the interpretation: 
REDUNDANCY ANALYSIS TUTORIAL: Landscape Genetics Paul Gugger
redundancy-analysis-for-landscape-genetics.pdf on mac


```
library(vegan)

GenData <- SE.MAF4
GenData$X.CLST <- NULL

ClimData <- read.csv(Climate.Data <- read.csv("SE.Full.Climate.test.csv", header=T)

GenData <- SE.MAF4
GenData$X.CLST <- NULL

##1. Run RDA with only climate data 
##H0: climate data does not affect genotype

RDA.SEfull <- rda(GenData~Avgtemp.90days+bio2+bio15+bio18, Climate.Data)
anova(RDA.SEfull)

> anova(RDA.SEfull)
Set of permutations < 'minperm'. Generating entire set.
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 999

Model: rda(formula = GenData ~ Avgtemp.90days + bio2 + bio15 + bio18, data = Climate.Data)
         Df Variance      F Pr(>F)  
Model     4   106.47 2.7022  0.027 *
Residual  2    19.70                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##So H0 rejected -> climate affects genotype
##to see which variables are most important, we can plot the results in a biplot

plot(RDA.SEfull)

summary(RDA.SEfull)
plot(RDA.SEfull)

##Partial out geog
pRDA.geog <- rda(GenData~bio2+bio15+bio18+Avgtemp.90days+ Condition(x + y), Climate.Data)
head(summary(pRDA.geog))

##Partial out climate
pRDA.climate <- rda(GenData~x+y + Condition(bio2+bio15+bio18+Avgtemp.90days), Climate.Data)
head(summary(pRDA.climate))


```
