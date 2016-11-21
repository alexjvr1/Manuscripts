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

Step 5. Remove indivs with >0.5 missing data

  - SFS and LD

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
plink --file --SE.s4.plink --recode --recodeA
plink --file SE.s4.plink --recode --recodeA
plink --file SE.s4.plink --freq --out SE.s4
```

```
##R

SE.s1.freq <- read.table("SE.s1.frq", header=T)
SE.s2.freq <- read.table("SE.s2.frq", header=T)
SE.s3.freq <- read.table("SE.s3.frq", header=T)
SE.s4.freq <- read.table("SE.s4.frq", header=T)

my.bin.width <- 0.05

par(mfrow=c(2,2))
hist(SE.s1.freq$MAF, main="SE s1 (max.miss0.8; 15294; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s2.freq$MAF, main="SE s2 (max.miss0.8; maf 0.05; 4519; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s3.freq$MAF, main="SE s3 (max.miss0.8; maf 0.05; thin; 2199; 193) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(SE.s4.freq$MAF, main="SE s4 (max.miss0.8; maf 0.05; thin; 2199; 171) SFS", breaks=seq(0,0.5, by=my.bin.width))

```


![alt_txt][SFS.s1-4]
[SFS.s1-4]:https://cloud.githubusercontent.com/assets/12142475/20431680/a805d262-ad9b-11e6-95c9-4390efaeb9c9.png


Missingness across individuals
```
alexjvr$ vcftools --vcf SE.s4.171.recode.vcf --missing-indv --out SE.s4
 
##R
library(ggplot2)
SE.s4.miss <- read.table("SE.s4.imiss", header=T)
pop.s4 <- read.table("SE.s4.popnames", header=F)
SE.s4.miss$pop <- pop.s4$V1
SE.s4.miss$pop.order <- pop.s4$V2

SE.s4.miss.sort <- SE.s4.miss[order(SE.s4.miss$pop.order),]

SE.s4.miss.sort <- SE.s4.miss[order(SE.s4.miss$pop.order),]

SE.s4.miss.sort$pop <- factor(SE.s4.miss.sort$pop, levels=SE.s4.miss.sort$pop)
```

![alt_txt][miss.s4]
[miss.s4]:https://cloud.githubusercontent.com/assets/12142475/20431948/394aada0-ad9d-11e6-9162-12229591b63e.png




##Subset data

into populations and regions


/Users/alexjvr/2016RADAnalysis/5_SE.MS1/SE.FullData.Analyses/sumstats/DataFiltering/SE.Full_Filtered_20161118
```
mkdir subset.data
for i in $(ls popnames.plink.folder/); do plink --file SE.s4.plink --keep popnames.plink.folder/$i --recode --recodeA --out subset.data/$i.plink; done
```

and calculate SFS and LD for each 

```
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


#And with only variable loci: 

DE.3pop.frq.var <- subset(DE.3pop.frq, MAF>0.001)
Sk.3pop.frq.var <- subset(Sk.3pop.frq, MAF>0.001)
Upp.3pop.frq.var <- subset(Upp.3pop.frq, MAF>0.001)
Umea.3pop.frq.var <- subset(Umea.3pop.frq, MAF>0.001)
Lulea.3pop.frq.var <- subset(Lulea.3pop.frq, MAF>0.001)
Kir.2pop.frq.var <- subset(Kir.2pop.frq, MAF>0.001)
FIN.frq.var <- subset(FIN.frq, MAF>0.001)

my.bin.width <- 0.01

par(mfrow=c(2, 4))
hist(DE.3pop.frq.var$MAF, main="DE.3pop.var (29, 1680, 23.6% fixed) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Sk.3pop.frq.var$MAF, main="Sk.3pop (48, 2091, 4.1% fixed) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Upp.3pop.frq.var$MAF, main="Upp.3pop (28, 1261, 42.7%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Umea.3pop.frq.var$MAF, main="Umea.3pop (19, 1146, 47.9%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Lulea.3pop.frq.var$MAF, main="Lulea.3pop (19, 1123, 48.9%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(Kir.2pop.frq.var$MAF, main="Kir.2pop (20, 1196, 45.6%) SFS", breaks=seq(0,0.5, by=my.bin.width))
hist(FIN.frq.var$MAF, main="FIN (7, 901, 59%) SFS", breaks=seq(0,0.5, by=my.bin.width))


###Looking only at fixed loci

#And without the fixed loci: (Frequency here is based on 1/(nx2) for MAC >1)

DE.3.pop.frq.fix <- subset(DE.3pop.frq, MAF<0.0018)
Sk.3pop.frq.fix <- subset(Sk.3pop.frq, MAF<0.00105)
Upp.3pop.frq.fix <- subset(Upp.3pop.frq, MAF<0.00179)
Umea.3pop.frq.fix <- subset(Umea.3pop.frq, MAF<0.0027)
Lulea.3pop.frq.fix <- subset(Lulea.3pop.frq, MAF<0.0027)
Kir.2pop.frq.fix <- subset(Kir.2pop.frq, MAF<0.0026)
FIN.frq.fix <- subset(FIN.frq, MAF<0.070)

par(mfrow=c(2, 4))

SE.region.frq.fix.table <- do.call(rbind, lapply(ls(pattern="fix$"), get))

SE.fix.region.table.keep <- data.frame(table(SE.region.frq.fix.table$SNP))
my.bin.width=1
hist(SE.fix.region.table.keep$Freq, xlab="Number of pops", ylab="Frequency", main="Frequency of fixed loci found in increasing number of populations", breaks=seq(0,7, by=my.bin.width))
```


![alt_txt][SFS.s4.regions]
[SFS.s4.regions]:https://cloud.githubusercontent.com/assets/12142475/20433442/df6d3152-ada3-11e6-9a02-4123ad1ff7f4.png

![alt_txt][var.only]
[var.only]:https://cloud.githubusercontent.com/assets/12142475/20435155/4d959fb4-adab-11e6-99d7-cc0fd54519f2.png

![alt_txt][freq.fixed]
[freq.fixed]:https://cloud.githubusercontent.com/assets/12142475/20435161/56ab279a-adab-11e6-8ad5-2a512cb2c689.png




Variable sites only
![alt_txt][SFS.var.allregions]
[SFS.var.allregions]:https://cloud.githubusercontent.com/assets/12142475/20478502/ad6c623e-afda-11e6-8607-c630767538ba.png

Adjusted bin width increments of ~MAC 1 
![alt_txt][SFS.var.allregions.adjbinwidth]
[SFS.var.allregions.adjbinwidth]:https://cloud.githubusercontent.com/assets/12142475/20478501/ad4d9f84-afda-11e6-9e3e-c936431fac6e.png







##R2 LD

Checking for linkage 
```
##R
write.table(DE.3pop.frq.fix$SNP, "DE.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Sk.3pop.frq.fix$SNP, "Sk.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Upp.3pop.frq.fix$SNP, "Upp.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Umea.3pop.frq.fix$SNP, "Umea.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Lulea.3pop.frq.fix$SNP, "Lulea.3pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(Kir.2pop.frq.fix$SNP, "Kir.2pop.fixedloci", row.names=F, col.names=F, quote=F)
write.table(FIN.frq.fix$SNP, "FIN.fixedloci", row.names=F, col.names=F, quote=F)


##linux
plink --file subset.data/DE.3pop.plink --exclude DE.3pop.fixedloci --recode --recodeA --out subset.data/DE.3pop.var.plink
plink --file subset.data/DE.3pop.var.plink --r2 --out subset.data/DE.3pop.var
plink --file subset.data/Sk.3pop.plink --exclude Sk.3pop.fixedloci --recode --recodeA --out subset.data/Sk.3pop.var.plink
plink --file subset.data/Sk.3pop.var.plink --r2 --out subset.data/Sk.3pop.var
plink --file subset.data/Upp.3pop.plink --exclude Upp.3pop.fixedloci --recode --recodeA --out subset.data/Upp.3pop.var.plink
plink --file subset.data/Upp.3pop.var.plink --r2 --out subset.data/Upp.3pop.var
plink --file subset.data/Umea.3pop.plink --exclude Umea.3pop.fixedloci --recode --recodeA --out subset.data/Um.3pop.var.plink
plink --file subset.data/Um.3pop.var.plink --r2 --out subset.data/Um.3pop.var
plink --file subset.data/Lulea.3pop.plink --exclude Lulea.3pop.fixedloci --recode --recodeA --out subset.data/Lulea.3pop.var.plink
plink --file subset.data/Lulea.3pop.var.plink --r2 --out subset.data/Lulea.3pop.var
plink --file subset.data/Kir.2pop.plink --exclude Kir.2pop.fixedloci --recode --recodeA --out subset.data/Kir.2pop.var.plink
plink --file subset.data/Kir.2pop.var.plink --r2 --out subset.data/Kir.2pop.var
plink --file subset.data/FIN.plink --exclude FIN.fixedloci --recode --recodeA --out subset.data/FIN.var.plink
plink --file subset.data/FIN.var.plink --r2 --out subset.data/FIN.var




##R
DE.3pop.var.ld <- read.table("subset.data/DE.3pop.var.ld", header=T)
Sk.3pop.var.ld <- read.table("subset.data/Sk.3pop.var.ld", header=T)
Upp.3pop.var.ld <- read.table("subset.data/Upp.3pop.var.ld", header=T)
Lulea.3pop.var.ld <- read.table("subset.data/Lulea.3pop.var.ld", header=T)
Um.3pop.var.ld <- read.table("subset.data/Um.3pop.var.ld", header=T)
Kir.2pop.var.ld <- read.table("subset.data/Kir.2pop.var.ld", header=T)
FIN.var.ld <- read.table("subset.data/FIN.var.ld", header=T)

my.bin.width=0.05

par(mfrow=c(2,4))
hist(DE.3pop.var.ld$R2, main="DE.3pop.var (29, 1680 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Sk.3pop.var.ld$R2, main="Sk.3pop.var (29, 2091 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Upp.3pop.var.ld$R2, main="Upp.3pop (28, 1261 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Um.3pop.var.ld$R2, main="Umea.3pop (19, 1146 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Lulea.3pop.var.ld$R2, main="Lulea.3pop (19, 1123 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
hist(Kir.2pop.var.ld$R2,main="Kir.2pop (20, 1196 loci) R2",  breaks=seq(0,1.0, by=my.bin.width))
hist(FIN.var.ld$R2, main="FIN (7, 901 loci) R2", breaks=seq(0,1.0, by=my.bin.width))
```

R2 of variable loci only

Starts at 0.2, because: 

Finland still has high R2, because: 

![alt_txt][R2.variableLoci]
[R2.variableLoci]:https://cloud.githubusercontent.com/assets/12142475/20479938/91400c40-afe0-11e6-9fff-e3b6b14da7f6.png



