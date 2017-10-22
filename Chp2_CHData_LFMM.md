#  LFMM analyses of CHdata

I want to identify loci associated with environment for the 6 CH datasets. 

1. CHall.932.7744.recode.vcf

2. CHN.229.5265.recode.vcf

3. CHS.275.6339.recode.vcf

4. CHS.VS.135.5835.recode.vcf

5. CHS.TI.140.5692.recode.vcf

6. CZ.404.7288.recode.vcf


#### Env Variables

Identified in https://github.com/alexjvr1/Manuscripts/blob/master/5.CHP2_CH_LandscapeGenomics.md

I've chosen a different 5 environmental parameters for which to run BayEnv2 and LFMM. Based on the site-specific parameters that were calculated by Josh, and the reduction based on corrolation of 0.8. 

In the order that they are in the lfmm.env file 

1. shadow.days

2. solar.rad.60d

3. pcpt.60d

4. day10cm

5. temp.laying.date



### Input files for LFMM

Working in this directory: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM/Oct2017

See tutorial for command line version: http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html



Input files

1. .env environmental file

2. Genotype file


#### Prepare input files

Make sure that R3.3.1 is launched. (Or later than R3.2)

LEA is a bioconductor package.

http://www.bioconductor.org/packages/release/bioc/html/LEA.html

In R:

```
#source("http://bioconductor.org/biocLite.R")
#biocLite("LEA")

library(LEA)
```


###### 1. .env file

Convert the environmental data to the specific .ENV format needed by lfmm. Remember that .csv environmental input file needs env data for 
all individuals. Make sure that all the empty cells are deleted in the .csv. Any NA's in the .env file will cause problems 
(mismatch between nr of indivs and env data). I got a strange error for this: Error: SNP 9 same across all indivs. 
But once I removed the NA's from the .env file the run started without problems. 

```
env <- read.csv("CHN.n5.env.txt", header=F) ###read in the environmental data
write.env(env, "CHN.n5.env")   ##convert to correct .env format
[1] "CHN.n5.env"
```


###### 2. .geno file

Convert the vcf file to lfmm format in R

```
library(LEA)

genotype = vcf2geno("name.vcf")
```


##### Pop structure


use snmf to determine how many latent factors best describe the population structure:
```
obj.snmf = snmf(genotype, K = 1:20, entropy = T, ploidy = 2, project="new")   ###K range to be tested should be increased as needed. 
plot(obj.snmf)  ###see at which K cross-entropy is the lowest

```


CHN: K4

![alt_txt][CHN.K4]
[CHN.K4]:https://cloud.githubusercontent.com/assets/12142475/19263984/339c8016-8f9f-11e6-9c78-6782777aa717.png



CHS.TI: K6

![alt_txt][CHS.TI.K6]
[CHS.TI.K6]:https://cloud.githubusercontent.com/assets/12142475/19264160/ed7b7992-8f9f-11e6-82af-74788205ae5b.png


CHS.VS: K4
![alt_txt][CHS.VS.K4]
[CHS.VS.K4]:https://cloud.githubusercontent.com/assets/12142475/19264069/94ee38c8-8f9f-11e6-9667-bb48bfef4c83.png


CHS: K10
![alt_txt][CHS.K10]
[CHS.K10]:https://cloud.githubusercontent.com/assets/12142475/19264518/487c5ab8-8fa1-11e6-9d1e-3690e948430c.png


CZ: K9

![alt_txt][CZ.K9]
[CZ.K9]:https://cloud.githubusercontent.com/assets/12142475/19264784/800d6fde-8fa2-11e6-8517-bcb3832d2619.png


CHall K19 or 20
![alt_txt][CHall.K19]
[CHall.K19]:https://cloud.githubusercontent.com/assets/12142475/19264778/7b663f9c-8fa2-11e6-8998-6faa85e504be.png



Take note of the K for each sub-group. ie. It makes sense that CHS.VS + CHS.TI = CHStotal



##### Run LFMM

And run lfmm at the chosen K. This should be run on the server. 
```
project = lfmm("SEsubset.Final2.geno", "SE.env", K = 2, repetitions = 5, project = "new")  ##run LFMM with this .env file

```

Runs can be exported for easier transfer from server to local computer

```
export.lfmmProject(obj.lfmmProject)   ##this creates a .zip file in the current directory, which can easily be transferred to other folders. 

import.lfmmProject(obj.lfmmProject) ##to import the project into R from the new folder

```



###### Postprocessing: 

Based on the lambda calculated from run 1, I have decided to run 3 K for each dataset. Lambda should be close to or just below 1. Mine were all just above 1. 

stba had to be removed from CHS.TI, because of missing data in the env. variables. So this is also removed from CHS. 



First calculate the lambda for each K to determine what the correct K is: 


##### CHN
snmf K = 4

This is the old way of doing it (but I keep this here since this is how these projects were run). The new way is to export and import LFMM: 

Copy everything to the mac. And remember to change all the paths in the .lfmmProject file.

delimiter can be anything. I'm using # here because I have / in my filenames. 
```
sed -i -e 's#old/path#new/path#g' file.lfmmProject
```


Open R and import project. This CHN project contains K3-5. 5 runs for each K. d1-5. I show only results from K=4 here: 
```
project=import.lfmm("CHN.229.8951.recode_CHN.lfmmProject")


zs.d1 <- z.scores(project, K=4, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.448757
    
lambda=median(zs.d1.median^2)/0.6
lambda
  #  [1] 1.098487

adj.p.values.d1 =pchisq(zs.d1.median^2/0.6, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=4, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1]  1.077905
lambda=median(zs.d2.median^2)/0.85
lambda
   # [1] 0.5769154
adj.p.values.d2 =pchisq(zs.d2.median^2/0.4549364, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=4, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.613051

lambda=median(zs.d3.median^2)/0.73
lambda
   # [1] 1.005254
adj.p.values.d3 =pchisq(zs.d3.median^2/0.73, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=4, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.851754
lambda=median(zs.d4.median^2)/0.8
lambda
 #   [1] 1.053038
adj.p.values.d4 =pchisq(zs.d4.median^2/0.8, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=4, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 1.508529
lambda=median(zs.d5.median^2)/0.70
lambda
   # [1] 0.9804067
adj.p.values.d5 =pchisq(zs.d5.median^2/0.7, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d5)
w = which(sort(adj.p.values.d5) < q * (1:L) / L)
candidates.d5.k5 = order(adj.p.values.d5)[w]


pdf(file="CHN.LFMM.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
hist(adj.p.values.d3)
hist(adj.p.values.d4)
hist(adj.p.values.d5)
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




candidates.d3.k5 <- as.character(candidates.d3.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d3.k5.names <- locus.names[locus.names$ID %in% candidates.d3.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d3.k5.names <- paste("X", candidates.d3.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d4.k5 <- as.character(candidates.d4.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d4.k5.names <- locus.names[locus.names$ID %in% candidates.d4.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d4.k5.names <- paste("X", candidates.d4.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d5.k5 <- as.character(candidates.d5.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d5.k5.names <- locus.names[locus.names$ID %in% candidates.d5.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d5.k5.names <- paste("X", candidates.d5.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel
```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)
d3 <- length(candidates.d3.k5.names)
d4 <- length(candidates.d4.k5.names)
d5 <- length(candidates.d5.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))
d13 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names)))
d14 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names)))
d15 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d5.k5.names)))
d23 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names)))
d24 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names)))
d25 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d5.k5.names)))
d34 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names)))
d35 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d5.k5.names)))
d45 <- length(Reduce(intersect, list(candidates.d4.k5.names, candidates.d5.k5.names)))

d123 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d3.k5.names)))
d124 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d4.k5.names)))
d125 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d5.k5.names)))
d234 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d134 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d135 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d145 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d235 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))

d1234 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d4.k5.names)))
d1235 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d5.k5.names)))
d2345 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d12345 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, 
candidates.d4.k5.names, candidates.d5.k5.names)))

pdf(file="Venn.CHN.n5.LFMMonly_20170310.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("shadow.days", "sol.rad.60d", "pcpt.60d", "day10cm", "temp.laying.date"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```

![alt_txt][CHN.Venn]

[CHN.Venn]:https://user-images.githubusercontent.com/12142475/31121243-52d11a88-a82f-11e7-81b8-6cce216ca8e9.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
d2.names <- as.data.frame(candidates.d2.k5.names)
d3.names <- as.data.frame(candidates.d3.k5.names)
d4.names <- as.data.frame(candidates.d4.k5.names)
d5.names <- as.data.frame(candidates.d5.k5.names)

colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12345.names <- as.data.frame(d12345.names)

write.table(d12345.names$d12345.names, "CHN.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

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

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM/CHS

Open R and import project. This CHN project contains K3-5. 5 runs for each K. d1-5. I show only results from K=4 here: 
```
CHS.lfmm=import.lfmmProject("CHS.275.6339.recode_CHS.275_lfmmProject.zip")

project <- CHS.lfmm

zs.d1 <- z.scores(project, K=10, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1]  1.697623
    
lambda=median(zs.d1.median^2)/0.75
lambda
  #  [1]  1.029747

adj.p.values.d1 =pchisq(zs.d1.median^2/0.75, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=10, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.537919
lambda=median(zs.d2.median^2)/0.65
lambda
   # [1] 1.076393
adj.p.values.d2 =pchisq(zs.d2.median^2/0.65, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=10, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.770662

lambda=median(zs.d3.median^2)/0.75
lambda
   # [1] 1.074051
adj.p.values.d3 =pchisq(zs.d3.median^2/0.75, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=10, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.670943
lambda=median(zs.d4.median^2)/0.75
lambda
 #   [1] 1.013564
adj.p.values.d4 =pchisq(zs.d4.median^2/0.75, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=10, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 1.755527
lambda=median(zs.d5.median^2)/0.8
lambda
   # [1] 0.9983165
adj.p.values.d5 =pchisq(zs.d5.median^2/0.8, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d5)
w = which(sort(adj.p.values.d5) < q * (1:L) / L)
candidates.d5.k5 = order(adj.p.values.d5)[w]


pdf(file="CHS.LFMM.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
hist(adj.p.values.d3)
hist(adj.p.values.d4)
hist(adj.p.values.d5)
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




candidates.d3.k5 <- as.character(candidates.d3.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d3.k5.names <- locus.names[locus.names$ID %in% candidates.d3.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d3.k5.names <- paste("X", candidates.d3.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d4.k5 <- as.character(candidates.d4.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d4.k5.names <- locus.names[locus.names$ID %in% candidates.d4.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d4.k5.names <- paste("X", candidates.d4.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d5.k5 <- as.character(candidates.d5.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d5.k5.names <- locus.names[locus.names$ID %in% candidates.d5.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d5.k5.names <- paste("X", candidates.d5.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel
```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)
d3 <- length(candidates.d3.k5.names)
d4 <- length(candidates.d4.k5.names)
d5 <- length(candidates.d5.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))
d13 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names)))
d14 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names)))
d15 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d5.k5.names)))
d23 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names)))
d24 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names)))
d25 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d5.k5.names)))
d34 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names)))
d35 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d5.k5.names)))
d45 <- length(Reduce(intersect, list(candidates.d4.k5.names, candidates.d5.k5.names)))

d123 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d3.k5.names)))
d124 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d4.k5.names)))
d125 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d5.k5.names)))
d234 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d134 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d135 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d145 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d235 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))

d1234 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d4.k5.names)))
d1235 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d5.k5.names)))
d2345 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d12345 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, 
candidates.d4.k5.names, candidates.d5.k5.names)))

pdf(file="Venn.CHS.n5.LFMMonly_20170310.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("shadow.days", "sol.rad.60d", "pcpt.60d", "day10cm", "temp.laying.date"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```

![alt_txt][CHS.Venn]

[CHS.Venn]:https://user-images.githubusercontent.com/12142475/31120683-408a8dfc-a82d-11e7-8de8-7a78297aceac.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
d2.names <- as.data.frame(candidates.d2.k5.names)
d3.names <- as.data.frame(candidates.d3.k5.names)
d4.names <- as.data.frame(candidates.d4.k5.names)
d5.names <- as.data.frame(candidates.d5.k5.names)

colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists
d12345.names <- as.data.frame(d12345.names)
write.table(d12345.names$d12345.names, "CHS.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

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

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM/CHS.VS

Open R and import project.  
```
project=import.lfmm("CHS.275.9315.n5.recode.lfmm")


zs.d1 <- z.scores(project, K=4, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.559495
    
lambda=median(zs.d1.median^2)/0.65
lambda
  #  [1] 1.091494

adj.p.values.d1 =pchisq(zs.d1.median^2/0.65, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=4, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.946384
lambda=median(zs.d2.median^2)/0.9
lambda
   # [1]  0.9838678
adj.p.values.d2 =pchisq(zs.d2.median^2/0.9, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=4, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.533234

lambda=median(zs.d3.median^2)/0.65
lambda
   # [1] 1.073114
adj.p.values.d3 =pchisq(zs.d3.median^2/0.65, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=4, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.429176
lambda=median(zs.d4.median^2)/0.65
lambda
 #   [1] 1.000283
adj.p.values.d4 =pchisq(zs.d4.median^2/0.65, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=4, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 1.596478
lambda=median(zs.d5.median^2)/0.65
lambda
   # [1] 1.117378
adj.p.values.d5 =pchisq(zs.d5.median^2/0.65, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d5)
w = which(sort(adj.p.values.d5) < q * (1:L) / L)
candidates.d5.k5 = order(adj.p.values.d5)[w]


pdf(file="CHS.VS.LFMM.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
hist(adj.p.values.d3)
hist(adj.p.values.d4)
hist(adj.p.values.d5)
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




candidates.d3.k5 <- as.character(candidates.d3.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d3.k5.names <- locus.names[locus.names$ID %in% candidates.d3.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d3.k5.names <- paste("X", candidates.d3.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d4.k5 <- as.character(candidates.d4.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d4.k5.names <- locus.names[locus.names$ID %in% candidates.d4.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d4.k5.names <- paste("X", candidates.d4.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d5.k5 <- as.character(candidates.d5.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d5.k5.names <- locus.names[locus.names$ID %in% candidates.d5.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d5.k5.names <- paste("X", candidates.d5.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel
```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)
d3 <- length(candidates.d3.k5.names)
d4 <- length(candidates.d4.k5.names)
d5 <- length(candidates.d5.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))
d13 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names)))
d14 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names)))
d15 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d5.k5.names)))
d23 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names)))
d24 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names)))
d25 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d5.k5.names)))
d34 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names)))
d35 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d5.k5.names)))
d45 <- length(Reduce(intersect, list(candidates.d4.k5.names, candidates.d5.k5.names)))

d123 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d3.k5.names)))
d124 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d4.k5.names)))
d125 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d5.k5.names)))
d234 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d134 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d135 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d145 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d235 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))

d1234 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d4.k5.names)))
d1235 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d5.k5.names)))
d2345 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d12345 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, 
candidates.d4.k5.names, candidates.d5.k5.names)))

pdf(file="Venn.CHS.VS.n5.LFMMonly_20170310.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("shadow.days", "sol.rad.60d", "pcpt.60d", "day10cm", "temp.laying.date"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```

![alt_txt][CHS.VS.Venn]

[CHS.VS.Venn]:https://user-images.githubusercontent.com/12142475/31120999-69f379dc-a82e-11e7-855f-97d8d828f30e.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
d2.names <- as.data.frame(candidates.d2.k5.names)
d3.names <- as.data.frame(candidates.d3.k5.names)
d4.names <- as.data.frame(candidates.d4.k5.names)
d5.names <- as.data.frame(candidates.d5.k5.names)

colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists
d12345.names <- as.data.frame(d12345.names)
write.table(d12345.names$d12345.names, "CHS.VS.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

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

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM/Oct2017/CHS.TI

Open R and import project.  
```
project=import.lfmmProject("CHS.TI.140.5692.recode_CHS.TI_lfmmProject.zip")


zs.d1 <- z.scores(project, K=6, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 1.396319
    
lambda=median(zs.d1.median^2)/0.6
lambda
  #  [1] 1.058727

adj.p.values.d1 =pchisq(zs.d1.median^2/0.6, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=6, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.302821
lambda=median(zs.d2.median^2)/0.60
lambda
   # [1]  0.9878343
adj.p.values.d2 =pchisq(zs.d2.median^2/0.6, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=6, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.184177

lambda=median(zs.d3.median^2)/0.52
lambda
   # [1] 1.03601
adj.p.values.d3 =pchisq(zs.d3.median^2/0.52, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=6, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.871516
lambda=median(zs.d4.median^2)/0.82
lambda
 #   [1] 1.038318
adj.p.values.d4 =pchisq(zs.d4.median^2/0.82, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=6, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 1.43036
lambda=median(zs.d5.median^2)/0.63
lambda
   # [1] 1.032894
adj.p.values.d5 =pchisq(zs.d5.median^2/0.63, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d5)
w = which(sort(adj.p.values.d5) < q * (1:L) / L)
candidates.d5.k5 = order(adj.p.values.d5)[w]


pdf(file="CHS.TI.LFMM.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
hist(adj.p.values.d3)
hist(adj.p.values.d4)
hist(adj.p.values.d5)
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




candidates.d3.k5 <- as.character(candidates.d3.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d3.k5.names <- locus.names[locus.names$ID %in% candidates.d3.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d3.k5.names <- paste("X", candidates.d3.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d4.k5 <- as.character(candidates.d4.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d4.k5.names <- locus.names[locus.names$ID %in% candidates.d4.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d4.k5.names <- paste("X", candidates.d4.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d5.k5 <- as.character(candidates.d5.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d5.k5.names <- locus.names[locus.names$ID %in% candidates.d5.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d5.k5.names <- paste("X", candidates.d5.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel
```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)
d3 <- length(candidates.d3.k5.names)
d4 <- length(candidates.d4.k5.names)
d5 <- length(candidates.d5.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))
d13 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names)))
d14 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names)))
d15 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d5.k5.names)))
d23 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names)))
d24 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names)))
d25 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d5.k5.names)))
d34 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names)))
d35 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d5.k5.names)))
d45 <- length(Reduce(intersect, list(candidates.d4.k5.names, candidates.d5.k5.names)))

d123 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d3.k5.names)))
d124 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d4.k5.names)))
d125 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d5.k5.names)))
d234 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d134 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d135 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d145 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d235 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))

d1234 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d4.k5.names)))
d1235 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d5.k5.names)))
d2345 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d12345 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, 
candidates.d4.k5.names, candidates.d5.k5.names)))

pdf(file="Venn.CHS.TI.n5.LFMMonly_20170310.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("shadow.days", "sol.rad.60d", "pcpt.60d", "day10cm", "temp.laying.date"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```

![alt_txt][CHS.TI.Venn]

[CHS.TI.Venn]:https://user-images.githubusercontent.com/12142475/31121490-47e83ea2-a830-11e7-8269-b799ae35fbcf.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
d2.names <- as.data.frame(candidates.d2.k5.names)
d3.names <- as.data.frame(candidates.d3.k5.names)
d4.names <- as.data.frame(candidates.d4.k5.names)
d5.names <- as.data.frame(candidates.d5.k5.names)

colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12345.names <- as.data.frame(d12345.names)

write.table(d12345.names$d12345.names, "CHS.TI.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

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

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM/CZ

Open R and import project.  
```
project=import.lfmmProject("CZ.404.7288.recode_CZ_lfmmProject.zip")


zs.d1 <- z.scores(project, K=9, d=1)
zs.d1.median =apply(zs.d1, MARGIN=1, median)


lambda=median(zs.d1.median^2)/0.4549364
lambda
   # [1] 0.9234408
    
lambda=median(zs.d1.median^2)/0.4
lambda
  #  [1] 1.050267

adj.p.values.d1 =pchisq(zs.d1.median^2/0.4, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=9, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 1.386879
lambda=median(zs.d2.median^2)/0.60
lambda
   # [1]  1.05157
adj.p.values.d2 =pchisq(zs.d2.median^2/0.6, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=9, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.453828

lambda=median(zs.d3.median^2)/0.62
lambda
   # [1] 1.066773
adj.p.values.d3 =pchisq(zs.d3.median^2/0.62, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=9, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.45827
lambda=median(zs.d4.median^2)/0.64
lambda
 #   [1] 1.036594
adj.p.values.d4 =pchisq(zs.d4.median^2/0.64, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=9, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 1.306648
lambda=median(zs.d5.median^2)/0.58
lambda
   # [1] 1.0249
adj.p.values.d5 =pchisq(zs.d5.median^2/0.58, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d5)
w = which(sort(adj.p.values.d5) < q * (1:L) / L)
candidates.d5.k5 = order(adj.p.values.d5)[w]


pdf(file="CZ.LFMM.hist.pdf")
par(mfrow=c(3,2))
hist(adj.p.values.d1)
hist(adj.p.values.d2)
hist(adj.p.values.d3)
hist(adj.p.values.d4)
hist(adj.p.values.d5)
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




candidates.d3.k5 <- as.character(candidates.d3.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d3.k5.names <- locus.names[locus.names$ID %in% candidates.d3.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d3.k5.names <- paste("X", candidates.d3.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d4.k5 <- as.character(candidates.d4.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d4.k5.names <- locus.names[locus.names$ID %in% candidates.d4.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d4.k5.names <- paste("X", candidates.d4.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel



candidates.d5.k5 <- as.character(candidates.d5.k5)  ##change the list of candidates from LFMM output to a list of characters
candidates.d5.k5.names <- locus.names[locus.names$ID %in% candidates.d5.k5,]  ##select from locus.names$ID the rows that match candidates vector

candidates.d5.k5.names <- paste("X", candidates.d5.k5.names$SNP, sep=".") #rename the SNPs so that they don't get renamed in excel
```

VennDiagram of data

```
library(VennDiagram)

d1 <- length(candidates.d1.k5.names)
d2 <- length(candidates.d2.k5.names)
d3 <- length(candidates.d3.k5.names)
d4 <- length(candidates.d4.k5.names)
d5 <- length(candidates.d5.k5.names)

d12 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names)))
d13 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names)))
d14 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names)))
d15 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d5.k5.names)))
d23 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names)))
d24 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names)))
d25 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d5.k5.names)))
d34 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names)))
d35 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d5.k5.names)))
d45 <- length(Reduce(intersect, list(candidates.d4.k5.names, candidates.d5.k5.names)))

d123 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d3.k5.names)))
d124 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d4.k5.names)))
d125 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names,candidates.d5.k5.names)))
d234 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d134 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d4.k5.names)))
d135 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d145 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d235 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,candidates.d5.k5.names)))
d245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))
d345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d4.k5.names,candidates.d5.k5.names)))

d1234 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d4.k5.names)))
d1235 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, candidates.d5.k5.names)))
d2345 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d3.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1245 <- length(Reduce(intersect, list(candidates.d2.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d1345 <- length(Reduce(intersect, list(candidates.d3.k5.names, candidates.d1.k5.names,
candidates.d4.k5.names, candidates.d5.k5.names)))
d12345 <- length(Reduce(intersect, list(candidates.d1.k5.names, candidates.d2.k5.names, candidates.d3.k5.names, 
candidates.d4.k5.names, candidates.d5.k5.names)))

pdf(file="Venn.CZ.n5.LFMMonly_20170310.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("shadow.days", "sol.rad.60d", "pcpt.60d", "day10cm", "temp.laying.date"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```

![alt_txt][CZ.Venn]

[CZ.Venn]:https://user-images.githubusercontent.com/12142475/31121744-612b5f2e-a831-11e7-9af5-bd66a3c3c470.png


Prepare the data to incorporate with the Fst outlier graph
```
d1.names <- as.data.frame(candidates.d1.k5.names)
d2.names <- as.data.frame(candidates.d2.k5.names)
d3.names <- as.data.frame(candidates.d3.k5.names)
d4.names <- as.data.frame(candidates.d4.k5.names)
d5.names <- as.data.frame(candidates.d5.k5.names)

colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

d12345.names <- as.data.frame(d12345.names)

write.table(d12345.names$d12345.names, "CZ.LFMM.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

CZ.LFMM.outliers <- read.table("CZ.LFMM.alloutliers", header=F)
colnames(CZ.LFMM.outliers) <- ("loci")
CZ.LFMM.outliers <- as.character(CZ.LFMM.outliers$loci)
```



##############################
##############################
# Assume that a fraction of loci (frac) has genomic effects drawn from a different
# distribution (Chisq with 3 df, median z.2 = 2.38).
# This means that some loci have truly larger effects than expected due to chance.

frac <- 0.1              # fraction of loci sampled from the different distribution
lambda.adjust <- "no"    # should we correct for bias using the genomic inflation factor?
loci <- 10000            # nr of loci
set.seed(2)

z.2 <- rchisq(n = loci, df = 1, ncp = 0)               # generate z.2 values from the chi-sq distribution
if(frac > 0) z.2[1:(round(frac*loci))] <- rchisq(n = round(frac*loci), df = 3, ncp = 0)
hist(z.2, breaks = 19)
lambda <- median(z.2) / 0.456
par(mfrow = c(2, 1))
p.values <- pchisq(z.2, df = 1, lower.tail = FALSE)
hist(p.values, breaks = 19, main = "Initial distribution of p-values")
abline(a=(loci/20), b = 0, lty = 3, col = "red")
if(lambda.adjust == "yes") p.values = pchisq(z.2/lambda, df = 1, lower = FALSE)
h <- hist(p.values, col = "green", main = "Distribution after lambda adjustment, if requested")$counts
abline(a=(loci/20), b = 0, lty = 3, col = "red")


# This is from the appendix of Francois et al. (2016, Mol. Ecol.)
q.level = 0.1         # set the false discovery rate
L = length(p.values)
       # Benjamini-Hochberg algorithm
w = which(sort(p.values) < q.level * (1:L) / L)
candidates = order(p.values)[w]
candidates
adjusted.p <- p.adjust(p.values, method = "BH", n = length(p))    # this returns adjusted p-values, I think under an FDR = 0.05
length(adjusted.p[adjusted.p < 0.05])                             # How many are still significant?


# From Story & Tibshirani (2003, PNAS).
# Flat part of p-value distribution indicates expected number of false positives.
# Tuning parameter indicates the p-value beyond which the distribution appears flat. Most of the
# p-values in this area are presumably truly null.
tuning.parameter <- 0.5
nr.null <- 20 * sum( h[(1+tuning.parameter*20):20] ) / length((1+tuning.parameter*20):20)
abline(a=(nr.null/20), b = 0, lty = 3, col = "blue")
nr.null / loci                                                          # this is the estimated proportion of loci truly null.
0.05 * loci * (nr.null / loci) / (loci - sum( h[ (1+0.05*20):20 ] ))    # false discovery rate, pg. 9443
#
# Red dotted line is expected if all loci were null.
# Blue dotted line is our estimate of the true number of null loci.



for i in $(ls *.fq.trim.gz.sai.bam.sorted.bam)
do samtools view -q 30 -b ${i} > ${i}.q30.bam
done

plot(fstall.hist, xlim=c(-0.15,1.0))
plot(EAA.hist, xlim=c(-0.15,1.0), col=rgb(0,0,1,1/4), add=T)
plot(Outlier.hist, xlim=c(-0.15,1.0), col=rgb(1,0,0,1/4), add=T)



