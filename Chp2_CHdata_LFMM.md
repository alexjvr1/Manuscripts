#LFMM analyses of CHdata

I want to identify loci associated with environment for the 6 CH datasets

1. CHall.940.9608

2. CHN.229.8951 (lfmm removes all non-variant loci)

3. CHS.283.9608

4. CHS.VS.135.9608

5. CHS.TI.148.9608

6. CZ.404.9608


####Env Variables

Identified in https://github.com/alexjvr1/Manuscripts/blob/master/5.CHP2_CH_LandscapeGenomics.md

I've chosen a different 5 environmental parameters for which to run BayEnv2 and LFMM. Based on the site-specific parameters that were calculated by Josh, and the reduction based on corrolation of 0.8. 

In the order that they are in the lfmm.env file 

1. shadow.days

2. solar.rad.60d

3. pcpt.60d

4. day10cm

5. temp.laying.date

For CHS.VS temp.laying.date was mistakenly repeated as d=1. So there are 6 env variables; use d=2-6



###Input files for LFMM

Working in this directory: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/LFMM

Note that there are 2 runs. Run_2 for all datasets uses the 5 environmental vaiables I want to use for the final analysis


See tutorial for command line version: http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LEA_1.html

/Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/LFMM/SE132.FINAL

Input files

1. .env environmental file

2. Genotype file


####Prepare input files

Make sure that R3.3.1 is launched. (Or later than R3.2)

LEA is a bioconductor package.

http://www.bioconductor.org/packages/release/bioc/html/LEA.html

In R:

```
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

library(LEA)
```


######1. .env file

Convert the environmental data to the specific .ENV format needed by lfmm. Remember that .csv environmental input file needs env data for 
all individuals. Make sure that all the empty cells are deleted in the .csv. Any NA's in the .env file will cause problems 
(mismatch between nr of indivs and env data). I got a strange error for this: Error: SNP 9 same across all indivs. 
But once I removed the NA's from the .env file the run started without problems. 

```
env <- read.csv("CHN.n5.env.txt", header=F) ###read in the environmental data
write.env(env, "CHN.n5.env")   ##convert to correct .env format
[1] "CHN.n5.env"
```


######2. .geno file

Convert the vcf file to lfmm format in R

```
library(LEA)

genotype = vcf2geno("name.vcf")
```


#####Pop structure


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


CHall K19 or 20  - but what about overfitting???
![alt_txt][CHall.K19]
[CHall.K19]:https://cloud.githubusercontent.com/assets/12142475/19264778/7b663f9c-8fa2-11e6-8998-6faa85e504be.png



Take note of the K for each sub-group. ie. It makes sense that CHS.VS + CHS.TI = CHStotal



#####Run LFMM

And run lfmm at the chosen K. This should be run on the server. 
```
project = lfmm("SEsubset.Final2.geno", "SE.env", K = 2, repetitions = 5, project = "new")  ##run LFMM with this .env file

```

Runs can be exported for easier transfer from server to local computer

```
export.lfmmProject(obj.lfmmProject)   ##this creates a .zip file in the current directory, which can easily be transferred to other folders. 

import.lfmmProject(obj.lfmmProject) ##to import the project into R from the new folder

```



######Postprocessing: 

Based on the lambda calculated from run 1, I have decided to run 3 K for each dataset. Lambda should be close to or just below 1. Mine were all just above 1. 

stba had to be removed from CHS.TI, because of missing data in the env. variables. So this is also removed from CHS. 



First calculate the lambda for each K to determine what the correct K is: 

#####CHN
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
   # [1] 1.330662
    
lambda=median(zs.d1.median^2)/0.85
lambda
  #  [1] 0.7121958

adj.p.values.d1 =pchisq(zs.d1.median^2/0.4549364, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d1)
w = which(sort(adj.p.values.d1) < q * (1:L) / L)
candidates.d1.k5 = order(adj.p.values.d1)[w]



zs.d2 <- z.scores(project, K=5, d=2)
zs.d2.median =apply(zs.d2, MARGIN=1, median)
lambda=median(zs.d2.median^2)/0.4549364
lambda
 #   [1] 0.9810137
lambda=median(zs.d2.median^2)/0.85
lambda
   # [1] 0.5250574
adj.p.values.d2 =pchisq(zs.d2.median^2/0.4549364, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d2)
w = which(sort(adj.p.values.d2) < q * (1:L) / L)
candidates.d2.k5 = order(adj.p.values.d2)[w]




zs.d3 <- z.scores(project, K=5, d=3)
zs.d3.median =apply(zs.d3, MARGIN=1, median)
lambda=median(zs.d3.median^2)/0.4549364
lambda
   # [1] 1.377389

lambda=median(zs.d3.median^2)/0.85
lambda
   # [1] 0.737205
adj.p.values.d3 =pchisq(zs.d3.median^2/0.4549364, df=1, lower=F)
q = 0.05
L = length(adj.p.values.d3)
w = which(sort(adj.p.values.d3) < q * (1:L) / L)
candidates.d3.k5 = order(adj.p.values.d3)[w]




zs.d4 <- z.scores(project, K=5, d=4)
zs.d4.median =apply(zs.d4, MARGIN=1, median)
lambda=median(zs.d4.median^2)/0.4549364
lambda
 #   [1] 1.569119
lambda=median(zs.d4.median^2)/0.85
lambda
 #   [1] 0.8398226
adj.p.values.d4 =pchisq(zs.d4.median^2/0.85, df=1, lower=F)

q = 0.05
L = length(adj.p.values.d4)
w = which(sort(adj.p.values.d4) < q * (1:L) / L)
candidates.d4.k5 = order(adj.p.values.d4)[w]



zs.d5 <- z.scores(project, K=5, d=5)
zs.d5.median =apply(zs.d5, MARGIN=1, median)
lambda=median(zs.d5.median^2)/0.4549364
lambda
   # [1] 2.888479
lambda=median(zs.d5.median^2)/0.85
lambda
   # [1] 1.034839
adj.p.values.d5 =pchisq(zs.d5.median^2/0.4549364, df=1, lower=F)

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
[CHN.hist]:https://cloud.githubusercontent.com/assets/12142475/23806921/e38370da-05c3-11e7-8ae0-90025261e203.png



Select all the candidates from the full list of SNPs

```
#read the SNP names into R

locus.names <- read.table("CHN/CHN.229.8951.plink.map", header=F)
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


Prepare the data to incorporate with the Fst outlier graph
```
colnames(d1.names) <- "names"
colnames(d2.names) <- "names"
colnames(d3.names) <- "names"
colnames(d4.names) <- "names"
colnames(d5.names) <- "names"

d12345.names <- rbind(d1.names, d2.names, d3.names, d4.names, d5.names)  ##Join all data.frames by "name" column. This only works of colnames are the same (at least one column name)

d12345.names <- lapply(d12345.names, unique)  #select only the unique rows. 

d12345.names <- sub(":", ".", d12345.names$names) ##replace the ":" in the locus names so that they're in the same format as the Fst and RDA lists

write.table(d12345.names$d12345.names, "LFMM.alloutliers", col.names=F, row.names=F, quote=F)

##linux.
##copy the list over to /Users/alexjvr/2016RADAnalysis/5_SE.MS1/DEC2016_SEonly/SumStats

#Read into R

colnames(LFMM.outliers) <- ("loci")
LFMM.outliers <- as.character(LFMM.outliers$loci)
colnames(LFMM.outliers) <- ("loci")
LFMM.outliers <- as.character(LFMM.outliers$loci)
```
