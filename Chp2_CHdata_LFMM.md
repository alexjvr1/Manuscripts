#LFMM analyses of CHdata

I want to identify loci associated with environment for the 6 CH datasets

1. CHall.940.9608

2. CHN.229.9608

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
env <- read.csv("BIOclim_allindivs_SEsubset.csv.env", header=F) ###read in the environmental data
write.env(env, "SE.env")   ##convert to correct .env format
[1] "SE.env"
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




And run lfmm at the chosen K. This should be run on the server. 
```

project = lfmm("SEsubset.Final2.geno", "SE.env", K = 2, repetitions = 5, project = "new")  ##run LFMM with this .env file

```



######Postprocessing: 

Based on the lambda calculated from run 1, I have decided to run 3 K for each dataset. Lambda should be close to or just below 1. Mine were all just above 1. 

stba had to be removed from CHS.TI, because of missing data in the env. variables. So this is also removed from CHS. 



First calculate the lambda for each K to determine what the correct K is: 

#####CHN
snmf K = 4

Copy everything to the mac. And remember to change all the paths in the .lfmmProject file

delimiter can be anything. I'm using # here because I have / in my filenames. 
```
sed -i -e 's#old/path#new/path#g' file.lfmmProject
```


```
zs3.1 = z.scores(project, K = 3, d=1)
zs.stouffer3.1 = apply(zs3.1, MARGIN = 1, median)
lambda3.1 = median(zs3.1^2)/.456  
lambda3.1  ##note down lambda for the different K

# calculate adjusted p-values
cp.values3.1 = pchisq(zs.stouffer3.1^2/lambda3.1, df = 1, lower = FALSE)

zs3.2 = z.scores(project, K = 3, d=2)
zs.stouffer3.2 = apply(zs3.2, MARGIN = 1, median)
lambda3.2 = median(zs3.2^2)/.456  
lambda3.2  ##note down lambda for the different K

# calculate adjusted p-values
cp.values3.2 = pchisq(zs.stouffer3.2^2/lambda3.2, df = 1, lower = FALSE)

zs3.3 = z.scores(project, K = 3, d=3)
zs.stouffer3.3 = apply(zs3.3, MARGIN = 1, median)
lambda3.3 = median(zs3.3^2)/.456  
lambda3.3  ##note down lambda for the different K

# calculate adjusted p-values
cp.values3.3 = pchisq(zs.stouffer3.3^2/lambda3.3, df = 1, lower = FALSE)

zs3.4 = z.scores(project, K = 3, d=4)
zs.stouffer3.4 = apply(zs3.4, MARGIN = 1, median)
lambda3.4 = median(zs3.4^2)/.456  
lambda3.4  ##note down lambda for the different K
# calculate adjusted p-values
cp.values3.4 = pchisq(zs.stouffer3.4^2/lambda3.4, df = 1, lower = FALSE)

zs3.5 = z.scores(project, K = 3, d=5)
zs.stouffer3.5 = apply(zs3.5, MARGIN = 1, median)
lambda3.5 = median(zs3.5^2)/.456  
lambda3.5  ##note down lambda for the different K
# calculate adjusted p-values
cp.values3.5 = pchisq(zs.stouffer3.5^2/lambda3.5, df = 1, lower = FALSE)


zs4.1 = z.scores(project, K = 4, d=1)
zs.stouffer4.1 = apply(zs4.1, MARGIN = 1, median)
lambda4.1 = median(zs4.1^2)/.456  
lambda4.1  ##note down lambda for the different K

# calculate adjusted p-values
cp.values4.1 = pchisq(zs.stouffer4.1^2/lambda4.1, df = 1, lower = FALSE)


zs4.2 = z.scores(project, K = 4, d=2)
zs.stouffer4.2 = apply(zs4.2, MARGIN = 1, median)
lambda4.2 = median(zs4.2^2)/.456  
lambda4.2  ##note down lambda for the different K2

# calculate adjusted p-values
cp.values4.2 = pchisq(zs.stouffer4.2^2/lambda4.2, df = 1, lower = FALSE)

zs4.3 = z.scores(project, K = 4, d=3)
zs.stouffer4.3 = apply(zs4.3, MARGIN = 1, median)
lambda4.3 = median(zs4.3^2)/.456  
lambda4.3  ##note down lambda for the different K

# calculate adjusted p-values
cp.values4.3 = pchisq(zs.stouffer4.3^2/lambda4.3, df = 1, lower = FALSE)

zs4.4 = z.scores(project, K = 4, d=4)
zs.stouffer4.4 = apply(zs4.4, MARGIN = 1, median)
lambda4.4 = median(zs4.4^2)/.456  
lambda4.4  ##note down lambda for the different K

# calculate adjusted p-values
cp.values4.4 = pchisq(zs.stouffer4.4^2/lambda4.4, df = 1, lower = FALSE)

zs4.5 = z.scores(project, K = 4, d=5)
zs.stouffer4.5 = apply(zs4.5, MARGIN = 1, median)
lambda4.5 = median(zs4.5^2)/.456  
lambda4.5  ##note down lambda for the different K

# calculate adjusted p-values
cp.values4.5 = pchisq(zs.stouffer4.5^2/lambda4.5, df = 1, lower = FALSE)

zs5.1 = z.scores(project, K = 5, d=1)
zs.stouffer5.1 = apply(zs5.1, MARGIN = 1, median)
lambda5.1 = median(zs5.1^2)/.456  
lambda5.1  ##note down lambda for the different K

# calculate adjusted p-values
cp.values5.1 = pchisq(zs.stouffer5.1^2/lambda5.1, df = 1, lower = FALSE)

zs5.2 = z.scores(project, K = 5, d=2)
zs.stouffer5.2 = apply(zs5.2, MARGIN = 1, median)
lambda5.2 = median(zs5.2^2)/.456  
lambda5.2  ##note down lambda for the different K

# calculate adjusted p-values
cp.values5.2 = pchisq(zs.stouffer5.2^2/lambda5.2, df = 1, lower = FALSE)

zs5.3 = z.scores(project, K = 5, d=3)
zs.stouffer5.3 = apply(zs5.3, MARGIN = 1, median)
lambda5.3 = median(zs5.3^2)/.456  
lambda5.3  ##note down lambda for the different K

# calculate adjusted p-values
cp.values5.3 = pchisq(zs.stouffer5.3^2/lambda5.3, df = 1, lower = FALSE)

zs5.4 = z.scores(project, K = 5, d=4)
zs.stouffer5.4 = apply(zs5.4, MARGIN = 1, median)
lambda5.4 = median(zs5.4^2)/.456  
lambda5.4  ##note down lambda for the different K

# calculate adjusted p-values
cp.values5.4 = pchisq(zs.stouffer5.4^2/lambda5.4, df = 1, lower = FALSE)

zs5.5 = z.scores(project, K = 5, d=5)
zs.stouffer5.5 = apply(zs5.5, MARGIN = 1, median)
lambda5.5 = median(zs5.5^2)/.456  
lambda5.5  ##note down lambda for the different K

# calculate adjusted p-values
cp.values5.5 = pchisq(zs.stouffer5.5^2/lambda5.5, df = 1, lower = FALSE)

zs6.1 = z.scores(project, K = 6, d=1)
zs.stouffer6.1 = apply(zs6.1, MARGIN = 1, median)
lambda6.1 = median(zs6.1^2)/.456  
lambda  ##note down lambda for the different K

# calculate adjusted p-values
cp.values6.1 = pchisq(zs.stouffer6.1^2/lambda6.1, df = 1, lower = FALSE)

zs6.2 = z.scores(project, K = 6, d=2)
zs.stouffer6.2 = apply(zs6.2, MARGIN = 1, median)
lambda6.2 = median(zs6.2^2)/.456  
lambda  ##note down lambda for the different K

# calculate adjusted p-values
cp.values6.2 = pchisq(zs.stouffer6.2^2/lambda6.2, df = 1, lower = FALSE)

zs6.3 = z.scores(project, K = 6, d=3)
zs.stouffer6.3 = apply(zs6.3, MARGIN = 1, median)
lambda6.3 = median(zs6.3^2)/.456  
lambda  ##note down lambda for the different K

# calculate adjusted p-values
cp.values6.3 = pchisq(zs.stouffer6.3^2/lambda6.3, df = 1, lower = FALSE)

zs6.4 = z.scores(project, K = 6, d=4)
zs.stouffer6.4 = apply(zs6.4, MARGIN = 1, median)
lambda6.4 = median(zs6.4^2)/.456  
lambda  ##note down lambda for the different K
# calculate adjusted p-values
cp.values6.4 = pchisq(zs.stouffer6.4^2/lambda6.4, df = 1, lower = FALSE)

zs6.5 = z.scores(project, K = 6, d=5)
zs.stouffer6.5 = apply(zs6.5, MARGIN = 1, median)
lambda6.5 = median(zs6.5^2)/.456  
lambda  ##note down lambda for the different K

# calculate adjusted p-values
cp.values6.5 = pchisq(zs.stouffer6.5^2/lambda6.5, df = 1, lower = FALSE)



par(mfrow=c(5,3))

hist(cp.values4.1)
hist(cp.values4.2)
hist(cp.values4.3)
hist(cp.values4.4)
hist(cp.values4.5)
hist(cp.values5.1)
hist(cp.values5.2)
hist(cp.values5.3)
hist(cp.values5.4)
hist(cp.values5.5)
hist(cp.values6.1)
hist(cp.values6.2)
hist(cp.values6.3)
hist(cp.values6.4)
hist(cp.values6.5)

#FDR correction
q = 0.05 
L = length(cp.values4.1) 
w = which(sort(cp.values4.1) < q * (1:L) / L) 
candidates4.1 = order(cp.values4.1)[w]

##print out candidates

```
LFMM lists the loci 1-nLoci

To get the real name of the locus, I need to look at the .loci file output from pgdspider when the BayENV2 input is created.

Copy the .loci file to the LFMM subfolder. 

Number the lines in this file: 
```
nl CHS.VS.pgdspider_loci.txt > CHS.VS.8779.loci
```

And select the specific lines that were identified as outliers
```
sed -n "$(sed 's/$/p/' numbers.n1)" CHS.VS.8779.loci
```


