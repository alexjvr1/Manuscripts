#Bayenv2: CHall data


I want to identify loci associated with environment for the 6 CH datasets. stba is removed from all datasets (due to missing env data),
hence 8 less individuals in CHall, CHS, and CHS.TI.

1. CHall.932.9608

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

#Input files:

https://bitbucket.org/tguenther/bayenv2_public/src/8e4039f64d61?at=default

1. Input file with subset of 1000 SNPs for co-variance matrix (pop structure)

2. Input file with all SNPs for the association analysis

3. ENV input file with normalised environmental parameters. I will first use only temperature.

NB: population order in the input files should all be the same.


##1. Covariance matrix estimation

Use a random set of 1000SNPs for matrix estimation. 

I estimated Fst distributions per locus for each of the six datasets. To make sure I'm excluding loci under selection for the estimation of the neutral covariance matrix, I will remove loci in the upper Fst tail from each dataset:

1. Identify neutral loci (ie. exclude upper Fst tail) and select 1000 at random

2. Select 1000 loci from vcf file

3. Convert vcf to bayenv2 format

####CHall

1. Fst<0.6
```
CHall.FstNeutral <- subset(stats.hier.CHall940.perlocus, Fst<0.6)
summary(CHall.FstNeutral)
CHall940.FstNeutral.names <- row.names(CHall.FstNeutral)  ##get the locus names
CHall940.FstNeutral.names <- gsub("X", "", CHall940.FstNeutral.names)
CHall940.FstNeutral.names <- gsub("\\.", ":", CHall940.FstNeutral.names)
CHall940.FstNeutral.names 
CHall940.FstNeutral.names <- as.data.frame(CHall940.FstNeutral.names )
CHall940.1000loci.names <- CHall940.FstNeutral.names[sample(nrow(CHall940.FstNeutral.names), 1000), ]
write.table(CHall940.1000loci.names, "CHall.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHall.940.9608.newnames.vcf --snps CHall.1000loci.names  --recode --recode-INFO-all --out CHall940.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHall940.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/CHall.spidpopfile

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHall/CHall940.1000



####CHN

1. Fst<0.4
```
CHN.FstNeutral <- subset(stats.hier.CHN229.perlocus, Fst<0.4)
summary(CHN.FstNeutral)
CHN229.FstNeutral.names <- row.names(CHN.FstNeutral)  ##get the locus names
CHN229.FstNeutral.names <- gsub("X", "", CHN229.FstNeutral.names)
CHN229.FstNeutral.names <- gsub("\\.", ":", CHN229.FstNeutral.names)
CHN229.FstNeutral.names 
CHN229.FstNeutral.names <- as.data.frame(CHN229.FstNeutral.names )
CHN229.1000loci.names <- CHN229.FstNeutral.names[sample(nrow(CHN229.FstNeutral.names), 1000), ]
write.table(CHN229.1000loci.names, "CHN229.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHN.229.9608.recode.vcf --snps CHN229.1000loci.names --recode --recode-INFO-all --out CHN229.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHN229.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHN/CHN229.1000


####CHS

1. Fst<0.5
```
CHS.FstNeutral <- subset(stats.hier.CHS283.perlocus, Fst<0.5)
summary(CHS.FstNeutral)
CHS283.FstNeutral.names <- row.names(CHS.FstNeutral)  ##get the locus names
CHS283.FstNeutral.names <- gsub("X", "", CHS283.FstNeutral.names)
CHS283.FstNeutral.names <- gsub("\\.", ":", CHS283.FstNeutral.names)
CHS283.FstNeutral.names 
CHS283.FstNeutral.names <- as.data.frame(CHS283.FstNeutral.names )
CHS283.1000loci.names <- CHS283.FstNeutral.names[sample(nrow(CHS283.FstNeutral.names), 1000), ]
write.table(CHS283.1000loci.names, "CHS283.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.283.9608.recode.vcf --snps CHS283.1000loci.names --recode --recode-INFO-all --out CHS283.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS283.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS/CHS283.1000


####CHS.TI

1. Fst<0.5
```
CHS.TI.FstNeutral <- subset(stats.hier.CHS.TI148.perlocus, Fst<0.5)
summary(CHS.FstNeutral)
CHS.TI148.FstNeutral.names <- row.names(CHS.TI.FstNeutral)  ##get the locus names
CHS.TI148.FstNeutral.names <- gsub("X", "", CHS.TI148.FstNeutral.names)
CHS.TI148.FstNeutral.names <- gsub("\\.", ":", CHS.TI148.FstNeutral.names)
CHS.TI148.FstNeutral.names 
CHS.TI148.FstNeutral.names <- as.data.frame(CHS.TI148.FstNeutral.names )
CHS.TI148.1000loci.names <- CHS.TI148.FstNeutral.names[sample(nrow(CHS.TI148.FstNeutral.names), 1000), ]
write.table(CHS.TI148.1000loci.names, "CHS.TI148.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.TI.148.9608.recode.vcf --snps CHS.TI148.1000loci.names --recode --recode-INFO-all --out CHS.TI148.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS.TI148.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.TI/CHS.TI148.1000


####CHS.VS

1. Fst<0.5
```
CHS.VS.FstNeutral <- subset(stats.hier.CHS.VS135.perlocus, Fst<0.5)
summary(CHS.VS.FstNeutral)
CHS.VS135.FstNeutral.names <- row.names(CHS.TI.FstNeutral)  ##get the locus names
CHS.VS135.FstNeutral.names <- gsub("X", "", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names <- gsub("\\.", ":", CHS.VS135.FstNeutral.names)
CHS.VS135.FstNeutral.names 
CHS.VS135.FstNeutral.names <- as.data.frame(CHS.VS135.FstNeutral.names )
CHS.VS135.1000loci.names <- CHS.VS135.FstNeutral.names[sample(nrow(CHS.VS135.FstNeutral.names), 1000), ]
write.table(CHS.VS135.1000loci.names, "CHS.VS135.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CHS.VS.135.9608.recode.vcf --snps CHS.VS135.1000loci.names --recode --recode-INFO-all --out CHS.VS135.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CHS.VS135.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CHS.VS/CHS.VS135.1000



####CZ

1. Fst<0.4
```
CZ404.FstNeutral <- subset(stats.hier.CZ404.perlocus, Fst<0.4)
summary(CZ404.FstNeutral)
CZ404.FstNeutral.names <- row.names(CZ404.FstNeutral)  ##get the locus names
CZ404.FstNeutral.names <- gsub("X", "", CZ404.FstNeutral.names)
CZ404.FstNeutral.names <- gsub("\\.", ":", CZ404.FstNeutral.names)
CZ404.FstNeutral.names 
CZ404.FstNeutral.names <- as.data.frame(CHS.VS135.FstNeutral.names )
CZ404.1000loci.names <- CZ404.FstNeutral.names[sample(nrow(CZ404.FstNeutral.names), 1000), ]
write.table(CZ404.1000loci.names, "CZ404.1000loci.names", col.names=F, quote=F, row.names=F, sep=" ")
```

2. Subset vcf file
```
vcftools --vcf CZ.404.9608.recode.vcf --snps CZ404.1000loci.names --recode --recode-INFO-all --out CZ404.1000.Neutral
```

3. Convert to bayenv2 using pgdspider

vcf input: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/SumStats/CZ404.1000.Neutral.recode.vcf

bayenv output: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000Neutral.bayenv.txt

spidfile: /Users/alexjvr/2016RADAnalysis/spidfiles/vcf2bayenv1000.spid

pop def file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/input.files_subsets/

bayenv sample file: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000.pgdspider.popfile

bayenv loci names: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2/CZ/CZ404.1000


