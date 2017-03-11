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





