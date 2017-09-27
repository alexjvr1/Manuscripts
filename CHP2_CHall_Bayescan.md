# Bayescan

Since the EAA analyses didn't yield very convincing results, I've decided to identify outlier loci using Bayescan. 

Website: http://cmpg.unibe.ch/software/BayeScan/

This is available on the Sork computers:

soney@sorklab3.eeb.ucla.edu

signem123

Check with Sorel before running anything: sorel@ucla.edu

more programs in ~/data


## Input files

working directory: /Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/Bayescan

Input files are allele frequencies. They can be produced using pgdspider. 

I need allele frequencies per population since this is the level at which my environmental data is sampled. 

I'll convert vcf files to bayescan input files. 

First I need to make sure I exclude low MAF, as this could bias the results. (see the instruction manual). I'll filter for MAF of 5%


```
vcftools --vcf input.files_subsets/CHall.940.9608.10March.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Bayescan/CHall.940.9078
vcftools --vcf input.files_subsets/CHN.229.9608.recode.vcf	--maf 0.05 --recode --recode-INFO-all --out Bayescan/CHN.229.5960
vcftools --vcf input.files_subsets/CHS.283.9608.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Bayescan/CHS.283.7099
vcftools --vcf input.files_subsets/CHS.TI.148.9608.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Bayescan/CHS.TI.148.6298
vcftools --vcf input.files_subsets/CHS.VS.135.9608.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Bayescan/CHS.VS.135.6593
vcftools --vcf input.files_subsets/CZ.404.9608.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Bayescan/CZ.404.8115
```

Okay. So this means all the datasets have fewer loci than the 9608 loci I used for LFMM and BayEnv2. These low MAF loci probably affect 
the output from these analyses. So I need to check whether I need to rerun these analyses... 

From the LFMM pdf: http://membres-timc.imag.fr/Olivier.Francois/lfmm/faq.htm
    
    I need to filter for MAF of 5%

From the BayEnv2 pdf: Can't find anything specific, but it's sensible to filter for MAF 5%. 

Okay. So I need to refilter these datasets and rerun all the analyses...
