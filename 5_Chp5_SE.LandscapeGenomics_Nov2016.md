#SE Landscape genomics



SNP filters: 

1. Filter for 50% missingness and MAF of 0.05 (as recommended my Roesti et al). 

```
vcftools --vcf SEsubset.vcf --max-missing 0.5 --maf 0.05 --recode --recode-INFO-all --out subset.Filter_14Nov2016/SEsubset.Nov.vcf.s2

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SEsubset.vcf
	--recode-INFO-all
	--maf 0.05
	--max-missing 0.5
	--out subset.Filter_14Nov2016/SEsubset.Nov.vcf.s2
	--recode

Eighth Header entry should be INFO: INFO    
After filtering, kept 177 out of 177 Individuals
Outputting VCF file...
After filtering, kept 10326 out of a possible 410602 Sites
Run Time = 25.00 seconds

```

Check the missingness: 

```
vcftools --vcf SEsubset.Nov.vcf.s2.recode.vcf --missing-indv

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

```

![alt_txt][Fig1]
[Fig1]:https://cloud.githubusercontent.com/assets/12142475/20273129/6a5aba12-aa90-11e6-83fa-015ec87310ac.png


rename
```
bcftools query -l SEsubset.Nov.s2.names.vcf
```


2. site-pi per window

Check the frequency of SNPs in the sequence: 

Count the SNPs in windows. 
```
vcftools --vcf SEsubset.Nov.vcf.s2.recode.vcf --window-pi 10
```

Import into R and plot: 
```
SE.10pi <- read.table("out.windowed.pi", header=T)
SE.10pisort <- SE.10pi[order(BIN_START),]  #sort by bin
SE.10pisort$BIN_START <- factor(SE.10pisort$BIN_START, levels=SE.10pisort$BIN_START)

qplot(BIN_START, N_VARIANTS, data=SE.10pisort, geom=c("boxplot", "jitter"), main="Nr SNPs in 10bp sliding window: SE dataset")
qplot(BIN_START, PI, data=SE.10pisort, geom=c("boxplot", "jitter"), main="Nucleotide diversity SNPs in 10bp sliding window: SE dataset")
```

![alt_txt][Fig2]
[Fig2]:https://cloud.githubusercontent.com/assets/12142475/20315078/ca502278-ab5c-11e6-950b-c76421030338.png

test for significant difference: 

My data have a poisson distribution, and are count data. There is also not a significant difference in the mean and SD between classes (with a superficial look at the summary stats)
```
start <- factor(BIN_START)

with(SE.10pisort, tapply(N_VARIANTS, start, function(x) {sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))}))

                     1                     11                     21 
"M (SD) = 1.16 (0.40)" "M (SD) = 1.14 (0.40)" "M (SD) = 1.12 (0.34)" 
                    31                     41                     51 
"M (SD) = 1.11 (0.33)" "M (SD) = 1.13 (0.36)" "M (SD) = 1.14 (0.39)" 
                    61                     71                     81 
"M (SD) = 1.14 (0.41)" "M (SD) = 1.14 (0.39)" "M (SD) = 1.13 (0.37)" 
                    91                    101                    111 
"M (SD) = 1.12 (0.37)" "M (SD) = 1.16 (0.43)" "M (SD) = 1.08 (0.27)" 

ggplot(SE.10pisort, aes(N_VARIANTS, fill=start))+ geom_histogram(binwidth=0.5, position="dodge")

```

![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/20316362/c72a86f6-ab61-11e6-8496-86a974ab15b8.png



Useful website: http://www.ats.ucla.edu/stat/r/dae/poissonreg.htm


2. Filter loci with Linkage of r2 > 0.8 in >10 pops (>60% of populations)
```
vcftools --vcf SEsubset.Nov.vcf.s2.recode.vcf --plink --out SEsubset.Nov.plink

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SEsubset.Nov.vcf.s2.recode.vcf
	--out SEsubset.Nov.plink
	--plink

After filtering, kept 177 out of 177 Individuals
Writing PLINK PED and MAP files ... 
	PLINK: Only outputting biallelic loci.
Done.
After filtering, kept 10326 out of a possible 10326 Sites
Run Time = 1.00 seconds

```




Calculate r2 per population

1. subset the plink file to all the populations. the --keep file should have two columns, duplicated in this case (FamID, IndID)
```
plink --file SE.subset.Nov.plink --keep Upp_O --recode --recodeA --out Upp_O.plink

```

Calculate r2 in all populations
```
plink --file Upp_O.plink --r2 --out Upp_O.ld
```

Read into R. Combine all loci with r2 >0.8. If locus occurs > 10 times, write the name to a new file. 
```

```




