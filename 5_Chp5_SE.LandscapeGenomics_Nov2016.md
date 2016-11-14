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





2. Filter loci with Linkage of r2 > 0.8 in >10 pops (>60% of populations)
```


```
