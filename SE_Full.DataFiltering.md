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


```









