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




