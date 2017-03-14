# SE: Adaptation to temperature across latitude

Population structure and adaptation to temperature across latitude in Scandinavia. 

Main aim: 

Genomic patterns across environment. 

1. Population structure

2. Find loci that are likely to be under selection (with EAA)

3. Random Forest: allele frequency turnover across the environmental gradient. 
		
		- is there a geographic region that is particularly important?
	
		- Environmental parameter that is particularly important?


# Data

193 samples were de novo assembed using pyRAD

clustering threshold = 94%

minDP = 6x within individual

minDP = 4x between individuals

max sites with Quality <20 = 4


## Filtering: 


#### Step 1: 

Remove loci with >0.5 missingness across all individuals: 

```
vcftools --vcf 

```





Remove all loci with excess Heterozygosity
```
vcftools --vcf SEFinalc94d6m4p3.vcf --plink --out SEFinal.plink

plink --file SEFinal.plink --hardy
```

Read into R. Sort and subset based on markers >0.6 Heterozygosity. 
```
SE.hwe <- read.table("plink.hwe", header=T)

SE.remove <- subset(SE.hwe, O.HET>0.6)

SE.remove.SNPs <- SE.remove$SNP

write.table(SE.remove.SNPs, "SNPs.remove.Het0.6")
```
and remove all quotes using sed. And split the file into 2 columns. (V1=locus, V2=position)

```
sed 's/"//g' SNPs.remove.Het0.6 > SNPs2.remove.Het0.6

awk -F":" '$1=$1' OFS="\t" SNPs2.remove.Het0.6 > SNP3  ##can't use this method. vcf doesn't recognise these names
```

Remove these SNPs from the plink file

```
plink --file SEFinal.plink --exclude SNPs2.remove.Het0.6 --recode --recodeA --out SEFinal.s1


Web-based version check ( --noweb to skip )
Recent cached web-check found... OK, v1.07 is current

+++ PLINK 1.9 is now available! See above website for details +++ 

Writing this text to log file [ SEFinal.s1.log ]
Analysis started: Fri Nov  4 17:11:32 2016

Options in effect:
	--file SEFinal.plink
	--exclude SNPs2.remove.Het0.6
	--recode
	--recodeA
	--out SEFinal.s1

** For gPLINK compatibility, do not use '.' in --out **
1263471 (of 1263471) markers to be included from [ SEFinal.plink.map ]
Warning, found 193 individuals with ambiguous sex codes
Writing list of these individuals to [ SEFinal.s1.nosex ]
193 individuals read from [ SEFinal.plink.ped ] 
0 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
0 cases, 0 controls and 193 missing
0 males, 0 females, and 193 of unspecified sex
Reading list of SNPs to exclude [ SNPs2.remove.Het0.6 ] ... 20897 read
Before frequency and genotyping pruning, there are 1242574 SNPs
193 founders and 0 non-founders found
16 heterozygous haploid genotypes; set to missing
Writing list of heterozygous haploid genotypes to [ SEFinal.s1.hh ]
1 SNPs with no founder genotypes observed
Warning, MAF set to 0 for these SNPs (see --nonfounders)
Writing list of these SNPs to [ SEFinal.s1.nof ]
Total genotyping rate in remaining individuals is 0.178094
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 1242574 SNPs
After filtering, 0 cases, 0 controls and 193 missing
After filtering, 0 males, 0 females, and 193 of unspecified sex
Writing recoded ped file to [ SEFinal.s1.ped ] 
Writing new map file to [ SEFinal.s1.map ] 

```


And convert to vcf using pgdspider
```
java -Xmx2024m -Xms512m -jar /usr/local/PGDSpider_2.0.7.2/PGDSpider2-cli.jar -inputfile SEFinal.s1.ped -inputformat PED -outputfile SEFinal.s1.vcf -outputformat VCF -spid ped2vcf.spid 
```
I had to increase the memory allowance here. And remember to include the path to file.spid in the spid file. Even if it is in the working directory. 

The conversion ran for >1 hour. (not sure how long. I left it running overnight)


#### step 2: rename

```


```




#### Step 3: Missingness

Check the missingness across individuals
```
--vcf SEFinalc94d6m4p3.vcf --exclude SNPs2.remove.Het0.6 --recode --recode-INFO-all --out SEFinal.het0.6.s1.recode.vcf

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





## Dataset1

This is the dataset filtered for missingnes with which I can calculate summary statistics. 

1. Average heterozygosity

