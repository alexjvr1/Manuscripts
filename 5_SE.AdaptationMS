#SE: Adaptation to temperature across latitude

Population structure and adaptation to temperature across latitude in Scandinavia. 

Main aim: 

Genomic patterns across environment. 

1. Population structure

2. Find loci that are likely to be under selection 

3. Random Forest: allele frequency turnover across the environmental gradient. 


#Data

193 samples were de novo assembed using pyRAD

clustering threshold = 94%

minDP = 6x within individual

minDP = 4x between individuals

max sites with Quality <20 = 4


##Filtering: 


Step 1: 

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
and remove all quotes using sed

```
sed 's/"//g' SNPs.remove.Het0.6 > SNPs2.remove.Het0.6
```


Step 2: 
```

```





##Dataset1

This is the dataset filtered for missingnes with which I can calculate summary statistics. 

1. Average heterozygosity

