# SE Landscape genomics



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

![Fig1](https://cloud.githubusercontent.com/assets/12142475/20273129/6a5aba12-aa90-11e6-83fa-015ec87310ac.png)


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

![Fig2](https://cloud.githubusercontent.com/assets/12142475/20315078/ca502278-ab5c-11e6-950b-c76421030338.png)

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

![Fig3](https://cloud.githubusercontent.com/assets/12142475/20316362/c72a86f6-ab61-11e6-8496-86a974ab15b8.png)

```
summary(m1 <- glm(N_VARIANTS ~start, family="poisson", data=SE.10pisort))

Call:
glm(formula = N_VARIANTS ~ start, family = "poisson", data = SE.10pisort)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-0.1517  -0.1373  -0.1272  -0.1135   2.1052  

Coefficients:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept)  0.129844   0.032827   3.955 7.64e-05 ***
start11     -0.017893   0.048053  -0.372    0.710    
start21      0.008153   0.046012   0.177    0.859    
start31     -0.005337   0.046601  -0.115    0.909    
start41     -0.000453   0.045988  -0.010    0.992    
start51     -0.010669   0.046797  -0.228    0.820    
start61     -0.024848   0.046614  -0.533    0.594    
start71      0.018140   0.046576   0.389    0.697    
start81      0.007686   0.046850   0.164    0.870    
start91     -0.018483   0.047149  -0.392    0.695    
start101     0.004406   0.045636   0.097    0.923    
start111    -0.030998   0.070705  -0.438    0.661    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 922.71  on 9104  degrees of freedom
Residual deviance: 920.93  on 9093  degrees of freedom
AIC: 19874

Number of Fisher Scoring iterations: 4

```

None of the levels are significant. 

```
cov.m1 <- vcovHC(m1, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(Estimate= coef(m1), "Robust SE" = std.err,
"Pr(>|z|)" = 2 * pnorm(abs(coef(m1)/std.err), lower.tail=FALSE),
LL = coef(m1) - 1.96 * std.err,
UL = coef(m1) + 1.96 * std.err)

r.est

                 Estimate  Robust SE     Pr(>|z|)          LL          UL
(Intercept)  0.1298436195 0.01186972 7.497351e-28  0.10657896 0.153108279
start11     -0.0178932942 0.01682772 2.876357e-01 -0.05087563 0.015089038
start21      0.0081527306 0.01689110 6.293341e-01 -0.02495382 0.041259283
start31     -0.0053367164 0.01675469 7.500900e-01 -0.03817591 0.027502480
start41     -0.0004530196 0.01639421 9.779549e-01 -0.03258568 0.031679638
start51     -0.0106691825 0.01694585 5.289534e-01 -0.04388305 0.022544686
start61     -0.0248481340 0.01611111 1.230011e-01 -0.05642590 0.006729634
start71      0.0181397997 0.01780911 3.084078e-01 -0.01676606 0.053045655
start81      0.0076861790 0.01697815 6.507576e-01 -0.02559100 0.040963356
start91     -0.0184832135 0.01630198 2.568777e-01 -0.05043509 0.013468666
start101     0.0044059060 0.01659129 7.905817e-01 -0.02811303 0.036924842
start111    -0.0309977849 0.02308998 1.794411e-01 -0.07625414 0.014258572
```

And do a ChiSq test ti see if the model fit. 
We can use the residual deviance to perform a goodness of fit test for the overall model. The residual deviance is the difference between the deviance of the current model and the maximum deviance of the ideal model where the predicted values are identical to the observed. Therefore, if the residual difference is small enough, the goodness of fit test will not be significant, indicating that the model fits the data. 
H0 = model is close to ideal model
```
with(m1, cbind(res.deviance = deviance, df = df.residual,
+   p = pchisq(deviance, df.residual, lower.tail=FALSE)))
     res.deviance   df p
[1,]     920.9293 9093 1
```
Not rejected. i.e. model fits. 



Useful website: http://www.ats.ucla.edu/stat/r/dae/poissonreg.htm






4. HWE

Convert to plink and recodeA to get all the right files
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
Subset the plink file to all the populations. the --keep file should have two columns, duplicated in this case (FamID, IndID)
```
plink --file SE.subset.Nov.plink --keep Upp_O --recode --recodeA --out Upp_O.plink
```

And keep only SNPs that are variable in the population (because otherwise the A1 & A2 columns have missing data, and it is impossible to read it into R)
```
plink --file DE_B.plink --maf 0.05 --recode --recodeA --out DE.B.SNPonly
plink --file DE_W.plink --maf 0.05 --recode --recodeA --out DE.W.SNPonly
plink --file DE_K.plink --maf 0.05 --recode --recodeA --out DE.K.SNPonly

plink --file Sk_Ho.plink --maf 0.05 --recode --recodeA --out Sk_Ho.SNPonly
plink --file Sk_SF.plink --maf 0.05 --recode --recodeA --out Sk_SF.SNPonly
plink --file Sk_SL.plink --maf 0.05 --recode --recodeA --out Sk_SL.SNPonly

plink --file Upp_Gra.plink --maf 0.05 --recode --recodeA --out Upp_Gra.SNPonly
plink --file Upp_K.plink --maf 0.05 --recode --recodeA --out Upp_K.SNPonly
plink --file Upp_O.plink --maf 0.05 --recode --recodeA --out Upp_O.SNPonly

plink --file Um_UT3.plink --maf 0.05 --recode --recodeA --out Um_UT3.SNPonly
plink --file Um_Taf.plink --maf 0.05 --recode --recodeA --out Um_Taf.SNPonly
plink --file Um_Gr.plink --maf 0.05 --recode --recodeA --out Um_Gr.SNPonly

plink --file LT1.plink --maf 0.05 --recode --recodeA --out LT1.SNPonly
plink --file LT2.plink --maf 0.05 --recode --recodeA --out LT2.SNPonly
plink --file LT3.plink --maf 0.05 --recode --recodeA --out LT3.SNPonly

plink --file Kir_L.plink --maf 0.05 --recode --recodeA --out Kir_L.SNPonly
plink --file Kir_G.plink --maf 0.05 --recode --recodeA --out Kir_G.SNPonly

plink --file FIN.plink --maf 0.05 --recode --recodeA --out FIN.SNPonly
```


Calculate HWE within each pop
```
plink --file DE.B.SNPonly --hardy --out DE.B
plink --file DE.W.SNPonly --hardy --out DE.W
plink --file DE.K.SNPonly --hardy --out DE.K

plink --file Sk_Ho.SNPonly --hardy --out Sk_Ho
plink --file Sk_SL.SNPonly --hardy --out Sk_SL
plink --file Sk_SF.SNPonly --hardy --out Sk_SF

plink --file Upp_Gra.SNPonly --hardy --out Upp_Gra
plink --file Upp_K.SNPonly --hardy --out Upp_K
plink --file Upp_O.SNPonly --hardy --out Upp_O

plink --file Um_UT3.SNPonly --hardy --out Um_UT3
plink --file Um_Taf.SNPonly --hardy --out Um_Taf
plink --file Um_Gr.SNPonly --hardy --out Um_Gr

plink --file LT1.SNPonly --hardy --out LT1
plink --file LT2.SNPonly --hardy --out LT2
plink --file LT3.SNPonly --hardy --out LT3

plink --file Kir_G.SNPonly --hardy --out Kir_G
plink --file Kir_L.SNPonly --hardy --out Kir_L

plink --file FIN.SNPonly --hardy --out FIN
```

Read into R and 
```
DE.B.hwe <- read.table("DE.B.hwe", header=T)
DE.K.hwe <- read.table("DE.K.hwe", header=T)
DE.W.hwe <- read.table("DE.W.hwe", header=T)

Sk.Ho.hwe <- read.table("Sk_Ho.hwe", header=T)
Sk.SF.hwe <- read.table("Sk_SF.hwe", header=T)
Sk.SL.hwe <- read.table("Sk_SL.hwe", header=T)

Upp.Gra.hwe <- read.table("Upp_Gra.hwe", header=T)
Upp.K.hwe <- read.table("Upp_K.hwe", header=T)
Upp.O.hwe <- read.table("Upp_O.hwe", header=T)

Um.Gr.hwe <- read.table("Um_Gr.hwe", header=T)
Um.UT3.hwe <- read.table("Um_UT3.hwe", header=T)
Um.Taf.hwe <- read.table("Um_Taf.hwe", header=T)

LT1.hwe <- read.table("LT1.hwe", header=T)
LT2.hwe <- read.table("LT2.hwe", header=T)
LT3.hwe <- read.table("LT3.hwe", header=T)

Kir.G.hwe <- read.table("Kir_G.hwe", header=T)
Kir.L.hwe <- read.table("Kir_L.hwe", header=T)

FIN.hwe <- read.table("FIN.hwe", header=T)
```

Keep only the Test=All columns
```
DE.B.hwe <- subset(DE.B.hwe, TEST=="ALL")
DE.W.hwe <- subset(DE.W.hwe, TEST=="ALL")
DE.K.hwe <- subset(DE.K.hwe, TEST=="ALL")
Sk.Ho.hwe <- subset(Sk.Ho.hwe, TEST=="ALL")
Sk.SF.hwe <- subset(Sk.SF.hwe, TEST=="ALL")
Sk.SL.hwe <- subset(Sk.SL.hwe, TEST=="ALL")
Upp.Gra.hwe <- subset(Upp.Gra.hwe, TEST=="ALL")
Upp.K.hwe <- subset(Upp.K.hwe, TEST=="ALL")
Upp.O.hwe <- subset(Upp.O.hwe, TEST=="ALL")
Um.Gr.hwe <- subset(Um.Gr.hwe, TEST=="ALL")
Um.UT3.hwe <- subset(Um.UT3.hwe, TEST=="ALL")
Um.Taf.hwe <- subset(Um.Taf.hwe, TEST=="ALL")
LT1.hwe <- subset(LT1.hwe, TEST=="ALL")
LT2.hwe <- subset(LT2.hwe, TEST=="ALL")
LT3.hwe <- subset(LT3.hwe, TEST=="ALL")
Kir.G.hwe <- subset(Kir..hwe, TEST=="ALL")
Kir.G.hwe <- subset(Kir.G.hwe, TEST=="ALL")
Kir.L.hwe <- subset(Kir.L.hwe, TEST=="ALL")
FIN.hwe <- subset(FIN.hwe, TEST=="ALL")
```

Split into loci out of HWE + O.Het > 0.5, and those to keep
```
DE.B.hwe.remove <- subset(DE.B.hwe, P<0.050001)
DE.B.hwe.remove <- rbind(subset(DE.B.hwe, O.HET.>0.5))
DE.B.hwe.keep <- subset(DE.B.hwe, O.HET. <0.51)
DE.B.hwe.keep <- subset(DE.B.hwe.keep, P>0.05)

DE.W.hwe.remove <- subset(DE.W.hwe, P<0.050001)
DE.W.hwe.remove <- rbind(subset(DE.W.hwe, O.HET.>0.5))
DE.W.hwe.keep <- subset(DE.W.hwe, O.HET. <0.51)
DE.W.hwe.keep <- subset(DE.W.hwe.keep, P>0.05)

DE.K.hwe.remove <- subset(DE.K.hwe, P<0.050001)
DE.K.hwe.remove <- rbind(subset(DE.K.hwe, O.HET.>0.5))
DE.K.hwe.keep <- subset(DE.K.hwe, O.HET. <0.51)
DE.K.hwe.keep <- subset(DE.K.hwe.keep, P>0.05)

Sk.Ho.hwe.remove <- subset(Sk.Ho.hwe, P<0.050001)
Sk.Ho.hwe.remove <- rbind(subset(Sk.Ho.hwe, O.HET.>0.5))
Sk.Ho.hwe.keep <- subset(Sk.Ho.hwe, O.HET. <0.51)
Sk.Ho.hwe.keep <- subset(Sk.Ho.hwe.keep, P>0.05)

Sk.SF.hwe.remove <- subset(Sk.SF.hwe, P<0.050001)
Sk.SF.hwe.remove <- rbind(subset(Sk.SF.hwe, O.HET.>0.5))
Sk.SF.hwe.keep <- subset(Sk.SF.hwe, O.HET. <0.51)
Sk.SF.hwe.keep <- subset(Sk.SF.hwe.keep, P>0.05)

Sk.SL.hwe.remove <- subset(Sk.SL.hwe, P<0.050001)
Sk.SL.hwe.remove <- rbind(subset(Sk.SL.hwe, O.HET.>0.5))
Sk.SL.hwe.keep <- subset(Sk.SL.hwe, O.HET. <0.51)
Sk.SL.hwe.keep <- subset(Sk.SL.hwe.keep, P>0.05)

Upp.Gra.hwe.remove <- subset(Upp.Gra.hwe, P<0.050001)
Upp.Gra.hwe.remove <- rbind(subset(Upp.Gra.hwe, O.HET.>0.5))
Upp.Gra.hwe.keep <- subset(Upp.Gra.hwe, O.HET. <0.51)
Upp.Gra.hwe.keep <- subset(Upp.Gra.hwe.keep, P>0.05)

Upp.K.hwe.remove <- subset(Upp.K.hwe, P<0.050001)
Upp.K.hwe.remove <- rbind(subset(Upp.K.hwe, O.HET.>0.5))
Upp.K.hwe.keep <- subset(Upp.K.hwe, O.HET. <0.51)
Upp.K.hwe.keep <- subset(Upp.K.hwe.keep, P>0.05)

Upp.O.hwe.remove <- subset(Upp.O.hwe, P<0.050001)
Upp.O.hwe.remove <- rbind(subset(Upp.O.hwe, O.HET.>0.5))
Upp.O.hwe.keep <- subset(Upp.O.hwe, O.HET. <0.51)
Upp.O.hwe.keep <- subset(Upp.O.hwe.keep, P>0.05)

Um.Gr.hwe.remove <- subset(Um.Gr.hwe, P<0.050001)
Um.Gr.hwe.remove <- rbind(subset(Um.Gr.hwe, O.HET.>0.5))
Um.Gr.hwe.keep <- subset(Um.Gr.hwe, O.HET. <0.51)
Um.Gr.hwe.keep <- subset(Um.Gr.hwe.keep, P>0.05)

Um.Taf.hwe.remove <- subset(Um.Taf.hwe, P<0.050001)
Um.Taf.hwe.remove <- rbind(subset(Um.Taf.hwe, O.HET.>0.5))
Um.Taf.hwe.keep <- subset(Um.Taf.hwe, O.HET. <0.51)
Um.Taf.hwe.keep <- subset(Um.Taf.hwe.keep, P>0.05)

Um.UT3.hwe.remove <- subset(Um.UT3.hwe, P<0.050001)
Um.UT3.hwe.remove <- rbind(subset(Um.UT3.hwe, O.HET.>0.5))
Um.UT3.hwe.keep <- subset(Um.UT3.hwe, O.HET. <0.51)
Um.UT3.hwe.keep <- subset(Um.UT3.hwe.keep, P>0.05)

LT1.hwe.remove <- subset(LT1.hwe, P<0.050001)
LT1.hwe.remove <- rbind(subset(LT1.hwe, O.HET.>0.5))
LT1.hwe.keep <- subset(LT1.hwe, O.HET. <0.51)
LT1.hwe.keep <- subset(LT1.hwe.keep, P>0.05)

LT2.hwe.remove <- subset(LT2.hwe, P<0.050001)
LT2.hwe.remove <- rbind(subset(LT2.hwe, O.HET.>0.5))
LT2.hwe.keep <- subset(LT2.hwe, O.HET. <0.51)
LT2.hwe.keep <- subset(LT2.hwe.keep, P>0.05)

LT3.hwe.remove <- subset(LT3.hwe, P<0.050001)
LT3.hwe.remove <- rbind(subset(LT3.hwe, O.HET.>0.5))
LT3.hwe.keep <- subset(LT3.hwe, O.HET. <0.51)
LT3.hwe.keep <- subset(LT3.hwe.keep, P>0.05)

Kir.L.hwe.remove <- subset(Kir.L.hwe, P<0.050001)
Kir.L.hwe.remove <- rbind(subset(Kir.L.hwe, O.HET.>0.5))
Kir.L.hwe.keep <- subset(Kir.L.hwe, O.HET. <0.51)
Kir.L.hwe.keep <- subset(Kir.L.hwe.keep, P>0.05)

Kir.G.hwe.remove <- subset(Kir.G.hwe, P<0.050001)
Kir.G.hwe.remove <- rbind(subset(Kir.G.hwe, O.HET.>0.5))
Kir.G.hwe.keep <- subset(Kir.G.hwe, O.HET. <0.51)
Kir.G.hwe.keep <- subset(Kir.G.hwe.keep, P>0.05)

FIN.hwe.remove <- subset(FIN.hwe, P<0.050001)
FIN.hwe.remove <- rbind(subset(FIN.hwe, O.HET.>0.5))
FIN.hwe.keep <- subset(FIN.hwe, O.HET. <0.51)
FIN.hwe.keep <- subset(FIN.hwe.keep, P>0.05)
```

Bind all files together and plot the frequency of SNP names
```
SE.hwe.remove.file <- do.call(rbind, lapply(ls(pattern="remove$"), get))

test.table <- data.frame(table(SE.hwe.remove.file$SNP))
hist(test.table$Freq, xlab="Number of pops in which locus occurs", ylab="Frequency", main="Frequency of SNPs deviating from HWE and H.Obs>0.5")
```

![HWE.SNPs](https://cloud.githubusercontent.com/assets/12142475/20355442/a0eedf80-ac21-11e6-8966-a6368d69f7e9.png)


Sturgeon paper removes loci that deviate from HWE in >60% of the population. I will us 5pops as the cutoff
```
HWE.loci.remove.names <- subset(test.table, Freq>4)
write.table(HWE.loci.remove.names$Var1, quote=F, col.names=F, file="HWE.loci.remove.names", row.names=F)
```

And remove the 425 loci using plink 
```
 plink --file SEsubset.Nov.plink --exclude HWE.loci.remove.names --recode --recodeA --out SEsubset.Nov.s3

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Web-based version check ( --noweb to skip )
Recent cached web-check found... OK, v1.07 is current

+++ PLINK 1.9 is now available! See above website for details +++ 

Writing this text to log file [ SEsubset.Nov.s3.log ]
Analysis started: Wed Nov 16 17:49:55 2016

Options in effect:
	--file SEsubset.Nov.plink
	--exclude HWE.loci.remove.names
	--recode
	--recodeA
	--out SEsubset.Nov.s3

** For gPLINK compatibility, do not use '.' in --out **
10227 (of 10227) markers to be included from [ SEsubset.Nov.plink.map ]
Warning, found 177 individuals with ambiguous sex codes
Writing list of these individuals to [ SEsubset.Nov.s3.nosex ]
177 individuals read from [ SEsubset.Nov.plink.ped ] 
0 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
0 cases, 0 controls and 177 missing
0 males, 0 females, and 177 of unspecified sex
Reading list of SNPs to exclude [ HWE.loci.remove.names ] ... 425 read
Before frequency and genotyping pruning, there are 9802 SNPs
177 founders and 0 non-founders found
12 heterozygous haploid genotypes; set to missing
Writing list of heterozygous haploid genotypes to [ SEsubset.Nov.s3.hh ]
Total genotyping rate in remaining individuals is 0.683502
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 9802 SNPs
After filtering, 0 cases, 0 controls and 177 missing
After filtering, 0 males, 0 females, and 177 of unspecified sex
Writing recoded ped file to [ SEsubset.Nov.s3.ped ] 
Writing new map file to [ SEsubset.Nov.s3.map ]  
```


See how many SNPs are found in all pops:
Bind all files together and plot the frequency of SNP names
```
SE.hwe.keep.file <- do.call(rbind, lapply(ls(pattern="keep$"), get))

test.table.keep <- data.frame(table(SE.hwe.keep.file$SNP))
hist(test.table.keep$Freq, xlab="Number of pops", ylab="Frequency", main="Frequency of loci found in increasing number of populations")
```

![HWE.SNPtokeep](https://cloud.githubusercontent.com/assets/12142475/20355441/a0ed449a-ac21-11e6-8d28-a99956ed553e.png)





## 5. Filter loci with Linkage of r2 > 0.8 in >10 pops (>60% of populations)


Subset the plink file with the HWE s3 filter and look only at variable loci. 
```
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/DE_B --maf 0.00005 --recode --recodeA --out DE_B.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/DE_W --maf 0.00005 --recode --recodeA --out DE_W.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/DE_K --maf 0.00005 --recode --recodeA --out DE_K.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Sk_SF --maf 0.00005 --recode --recodeA --out Sk_SF.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Sk_SL --maf 0.00005 --recode --recodeA --out Sk_SL.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Sk_Ho --maf 0.00005 --recode --recodeA --out Sk_Ho.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Upp_O --maf 0.00005 --recode --recodeA --out Upp_O.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Upp_Gra --maf 0.00005 --recode --recodeA --out Upp_Gra.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Upp_K --maf 0.00005 --recode --recodeA --out Upp_K.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/LT1 --maf 0.00005 --recode --recodeA --out LT1.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/LT2 --maf 0.00005 --recode --recodeA --out LT2.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/LT3 --maf 0.00005 --recode --recodeA --out LT3.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Um_UT3 --maf 0.00005 --recode --recodeA --out Um_UT3.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Um_Taf --maf 0.00005 --recode --recodeA --out Um_Taf.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Um_Gr --maf 0.00005 --recode --recodeA --out Um_Gr.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Kir_L --maf 0.00005 --recode --recodeA --out Kir_L.s3.plink
plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/Kir_G --maf 0.00005 --recode --recodeA --out Kir_G.s3.plink

plink --file SEsubset.Nov.s3 --keep popnames.plink.folder/FIN --maf 0.00005 --recode --recodeA --out FIN.s3.plink
```



Calculate r2 per population for variable loci only

```
plink --file DE_B.s3.plink --r2 --out DE_B.ld
plink --file DE_K.s3.plink --r2 --out DE_W.ld
plink --file DE_W.s3.plink --r2 --out DE_K.ld

plink --file Sk_SF.s3.plink --r2 --out Sk_SF.ld
plink --file Sk_SL.s3.plink --r2 --out Sk_SL.ld
plink --file Sk_Ho.s3.plink --r2 --out Sk_Ho.ld

plink --file Upp_Gra.s3.plink --r2 --out Upp_Gr.ld
plink --file Upp_K.s3.plink --r2 --out Upp_K.ld
plink --file Upp_O.s3.plink --r2 --out Upp_O.ld

plink --file Um_Gr.s3.plink --r2 --out Um_Gr.ld
plink --file Um_UT3.s3.plink --r2 --out Um_UT3.ld
plink --file Um_Taf.s3.plink --r2 --out Um_Taf.ld

plink --file LT1.s3.plink --r2 --out LT1.ld
plink --file LT2.s3.plink --r2 --out LT2.ld
plink --file LT3.s3.plink --r2 --out LT3.ld

plink --file Kir_L.s3.plink --r2 --out Kir_L.ld
plink --file Kir_G.s3.plink --r2 --out Kir_G.ld

plink --file FIN.s3.plink --r2 --out FIN.ld
```

Read into R. 
```
DE.B.ld <- read.table("DE_B.ld.ld", header=T)
DE.K.ld <- read.table("DE_K.ld.ld", header=T)
DE.W.ld <- read.table("DE_W.ld.ld", header=T)

Sk.Ho.ld <- read.table("Sk_Ho.ld.ld", header=T)
Sk.SF.ld <- read.table("Sk_SF.ld.ld", header=T)
Sk.SL.ld <- read.table("Sk_SL.ld.ld", header=T)

Upp.Gra.ld <- read.table("Upp_Gr.ld.ld", header=T)
Upp.K.ld <- read.table("Upp_K.ld.ld", header=T)
Upp.O.ld <- read.table("Upp_O.ld.ld", header=T)

Um.Gr.ld <- read.table("Um_Gr.ld.ld", header=T)
Um.UT3.ld <- read.table("Um_UT3.ld.ld", header=T)
Um.Taf.ld <- read.table("Um_Taf.ld.ld", header=T)

LT1.ld <- read.table("LT1.ld.ld", header=T)
LT2.ld <- read.table("LT2.ld.ld", header=T)
LT3.ld <- read.table("LT3.ld.ld", header=T)

Kir.G.ld <- read.table("Kir_G.ld.ld", header=T)
Kir.L.ld <- read.table("Kir_L.ld.ld", header=T)

FIN.ld <- read.table("FIN.ld.ld", header=T)
```



Combine all loci with r2 >0.8. If locus occurs > 10 times, write the name to a new file. 
```
DE.B.ld.remove <- subset(DE.B.ld, R2>0.79999)
DE.B.ld.keep <- subset(DE.B.ld, R2<0.8)

DE.W.ld.remove <- subset(DE.W.ld, R2>0.79999)
DE.W.ld.keep <- subset(DE.W.ld, R2<0.8)

DE.K.ld.remove <- subset(DE.K.ld, R2>0.79999)
DE.K.ld.keep <- subset(DE.K.ld, R2<0.8)

Sk.Ho.ld.remove <- subset(Sk.Ho.ld, R2>0.79999)
Sk.Ho.ld.keep <- subset(Sk.Ho.ld, R2<0.8)

Sk.SF.ld.remove <- subset(Sk.SF.ld, R2>0.79999)
Sk.SF.ld.keep <- subset(Sk.SF.ld, R2<0.8)

Sk.SL.ld.remove <- subset(Sk.SL.ld, R2>0.79999)
Sk.SL.ld.keep <- subset(Sk.SL.ld, R2<0.8)

Upp.Gra.ld.remove <- subset(Upp.Gra.ld, R2>0.79999)
Upp.Gra.ld.keep <- subset(Upp.Gra.ld, R2<0.8)

Upp.K.ld.remove <- subset(Upp.K.ld, R2>0.79999)
Upp.K.ld.keep <- subset(Upp.K.ld, R2<0.8)

Upp.O.ld.remove <- subset(Upp.O.ld, R2>0.79999)
Upp.O.ld.keep <- subset(Upp.O.ld, R2<0.8)

LT1.ld.remove <- subset(LT1.ld, R2>0.79999)
LT1.ld.keep <- subset(LT1.ld, R2<0.8)

LT2.ld.remove <- subset(LT2.ld, R2>0.79999)
LT2.ld.keep <- subset(LT2.ld, R2<0.8)

LT3.ld.remove <- subset(LT3.ld, R2>0.79999)
LT3.ld.keep <- subset(LT3.ld, R2<0.8)

Um_UT3.ld.remove <- subset(Um.UT3.ld, R2>0.79999)
Um_UT3.ld.keep <- subset(Um.UT3.ld, R2<0.8)

Um_Taf.ld.remove <- subset(Um.Taf.ld, R2>0.79999)
Um_Taf.ld.keep <- subset(Um.Taf.ld, R2<0.8)

Um_Gr.ld.remove <- subset(Um.Gr.ld, R2>0.79999)
Um_Gr.ld.keep <- subset(Um.Gr.ld, R2<0.8)

Kir.G.ld.remove <- subset(Kir.G.ld, R2>0.79999)
Kir.G.ld.keep <- subset(Kir.G.ld, R2<0.8)

Kir.L.ld.remove <- subset(Kir.L.ld, R2>0.79999)
Kir.L.ld.keep <- subset(Kir.L.ld, R2<0.8)

FIN.ld.remove <- subset(FIN.ld, R2>0.79999)
FIN.ld.keep <- subset(FIN.ld, R2<0.8)
```

Bind all files together and plot the frequency of SNP names
```
SE.ld.remove.file <- do.call(rbind, lapply(ls(pattern="ld.remove$"), get))

test.table.ld <- data.frame(table(SE.ld.remove.file$SNP_A, SE.ld.remove.file$SNP_B))
test.table.ld.subset <- subset(test.table.ld, Freq>0)
hist(test.table.ld.subset$Freq, xlab="Number of pops in which locus occurs", ylab="Frequency", main="Frequency of SNP pairs with R2>0.8")

###and per pop
par(mfrow=c(3, 6))
hist(DE.B.ld$R2, main="DE.B")
hist(DE.K.ld$R2, main="DE.K")
hist(DE.W.ld$R2, main="DE.W")

hist(Sk.Ho.ld$R2, main="Sk.Ho")
hist(Sk.SF.ld$R2, main="Sk.SF")
hist(Sk.SL.ld$R2, main="Sk.SL")

hist(Upp.Gra.ld$R2, main="Upp.Gra")
hist(Upp.K.ld$R2, main="Upp.K")
hist(Upp.O.ld$R2, main="Upp.O")

hist(Um.Gr.ld$R2, main="Um.Gr")
hist(Um.Taf.ld$R2, main="Um.Taf")
hist(Um.UT3.ld$R2, main="Um.UT3")

hist(LT1.ld$R2, main="LT1")
hist(LT2.ld$R2, main="LT2")
hist(LT3.ld$R2, main="LT3")

hist(Kir.G.ld$R2, main="Kir.G")
hist(Kir.L.ld$R2, main="Kir.L")
hist(FIN.ld$R2, main="FIN")

```

![SE.R2](https://cloud.githubusercontent.com/assets/12142475/20363125/622b4168-ac3e-11e6-9b87-d74cb1e74738.png)


![SE.R2.2](https://cloud.githubusercontent.com/assets/12142475/20384698/0faffbaa-acb5-11e6-800b-dcb3eb42802f.png)


### 6. Remove individuals with too much missing data

```
vcftools --vcf SE.Final.3963.vcf --missing-indv
nano pop.imiss

###IN R
SE.imiss <- read.table("out.imiss", header=T)
pop.imiss <- read.table("pop.imiss", header=F)
SE.imiss$pop <- pop.imiss$V1
attach(SE.imiss)
dat.sort <- SE.imiss[order(pop),]
dat.sort$pop <- factor(dat.sort$pop, levels=dat.sort$pop)
indiv.remove.names <- subset(SE.imiss, F_MISS>0.55)
q <- qplot(pop, F_MISS, data=dat.sort, geom=c("boxplot", "jitter"))

##in linux
vcftools --vcf SE.Final.3963.vcf --remove remove.imiss.names --recode --recode-INFO-all --out SE.Final.3963.171
```

![SE.missing](https://cloud.githubusercontent.com/assets/12142475/20384838/82a04174-acb5-11e6-967f-2019a234b18f.png)




Final dataset: 

3963 loci

171 individuals

18 populations

Genotyping rate: 0.701687





### Final checks: 

convert the final dataset to plink and split into populations: 
```
vcftools --vcf SE.Final.3963.171.names.vcf --plink --out SE.Final.3963.171.plink
plink --file SE.Final.3963.171.plink --recode --recodeA

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/DE_B  --recode --recodeA --out FinalData.SE/DE_B.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/DE_W  --recode --recodeA --out FinalData.SE/DE_W.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/DE_K  --recode --recodeA --out FinalData.SE/DE_K.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Sk_SF  --recode --recodeA --out FinalData.SE/Sk_SF.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Sk_SL  --recode --recodeA --out FinalData.SE/Sk_SL.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Sk_Ho  --recode --recodeA --out FinalData.SE/Sk_Ho.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Upp_O  --recode --recodeA --out FinalData.SE/Upp_O.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Upp_Gra  --recode --recodeA --out FinalData.SE/Upp_Gra.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Upp_K  --recode --recodeA --out FinalData.SE/Upp_K.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/LT1  --recode --recodeA --out FinalData.SE/LT1.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/LT2  --recode --recodeA --out FinalData.SE/LT2.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/LT3  --recode --recodeA --out FinalData.SE/LT3.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Um_UT3  --recode --recodeA --out FinalData.SE/Um_UT3.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Um_Taf  --recode --recodeA --out FinalData.SE/Um_Taf.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Um_Gr  --recode --recodeA --out FinalData.SE/Um_Gr.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Kir_L  --recode --recodeA --out FinalData.SE/Kir_L.plink
plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/Kir_G  --recode --recodeA --out FinalData.SE/Kir_G.plink

plink --file SE.Final.3963.171.plink --keep popnames.plink.folder/FIN  --recode --recodeA --out FinalData.SE/FIN.plink
```

Check missingness overall, and MAF spectrum per pop, Nr of variable loci per pop
```
##linux
plink --file DE_K.plink --freq --out DE.K
plink --file DE_W.plink --freq --out DE.W
plink --file DE_B.plink --freq --out DE.B

plink --file Sk_Ho.plink --freq --out Sk_Ho
plink --file Sk_SF.plink --freq --out Sk_SL
plink --file Sk_SL.plink --freq --out Sk_SF

plink --file Upp_Gra.plink --freq --out Upp_Gra
plink --file Upp_K.plink --freq --out Upp_K
plink --file Upp_O.plink --freq --out Upp_O

plink --file Um_Taf.plink --freq --out Um_Taf
plink --file Um_Gr.plink --freq --out Um_Gr

plink --file LT1.plink --freq --out LT1
plink --file LT2.plink --freq --out LT2
plink --file LT3.plink --freq --out LT3

plink --file Kir_G.plink --freq --out Kir_G
plink --file Kir_L.plink --freq --out Kir_L

plink --file FIN.plink --freq --out FIN
```

Read into R and plot
```
DE.B.frq <- read.table("DE.B.frq", header=T)
DE.K.frq <- read.table("DE.K.frq", header=T)
DE.W.frq <- read.table("DE.W.frq", header=T)

Sk.Ho.frq <- read.table("Sk_Ho.frq", header=T)
Sk.SF.frq <- read.table("Sk_SF.frq", header=T)
Sk.SL.frq <- read.table("Sk_SL.frq", header=T)

Upp.Gra.frq <- read.table("Upp_Gra.frq", header=T)
Upp.K.frq <- read.table("Upp_K.frq", header=T)
Upp.O.frq <- read.table("Upp_O.frq", header=T)

Um.Gr.frq <- read.table("Um_Gr.frq", header=T)
Um.Taf.frq <- read.table("Um_Taf.frq", header=T)

LT1.frq <- read.table("LT1.frq", header=T)
LT2.frq <- read.table("LT2.frq", header=T)
LT3.frq <- read.table("LT3.frq", header=T)

Kir.G.frq <- read.table("Kir_G.frq", header=T)
Kir.L.frq <- read.table("Kir_L.frq", header=T)

FIN.frq <- read.table("FIN.frq", header=T)
```

```
par(mfrow=c(3,6))

hist(DE.B.frq$MAF, main="DE.B")
hist(DE.K.frq$MAF, main="DE.K")
hist(DE.W.frq$MAF, main="DE.W")

hist(Sk.Ho.frq$MAF, main="Sk.Ho")
hist(Sk.SF.frq$MAF, main="Sk.SF")
hist(Sk.SL.frq$MAF, main="Sk.SL")

hist(Upp.Gra.frq$MAF, main="Upp.Gra")
hist(Upp.K.frq$MAF, main="Upp.K")
hist(Upp.O.frq$MAF, main="Upp.O")

hist(Um.Gr.frq$MAF, main="Um.Gr")
hist(Um.Taf.frq$MAF, main="Um.Taf")
hist(Um.UT3.frq$MAF, main="Um.UT3")

hist(LT1.frq$MAF, main="LT1")
hist(LT2.frq$MAF, main="LT2")
hist(LT3.frq$MAF, main="LT3")

hist(Kir.G.frq$MAF, main="Kir.G")
hist(Kir.L.frq$MAF, main="Kir.L")
hist(FIN.frq$MAF, main="FIN")

###Add pop column to each .frq data.frame: 

DE.B.frq$pop <- c(rep_len("DE.B", length.out=3963))
DE.K.frq$pop <- c(rep_len("DE.K", length.out=3963))
DE.W.frq$pop <- c(rep_len("DE.W", length.out=3963))
Sk.Ho.frq$pop <- c(rep_len("Sk.Ho", length.out=3963))
Sk.SL.frq$pop <- c(rep_len("Sk.SL", length.out=3963))
Sk.SF.frq$pop <- c(rep_len("Sk.SF", length.out=3963))
Upp.Gra.frq$pop <- c(rep_len("Up.Gra", length.out=3963))
Upp.K.frq$pop <- c(rep_len("Up.K", length.out=3963))
Upp.O.frq$pop <- c(rep_len("Up.O", length.out=3963))
Um.Taf.frq$pop <- c(rep_len("Um.Taf", length.out=3963))
Um.Gr.frq$pop <- c(rep_len("Um.Gr", length.out=3963))
LT1.frq$pop <- c(rep_len("LT1", length.out=3963))
LT2.frq$pop <- c(rep_len("LT2", length.out=3963))
LT3.frq$pop <- c(rep_len("LT3", length.out=3963))
Kir.G.frq$pop <- c(rep_len("Kir.G", length.out=3963))
Kir.L.frq$pop <- c(rep_len("Kir.L", length.out=3963))
FIN.frq$pop <- c(rep_len("FIN", length.out=3963))

##bind all files together
SE.MAF.file <- do.call(rbind, lapply(ls(pattern=“frq$"), get))
pop <- factor(pop)
ggplot(SE.MAF.file, aes(MAF, fill=pop))+ geom_histogram(binwidth=0.05, position="dodge")
```


![maf.perpop.1](https://cloud.githubusercontent.com/assets/12142475/20389651/01358e94-accb-11e6-8872-08406c2bd1d8.png)

![maf.perpop.2](https://cloud.githubusercontent.com/assets/12142475/20389661/06b313a0-accb-11e6-8abb-47c8279af9db.png)


Nr of loci variable in >1 pop
```
SE.MAF.file.keep <- subset(SE.MAF.file, MAF>0.0001)   ##keep variable loci only
test.table.maf <- data.frame(table(SE.MAF.file.keep$SNP)) ##write table of frequencies of SNPs
hist(test.table.maf$Freq, main="Freq of loci variable in x pops")

SE.SNPs.keep.var10pops <- subset(test.table.maf, Freq>9)   ##Loci variable in >9pops
summary(SE.SNPs.keep.var10pops)
         Var1           Freq      
 100043:88 :   1   Min.   :10.00  
 100140:47 :   1   1st Qu.:11.00  
 100233:35 :   1   Median :12.00  
 100277:21 :   1   Mean   :12.53  
 100431:36 :   1   3rd Qu.:14.00  
 100440:100:   1   Max.   :17.00  
 (Other)   :1347                  
SE.SNPs.keep.var5pops <- subset(test.table.maf, Freq>4)
summary(SE.SNPs.keep.var5pops)
        Var1           Freq       
 100043:88:   1   Min.   : 5.000  
 100124:1 :   1   1st Qu.: 7.000  
 100140:47:   1   Median : 9.000  
 100233:35:   1   Mean   : 9.399  
 100261:84:   1   3rd Qu.:12.000  
 100277:21:   1   Max.   :17.000  
 (Other)  :3001   


```

![Variable.loci](https://cloud.githubusercontent.com/assets/12142475/20389920/6e275b30-accc-11e6-9c57-733562108857.png)



```
#Create the datasets

write.table(SE.SNPs.keep.var5pops$Var1, row.names=F, col.names=F, quote=F, "SE.SNPs.keep.var5pops")
write.table(SE.SNPs.keep.var10pops$Var1, row.names=F, col.names=F, quote=F, "SE.SNPs.keep.var10pops")

#linux

vcftools --vcf SE.Final.3963.171.names.vcf --snps SE.SNPs.keep.var5pops --recode --recode-INFO-all --out SE.var5pops.3007.171
vcftools --vcf SE.Final.3963.171.names.vcf --snps SE.SNPs.keep.var10pops --recode --recode-INFO-all --out SE.var10pops.1353.171
```


## RDA

Redundancy Analysis. 

Two input files: 

1. Genotype: Allele frequences for all loci by population

2. Env file: 

WD: /Users/alexjvr/2016RADAnalysis/5_SE.MS1/Analyses_old/summstats/subset.Filter_14Nov2016/SEsubset.version2/FinalData.SE

Genotype file
```
#convert to plink

vcftools --vcf SE.var5pops.3007.171.recode.vcf --plink --out SE.var5pops.3007.171.plink
plink --file SE.var5pops.3007.171.plink --recode --recodeA --out SE.var5pops.3007.171.plink

vcftools --vcf SE.var10pops.1353.171.recode.vcf --plink --out SE.var10pops.1353.171.plink
plink --file SE.var10pops.1353.171.plink --recode --recodeA --out SE.var10pops.1353.171.plink
```

Look at the *nosex* file and create a file with an third column specifying the cluster level (here pop): 


![cluster](https://cloud.githubusercontent.com/assets/12142475/20393131/022f5b26-acdb-11e6-87d4-7b8be12f5b14.png)


Calculate allele frequency within each cluster level for both plink files
```
plink --file SE.var5pops.3007.171.plink --within SE171.Clusters --freq SE.3007.171

plink --file SE.var10pops.1353.171.plink --within SE171.Clusters --freq --out SE.1353.171

```


Import into R to reformat the output - by population and loci as columns

```
SE.var5.maf <- read.table("SE.3007.171.frq.strat", header = T)
head(SE.var5.maf)

SE.var10.maf <- read.table("SE.1353.171.frq.strat", header = T)
head(SE.var10.maf)


SE.var5.maf$CHR = NULL
SE.var5.maf$A1 = NULL
SE.var5.maf$A2 = NULL
SE.var5.maf$MAC = NULL
SE.var5.maf$NCHROBS = NULL
head(SE.var5.maf)

SE.var5.maf <- SE.var5.maf[,c(2,1,3)]


SE.var10.maf$CHR = NULL
SE.var10.maf$A1 = NULL
SE.var10.maf$A2 = NULL
SE.var10.maf$MAC = NULL
SE.var10.maf$NCHROBS = NULL
head(SE.var10.maf)

SE.var10.maf <- SE.var10.maf[,c(2,1,3)]


library("ggplot2")
library("reshape2")

SE.var5.maf2 <- melt(SE.var5.maf, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.var5.maf2)
head(SE.var5.maf2)


SE.var5.maf3 <- dcast(SE.var5.maf2, formula= CLST ~ SNP)
head(SE.var5.maf3)
write.csv(SE.var5.maf3, file="SE.var5.maf.csv")


SE.var10.maf2 <- melt(SE.var10.maf, id.vars = c("CLST", "SNP"), variable_name = c("MAF"))
str(SE.var10.maf2)
head(SE.var10.maf2)


SE.var10.maf3 <- dcast(SE.var10.maf2, formula= CLST ~ SNP)
head(SE.var10.maf3)
write.csv(SE.var10.maf3, file="SE.var10.maf.csv")
```

Open in Excel and add "X" in front of all the locus names. 

```
##RDA in R
library(vegan)


```
