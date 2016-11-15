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




