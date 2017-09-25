# Env Comparison

Aim: to test whether the environmental gradients are significatntly different between the different regions in Switzerland. 

1. CHN, CHS, CZ

2. CHS.TI, CHS.VS


I will use a linear model to test the whether region or elevation have a significant effect on each environmental variable. 

I'm running both a straighforward linear model, and a model including the squared term of elevation for possible quadratic (i.e. non-linear) 
relationships between env variables and elevation. 



# CHN, CHS, and CZ


```
CHN.CHS.CZ.EnvData <- read.table("CHall.CHN.CHS.CZ.EnvData", header=T) ##read in all the EnvData for CHS.VS and CHS.TI

CHN.CHS.CZ.EnvData$elev.c <- CHN.CHS.CZ.EnvData$elev - mean(CHN.CHS.CZ.EnvData$elev) ##center the data
CHN.CHS.CZ.EnvData$elev.c <- CHN.CHS.CZ.EnvData$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

CHN.CHS.CZonly.EnvData <- CHN.CHS.CZ.EnvData[!(CHN.CHS.CZ.EnvData$Transect %in% ("CHall")),] ##remove CHS, so that we include only VS and TI in the model

CHN.CHS.CZonly.EnvData$elev.c.squared <- CHN.CHS.CZonly.EnvData$elev.c^2
```


```
m1.CHall.sol.rad <- lm(CHN.CHS.CZonly.EnvData$sol.rad.60d~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

m1.CHall.temp.laying.date <- lm(CHN.CHS.CZonly.EnvData$temp.laying.date ~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

m1.CHall.pcpt.60d <- lm(CHN.CHS.CZonly.EnvData$pcpt.60d ~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

m1.CHall.shadow.days <- lm(CHN.CHS.CZonly.EnvData$shadow.days ~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

m1.CHall.day10cm <- lm(CHN.CHS.CZonly.EnvData$day10cm ~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

m2.CHall.sol.rad.60d <- lm(CHN.CHS.CZonly.EnvData$sol.rad.60d~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

m2.CHall.temp.laying.date <- lm(CHN.CHS.CZonly.EnvData$temp.laying.date~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

m2.CHall.pcpt.60d <- lm(CHN.CHS.CZonly.EnvData$pcpt.60d~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

m2.CHall.shadow.days <- lm(CHN.CHS.CZonly.EnvData$shadow.day~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

m2.CHall.day10cm <- lm(CHN.CHS.CZonly.EnvData$day10cm~CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


```

First I can look at the distribution of the residuals from both sets of plots: 

m1 

![alt_txt][CHall.residuals.m1]

[CHall.residuals.m1]:https://user-images.githubusercontent.com/12142475/30804231-86147fa8-a1e4-11e7-8f1c-b540b64b8615.png


m2

![alt_txt][CHall.residuals.m2]

[CHall.residuals.m2]:https://user-images.githubusercontent.com/12142475/30804235-8629de3e-a1e4-11e7-9802-87bbd775f398.png


Then I can assess the model fit for m1 and m2

![alt_txt][CHall.sol.rad]

[CHall.sol.rad]:https://user-images.githubusercontent.com/12142475/30804230-86140492-a1e4-11e7-816e-12377ea7a590.png


![alt_txt][CHall.temp.laying.date]

[CHall.temp.laying.date]:https://user-images.githubusercontent.com/12142475/30804232-8616adc8-a1e4-11e7-8ac7-e909ca4ef037.png


![alt_txt][CHall.pcpt.60d]

[CHall.pcpt.60d]:https://user-images.githubusercontent.com/12142475/30804233-8616a60c-a1e4-11e7-9407-3e569a1f5c02.png


![alt_txt][CHall.shadow.days]

[CHall.shadow.days]:https://user-images.githubusercontent.com/12142475/30804229-8612c35c-a1e4-11e7-9452-d751d13ccf0b.png


![alt_txt][CHall.day10cm]

[CHall.day10cm]:https://user-images.githubusercontent.com/12142475/30804234-861d4192-a1e4-11e7-86dd-09118e848365.png



# 2. CHS.TI and CHS.VS

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/EnvData.Differences

To test whether elevation or region are important predictors of each env. variable, I ran a linear model: 


```
CHS.VS.TI.EnvData <- read.table("CHS.TI.VS.EnvData", header=T) ##read in all the EnvData for CHS.VS and CHS.TI

CHS.VS.TI.EnvData$elev.c <- CHS.VS.TI.EnvData$elev - mean(CHS.VS.TI.EnvData$elev) ##center the data
CHS.VS.TI.EnvData$elev.c <- CHS.VS.TI.EnvData$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

CHS.VS.TIonly.EnvData <- CHS.VS.TI.EnvData[!(CHS.VS.TI.EnvData$Transect %in% ("CHS")),] ##remove CHS, so that we include only VS and TI in the model

m1.sol.rad <- lm(CHS.VS.TIonly.EnvData$sol.rad.60d~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

m1.temp.laying.date <- lm(CHS.VS.TIonly.EnvData$temp.laying.date ~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

m1.pcpt.60d <- lm(CHS.VS.TIonly.EnvData$pcpt.60d ~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

m1.shadow.days <- lm(CHS.VS.TIonly.EnvData$shadow.days ~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

m1.day10cm <- lm(CHS.VS.TIonly.EnvData$day10cm ~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

```


To see whether the model fits, I can plot the distribution of the residuals for each model. Residuals are the difference between the model 
prediction and the actual values. We want a normal distribution around the mean. 

```
par(mfrow=c(3,2))
hist(m1.sol.rad$residuals)
hist(m1.pcpt.60d$residuals)
hist(m1.shadow.days$residuals)
hist(m1.temp.laying.date$residuals)
hist(m1.day10cm$residuals)
```

![alt_txt][CHS.VS.TI.Fig1]

[CHS.VS.TI.Fig1]:https://user-images.githubusercontent.com/12142475/30802217-aa282464-a1dd-11e7-9752-71c04b31a456.png


Then we can interpret the models further like this: 

https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R

OR 

https://rpubs.com/FelipeRego/MultipleLinearRegressionInRFirstSteps


#### sol.rad.60d
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$sol.rad.60d ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

Residuals:
   Min     1Q Median     3Q    Max 
-46010 -13727  -2863  15217  41472 

Coefficients:
                                                                  Estimate
(Intercept)                                                       92170.43
CHS.VS.TIonly.EnvData$elev.c                                      21679.91
CHS.VS.TIonly.EnvData$TransectCHS.VS                              17277.76
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS    81.54
                                                                  Std. Error
(Intercept)                                                          6813.60
CHS.VS.TIonly.EnvData$elev.c                                         9003.41
CHS.VS.TIonly.EnvData$TransectCHS.VS                                10802.23
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   14539.52
                                                                  t value
(Intercept)                                                        13.527
CHS.VS.TIonly.EnvData$elev.c                                        2.408
CHS.VS.TIonly.EnvData$TransectCHS.VS                                1.599
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.006
                                                                  Pr(>|t|)    
(Intercept)                                                       1.59e-11 ***
CHS.VS.TIonly.EnvData$elev.c                                        0.0258 *  
CHS.VS.TIonly.EnvData$TransectCHS.VS                                0.1254    
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.9956    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 24650 on 20 degrees of freedom
Multiple R-squared:  0.4535,	Adjusted R-squared:  0.3715 
F-statistic: 5.532 on 3 and 20 DF,  p-value: 0.006232
```

The coefficient - estimate of the intercept is where the predictor variable is 0. I.e. for elevation this would be sol radiation at sea level. But since we've centered and scaled our data, this is actually the sol rad at the median elevation for this dataset. 

```
median(CHS.VS.TI.EnvData$elev)

#1618.5m
```




#### temp.laying.date
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$temp.laying.date ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.9136 -0.8514 -0.1457  1.0332  2.1894 

Coefficients:
                                                                  Estimate
(Intercept)                                                         6.4204
CHS.VS.TIonly.EnvData$elev.c                                       -0.0854
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -1.0874
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   1.1649
                                                                  Std. Error
(Intercept)                                                           0.3635
CHS.VS.TIonly.EnvData$elev.c                                          0.4803
CHS.VS.TIonly.EnvData$TransectCHS.VS                                  0.5763
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS     0.7757
                                                                  t value
(Intercept)                                                        17.662
CHS.VS.TIonly.EnvData$elev.c                                       -0.178
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -1.887
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   1.502
                                                                  Pr(>|t|)    
(Intercept)                                                       1.14e-13 ***
CHS.VS.TIonly.EnvData$elev.c                                        0.8607    
CHS.VS.TIonly.EnvData$TransectCHS.VS                                0.0738 .  
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.1488    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.315 on 20 degrees of freedom
Multiple R-squared:  0.2126,	Adjusted R-squared:  0.09448 
F-statistic:   1.8 on 3 and 20 DF,  p-value: 0.1797
```


#### pcpt.60d

```
Call:
lm(formula = CHS.VS.TIonly.EnvData$pcpt.60d ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-6.6579 -3.3380 -0.6279  3.1537  8.9068 

Coefficients:
                                                                  Estimate
(Intercept)                                                         35.729
CHS.VS.TIonly.EnvData$elev.c                                        -1.574
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -16.093
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS    7.668
                                                                  Std. Error
(Intercept)                                                            1.336
CHS.VS.TIonly.EnvData$elev.c                                           1.766
CHS.VS.TIonly.EnvData$TransectCHS.VS                                   2.118
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS      2.851
                                                                  t value
(Intercept)                                                        26.738
CHS.VS.TIonly.EnvData$elev.c                                       -0.892
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -7.596
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   2.689
                                                                  Pr(>|t|)    
(Intercept)                                                        < 2e-16 ***
CHS.VS.TIonly.EnvData$elev.c                                        0.3832    
CHS.VS.TIonly.EnvData$TransectCHS.VS                              2.57e-07 ***
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.0141 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.835 on 20 degrees of freedom
Multiple R-squared:  0.7576,	Adjusted R-squared:  0.7212 
F-statistic: 20.84 on 3 and 20 DF,  p-value: 2.288e-06

```


#### shadow.days
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$shadow.days ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

Residuals:
   Min     1Q Median     3Q    Max 
-7.705 -3.710 -1.469  2.418 14.644 

Coefficients:
                                                                  Estimate
(Intercept)                                                        37.5892
CHS.VS.TIonly.EnvData$elev.c                                        1.5577
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -6.1537
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.1297
                                                                  Std. Error
(Intercept)                                                           1.7241
CHS.VS.TIonly.EnvData$elev.c                                          2.2782
CHS.VS.TIonly.EnvData$TransectCHS.VS                                  2.7334
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS     3.6791
                                                                  t value
(Intercept)                                                        21.802
CHS.VS.TIonly.EnvData$elev.c                                        0.684
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -2.251
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.035
                                                                  Pr(>|t|)    
(Intercept)                                                       2.07e-15 ***
CHS.VS.TIonly.EnvData$elev.c                                        0.5020    
CHS.VS.TIonly.EnvData$TransectCHS.VS                                0.0358 *  
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.9722    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 6.239 on 20 degrees of freedom
Multiple R-squared:  0.2056,	Adjusted R-squared:  0.08646 
F-statistic: 1.726 on 3 and 20 DF,  p-value: 0.1939
```


#### day10cm
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$day10cm ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-29.971  -7.752   0.251   7.025  36.700 

Coefficients:
                                                                  Estimate
(Intercept)                                                        121.534
CHS.VS.TIonly.EnvData$elev.c                                        89.175
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -22.244
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS    2.452
                                                                  Std. Error
(Intercept)                                                            4.840
CHS.VS.TIonly.EnvData$elev.c                                           6.395
CHS.VS.TIonly.EnvData$TransectCHS.VS                                   7.673
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS     10.328
                                                                  t value
(Intercept)                                                        25.111
CHS.VS.TIonly.EnvData$elev.c                                       13.944
CHS.VS.TIonly.EnvData$TransectCHS.VS                               -2.899
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.237
                                                                  Pr(>|t|)    
(Intercept)                                                        < 2e-16 ***
CHS.VS.TIonly.EnvData$elev.c                                      9.18e-12 ***
CHS.VS.TIonly.EnvData$TransectCHS.VS                               0.00887 ** 
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS  0.81473    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 17.51 on 20 degrees of freedom
Multiple R-squared:  0.9428,	Adjusted R-squared:  0.9342 
F-statistic: 109.8 on 3 and 20 DF,  p-value: 1.359e-12
```




The second line (coefficient - estimate) gives the slope, ie. the effect that elev has on sol.rad. For every km increase in elevation, sol.rad changes by x. 

coeffiecient - std error gives the amount that the coefficient varies from the actual average response variable. 

t-value is the number of std deviations the coefficient is away from 0. Far away from 0 = reject null hypothesis that predictor variable has no effect on the response variable. 



And then we can try to improve the model by including square terms of elevation: 
```
m2.sol.rad.60d <- lm(CHS.VS.TIonly.EnvData$sol.rad.60d~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

m2.temp.laying.date <- lm(CHS.VS.TIonly.EnvData$temp.laying.date~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

m2.pcpt.60d <- lm(CHS.VS.TIonly.EnvData$pcpt.60d~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

m2.shadow.days <- lm(CHS.VS.TIonly.EnvData$shadow.day~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

m2.day10cm <- lm(CHS.VS.TIonly.EnvData$day10cm~CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

```

