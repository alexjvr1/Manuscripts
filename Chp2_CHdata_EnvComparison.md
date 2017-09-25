# Env Comparison

Aim: to test whether the environmental gradients are significatntly different between the different regions in Switzerland. 

1. CHN, CHS, CZ

2. CHS.TI, CHS.VS


I will use a linear model to test the whether region or elevation have a significant effect on each environmental variable. 

I'm running both a straighforward linear model, and a model including the squared term of elevation for possible quadratic (i.e. non-linear) 
relationships between env variables and elevation. 


# 2. CHS.TI and CHS.VS

/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/EnvData.Differences

To test whether elevation or region are important predictors of each env. variable, I ran a linear model: 


```
CHS.VS.TI.EnvData <- read.table("CHS.TI.VS.EnvData", header=T) ##read in all the EnvData for CHS.VS and CHS.TI

CHS.VS.TI.EnvData$elev.c <- CHS.VS.TI.EnvData$elev - mean(CHS.VS.TI.EnvData$elev) ##center the data
CHS.VS.TI.EnvData$elev.c <- CHS.VS.TI.EnvData$elev.c/1000 ##divide by 1000 to work elevation in km - this avoids the values becoming too large in the model

m1.sol.rad <- lm(CHS.VS.TI.EnvData$sol.rad.60d~CHS.VS.TI.EnvData$elev.c + CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

m1.temp.laying.date <- lm(CHS.VS.TI.EnvData$temp.laying.date ~CHS.VS.TI.EnvData$elev.c + CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

m1.pcpt.60d <- lm(CHS.VS.TI.EnvData$pcpt.60d ~CHS.VS.TI.EnvData$elev.c + CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

m1.shadow.days <- lm(CHS.VS.TI.EnvData$shadow.days ~CHS.VS.TI.EnvData$elev.c + CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

m1.day10cm <- lm(CHS.VS.TI.EnvData$day10cm ~CHS.VS.TI.EnvData$elev.c + CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

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

[CHS.VS.TI.Fig1]:https://user-images.githubusercontent.com/12142475/30785632-c3c04428-a161-11e7-95c4-da2de253d356.png


Then we can interpret the models further like this: 

https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R
```
summary(m1.sol.rad)

Call:
lm(formula = CHS.VS.TI.EnvData$sol.rad.60d ~ CHS.VS.TI.EnvData$elev.c + 
    CHS.VS.TI.EnvData$Transect + CHS.VS.TI.EnvData$Transect:CHS.VS.TI.EnvData$elev.c)

Residuals:
   Min     1Q Median     3Q    Max 
-46010 -18365  -2422  15000  49656 

Coefficients:
                                                          Estimate Std. Error
(Intercept)                                                  99379       5068
CHS.VS.TI.EnvData$elev.c                                     25188       6781
CHS.VS.TI.EnvData$TransectCHS.TI                             -7208       8531
CHS.VS.TI.EnvData$TransectCHS.VS                             10070       9847
CHS.VS.TI.EnvData$elev.c:CHS.VS.TI.EnvData$TransectCHS.TI    -3508      11323
CHS.VS.TI.EnvData$elev.c:CHS.VS.TI.EnvData$TransectCHS.VS    -3426      13348
                                                          t value Pr(>|t|)    
(Intercept)                                                19.607  < 2e-16 ***
CHS.VS.TI.EnvData$elev.c                                    3.715 0.000595 ***
CHS.VS.TI.EnvData$TransectCHS.TI                           -0.845 0.402930    
CHS.VS.TI.EnvData$TransectCHS.VS                            1.023 0.312337    
CHS.VS.TI.EnvData$elev.c:CHS.VS.TI.EnvData$TransectCHS.TI  -0.310 0.758253    
CHS.VS.TI.EnvData$elev.c:CHS.VS.TI.EnvData$TransectCHS.VS  -0.257 0.798690    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 24830 on 42 degrees of freedom
Multiple R-squared:  0.418,	Adjusted R-squared:  0.3487 
F-statistic: 6.032 on 5 and 42 DF,  p-value: 0.0002714
```

The coefficient - estimate of the intercept is where the predictor variable is 0. I.e. for elevation this would be sol radiation at sea level. But since we've centered and scaled our data, this is actually the sol rad at the median elevation for this dataset. 

```
median(CHS.VS.TI.EnvData$elev)

#1618.5m
```

The second line (coefficient - estimate) gives the slope, ie. the effect that elev has on sol.rad. For every km increase in elevation, sol.rad changes by x. 

coeffiecient - std error gives the amount that the coefficient varies from the actual average response variable. 

t-value is the number of std deviations the coefficient is away from 0. Far away from 0 = reject null hypothesis that predictor variable has no effect on the response variable. 

