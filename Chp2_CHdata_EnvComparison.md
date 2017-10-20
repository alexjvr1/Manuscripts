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

```
par(mfrow=c(2,4))
plot(m1.CHall.sol.rad)
plot(m2.CHall.sol.rad.60d)
par(mfrow=c(2,4))
plot(m1.CHall.temp.laying.date)
plot(m2.CHall.temp.laying.date)
par(mfrow=c(2,4))
plot(m1.CHall.pcpt.60d)
plot(m2.CHall.pcpt.60d)
par(mfrow=c(2,4))
plot(m1.CHall.shadow.days)
plot(m2.CHall.shadow.days)
par(mfrow=c(2,4))
plot(m1.CHall.day10cm)
plot(m2.CHall.day10cm)
```

sol.rad

![alt_txt][CHall.sol.rad]

[CHall.sol.rad]:https://user-images.githubusercontent.com/12142475/30805298-b6a79714-a1e8-11e7-8705-0c6aec9c1fd3.png



temp.laying.date
![alt_txt][CHall.temp.laying.date]

[CHall.temp.laying.date]:https://user-images.githubusercontent.com/12142475/30805300-b6acbbe0-a1e8-11e7-8d1d-c4741e2e3b71.png


pcpt.60d
![alt_txt][CHall.pcpt.60d]

[CHall.pcpt.60d]:https://user-images.githubusercontent.com/12142475/30805301-b6af554e-a1e8-11e7-9073-420fc16307f0.png


shadow.days
![alt_txt][CHall.shadow.days]

[CHall.shadow.days]:https://user-images.githubusercontent.com/12142475/30805297-b6a74f98-a1e8-11e7-9e0b-e2f6fd1be796.png


day10cm

![alt_txt][CHall.day10cm]

[CHall.day10cm]:https://user-images.githubusercontent.com/12142475/30805299-b6a806ea-a1e8-11e7-8fa1-25b8d6304dae.png


And then assess the results from each model: 

##### sol.rad

m1
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$sol.rad.60d ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

Residuals:
   Min     1Q Median     3Q    Max 
-48507 -21517  -1751  20954  50213 

Coefficients:
                                                                 Estimate
(Intercept)                                                         93173
CHN.CHS.CZonly.EnvData$elev.c                                       18035
CHN.CHS.CZonly.EnvData$TransectCHS                                   5668
CHN.CHS.CZonly.EnvData$TransectCZ                                   12073
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS     7153
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ     13790
                                                                 Std. Error
(Intercept)                                                            6381
CHN.CHS.CZonly.EnvData$elev.c                                         11259
CHN.CHS.CZonly.EnvData$TransectCHS                                     8487
CHN.CHS.CZonly.EnvData$TransectCZ                                      7780
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS      13519
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ       12690
                                                                 t value
(Intercept)                                                       14.602
CHN.CHS.CZonly.EnvData$elev.c                                      1.602
CHN.CHS.CZonly.EnvData$TransectCHS                                 0.668
CHN.CHS.CZonly.EnvData$TransectCZ                                  1.552
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS   0.529
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ    1.087
                                                                 Pr(>|t|)    
(Intercept)                                                        <2e-16 ***
CHN.CHS.CZonly.EnvData$elev.c                                       0.113    
CHN.CHS.CZonly.EnvData$TransectCHS                                  0.506    
CHN.CHS.CZonly.EnvData$TransectCZ                                   0.125    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS    0.598    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ     0.281    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 27400 on 75 degrees of freedom
Multiple R-squared:  0.3867,	Adjusted R-squared:  0.3458 
F-statistic: 9.459 on 5 and 75 DF,  p-value: 5.05e-07
```

m2

```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$sol.rad.60d ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

Residuals:
   Min     1Q Median     3Q    Max 
-51944 -20340    214  16986  51973 

Coefficients:
                                                                         Estimate
(Intercept)                                                                 80539
CHN.CHS.CZonly.EnvData$elev.c                                               27205
CHN.CHS.CZonly.EnvData$elev.c.squared                                       42102
CHN.CHS.CZonly.EnvData$TransectCHS                                          33528
CHN.CHS.CZonly.EnvData$TransectCZ                                           20366
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS            -6542
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ              4748
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   -69159
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ    -34601
                                                                         Std. Error
(Intercept)                                                                    8906
CHN.CHS.CZonly.EnvData$elev.c                                                 11868
CHN.CHS.CZonly.EnvData$elev.c.squared                                         21341
CHN.CHS.CZonly.EnvData$TransectCHS                                            13251
CHN.CHS.CZonly.EnvData$TransectCZ                                             12144
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS              14122
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ               13157
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS      25816
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ       24565
                                                                         t value
(Intercept)                                                                9.044
CHN.CHS.CZonly.EnvData$elev.c                                              2.292
CHN.CHS.CZonly.EnvData$elev.c.squared                                      1.973
CHN.CHS.CZonly.EnvData$TransectCHS                                         2.530
CHN.CHS.CZonly.EnvData$TransectCZ                                          1.677
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -0.463
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ            0.361
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS  -2.679
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   -1.409
                                                                         Pr(>|t|)
(Intercept)                                                              1.73e-13***
CHN.CHS.CZonly.EnvData$elev.c                                             0.02481*
CHN.CHS.CZonly.EnvData$elev.c.squared                                     0.05235.
CHN.CHS.CZonly.EnvData$TransectCHS                                        0.01359*
CHN.CHS.CZonly.EnvData$TransectCZ                                         0.09788.
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          0.64456
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           0.71926
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS  0.00915**
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   0.16327
                                                                            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 26580 on 72 degrees of freedom
Multiple R-squared:  0.4463,	Adjusted R-squared:  0.3847 
F-statistic: 7.253 on 8 and 72 DF,  p-value: 5.162e-07

```


##### temp.laying.date

m1
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$temp.laying.date ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-2.3837 -0.9034 -0.1189  0.8848  3.5392 

Coefficients:
                                                                 Estimate
(Intercept)                                                        5.3995
CHN.CHS.CZonly.EnvData$elev.c                                      1.3849
CHN.CHS.CZonly.EnvData$TransectCHS                                 0.6952
CHN.CHS.CZonly.EnvData$TransectCZ                                 -0.7491
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  -1.2206
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -1.6719
                                                                 Std. Error
(Intercept)                                                          0.3037
CHN.CHS.CZonly.EnvData$elev.c                                        0.5358
CHN.CHS.CZonly.EnvData$TransectCHS                                   0.4039
CHN.CHS.CZonly.EnvData$TransectCZ                                    0.3702
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS     0.6434
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ      0.6039
                                                                 t value
(Intercept)                                                       17.780
CHN.CHS.CZonly.EnvData$elev.c                                      2.585
CHN.CHS.CZonly.EnvData$TransectCHS                                 1.721
CHN.CHS.CZonly.EnvData$TransectCZ                                 -2.023
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  -1.897
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -2.768
                                                                 Pr(>|t|)    
(Intercept)                                                       < 2e-16 ***
CHN.CHS.CZonly.EnvData$elev.c                                     0.01169 *  
CHN.CHS.CZonly.EnvData$TransectCHS                                0.08934 .  
CHN.CHS.CZonly.EnvData$TransectCZ                                 0.04661 *  
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  0.06166 .  
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   0.00709 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.304 on 75 degrees of freedom
Multiple R-squared:  0.2601,	Adjusted R-squared:  0.2107 
F-statistic: 5.272 on 5 and 75 DF,  p-value: 0.0003337
```


m2
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$temp.laying.date ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

Residuals:
    Min      1Q  Median      3Q     Max 
-2.2104 -0.8604 -0.1103  0.6637  3.5441 

Coefficients:
                                                                         Estimate
(Intercept)                                                               5.27907
CHN.CHS.CZonly.EnvData$elev.c                                             1.47226
CHN.CHS.CZonly.EnvData$elev.c.squared                                     0.40115
CHN.CHS.CZonly.EnvData$TransectCHS                                       -0.09071
CHN.CHS.CZonly.EnvData$TransectCZ                                        -0.23416
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS         -1.03869
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ          -1.77093
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS  1.20941
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ  -1.08304
                                                                         Std. Error
(Intercept)                                                                 0.42609
CHN.CHS.CZonly.EnvData$elev.c                                               0.56780
CHN.CHS.CZonly.EnvData$elev.c.squared                                       1.02106
CHN.CHS.CZonly.EnvData$TransectCHS                                          0.63400
CHN.CHS.CZonly.EnvData$TransectCZ                                           0.58106
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS            0.67564
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ             0.62949
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS    1.23517
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ     1.17532
                                                                         t value
(Intercept)                                                               12.390
CHN.CHS.CZonly.EnvData$elev.c                                              2.593
CHN.CHS.CZonly.EnvData$elev.c.squared                                      0.393
CHN.CHS.CZonly.EnvData$TransectCHS                                        -0.143
CHN.CHS.CZonly.EnvData$TransectCZ                                         -0.403
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -1.537
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -2.813
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   0.979
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   -0.921
                                                                         Pr(>|t|)
(Intercept)                                                               < 2e-16***
CHN.CHS.CZonly.EnvData$elev.c                                             0.01152*
CHN.CHS.CZonly.EnvData$elev.c.squared                                     0.69557
CHN.CHS.CZonly.EnvData$TransectCHS                                        0.88663
CHN.CHS.CZonly.EnvData$TransectCZ                                         0.68815
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          0.12860
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           0.00632**
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS  0.33079
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   0.35988
                                                                            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.272 on 72 degrees of freedom
Multiple R-squared:  0.3247,	Adjusted R-squared:  0.2497 
F-statistic: 4.328 on 8 and 72 DF,  p-value: 0.00027
```


##### pcpt.60d

m1
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$pcpt.60d ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

Residuals:
     Min       1Q   Median       3Q      Max 
-19.3829  -3.6412  -0.6273   4.0841  17.3952 

Coefficients:
                                                                 Estimate
(Intercept)                                                        28.366
CHN.CHS.CZonly.EnvData$elev.c                                      12.063
CHN.CHS.CZonly.EnvData$TransectCHS                                  1.556
CHN.CHS.CZonly.EnvData$TransectCZ                                  -2.368
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  -13.791
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ    -3.099
                                                                 Std. Error
(Intercept)                                                           1.641
CHN.CHS.CZonly.EnvData$elev.c                                         2.895
CHN.CHS.CZonly.EnvData$TransectCHS                                    2.182
CHN.CHS.CZonly.EnvData$TransectCZ                                     2.000
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS      3.476
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ       3.262
                                                                 t value
(Intercept)                                                       17.291
CHN.CHS.CZonly.EnvData$elev.c                                      4.167
CHN.CHS.CZonly.EnvData$TransectCHS                                 0.713
CHN.CHS.CZonly.EnvData$TransectCZ                                 -1.184
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  -3.968
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -0.950
                                                                 Pr(>|t|)    
(Intercept)                                                       < 2e-16 ***
CHN.CHS.CZonly.EnvData$elev.c                                    8.17e-05 ***
CHN.CHS.CZonly.EnvData$TransectCHS                               0.477981    
CHN.CHS.CZonly.EnvData$TransectCZ                                0.240178    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS 0.000165 ***
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ  0.345174    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7.045 on 75 degrees of freedom
Multiple R-squared:  0.4341,	Adjusted R-squared:  0.3964 
F-statistic: 11.51 on 5 and 75 DF,  p-value: 2.917e-08
```

m2
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$pcpt.60d ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

Residuals:
     Min       1Q   Median       3Q      Max 
-19.7214  -3.0089   0.2209   3.1759  18.1329 

Coefficients:
                                                                         Estimate
(Intercept)                                                                28.189
CHN.CHS.CZonly.EnvData$elev.c                                              12.191
CHN.CHS.CZonly.EnvData$elev.c.squared                                       0.588
CHN.CHS.CZonly.EnvData$TransectCHS                                          3.046
CHN.CHS.CZonly.EnvData$TransectCZ                                           2.872
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -14.309
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ            -3.377
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   -2.922
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ    -9.339
                                                                         Std. Error
(Intercept)                                                                   2.282
CHN.CHS.CZonly.EnvData$elev.c                                                 3.041
CHN.CHS.CZonly.EnvData$elev.c.squared                                         5.468
CHN.CHS.CZonly.EnvData$TransectCHS                                            3.395
CHN.CHS.CZonly.EnvData$TransectCZ                                             3.112
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS              3.618
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ               3.371
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS      6.614
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ       6.294
                                                                         t value
(Intercept)                                                               12.354
CHN.CHS.CZonly.EnvData$elev.c                                              4.009
CHN.CHS.CZonly.EnvData$elev.c.squared                                      0.108
CHN.CHS.CZonly.EnvData$TransectCHS                                         0.897
CHN.CHS.CZonly.EnvData$TransectCZ                                          0.923
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -3.955
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -1.002
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS  -0.442
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   -1.484
                                                                         Pr(>|t|)
(Intercept)                                                               < 2e-16***
CHN.CHS.CZonly.EnvData$elev.c                                            0.000147***
CHN.CHS.CZonly.EnvData$elev.c.squared                                    0.914659
CHN.CHS.CZonly.EnvData$TransectCHS                                       0.372643
CHN.CHS.CZonly.EnvData$TransectCZ                                        0.359020
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS         0.000177***
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ          0.319808
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS 0.660002
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ  0.142204
                                                                            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 6.809 on 72 degrees of freedom
Multiple R-squared:  0.4925,	Adjusted R-squared:  0.4361 
F-statistic: 8.734 on 8 and 72 DF,  p-value: 2.954e-08
```



##### shadow.days

m1
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$shadow.days ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-11.581  -4.395  -3.477   0.978  25.783 

Coefficients:
                                                                 Estimate
(Intercept)                                                        39.311
CHN.CHS.CZonly.EnvData$elev.c                                       4.727
CHN.CHS.CZonly.EnvData$TransectCHS                                 -4.279
CHN.CHS.CZonly.EnvData$TransectCZ                                  -4.201
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS   -4.355
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ    -5.616
                                                                 Std. Error
(Intercept)                                                           1.977
CHN.CHS.CZonly.EnvData$elev.c                                         3.488
CHN.CHS.CZonly.EnvData$TransectCHS                                    2.630
CHN.CHS.CZonly.EnvData$TransectCZ                                     2.410
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS      4.189
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ       3.932
                                                                 t value
(Intercept)                                                       19.884
CHN.CHS.CZonly.EnvData$elev.c                                      1.355
CHN.CHS.CZonly.EnvData$TransectCHS                                -1.627
CHN.CHS.CZonly.EnvData$TransectCZ                                 -1.743
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS  -1.040
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -1.428
                                                                 Pr(>|t|)    
(Intercept)                                                        <2e-16 ***
CHN.CHS.CZonly.EnvData$elev.c                                      0.1794    
CHN.CHS.CZonly.EnvData$TransectCHS                                 0.1079    
CHN.CHS.CZonly.EnvData$TransectCZ                                  0.0854 .  
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS   0.3018    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ    0.1573    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 8.49 on 75 degrees of freedom
Multiple R-squared:  0.06248,	Adjusted R-squared:  -2.21e-05 
F-statistic: 0.9996 on 5 and 75 DF,  p-value: 0.424

```

m2
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$shadow.day ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

Residuals:
    Min      1Q  Median      3Q     Max 
-10.717  -4.222  -3.708   2.577  25.081 

Coefficients:
                                                                         Estimate
(Intercept)                                                               40.0981
CHN.CHS.CZonly.EnvData$elev.c                                              4.1559
CHN.CHS.CZonly.EnvData$elev.c.squared                                     -2.6236
CHN.CHS.CZonly.EnvData$TransectCHS                                        -5.1704
CHN.CHS.CZonly.EnvData$TransectCZ                                         -2.9693
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -3.7528
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -5.1039
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   2.8085
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   -0.8664
                                                                         Std. Error
(Intercept)                                                                  2.8851
CHN.CHS.CZonly.EnvData$elev.c                                                3.8447
CHN.CHS.CZonly.EnvData$elev.c.squared                                        6.9138
CHN.CHS.CZonly.EnvData$TransectCHS                                           4.2929
CHN.CHS.CZonly.EnvData$TransectCZ                                            3.9344
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS             4.5749
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ              4.2624
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS     8.3636
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ      7.9583
                                                                         t value
(Intercept)                                                               13.898
CHN.CHS.CZonly.EnvData$elev.c                                              1.081
CHN.CHS.CZonly.EnvData$elev.c.squared                                     -0.379
CHN.CHS.CZonly.EnvData$TransectCHS                                        -1.204
CHN.CHS.CZonly.EnvData$TransectCZ                                         -0.755
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          -0.820
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -1.197
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   0.336
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ   -0.109
                                                                         Pr(>|t|)
(Intercept)                                                                <2e-16
CHN.CHS.CZonly.EnvData$elev.c                                               0.283
CHN.CHS.CZonly.EnvData$elev.c.squared                                       0.705
CHN.CHS.CZonly.EnvData$TransectCHS                                          0.232
CHN.CHS.CZonly.EnvData$TransectCZ                                           0.453
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS            0.415
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ             0.235
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS    0.738
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ     0.914
                                                                            
(Intercept)                                                              ***
CHN.CHS.CZonly.EnvData$elev.c                                               
CHN.CHS.CZonly.EnvData$elev.c.squared                                       
CHN.CHS.CZonly.EnvData$TransectCHS                                          
CHN.CHS.CZonly.EnvData$TransectCZ                                           
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS            
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ             
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS    
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 8.61 on 72 degrees of freedom
Multiple R-squared:  0.07443,	Adjusted R-squared:  -0.02841 
F-statistic: 0.7237 on 8 and 72 DF,  p-value: 0.67
```



##### day10cm

m1
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$day10cm ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect + CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c)

Residuals:
    Min      1Q  Median      3Q     Max 
-63.108 -17.022   0.388  10.250 102.174 

Coefficients:
                                                                 Estimate
(Intercept)                                                       127.300
CHN.CHS.CZonly.EnvData$elev.c                                      78.528
CHN.CHS.CZonly.EnvData$TransectCHS                                -16.586
CHN.CHS.CZonly.EnvData$TransectCZ                                  -4.294
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS    7.159
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -14.876
                                                                 Std. Error
(Intercept)                                                           6.641
CHN.CHS.CZonly.EnvData$elev.c                                        11.717
CHN.CHS.CZonly.EnvData$TransectCHS                                    8.833
CHN.CHS.CZonly.EnvData$TransectCZ                                     8.096
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS     14.069
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ      13.206
                                                                 t value
(Intercept)                                                       19.170
CHN.CHS.CZonly.EnvData$elev.c                                      6.702
CHN.CHS.CZonly.EnvData$TransectCHS                                -1.878
CHN.CHS.CZonly.EnvData$TransectCZ                                 -0.530
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS   0.509
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ   -1.126
                                                                 Pr(>|t|)    
(Intercept)                                                       < 2e-16 ***
CHN.CHS.CZonly.EnvData$elev.c                                    3.34e-09 ***
CHN.CHS.CZonly.EnvData$TransectCHS                                 0.0643 .  
CHN.CHS.CZonly.EnvData$TransectCZ                                  0.5974    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS   0.6123    
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ    0.2636    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 28.52 on 75 degrees of freedom
Multiple R-squared:  0.7876,	Adjusted R-squared:  0.7734 
F-statistic: 55.61 on 5 and 75 DF,  p-value: < 2.2e-16
```


m2
```
Call:
lm(formula = CHN.CHS.CZonly.EnvData$day10cm ~ CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)

Residuals:
    Min      1Q  Median      3Q     Max 
-55.530  -9.073  -1.184   2.247 109.012 

Coefficients:
                                                                         Estimate
(Intercept)                                                              138.8801
CHN.CHS.CZonly.EnvData$elev.c                                             70.1221
CHN.CHS.CZonly.EnvData$elev.c.squared                                    -38.5906
CHN.CHS.CZonly.EnvData$TransectCHS                                       -11.0760
CHN.CHS.CZonly.EnvData$TransectCZ                                          0.8811
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS          10.4862
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -6.9653
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   8.2206
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ    9.6332
                                                                         Std. Error
(Intercept)                                                                  8.9705
CHN.CHS.CZonly.EnvData$elev.c                                               11.9540
CHN.CHS.CZonly.EnvData$elev.c.squared                                       21.4966
CHN.CHS.CZonly.EnvData$TransectCHS                                          13.3476
CHN.CHS.CZonly.EnvData$TransectCZ                                           12.2331
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS            14.2245
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ             13.2528
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS    26.0043
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ     24.7442
                                                                         t value
(Intercept)                                                               15.482
CHN.CHS.CZonly.EnvData$elev.c                                              5.866
CHN.CHS.CZonly.EnvData$elev.c.squared                                     -1.795
CHN.CHS.CZonly.EnvData$TransectCHS                                        -0.830
CHN.CHS.CZonly.EnvData$TransectCZ                                          0.072
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS           0.737
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ           -0.526
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   0.316
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ    0.389
                                                                         Pr(>|t|)
(Intercept)                                                               < 2e-16***
CHN.CHS.CZonly.EnvData$elev.c                                            1.25e-07***
CHN.CHS.CZonly.EnvData$elev.c.squared                                      0.0768.
CHN.CHS.CZonly.EnvData$TransectCHS                                         0.4094
CHN.CHS.CZonly.EnvData$TransectCZ                                          0.9428
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCHS           0.4634
CHN.CHS.CZonly.EnvData$elev.c:CHN.CHS.CZonly.EnvData$TransectCZ            0.6008
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCHS   0.7528
CHN.CHS.CZonly.EnvData$elev.c.squared:CHN.CHS.CZonly.EnvData$TransectCZ    0.6982
                                                                              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 26.77 on 72 degrees of freedom
Multiple R-squared:  0.8203,	Adjusted R-squared:  0.8003 
F-statistic: 41.08 on 8 and 72 DF,  p-value: < 2.2e-16
```




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

par(mfrow=c(3,2))
hist(m2.sol.rad.60d$residuals)
hist(m2.pcpt.60d$residuals)
hist(m2.shadow.days$residuals)
hist(m2.temp.laying.date$residuals)
hist(m2.day10cm$residuals)
```

![alt_txt][CHS.VS.TI.Fig1]

[CHS.VS.TI.Fig1]:https://user-images.githubusercontent.com/12142475/30802217-aa282464-a1dd-11e7-9752-71c04b31a456.png


![alt_txt][CHS.VS.TI.m2.residuals]

[CHS.VS.TI.m2.residuals]:https://user-images.githubusercontent.com/12142475/30808648-e29be808-a1f6-11e7-8fa5-2822360d0a44.png


Then we can interpret the models further like this: 

https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R

OR 

https://rpubs.com/FelipeRego/MultipleLinearRegressionInRFirstSteps


The coefficient - estimate of the intercept is where the predictor variable is 0. I.e. for elevation this would be sol radiation at sea level. But since we've centered and scaled our data, this is actually the sol rad at the median elevation for this dataset. 

```
median(CHS.VS.TI.EnvData$elev)

#1618.5m
```

#### sol.rad.60d

m1
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

m2
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$sol.rad.60d ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

Residuals:
   Min     1Q Median     3Q    Max 
-36268  -9377  -3370   8692  39411 

Coefficients:
                                                                          Estimate
(Intercept)                                                                 102812
CHS.VS.TIonly.EnvData$elev.c                                                 12252
CHS.VS.TIonly.EnvData$elev.c.squared                                        -21752
CHS.VS.TIonly.EnvData$TransectCHS.VS                                         33977
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS            15576
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   -32000
                                                                          Std. Error
(Intercept)                                                                    11205
CHS.VS.TIonly.EnvData$elev.c                                                   11401
CHS.VS.TIonly.EnvData$elev.c.squared                                           19696
CHS.VS.TIonly.EnvData$TransectCHS.VS                                           16016
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS              15026
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS      26512
                                                                          t value
(Intercept)                                                                 9.176
CHS.VS.TIonly.EnvData$elev.c                                                1.075
CHS.VS.TIonly.EnvData$elev.c.squared                                       -1.104
CHS.VS.TIonly.EnvData$TransectCHS.VS                                        2.121
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           1.037
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS  -1.207
                                                                          Pr(>|t|)
(Intercept)                                                                3.3e-08***
CHS.VS.TIonly.EnvData$elev.c                                                 0.297
CHS.VS.TIonly.EnvData$elev.c.squared                                         0.284
CHS.VS.TIonly.EnvData$TransectCHS.VS                                         0.048*
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS            0.314
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS    0.243
                                                                              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 20690 on 18 degrees of freedom
Multiple R-squared:  0.6535,	Adjusted R-squared:  0.5573 
F-statistic: 6.791 on 5 and 18 DF,  p-value: 0.001014

```

And plot both models
```
par(mfrow=c(2,4))
plot(m1.sol.rad)
plot(m2.sol.rad.60d)
```


![alt_txt][CHSall.sol.rad.m1m2]

[CHSall.sol.rad.m1m2]:https://user-images.githubusercontent.com/12142475/30808778-4b5a76b6-a1f7-11e7-9edd-7bd2386d4b2c.png




#### temp.laying.date

m1
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

m2
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$temp.laying.date ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.7693 -0.7813 -0.2319  0.7575  2.2666 

Coefficients:
                                                                          Estimate
(Intercept)                                                                5.28941
CHS.VS.TIonly.EnvData$elev.c                                               0.91653
CHS.VS.TIonly.EnvData$elev.c.squared                                       2.31167
CHS.VS.TIonly.EnvData$TransectCHS.VS                                      -0.46974
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS          0.04905
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS -1.30249
                                                                          Std. Error
(Intercept)                                                                  0.66712
CHS.VS.TIonly.EnvData$elev.c                                                 0.67877
CHS.VS.TIonly.EnvData$elev.c.squared                                         1.17266
CHS.VS.TIonly.EnvData$TransectCHS.VS                                         0.95359
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS            0.89464
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS    1.57847
                                                                          t value
(Intercept)                                                                 7.929
CHS.VS.TIonly.EnvData$elev.c                                                1.350
CHS.VS.TIonly.EnvData$elev.c.squared                                        1.971
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -0.493
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           0.055
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS  -0.825
                                                                          Pr(>|t|)
(Intercept)                                                               2.78e-07***
CHS.VS.TIonly.EnvData$elev.c                                                0.1937
CHS.VS.TIonly.EnvData$elev.c.squared                                        0.0643.
CHS.VS.TIonly.EnvData$TransectCHS.VS                                        0.6282
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           0.9569
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.4201
                                                                                
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.232 on 18 degrees of freedom
Multiple R-squared:  0.3783,	Adjusted R-squared:  0.2056 
F-statistic: 2.191 on 5 and 18 DF,  p-value: 0.1006
```


Model comparison
```
par(mfrow=c(2,4))
plot(m1.temp.laying.date)
plot(m2.temp.laying.date)
```

![alt_txt][CHSall.temp.layind.date.m1m2]

[CHSall.temp.layind.date.m1m2]:https://user-images.githubusercontent.com/12142475/30808949-022035ca-a1f8-11e7-9c87-a3fd26866e71.png


#### pcpt.60d

m1
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

m2
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$pcpt.60d ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

Residuals:
    Min      1Q  Median      3Q     Max 
-6.9030 -2.9233 -0.5687  3.1934  8.2335 

Coefficients:
                                                                          Estimate
(Intercept)                                                                 41.257
CHS.VS.TIonly.EnvData$elev.c                                                -6.472
CHS.VS.TIonly.EnvData$elev.c.squared                                       -11.300
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -20.986
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           12.707
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   10.052
                                                                          Std. Error
(Intercept)                                                                    2.299
CHS.VS.TIonly.EnvData$elev.c                                                   2.339
CHS.VS.TIonly.EnvData$elev.c.squared                                           4.042
CHS.VS.TIonly.EnvData$TransectCHS.VS                                           3.287
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS              3.083
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS      5.440
                                                                          t value
(Intercept)                                                                17.943
CHS.VS.TIonly.EnvData$elev.c                                               -2.766
CHS.VS.TIonly.EnvData$elev.c.squared                                       -2.796
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -6.385
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           4.121
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   1.848
                                                                          Pr(>|t|)
(Intercept)                                                               6.22e-13***
CHS.VS.TIonly.EnvData$elev.c                                              0.012719*
CHS.VS.TIonly.EnvData$elev.c.squared                                      0.011945*
CHS.VS.TIonly.EnvData$TransectCHS.VS                                      5.17e-06***
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS         0.000642***
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS 0.081148.
                                                                             
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4.246 on 18 degrees of freedom
Multiple R-squared:  0.8318,	Adjusted R-squared:  0.785 
F-statistic:  17.8 on 5 and 18 DF,  p-value: 2.083e-06
```


Model comparison
```
par(mfrow=c(2,4))
plot(m1.pcpt.60d)
plot(m2.pcpt.60d)

```

![alt_txt][CHSall.pcpt.60d.m1m2]

[CHSall.pcpt.60d.m1m2]:https://user-images.githubusercontent.com/12142475/30809015-34128a7e-a1f8-11e7-9d43-f19eebfd3ff2.png

#### shadow.days

m1
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

m2
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$shadow.day ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

Residuals:
   Min     1Q Median     3Q    Max 
-7.371 -4.282 -1.119  2.680 14.539 

Coefficients:
                                                                          Estimate
(Intercept)                                                                39.1274
CHS.VS.TIonly.EnvData$elev.c                                                0.1949
CHS.VS.TIonly.EnvData$elev.c.squared                                       -3.1443
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -9.9958
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           0.9813
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   7.6735
                                                                          Std. Error
(Intercept)                                                                   3.4712
CHS.VS.TIonly.EnvData$elev.c                                                  3.5318
CHS.VS.TIonly.EnvData$elev.c.squared                                          6.1016
CHS.VS.TIonly.EnvData$TransectCHS.VS                                          4.9617
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS             4.6550
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS     8.2131
                                                                          t value
(Intercept)                                                                11.272
CHS.VS.TIonly.EnvData$elev.c                                                0.055
CHS.VS.TIonly.EnvData$elev.c.squared                                       -0.515
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -2.015
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           0.211
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.934
                                                                          Pr(>|t|)
(Intercept)                                                               1.37e-09***
CHS.VS.TIonly.EnvData$elev.c                                                0.9566
CHS.VS.TIonly.EnvData$elev.c.squared                                        0.6126
CHS.VS.TIonly.EnvData$TransectCHS.VS                                        0.0591.
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           0.8354
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   0.3625
  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 6.41 on 18 degrees of freedom
Multiple R-squared:  0.2452,	Adjusted R-squared:  0.03555 
F-statistic:  1.17 on 5 and 18 DF,  p-value: 0.362
```


Model comparison
```
par(mfrow=c(2,4))
plot(m1.shadow.days)
plot(m2.shadow.days)
```

![alt_txt][CHSall.shadow.days.m1m2]

[CHSall.shadow.days.m1m2]:https://user-images.githubusercontent.com/12142475/30809073-661713f0-a1f8-11e7-9783-5804fc277811.png




#### day10cm

m1
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

m2
```
Call:
lm(formula = CHS.VS.TIonly.EnvData$day10cm ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)

Residuals:
     Min       1Q   Median       3Q      Max 
-22.4727  -2.7720  -0.8752   2.7838  19.7137 

Coefficients:
                                                                          Estimate
(Intercept)                                                                138.481
CHS.VS.TIonly.EnvData$elev.c                                                74.161
CHS.VS.TIonly.EnvData$elev.c.squared                                       -34.640
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -17.067
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           22.375
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS   -8.855
                                                                          Std. Error
(Intercept)                                                                    6.038
CHS.VS.TIonly.EnvData$elev.c                                                   6.143
CHS.VS.TIonly.EnvData$elev.c.squared                                          10.614
CHS.VS.TIonly.EnvData$TransectCHS.VS                                           8.631
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS              8.097
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS     14.286
                                                                          t value
(Intercept)                                                                22.935
CHS.VS.TIonly.EnvData$elev.c                                               12.072
CHS.VS.TIonly.EnvData$elev.c.squared                                       -3.264
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       -1.977
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS           2.763
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS  -0.620
                                                                          Pr(>|t|)
(Intercept)                                                               8.95e-15***
CHS.VS.TIonly.EnvData$elev.c                                              4.58e-10***
CHS.VS.TIonly.EnvData$elev.c.squared                                       0.00431**
CHS.VS.TIonly.EnvData$TransectCHS.VS                                       0.06351.
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$TransectCHS.VS          0.01280*
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$TransectCHS.VS  0.54316

---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 11.15 on 18 degrees of freedom
Multiple R-squared:  0.9791,	Adjusted R-squared:  0.9733 
F-statistic: 168.8 on 5 and 18 DF,  p-value: 1.81e-14
```


Model comparison
```
par(mfrow=c(2,4))
plot(m1.day10cm)
plot(m2.day10cm)

```

![alt_txt][CHSall.day10cm.m1m2]

[CHSall.day10cm.m1m2]:https://user-images.githubusercontent.com/12142475/30809125-99d4c87c-a1f8-11e7-925b-2d329a00a71f.png






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

## Model comparison

To see which model is the best: 

```
## CHS.VS and TI
AIC(m1.sol.rad)
[1] 559.1436
AIC(m2.sol.rad.60d)
[1] 552.2046
 
AIC(m1.temp.laying.date)
[1] 86.88899
AIC(m2.temp.laying.date)
[1] 85.21746

AIC(m1.pcpt.60d)
[1] 149.3764
AIC(m2.pcpt.60d)
[1] 144.6123
 
AIC(m1.shadow.days)
[1] 161.6098
AIC(m2.shadow.days)
[1] 164.3827
 
AIC(m1.day10cm)
[1] 211.1539
AIC(m2.day10cm)
[1] 190.9545

##CHS, CZ, and CHN

AIC(m1.CHall.sol.rad)
[1] 1893.023
AIC(m2.CHall.sol.rad.60d)
[1] 1890.752

AIC(m1.CHall.temp.laying.date)
[1] 280.659
AIC(m2.CHall.temp.laying.date)
[1] 279.2502
 
AIC(m1.CHall.pcpt.60d)
[1] 553.919
AIC(m2.CHall.pcpt.60d)
[1] 551.0943
 
AIC(m1.CHall.shadow.days)
[1] 584.143
AIC(m2.CHall.shadow.days)
[1] 589.1039
 
AIC(m1.CHall.day10cm)
[1] 780.4261
AIC(m2.CHall.day10cm)
[1] 772.8727
```


### Assess the best models

##### CHS.TI and CHS.VS
```
library(car)
Anova(m2.sol.rad.60d)

Anova Table (Type II tests)

Response: CHS.VS.TIonly.EnvData$sol.rad.60d
                                                                        Sum Sq
CHS.VS.TIonly.EnvData$elev.c                                        3494817539
CHS.VS.TIonly.EnvData$elev.c.squared                                3826335562
CHS.VS.TIonly.EnvData$Transect                                      1960869323
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$Transect          460077479
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$Transect  623763911
Residuals                                                           7706781316
                                                                    Df F value
CHS.VS.TIonly.EnvData$elev.c                                         1  8.1625
CHS.VS.TIonly.EnvData$elev.c.squared                                 1  8.9368
CHS.VS.TIonly.EnvData$Transect                                       1  4.5798
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$Transect          1  1.0746
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$Transect  1  1.4569
Residuals                                                           18        
                                                                      Pr(>F)   
CHS.VS.TIonly.EnvData$elev.c                                        0.010473 * 
CHS.VS.TIonly.EnvData$elev.c.squared                                0.007864 **
CHS.VS.TIonly.EnvData$Transect                                      0.046298 * 
CHS.VS.TIonly.EnvData$elev.c:CHS.VS.TIonly.EnvData$Transect         0.313644   
CHS.VS.TIonly.EnvData$elev.c.squared:CHS.VS.TIonly.EnvData$Transect 0.243053   
Residuals                                                                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


```


## Figures

I want to plot elevation x env variable for 1) global and 2) local comparisons.

I then need to add the linear model (lm) i.s.o spline to draw the figure.


CHS.TI.VS
```
##sol.rad.60d
sol.rad.60d.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$sol.rad.60d ~ CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)


fit.sol.rad.60d.CHS.VS.TI <- augment(sol.rad.60d.fit.CHS.VS.TI)
fit.sol.rad.60d.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.sol.rad.60d.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.sol.rad.60d.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.sol.rad.60d.CHS.VS.TI$colour)

p1.sol.rad.60d.CHS.TI.VS <- ggplot(fit.sol.rad.60d.CHS.VS.TI, aes(x=CHS.VS.TIonly.EnvData.elev, y=CHS.VS.TIonly.EnvData.sol.rad.60d, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.sol.rad.60d.CHS.VS.TI$colour) + geom_line(aes(x=CHS.VS.TIonly.EnvData.elev, y=.fitted, colour=CHS.VS.TIonly.EnvData.Transect, group=CHS.VS.TIonly.EnvData.Transect, color=fit.sol.rad.60d.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))


##temp.laying.date
temp.laying.date.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$temp.laying.date ~ CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)


fit.temp.laying.date.CHS.VS.TI <- augment(temp.laying.date.fit.CHS.VS.TI)
fit.temp.laying.date.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.temp.laying.date.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.temp.laying.date.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.temp.laying.date.CHS.VS.TI$colour)

p2.temp.laying.date.CHS.TI.VS <- ggplot(fit.temp.laying.date.CHS.VS.TI, aes(x=CHS.VS.TIonly.EnvData.elev, y=CHS.VS.TIonly.EnvData.temp.laying.date, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.temp.laying.date.CHS.VS.TI$colour) + geom_line(aes(x=CHS.VS.TIonly.EnvData.elev, y=.fitted, colour=CHS.VS.TIonly.EnvData.Transect, group=CHS.VS.TIonly.EnvData.Transect, color=fit.temp.laying.date.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))




##pcpt.60d
pcpt.60d.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$pcpt.60d ~ CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)


fit.pcpt.60d.CHS.VS.TI <- augment(pcpt.60d.fit.CHS.VS.TI)
fit.pcpt.60d.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.pcpt.60d.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.pcpt.60d.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.pcpt.60d.CHS.VS.TI$colour)

p3.pcpt.60d.CHS.TI.VS <- ggplot(fit.pcpt.60d.CHS.VS.TI, aes(x=CHS.VS.TIonly.EnvData.elev, y=CHS.VS.TIonly.EnvData.pcpt.60d, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.pcpt.60d.CHS.VS.TI$colour) + geom_line(aes(x=CHS.VS.TIonly.EnvData.elev, y=.fitted, colour=CHS.VS.TIonly.EnvData.Transect, group=CHS.VS.TIonly.EnvData.Transect, color=fit.pcpt.60d.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))






##shadow.days
shadow.days.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$shadow.days ~ CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)


fit.shadow.days.CHS.VS.TI <- augment(shadow.days.fit.CHS.VS.TI)
fit.shadow.days.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.shadow.days.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.shadow.days.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.shadow.days.CHS.VS.TI$colour)

p4.shadow.days.CHS.TI.VS <- ggplot(fit.shadow.days.CHS.VS.TI, aes(x=CHS.VS.TIonly.EnvData.elev, y=CHS.VS.TIonly.EnvData.shadow.days, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.shadow.days.CHS.VS.TI$colour) + geom_line(aes(x=CHS.VS.TIonly.EnvData.elev, y=.fitted, colour=CHS.VS.TIonly.EnvData.Transect, group=CHS.VS.TIonly.EnvData.Transect, color=fit.shadow.days.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))






##day10cm
day10cm.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$day10cm ~ CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)


fit.day10cm.CHS.VS.TI <- augment(day10cm.fit.CHS.VS.TI)
fit.day10cm.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.day10cm.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.day10cm.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.day10cm.CHS.VS.TI$colour)

p5.day10cm.CHS.TI.VS <- ggplot(fit.day10cm.CHS.VS.TI, aes(x=CHS.VS.TIonly.EnvData.elev, y=CHS.VS.TIonly.EnvData.day10cm, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.day10cm.CHS.VS.TI$colour) + geom_line(aes(x=CHS.VS.TIonly.EnvData.elev, y=.fitted, colour=CHS.VS.TIonly.EnvData.Transect, group=CHS.VS.TIonly.EnvData.Transect, color=fit.day10cm.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))

```


multiple plot
```
pdf("CHS.TI.VS.envModel.pdf")
multiplot(p1.sol.rad.60d.CHS.TI.VS, p2.temp.laying.date.CHS.TI.VS, p3.pcpt.60d.CHS.TI.VS, p4.shadow.days.CHS.TI.VS, p5.day10cm.CHS.TI.VS, cols=1)
dev.off()
```


CHN.CHS.CZ
```
CHN.CHS.CZonly.EnvData <- subset(CHN.CHS.CZ.EnvData, CHN.CHS.CZ.EnvData$Transect!="CHall")  
CHN.CHS.CZonly.EnvData$elev.c.squared <- CHN.CHS.CZonly.EnvData$elev.c^2
```

CHN.CHS.CZ
```
##sol.rad.60d
sol.rad.60d.fit.CHN.CHS.CZ <- lm(formula = CHN.CHS.CZonly.EnvData$sol.rad.60d ~ CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


fit.sol.rad.60d.CHN.CHS.CZ <- augment(sol.rad.60d.fit.CHN.CHS.CZ)
fit.sol.rad.60d.CHN.CHS.CZ$colour <- gsub("CHN", "palegreen3", fit.sol.rad.60d.CHN.CHS.CZ$CHN.CHS.CZonly.EnvData.Transect)
fit.sol.rad.60d.CHN.CHS.CZ$colour <- gsub("CHS", "mediumpurple4", fit.sol.rad.60d.CHN.CHS.CZ$colour)
fit.sol.rad.60d.CHN.CHS.CZ$colour <- gsub("CZ", "skyblue3", fit.sol.rad.60d.CHN.CHS.CZ$colour)

p6.sol.rad.60d.CHN.CHS.CZ <- ggplot(fit.sol.rad.60d.CHN.CHS.CZ, aes(x=CHN.CHS.CZonly.EnvData.elev, y=CHN.CHS.CZonly.EnvData.sol.rad.60d, group=CHN.CHS.CZonly.EnvData.Transect)) + geom_point(color=fit.sol.rad.60d.CHN.CHS.CZ$colour) + geom_line(aes(x=CHN.CHS.CZonly.EnvData.elev, y=.fitted, colour=CHN.CHS.CZonly.EnvData.Transect, group=CHN.CHS.CZonly.EnvData.Transect, color=fit.sol.rad.60d.CHN.CHS.CZ$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))


##temp.laying.date
temp.laying.date.fit.CHN.CHS.CZ <- lm(formula = CHN.CHS.CZonly.EnvData$temp.laying.date ~ CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


fit.temp.laying.date.CHN.CHS.CZ <- augment(temp.laying.date.fit.CHN.CHS.CZ)
fit.temp.laying.date.CHN.CHS.CZ$colour <- gsub("CHN", "palegreen3", fit.temp.laying.date.CHN.CHS.CZ$CHN.CHS.CZonly.EnvData.Transect)
fit.temp.laying.date.CHN.CHS.CZ$colour <- gsub("CHS", "mediumpurple4", fit.temp.laying.date.CHN.CHS.CZ$colour)
fit.temp.laying.date.CHN.CHS.CZ$colour <- gsub("CZ", "skyblue3", fit.temp.laying.date.CHN.CHS.CZ$colour)



p7.temp.laying.date.CHN.CHS.CZ <- ggplot(fit.temp.laying.date.CHN.CHS.CZ, aes(x=CHN.CHS.CZonly.EnvData.elev, y=CHN.CHS.CZonly.EnvData.temp.laying.date, group=CHN.CHS.CZonly.EnvData.Transect)) + geom_point(color=fit.temp.laying.date.CHN.CHS.CZ$colour) + geom_line(aes(x=CHN.CHS.CZonly.EnvData.elev, y=.fitted, colour=CHN.CHS.CZonly.EnvData.Transect, group=CHN.CHS.CZonly.EnvData.Transect, color=fit.temp.laying.date.CHN.CHS.CZ$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))




##pcpt.60d
pcpt.60d.fit.CHN.CHS.CZ <- lm(formula = CHN.CHS.CZonly.EnvData$pcpt.60d ~ CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


fit.pcpt.60d.CHN.CHS.CZ <- augment(pcpt.60d.fit.CHN.CHS.CZ)
fit.pcpt.60d.CHN.CHS.CZ$colour <- gsub("CHN", "palegreen3", fit.pcpt.60d.CHN.CHS.CZ$CHN.CHS.CZonly.EnvData.Transect)
fit.pcpt.60d.CHN.CHS.CZ$colour <- gsub("CHS", "mediumpurple4", fit.pcpt.60d.CHN.CHS.CZ$colour)
fit.pcpt.60d.CHN.CHS.CZ$colour <- gsub("CZ", "skyblue3", fit.pcpt.60d.CHN.CHS.CZ$colour)


p8.pcpt.60d.CHN.CHS.CZ <- ggplot(fit.pcpt.60d.CHN.CHS.CZ, aes(x=CHN.CHS.CZonly.EnvData.elev, y=CHN.CHS.CZonly.EnvData.pcpt.60d, group=CHN.CHS.CZonly.EnvData.Transect)) + geom_point(color=fit.pcpt.60d.CHN.CHS.CZ$colour) + geom_line(aes(x=CHN.CHS.CZonly.EnvData.elev, y=.fitted, colour=CHN.CHS.CZonly.EnvData.Transect, group=CHN.CHS.CZonly.EnvData.Transect, color=fit.pcpt.60d.CHN.CHS.CZ$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))






##shadow.days
shadow.days.fit.CHN.CHS.CZ <- lm(formula = CHN.CHS.CZonly.EnvData$shadow.days ~ CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


fit.shadow.days.CHN.CHS.CZ <- augment(shadow.days.fit.CHN.CHS.CZ)
fit.shadow.days.CHN.CHS.CZ$colour <- gsub("CHN", "palegreen3", fit.shadow.days.CHN.CHS.CZ$CHN.CHS.CZonly.EnvData.Transect)
fit.shadow.days.CHN.CHS.CZ$colour <- gsub("CHS", "mediumpurple4", fit.shadow.days.CHN.CHS.CZ$colour)
fit.shadow.days.CHN.CHS.CZ$colour <- gsub("CZ", "skyblue3", fit.shadow.days.CHN.CHS.CZ$colour)


p9.shadow.days.CHN.CHS.CZ <- ggplot(fit.shadow.days.CHN.CHS.CZ, aes(x=CHN.CHS.CZonly.EnvData.elev, y=CHN.CHS.CZonly.EnvData.shadow.days, group=CHN.CHS.CZonly.EnvData.Transect)) + geom_point(color=fit.shadow.days.CHN.CHS.CZ$colour) + geom_line(aes(x=CHN.CHS.CZonly.EnvData.elev, y=.fitted, colour=CHN.CHS.CZonly.EnvData.Transect, group=CHN.CHS.CZonly.EnvData.Transect, color=fit.shadow.days.CHN.CHS.CZ$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))






##day10cm
day10cm.fit.CHN.CHS.CZ <- lm(formula = CHN.CHS.CZonly.EnvData$day10cm ~ CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$elev.c.squared + CHN.CHS.CZonly.EnvData$Transect + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev + 
    CHN.CHS.CZonly.EnvData$Transect:CHN.CHS.CZonly.EnvData$elev.c.squared)


fit.day10cm.CHN.CHS.CZ <- augment(day10cm.fit.CHN.CHS.CZ)
fit.day10cm.CHN.CHS.CZ$colour <- gsub("CHN", "palegreen3", fit.day10cm.CHN.CHS.CZ$CHN.CHS.CZonly.EnvData.Transect)
fit.day10cm.CHN.CHS.CZ$colour <- gsub("CHS", "mediumpurple4", fit.day10cm.CHN.CHS.CZ$colour)
fit.day10cm.CHN.CHS.CZ$colour <- gsub("CZ", "skyblue3", fit.day10cm.CHN.CHS.CZ$colour)


p10.day10cm.CHN.CHS.CZ <- ggplot(fit.day10cm.CHN.CHS.CZ, aes(x=CHN.CHS.CZonly.EnvData.elev, y=CHN.CHS.CZonly.EnvData.day10cm, group=CHN.CHS.CZonly.EnvData.Transect)) + geom_point(color=fit.day10cm.CHN.CHS.CZ$colour) + geom_line(aes(x=CHN.CHS.CZonly.EnvData.elev, y=.fitted, colour=CHN.CHS.CZonly.EnvData.Transect, group=CHN.CHS.CZonly.EnvData.Transect, color=fit.day10cm.CHN.CHS.CZ$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))

```


multiple plot
```
pdf("CHN.CHS.CZ.envModel.pdf")
multiplot(p6.sol.rad.60d.CHN.CHS.CZ, p7.temp.laying.date.CHN.CHS.CZ, p8.pcpt.60d.CHN.CHS.CZ, p9.shadow.days.CHN.CHS.CZ, p10.day10cm.CHN.CHS.CZ, cols=1)
dev.off()
```


define multiplot
```
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

```

