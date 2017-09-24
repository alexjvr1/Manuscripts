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
