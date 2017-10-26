# Env Comparison


##SE
```
###Season

SE.variables <- read.csv("SE.Derived.EnvData_20171009.csv", header=T)
SE.variables$Lat.squared <- (SE.variables$Lat)^2
head(SE.variables)

library(ggplot2)
library(car)
library(broom)

SE.season.linearmodel <- lm(formula = SE.variables$days.6.degrees ~ SE.variables$Lat + SE.variables$Lat.squared)
Anova(SE.season.linearmodel)

fit.SE.season.linearmodel <- augment(SE.season.linearmodel)

p1.SE.season.linearmodel <- ggplot(fit.SE.season.linearmodel, aes(x=SE.variables.Lat, y=SE.variables.days.6.degrees)) + geom_point(color="steelblue3") + geom_line(aes(x=SE.variables.Lat, y=.fitted), size=2)


###Mean Temp

SE.variables <- read.csv("SE.Derived.EnvData_20171009.csv", header=T)
SE.variables$Lat.squared <- (SE.variables$Lat)^2
head(SE.variables)

library(ggplot2)
library(car)
library(broom)

SE.temp.linearmodel <- lm(formula = SE.variables$mean.temp.60.days ~ SE.variables$Lat + SE.variables$Lat.squared)
Anova(SE.temp.linearmodel)

fit.SE.temp.linearmodel <- augment(SE.temp.linearmodel)

p1.SE.temp.linearmodel <- ggplot(fit.SE.temp.linearmodel, aes(x=SE.variables.Lat, y=SE.variables.mean.temp.60.days)) + geom_point(color="steelblue3") + geom_line(aes(x=SE.variables.Lat, y=.fitted), size=2)

```


### CHS.TI & CHS.VS

```
CHS.VS.TIonly.EnvData <-  CHS.TI.VS.variables
CHS.VS.TIonly.EnvData <- CHS.VS.TIonly.EnvData[-15,] ##remove stba
CHS.VS.TIonly.EnvData$elev.c <- CHS.VS.TIonly.EnvData$elev - mean(CHS.VS.TIonly.EnvData$elev)
CHS.VS.TIonly.EnvData$elev.c.squared <- CHS.VS.TIonly.EnvData$elev.c^2

##Temp
Temp.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$mean.temp.60d ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)
Anova(Temp.fit.CHS.VS.TI)

fit.Temp.CHS.VS.TI <- augment(Temp.fit.CHS.VS.TI)
fit.Temp.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.Temp.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.Temp.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.Temp.CHS.VS.TI$colour)
fit.Temp.CHS.VS.TI$elev <- CHS.VS.TIonly.EnvData$elev

p1.mean.temp.60d.CHS.TI.VS <- ggplot(fit.Temp.CHS.VS.TI, aes(x=elev, y=CHS.VS.TIonly.EnvData.mean.temp.60d, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.Temp.CHS.VS.TI$colour) + geom_line(aes(x=elev, y=.fitted, group=CHS.VS.TIonly.EnvData.Transect, colour=CHS.VS.TIonly.EnvData.Transect,
color=fit.Temp.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))


##Season
Season.fit.CHS.VS.TI <- lm(formula = CHS.VS.TIonly.EnvData$days.above.6 ~ CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$elev.c.squared + CHS.VS.TIonly.EnvData$Transect + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c + 
    CHS.VS.TIonly.EnvData$Transect:CHS.VS.TIonly.EnvData$elev.c.squared)
Anova(Season.fit.CHS.VS.TI)

fit.Season.CHS.VS.TI <- augment(Season.fit.CHS.VS.TI)
fit.Season.CHS.VS.TI$colour <- gsub("CHS.TI", "thistle3", fit.Season.CHS.VS.TI$CHS.VS.TIonly.EnvData.Transect)
fit.Season.CHS.VS.TI$colour <- gsub("CHS.VS", "mediumpurple4", fit.Season.CHS.VS.TI$colour)
fit.Season.CHS.VS.TI$elev <- CHS.VS.TIonly.EnvData$elev

p2.Season.CHS.TI.VS <- ggplot(fit.Season.CHS.VS.TI, aes(x=elev, y=CHS.VS.TIonly.EnvData.days.above.6, group=CHS.VS.TIonly.EnvData.Transect)) + geom_point(color=fit.Season.CHS.VS.TI$colour) + geom_line(aes(x=elev, y=.fitted, group=CHS.VS.TIonly.EnvData.Transect, colour=CHS.VS.TIonly.EnvData.Transect,
color=fit.Season.CHS.VS.TI$colour), size=2) + scale_colour_manual(values= c("CHS.TI"="thistle3", "CHS.VS"="mediumpurple4"))

```






### CHall

```
CH.EnvData <-  CH.variables
CH.EnvData <- CH.EnvData[-34,] ##remove stba
CH.EnvData$elev.c <- CH.EnvData$elev - mean(CH.EnvData$elev)
CH.EnvData$elev.c.squared <- CH.EnvData$elev.c^2

##Temp
Temp.fit.CHall <- lm(formula = CH.EnvData$mean.temp.60d ~ CH.EnvData$elev.c + 
    CH.EnvData$elev.c.squared + CH.EnvData$Transect2 + 
    CH.EnvData$Transect2:CH.EnvData$elev.c + 
    CH.EnvData$Transect2:CH.EnvData$elev.c.squared)
Anova(Temp.fit.CHall)

fit.Temp.CHall <- augment(Temp.fit.CHall)
fit.Temp.CHall$colour <- gsub("CHN", "palegreen3", fit.Temp.CHall$CH.EnvData.Transect2)
fit.Temp.CHall$colour <- gsub("CHS", "mediumpurple4", fit.Temp.CHall$colour)
fit.Temp.CHall$colour <- gsub("CZ", "skyblue3", fit.Temp.CHall$colour)
fit.Temp.CHall$elev <- CH.EnvData$elev

p1.Temp.CHall <- ggplot(fit.Temp.CHall, aes(x=elev, y=CH.EnvData.mean.temp.60d, group=CH.EnvData.Transect2)) + geom_point(color=fit.Temp.CHall$colour) + geom_line(aes(x=elev, y=.fitted, group=CH.EnvData.Transect2, colour=CH.EnvData.Transect2,
color=fit.Temp.CHall$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))


##Season
Season.fit.CHall <- lm(formula = CH.EnvData$days.above.6 ~ CH.EnvData$elev.c + 
    CH.EnvData$elev.c.squared + CH.EnvData$Transect2 + 
    CH.EnvData$Transect2:CH.EnvData$elev.c + 
    CH.EnvData$Transect2:CH.EnvData$elev.c.squared)
Anova(Season.fit.CHall)

fit.Season.CHall <- augment(Season.fit.CHall)
fit.Season.CHall$colour <- gsub("CHN", "palegreen3", fit.Season.CHall$CH.EnvData.Transect2)
fit.Season.CHall$colour <- gsub("CHS", "mediumpurple4", fit.Season.CHall$colour)
fit.Season.CHall$colour <- gsub("CZ", "skyblue3", fit.Season.CHall$colour)
fit.Season.CHall$elev <- CH.EnvData$elev

p2.Season.CHall <- ggplot(fit.Season.CHall, aes(x=elev, y=CH.EnvData.days.above.6, group=CH.EnvData.Transect2)) + geom_point(color=fit.Season.CHall$colour) + geom_line(aes(x=elev, y=.fitted, group=CH.EnvData.Transect2, colour=CH.EnvData.Transect2,
color=fit.Season.CHall$colour), size=2) + scale_colour_manual(values= c("CHN"="palegreen3", "CHS"="mediumpurple4", CZ="skyblue3"))


```
