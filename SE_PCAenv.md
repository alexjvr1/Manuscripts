# SE PCA env

Aim: 

Determine if the ecological niche of R.temp differs between East and West Sweden. 

Also, determine if the samples collected for genomic analyses is similar to SE east. 


Data have been downloaded from gbif. 

Analyses done in R: 

```
#R version 3.5.0 (2018-04-23) -- "Joy in Playing"
getwd()

/Users/alexjvr/2016RADAnalysis/StudentENMproject

data <- read.csv("")


##add a new column with regions specified. 
data$region <- 'region'
head(data)
data$region[data$long > 15] <- 'East'
data$region[data$long < 15] <- 'West'
```

