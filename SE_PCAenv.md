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

data <- read.csv("gbifdf_SE.csv", header=T)


##add a new column with regions specified. 
data$region <- 'region'  #create a new column filled with the word region
head(data)
data$region[data$long > 15] <- 'East'  #rename the region column rows according to longitude
data$region[data$long < 15] <- 'West'
```



PCA of subset of variables (BIO2, 5, 6, 12, 18)
```
pca.data <- prcomp(data[,c(4,7,8,14,20)], center=T, scale.=T)

summary(pca.data)

summary(pca.data)
Importance of components:
                          PC1    PC2    PC3     PC4     PC5
Standard deviation     1.5261 1.3165 0.8785 0.34595 0.21543
Proportion of Variance 0.4658 0.3467 0.1543 0.02394 0.00928
Cumulative Proportion  0.4658 0.8124 0.9668 0.99072 1.00000

pca.data$rotation  ##shows the loadings of each variable on every axis.
             PC1         PC2       PC3          PC4         PC5
bio2  -0.3664481 -0.53103209 0.4812926 -0.214987790 -0.55304473
bio5   0.4688835 -0.06575763 0.7820183  0.111568894  0.38964588
bio6   0.4273046  0.55679423 0.1382850  0.003427836 -0.69875353
bio12 -0.4066736  0.54131371 0.2658444 -0.645798654  0.23209252
bio18 -0.5457943  0.33262668 0.2588607  0.724061448 -0.01393525
```



### Plot

```
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(pca.data)
```

first just the plain graph: 

