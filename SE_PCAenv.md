# SE PCA env

Aim: 

1. Determine if the ecological niche of R.temp differs between East and West Sweden. 

2. Determine if the samples collected for genomic analyses is similar to Sweden_East. 


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

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/49028289-f2c3ab00-f199-11e8-8e27-a08e432ed68c.png



Now change some features. See [here](https://stackoverflow.com/questions/40287943/ggbiplot-change-the-group-color-and-marker) for more info on this. 
```
data.region <- data$region  ##object with the region names in
ggbiplot(pca.data, group=data.region) + scale_color_manual(name="Region", values=c("goldenrod1", "cyan3")) + theme(legend.direction ="horizontal", 
      legend.position = "right")

#scale_shape_manual(name="Variety", values=c(17:19)) +
#geom_point(aes(colour=wine.class, shape=wine.class), size = 3)
```

Change the colours as you wish using the [R colours](http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)

![alt_txt][Fig2]

[Fig2]:https://user-images.githubusercontent.com/12142475/49029688-1f2cf680-f19d-11e8-9c79-e1e4f8eb754b.png



### Is the genomic data different? 

So the East vs West niche looks quite different. We'll quantify this in a moment. We probably also want to determine whether the samples collected for genomic analyses are a good representative of the eastern niche:

```
#read in the data
data.genom <- read.table("SEonly.BIOclim", header=T) 
head(data.genom)

## select only the variables of interest and create a df for the pca
data4pca <- rbind((data[,c(4,7,8,14,20)]),(data.genom[,c(4,7,8,14,20)])) 

#make a regions list for the new dataset
data.genom$region <- 'genom'  
data.region.genom <- data.genom$region #write the regions column to a variable 
data.region2 <- c(data.region, data.region.genom)  
summary(data.region2)  ##check that this is the right lenght and contains all the variables

#run PCA and plot as before

pca.data2 <- prcomp(data4pca, center=T, scale.=T)
ggbiplot(pca.data2, group=data.region2) + scale_color_manual(name="Region", values=c("goldenrod1", "black", "cyan3")) + theme(legend.direction ="horizontal", 
      legend.position = "right")
```


![alt_txt][Fig3]

[Fig3]:https://user-images.githubusercontent.com/12142475/49031456-e2173300-f1a1-11e8-9596-21959cf945c8.png


It's a bit unexpected that three of the locations where we collected genomic data from does not fall within the eastern ecoregion. The data look correct (I've checked the order of the data.region and data4pca files, and I've plotted a histogram of all the variables to check that they haven't been reordered or differently normalised: 
```
tail(data4pcr)

par(mfrow=c(3,2))
hist(data4pca$bio2)
hist(data4pca$bio5)
hist(data4pca$bio6)
hist(data4pca$bio12)
hist(data4pca$bio18)
```


