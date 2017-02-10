#Gradient Forest Analysis

Here I aim to: 

1)determine the environmental and spatial determinants of variation in allele frequency across the Swedish latitudinal gradient. 

2)map the current and future patterns of genomic variation in relation to climate

I'm using the package GradientForest and the method described by FitzPatrick & Keller 2013


##Data

I'm using 3 different datasets: 

1. Fst outliers as determined by BayEnv2 XtX, and RDA (RD1 & 2)

2. Loci associated with environment as determined by LFMM & BayEnv2 

3. "Neutral" loci (all - Data1+2)


##Presenting Results

As in the Fitzpatric paper.

1. Table1 = summary of the datasets

2. Fig1 = R2 heatmap of reletive importance of different predictor variables

3. Graph of the average allele turnover for each of the three datasets (on one Figure)

4. Graph of the SNPs for BIO5 for ENV & Fst datasets

5. PCA of env variables with sample sites on (PCA centered but not scaled)

6. Geographic map of PC space: colour according to PC space

7. Map difference between datasets (Neutral vs Fst, Neutral vs Env)

8. Pop-level vulnerability to climate change: Map mean genetic offset between current and future climates. 



###Table 1
```


```




###Fig1: R2 heatmap
```
#R

#Retrieve results from the individual outputs for the three GF models: 

R.sq.alldatasets <- (rowMeans(gf.ENVcandidates.model.SEtemp$imp.rsq, na.rm=T))  ##get the mean across all loci
R.sq.alldatasets <- as.data.frame(R.sq.alldatasets)
R.sq.alldatasets
colnames(R.sq.alldatasets) <- "ENV"
R.sq.Fst <- (rowMeans(gf.Fst.model.SEtemp$imp.rsq, na.rm=T))
R.sq.Fst <- as.data.frame(R.sq.Fst)
R.sq.Fst
R.sq.alldatasets$Fst <- R.sq.Fst$R.sq.Fst

R.sq.Neutral <- (rowMeans(gf.NEUTRAL.model.SEtemp$imp.rsq, na.rm=T))
R.sq.Neutral <- as.data.frame(R.sq.Neutral)
R.sq.alldatasets$Neutral <- R.sq.Neutral$R.sq.Neutral

R.sq.alldatasets <- as.matrix(R.sq.alldatasets)  ##turn into matrix for heatmap
R.sq.all.melt <- melt(R.sq.alldatasets)  ##melt for ggplot heatmap
R.sq.all.melt <- R.sq.all.melt[order(R.sq.all.melt$Var1),]  ##order by predictor variable

library(RColorBrewer)
#hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')  ##change the colour palette to red. Default is blue

ggplot(R.sq.all.melt, aes(x=Var2, y=Var1, fill=value)) + geom_tile() + coord_equal() +   ##specify x and y variable, coord_equal changes it to squares
scale_fill_gradient(name="R2") +   ##title of the legend
theme(axis.title.x=element_blank(), axis.title.y=element_blank())   ##remove names of axes
```

![alt_txt][Fig1]
[Fig1]:https://cloud.githubusercontent.com/assets/12142475/22826395/85140a42-ef92-11e6-8dea-5c885369b22b.png


