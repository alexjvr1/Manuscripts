## Venn Diagrams 

Calculating the overlap in loci identified for all the different datasets. 

getwd()
"/Users/alexjvr/2016RADAnalysis/3_CH.landscapeGenomics/subsets/BayENV2_results/CHP2.results"

I've saved the R project from fgcz47 where I was working: 

fgcz47:/srv/kenlab/alexjvr_p1795/CHcomplete/BayENV2/CHP2


I'm showing: 

1. Venn Diagram within analysis for 
  
  1. LFMM
  
  2. bayEnv2

2. overlap between LFMM, Bayenv2, PCadapt, XtX within each region (i.e. 6 plots)

3. Overlap between all candidate loci identified for CHall, CHN, CHS, CZ

4. Overlap between all candidate loci identified for CHS, CHS.TI, CHS.VS


### 1.1 Venn Diagram for loci within LFMM

These have already been drawn. 

See: https://github.com/alexjvr1/Manuscripts/blob/master/Chp2_CHdata_LFMM.md



Venn.CHN.n5.LFMMonly_20170310.pdf

Venn.CHS.n5.LFMMonly_20170310.pdf

Venn.CZ.n5.LFMMonly_20170310.pdf

Venn.CHS.VS.n5.LFMMonly_20170310.pdf

Venn.CHS.TI.n5.LFMMonly_20170310.pdf

### 1.2 Venn Diagram for loci within BayEnv2


##### CHall

```
load("CHP2.BayEnv2.fgcz.RData")

library(VennDiagram)

d1 <- length(CHall.rad.bayenv.candidates.names)
d2 <- length(CHall.shadow.days.bayenv.candidates.names)
d3 <- length(CHall.temp.bayenv.candidates.names)
d4 <- length(CHall.pcpt.bayenv.candidates.names)
d5 <- length(CHall.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CHall.temp.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CHall.temp.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CHall.pcpt.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.shadow.days.bayenv.candidates.names,CHall.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.shadow.days.bayenv.candidates.names,CHall.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.shadow.days.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,CHall.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,CHall.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CHall.temp.bayenv.candidates.names, CHall.pcpt.bayenv.candidates.names,CHall.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,
CHall.pcpt.bayenv.candidates.names, CHall.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,
CHall.day.10cm.bayenv.candidates.names, CHall.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names,
CHall.pcpt.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CHall.shadow.days.bayenv.candidates.names, CHall.rad.bayenv.candidates.names,
CHall.pcpt.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CHall.temp.bayenv.candidates.names, CHall.rad.bayenv.candidates.names,
CHall.pcpt.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CHall.rad.bayenv.candidates.names, CHall.shadow.days.bayenv.candidates.names, CHall.temp.bayenv.candidates.names, 
CHall.pcpt.bayenv.candidates.names, CHall.day.10cm.bayenv.candidates.names)))

pdf("CHall.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CHall.rad", "CHall.shadow.days", "CHall.pcpt", "CHall.temp", "CHall.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()
```


##### CHN

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHN.rad.bayenv.candidates.names
CHN.shadow.days.bayenv.candidates.names
CHN.temp.bayenv.candidates.names
CHN.pcpt.bayenv.candidates.names
CHN.day.10cm.bayenv.candidates.names


CHN.shadow.days.bayenv.candidates.names <- NULL
CHN.temp.bayenv.candidates.names <- NULL
CHN.pcpt.bayenv.candidates.names <- NULL
CHN.day.10cm.bayenv.candidates.names <- NULL


d1 <- length(CHN.rad.bayenv.candidates.names)
d2 <- length(CHN.shadow.days.bayenv.candidates.names)
d3 <- length(CHN.temp.bayenv.candidates.names)
d4 <- length(CHN.pcpt.bayenv.candidates.names)
d5 <- length(CHN.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CHN.temp.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CHN.temp.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CHN.pcpt.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.shadow.days.bayenv.candidates.names,CHN.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.shadow.days.bayenv.candidates.names,CHN.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.shadow.days.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,CHN.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,CHN.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CHN.temp.bayenv.candidates.names, CHN.pcpt.bayenv.candidates.names,CHN.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,
CHN.pcpt.bayenv.candidates.names, CHN.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,
CHN.day.10cm.bayenv.candidates.names, CHN.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names,
CHN.pcpt.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CHN.shadow.days.bayenv.candidates.names, CHN.rad.bayenv.candidates.names,
CHN.pcpt.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CHN.temp.bayenv.candidates.names, CHN.rad.bayenv.candidates.names,
CHN.pcpt.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CHN.rad.bayenv.candidates.names, CHN.shadow.days.bayenv.candidates.names, CHN.temp.bayenv.candidates.names, 
CHN.pcpt.bayenv.candidates.names, CHN.day.10cm.bayenv.candidates.names)))

pdf("CHN.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CHN.rad", "CHN.shadow.days", "CHN.pcpt", "CHN.temp", "CHN.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()


```



##### CHS

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.rad.bayenv.candidates.names
CHS.shadow.days.bayenv.candidates.names
CHS.temp.bayenv.candidates.names
CHS.pcpt.bayenv.candidates.names
CHS.day.10cm.bayenv.candidates.names



d1 <- length(CHS.rad.bayenv.candidates.names)
d2 <- length(CHS.shadow.days.bayenv.candidates.names)
d3 <- length(CHS.temp.bayenv.candidates.names)
d4 <- length(CHS.pcpt.bayenv.candidates.names)
d5 <- length(CHS.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CHS.temp.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CHS.temp.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CHS.pcpt.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.shadow.days.bayenv.candidates.names,CHS.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.shadow.days.bayenv.candidates.names,CHS.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.shadow.days.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,CHS.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,CHS.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CHS.temp.bayenv.candidates.names, CHS.pcpt.bayenv.candidates.names,CHS.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,
CHS.pcpt.bayenv.candidates.names, CHS.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,
CHS.day.10cm.bayenv.candidates.names, CHS.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names,
CHS.pcpt.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CHS.shadow.days.bayenv.candidates.names, CHS.rad.bayenv.candidates.names,
CHS.pcpt.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CHS.temp.bayenv.candidates.names, CHS.rad.bayenv.candidates.names,
CHS.pcpt.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CHS.rad.bayenv.candidates.names, CHS.shadow.days.bayenv.candidates.names, CHS.temp.bayenv.candidates.names, 
CHS.pcpt.bayenv.candidates.names, CHS.day.10cm.bayenv.candidates.names)))

pdf("CHS.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CHS.rad", "CHS.shadow.days", "CHS.pcpt", "CHS.temp", "CHS.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()


```



##### CZ
```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CZ.rad.bayenv.candidates.names
CZ.shadow.days.bayenv.candidates.names
CZ.temp.bayenv.candidates.names
CZ.pcpt.bayenv.candidates.names
CZ.day.10cm.bayenv.candidates.names



d1 <- length(CZ.rad.bayenv.candidates.names)
d2 <- length(CZ.shadow.days.bayenv.candidates.names)
d3 <- length(CZ.temp.bayenv.candidates.names)
d4 <- length(CZ.pcpt.bayenv.candidates.names)
d5 <- length(CZ.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CZ.temp.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CZ.temp.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CZ.pcpt.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.shadow.days.bayenv.candidates.names,CZ.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.shadow.days.bayenv.candidates.names,CZ.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.shadow.days.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,CZ.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,CZ.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CZ.temp.bayenv.candidates.names, CZ.pcpt.bayenv.candidates.names,CZ.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,
CZ.pcpt.bayenv.candidates.names, CZ.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,
CZ.day.10cm.bayenv.candidates.names, CZ.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names,
CZ.pcpt.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CZ.shadow.days.bayenv.candidates.names, CZ.rad.bayenv.candidates.names,
CZ.pcpt.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CZ.temp.bayenv.candidates.names, CZ.rad.bayenv.candidates.names,
CZ.pcpt.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CZ.rad.bayenv.candidates.names, CZ.shadow.days.bayenv.candidates.names, CZ.temp.bayenv.candidates.names, 
CZ.pcpt.bayenv.candidates.names, CZ.day.10cm.bayenv.candidates.names)))

pdf("CZ.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CZ.rad", "CZ.shadow.days", "CZ.pcpt", "CZ.temp", "CZ.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()

```



##### CHS.VS
```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.VS.rad.bayenv.candidates.names
CHS.VS.shadow.days.bayenv.candidates.names
CHS.VS.temp.bayenv.candidates.names
CHS.VS.pcpt.bayenv.candidates.names
CHS.VS.day.10cm.bayenv.candidates.names



d1 <- length(CHS.VS.rad.bayenv.candidates.names)
d2 <- length(CHS.VS.shadow.days.bayenv.candidates.names)
d3 <- length(CHS.VS.temp.bayenv.candidates.names)
d4 <- length(CHS.VS.pcpt.bayenv.candidates.names)
d5 <- length(CHS.VS.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CHS.VS.temp.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CHS.VS.temp.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.shadow.days.bayenv.candidates.names,CHS.VS.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.shadow.days.bayenv.candidates.names,CHS.VS.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.shadow.days.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,CHS.VS.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,CHS.VS.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CHS.VS.temp.bayenv.candidates.names, CHS.VS.pcpt.bayenv.candidates.names,CHS.VS.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,
CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,
CHS.VS.day.10cm.bayenv.candidates.names, CHS.VS.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names,
CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.rad.bayenv.candidates.names,
CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CHS.VS.temp.bayenv.candidates.names, CHS.VS.rad.bayenv.candidates.names,
CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CHS.VS.rad.bayenv.candidates.names, CHS.VS.shadow.days.bayenv.candidates.names, CHS.VS.temp.bayenv.candidates.names, 
CHS.VS.pcpt.bayenv.candidates.names, CHS.VS.day.10cm.bayenv.candidates.names)))

pdf("CHS.VS.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CHS.VS.rad", "CHS.VS.shadow.days", "CHS.VS.pcpt", "CHS.VS.temp", "CHS.VS.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()

```


##### CHS.TI

```
library(VennDiagram)


#Check that all the lists have names in them. Where no candidates were found, there’ll only be an “X.”
#These can be made to <- NULL so that they aren’t reflected in the Venn

CHS.TI.rad.bayenv.candidates.names
CHS.TI.shadow.days.bayenv.candidates.names
CHS.TI.temp.bayenv.candidates.names
CHS.TI.pcpt.bayenv.candidates.names
CHS.TI.day.10cm.bayenv.candidates.names



d1 <- length(CHS.TI.rad.bayenv.candidates.names)
d2 <- length(CHS.TI.shadow.days.bayenv.candidates.names)
d3 <- length(CHS.TI.temp.bayenv.candidates.names)
d4 <- length(CHS.TI.pcpt.bayenv.candidates.names)
d5 <- length(CHS.TI.day.10cm.bayenv.candidates.names)


d12 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.shadow.days.bayenv.candidates.names)))
d13 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names)))
d14 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names)))
d15 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d23 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names)))
d24 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names)))
d25 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d34 <- length(Reduce(intersect, list(CHS.TI.temp.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names)))
d35 <- length(Reduce(intersect, list(CHS.TI.temp.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d45 <- length(Reduce(intersect, list(CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))

d123 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.shadow.days.bayenv.candidates.names,CHS.TI.temp.bayenv.candidates.names)))
d124 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.shadow.days.bayenv.candidates.names,CHS.TI.pcpt.bayenv.candidates.names)))
d125 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.shadow.days.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))
d234 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,CHS.TI.pcpt.bayenv.candidates.names)))
d134 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,CHS.TI.pcpt.bayenv.candidates.names)))
d135 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))
d145 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))
d235 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))
d245 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))
d345 <- length(Reduce(intersect, list(CHS.TI.temp.bayenv.candidates.names, CHS.TI.pcpt.bayenv.candidates.names,CHS.TI.day.10cm.bayenv.candidates.names)))


d1234 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,
CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.rad.bayenv.candidates.names)))
d1235 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,
CHS.TI.day.10cm.bayenv.candidates.names, CHS.TI.rad.bayenv.candidates.names)))

d2345 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names,
CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d1245 <- length(Reduce(intersect, list(CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.rad.bayenv.candidates.names,
CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d1345 <- length(Reduce(intersect, list(CHS.TI.temp.bayenv.candidates.names, CHS.TI.rad.bayenv.candidates.names,
CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))
d12345 <- length(Reduce(intersect, list(CHS.TI.rad.bayenv.candidates.names, CHS.TI.shadow.days.bayenv.candidates.names, CHS.TI.temp.bayenv.candidates.names, 
CHS.TI.pcpt.bayenv.candidates.names, CHS.TI.day.10cm.bayenv.candidates.names)))

pdf("CHS.TI.bayenv2.VENN.pdf")
draw.quintuple.venn(area1=d1, area2=d2, area3=d3, area4=d4, area5=d5,
n12=d12, n13=d13, n14=d14, n15=d15, n23=d23, n24=d24, n25=d25, n34=d34, n35=d35, n45=d45,
n123=d123, n124=d124, n125=d125, n134=d134, n135=d135, n145=d145, n234=d234, n235=d235, n245=d245, n345=d345,
n1234=d1234, n1235=d1235, n1245=d1245, n1345=d1345, n2345=d2345, n12345=d12345, 
category=c("CHS.TI.rad", "CHS.TI.shadow.days", "CHS.TI.pcpt", "CHS.TI.temp", "CHS.TI.day.10cm"),
lty="blank", 
fill=c("yellow", "orange", "skyblue1", "skyblue3", "blue")
)
dev.off()


```


### 2 overlap between LFMM, Bayenv2, PCadapt, XtX within each region


##### CHall

##### CHN

##### CHS

##### CZ

##### CHS.VS

##### CHS.TI


### 3. Overlap between all candidate loci identified for CHall, CHN, CHS, CZ


### 4. Overlap between all candidate loci identified for CHS, CHS.TI, CHS.VS
