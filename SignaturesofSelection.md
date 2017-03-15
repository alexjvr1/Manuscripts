# Genomics of adaptation. 

## What to expect, and how to search for the signals. 

A Summary of the literature

## What to expect: 

Most of this is from theory (models)

Or tests from long-term selection experiments

Population genomics: looking for selection by looking for patterns of sweeps. The signature is reduced heterozygosity due to fixation of an allele. 

Hard sweep: sweep from a new mutation. (single, large effect). Theory predicts a strong pattern with a wide signature (due to hitchiking)
- Complete hard sweeps leave a signature ~1/10 the selection on the allele in recombination units. 
- Incomplete hard sweep: could be due to insufficient time, selection coefficients changing over time, and epistatic interactions influencing the sweep trajectory

Soft sweep: from standing variation, often multiple alleles. In asexual haplotypes, soft selection acts on multiple allele copies (i.e. in different strains). But these may or may not be Identical by Descent. 
i.e. if beneficial alleles enter the population independently, either through mutation of migration, these are considered multiple-origin soft sweeps.  
A single-origin soft sweep will have a narrower footprint of reduced heterogeneity; since the allele is older, recombination will break down the sweep pattern. 
However, multiple-origin soft sweeps will retain the heterozygosity around the beneficial site, since the genetic background is variable. This is almost indistinguishable from partial or incomplete soft sweeps. 

interesting case of multiple-origin soft-sweep: independent arising of the same mutation (allele) for pesticide resistance (Karasov et al 2010) - drosophila

1. sexual 

Reviewed in Burke (2012). Sexual recombination allows for the fixation of a beneficial mutation without reducing genome-wide variation. Thus good systems in which to detect hard selection (i.e. new mutations becoming fixed in the population). 
Theoretically both hard and soft sweeps can be studied in sexually reproducing systems. But practically, the population sizes and the number of generations for which they can be maintained, means that new mutations are rare. So most of the selection is acting on standing variation (i.e. soft sweeps). 


Something striking when thinking about sexual vs asexual populations, is that it seems extremely important to maintain genomic variation as much as possible. 
Through recombination, mutation, gene-flow and drift. So even if a locus sweeps to fixation, the neighbouring sites will regain heterzygosity through recombination (eventually). And it seems logical that selection on complex traits would rather be "flexible" in producing the same phenotype from different genomic determinents. 
In this way, selection is possible regardless of the standing genetic variation. And maximum genomic variation is maintained despite changing selective pressures. 



2. asexual 

Reviewed in Burke (2012). Asexual systems for experimental evolution provides a system in which to study hard selective sweeps. Theoretically, if there were a beneficial mutation in a haplotype, it would increasein frequency very rapidly (the whole haplotype, since there's no recombination). 
So in a new beneficial mutation would sweep to fixation. In practice, this seems more complicated. Where multiple mutations/haplotypes are maintained, and hard sweeps only lead to fixation in a fraction of cases. 
Barrick & Lenski 2009: Most of the time (~50% of cases) the beneficial allele will eventually be fixed. In other cases, the allele might be outcompeted by a competing haplotype, or kept at transient frequencies in the population (perhaps due to frequency-dependent selection).  
Lang et al. 2012 tested 2 selection advantages (s=0.6% and s=1.5%) with pop sizes Ne=10^5 and 10^6. Hard sweeps were infrequent and only seen in the smaller pops with higher selective advantage. Other cases: outcompeted, i.e. clonal interference, or remain at intermediate frequency. Or rise to some frequency, and then after a few generations move to fixation (i.e. complex dynamics). 

Conclusion: Even under the simplest constructed scenarios of biological selection, selection is complex, with some element of yet unpredictable dynamics. 


Or genome-enabled systems

1. Humans

2. Drosophila

3. Arabidopsis


4. Stickleback


5. Forest trees



What we know from Natural systems




Large effect loci vs quantitative traits

Effect of gene flow

Effect of selection strength

Standing mutations vs new mutations (i.e. mutation rate)

Population size (expansion/bottleneck/drift)

Recombination





##Population genomics/Fst outlier approaches

Population genomics aims to understand how adaptive and non-adaptive processes are shaping genomic variation: 
Review by Vitti et al. 2013



1. By looking at the effect of selection on loci

2. By describing the effect of population structure and demographic history on neutral variation (using ML- and ABC-based approaches - see below). 

Fst outlier approaches are looking for signatures of genomic sweeps. i.e. large effect loci (single or multiple) that are beneficial in one environment, but not another. 
i.e. it looks for regions with high Fst (genetic divergence), with reduced heterzygosity around the sweep region (sweep footprint). 
This footprint will change based on the type of sweep. - if it's a hard sweep, there will be a wide region around the beneficial locus that will show reduced heterozygosity. Soft sweeps (i.e. from standing variation) will show a narrower region of reduced variation (since the mutation is older). 

Soft-sweeps can be single-origin (i.e. the loci are all identical by descent), or multiple origin. For single-origin, there should be the same reduction in heterozygosity around the locus under selection, while for multiple origin, there will still be heterogeneity around the locus, since the genetic background is different. However, this is almost impossible to tell apart from incomplete soft-sweeps. 

However, current estimators vary massively in their results. - Although this could be explained by a variation in what they are actually measuring: 

1. Divergence-based methods (ie. Fst outlier/McDonald-Kreitman (MK) or HKA) count the effective fixation of many weakly selective mutation over longer evolutionary time. 

Drawbacks of these methods are that 1) populations need to be defined a priori, which is not always simple (weak structure, continuous pops, small sample sizes), 2) Fst will work better for hard sweeps, which is more likely to occur in large populations (time for new mutations to arise), 3) There is no ecological hypothesis to support the test (Pavladis et al. 2012 simulated datasets and performed outlier analyses. They showed that false positives showed significant enrichment of biologically interesting gene catagories)

Positive: genotype-first approach. I.e. no bias in deciding what the phenotype under selection might be. 


2. Polymorphism-based methods (e.g. SFS) are most impacted by the recent fixation of strongly beneficial mutations. 




![alt_txt][joost.table1]

[joost.table1]:https://cloud.githubusercontent.com/assets/12142475/15557079/f310ed42-2285-11e6-864e-468724ef4495.png


#### Demographic effects on loci: 

##### Gene-flow/selection balance

When gene flow is high, loci are unlikely to become fixed. But if gene-flow is low, a new selected allele is much more likely to become fixed in a geographically heterogeneous environment. 


####Pop genetics classic tests for selection/non-neutrality 

1. Lewontin-Krakauer (1973): Fst outlier. Tests the observed value of the inter-locus variance of Fst, assuming gene frequency at a locus for each subpopulation is randomly drawn from a given frequency distribution, and that Fst is the same for all loci. 

2. McDonald-Kreitman (1991): the ratio of fixed synonymous vs fixed non-synonymous mutations 



##### Population size changes

Large populations: more likey for new mutations to arise, thus hard sweeps. 





How do I identify sweeps under different demographic scenarios? 




## Landscape genomics/Environmental association analyses

These methods corrolate allele frequencies with environmental variables. They are more likely to identify ecologically relevant loci (Eckert et al 2010). But these methods do not take population structure or demographic history into account, which could lead to false positives. 

![alt_txt][Joost2]

[Joost2]:https://cloud.githubusercontent.com/assets/12142475/15559047/d54f6d84-2293-11e6-8b9d-af6c5d8d598e.png


Ways to correct for population structure: 

1. Use a control dataset to estimate the null distribution and correct the probability (p-value) at which the null distribution is rejected (Hancock et al. 2008). But a control dataset is not often available.

2. Incorporate a covariance matrix to account for spatial autocorrelation of individuals collected at the same place (Poncet et al. 2010). 

3. Incorporate a covariance matrix of Moran's Eigenvector Maps (MEMs) for unaccounted for environmental variables (Manel et al. 2010)

4. Covariance matrix of genetic structure (Hancock et al. 2008, Coop et al. 2010). - A control dataset is needed to estimate relatedness. Often the same dataset is used which leads to circularity and a loss of statistical power. 

5. LFMM (Frichot 2013) simultaneously fits a model of population structure and environmental effects. Population structure is modelled as K latent factors (independent, linear combinations of the genetic data estimated from joint distributions), and environmental covariates are modeled as fixed effects. (see also De Mita et al. 2013). This method seems to have reduced false positives and false negatives, so is better at detecting environmental associations with population structure. 









### Why are we detecting different loci with the different methods?




### Combining methods

We should be combining population genomic and landscape genomic methods when looking for adaptation/selection in natural populations. 

-> Explain why allele frequencies differ between populations + spatial genetic variation in relation to environmental adaptation, while correcting for demographic history and population structure. 

![alt_txt][comb]

[comb]:https://cloud.githubusercontent.com/assets/12142475/15559116/415d5d42-2294-11e6-9375-1a676905be4b.png


## Repeatability of selection


### How frequent is parallel patterns of adaptation? 




Transcriptomics

Gene pathways

Genes

Mutations


### Parallelism in quantitative traits?



## What are the appropriate methods to use and when? 



## Describing neutral demographic history of populations

Based on expected site frequency spectrum (ie. frequencies of the segregating sites), linkage disequilibrium (ie. association between segregating mutations), divergence data (i.e. fixed differences between populations/species). 

SFS - mostly for demographic estimation (Thorton & Andolfatto 2006)

LD - Genomic scans for adaptive sweeps/fixations (Pavlidis et al. 2010)

Recurrent hitchiking estimation (e.g. Jensen et al. 2008). 

### ML-based approaches


###Approximate Bayesian approaches

reviewed in Beaumont 2010
