#Review: NGS in Ecological studies


##Goal
This is a review stemming from the BAG workshop held in March 2015. 



 (Davey et al. 2010)

Probable B@G manuscript structure

Introduction
Conservation biology, best practices for genomic studies

Start from difficulties in using NGS data by conservation geneticists. This is well stated in the TREE paper. Emphasise that nowadays NGS can be cheaper than SSR and produce massive information (need refs, but Martin is very good at this, perhaps in his paper?). Our goal is to bridge the gap stated in the TREE paper and give an overview of the strategies available to address conservation issues.

1)	Project design
First step is to produce a molecular dataset, and this needs to be proportional to the questions one aims to address. If overall pop stats are required, sparse genome sampling in the order of few thousand of markers can be sufficient and GBS (+RAD) (as well as SNP chips for model systems) are nowadays possibly the most approachable methods (refer to upcoming paper from Martin showing shortcomings of SSR to infer pop stats). On the other hand, in case adaptive loci constitute the aim of the study, dense sampling is necessary in order to capture signal of differentiation islands whose size can be limited in wild pop because of short LD. For this, maximising number of markers in RAD or WGS would be suitable approaches. Pool-Seq could be mentioned and Martin’s paper cited.

However, it needs to be emphasised, that the use of adaptive genes in conservation is still highly controversial (TREE paper and Refs therein).

a.	What is the question?
Nr markers vs individuals – references for Sanger sequencing: Felsenstein 2006; Carling & Brumfield 2007. 
See mol ecol post on transcriptomics: samples vs sequences (http://www.molecularecologist.com/2015/05/next-generation-sequencing-more-replicates-or-more-sequence/?utm_source=feedburner&utm_medium=email&utm_campaign=Feed%3A+molecularecologist%2FfewY+%28The+Molecular+Ecologist%29) 

b.	Sampling strategy: indivs vs pops
c.	Genomic resources available?
d.	Density of markers are needed? 
Perhaps as a final comment: experimental design is key. Strategies to identify adaptive loci suffer from high false positive rate, and pops should be sampled in order to overcome this. Tricky thing is always to distinguish between local adaptation, drift, and adaptation to a more general factor (e.g. altitude). Possibly experimental design should include multiple pops, but perhaps this is not doable for endangered systems.

2)	Marker generation: 
a.	Wet lab approaches  
Screen genome-wide genetic diversity to infer the endangered populations or species. Probable approaches which can be used (showing advantage and disadvantages)
i.	RADseq
ii.	Whole Genome seq
iii.	RNAseq
iv.	Exomecapture
v.	Poolseq

b.	Post-sequencing processing & SNP calling
Beside the sequencing technology, producing a dataset includes SNP calling, and nowadays different strategies are available. Refer to Stacks originally developed for mapping pops, but potentially suitable for wild pops. Compare to novel approaches based on population sapling (vs. individual as Stacks), indel problem, and overall SNP calling accuracy (Jonathan’s paper). For WGS, refer to developments of the ANGSD package that requires minimal coverage per individuals as based on likelihood of SNPs inferred from the SFS. As alternative Pool-Seq is a good option, but with limitations (Schlotterer’s Nat Gen Rev).

i.	Pipelines 
ii.	Filtering
iii.	SNP calling
1.	Freebayes
2.	SAMtools
3.	GATK
Using probabilistic models instead of SNP or genotype calling: 
	ngsTools Fumagalli et al. 2014

Summaries of pipelines:
Singhal 2012 http://arxiv.org/pdf/1211.1737.pdf

iv.	Filter SNPs
Highlight the importance of filtering, and bring examples of the thorough procedure proposed by Jon. IDEA: if we involve him, I suspect he could be interested in putting his filtering strategy in a box, and this would be highly valuable for a wide audience. I don’t think he could publish that as stand-alone, and he could appreciate this opportunity! Besides, he could provide inputs on other parts as well. I’m in touch with him, I can ask him.

v.	Reference Mol Ecol 2015 (errors from RADseq..)

3)	Biological questions
Identification of the genetic diversity within an among populations to infer the survival and extinction risk

a.	Summary stats
Once the dataset is available, downstream analyses are tightly linked to the question. Pop stats on RAD directly output by Stacks, else more manual pipelines starting from a VCF file (e.g. Hierfstat, vcftools, etc.). Popoolation2 for pooled data. On Fst, highlight Sam’s notion to NOT average Fst within pop as this produces biased values. Would be possible to have a box on this showing an example. I’m sure this would be very valuable as it is virtually unknown to people!
i.	Within & between pop genetic diversity
ii.	Comparisons with other pops.. 
iii.	Considerations (e.g. a priori populations? Population summary stats?)

b.	Demographic inference
Making demographic inference (involve Daniel Wegmann) to infer whether the population/species went through a strong bottleneck or hybridisation event with alien or close relative species. Helps to infer whether the species is endangered because human activity or was forever at low density. Probable connectivity among isolated populations could also be inferred

Demography also relevant to infer pop history, as explained by Martin above. Perhaps emphasise that can be used to infer selection too. Daniel will write an encyclopaedia on this…



c.	Signals of selection
Identify the potential of local adapted population that should be treated separately. Identify the adaptive potential and the degree of local adaption with genome scans.

Identify the important environmental factors and genes involved in local adaptation using environmental association studies.

Once we have the outliers, what do we do with them? Interpretation is the trickiest thing! Annotations are very delicate tools to use and 95% of the time they are unreliable. Raise awareness on these points. GO terms is a strategy and define broader patterns and for conservation could perhaps be better suited than single-gene stories, but comment on papers showing that GO terms can be virtually obtained by chance! 

Associations are very valuable to find potential adaptation to abiotic and biotic factors, and could provide a straightforward evidence for conservation practice (e.g. if a plant is adapted to cold, won’t survive in the heat, etc.). However, stress the fact that this approach, compared to e.g. window-Fst outliers, is totally replying on the SNP calling, and this can lead to high rate of spurious associations. Refer to upcoming Rellstab’s review.




		
All these approaches should help to identify the best management strategy of endangered populations/species.


