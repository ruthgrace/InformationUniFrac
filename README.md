InformationUniFrac
==========

Calculate (UniFrac difference)[http://www.mothur.org/wiki/Weighted_UniFrac_algorithm] by weighting by Shannon entropy value per OTU, rather than proportional abundance. Preserves the distance property of UniFrac, while allowing for weighting.

##Information weighting method for UniFrac

####Shannon Entropy
The Shannon entropy is a value that describes, given a probability distribution (modelled by detected proportional abundance), how unlikely one is able to accurately determine the identity of a randomly selected bacterial cell from the sample. In other words, Shannon entropy describes the amount of uncertainty each OTU proportional abundance holds in terms of being able to predict further sampling. For example, if a sample is 100% the same bacteria, then the Shannon entropy is zero, because the identity of the bacteria can certainly be determined. If a sample is 25% bacteria A and 75% bacteria B, the entropy of A is gerater than that of B, because the certainty of sampling an A bacteria is less than a B bacteria. Shannon entropy is calculated for each OTU as follows:
```
(proportional abundance) * log base 2 (1 / proportional abundance)
```

####UniFrac
Unweighted UniFrac distance is calculated using a phylogenetic tree, accounting for which species are present in each sample, and how far apart they are in the tree. Clasically weighted UniFrac is weighted by the proportion of abundance of the species under each node in the tree.

####Information weighting
Information weighted UniFrac is weighted by the Shannon entropy value at each node of the phylogenetic tree, if all the read counts from all the OTUs under the node are treated as if they come from one OTU. 

##Information weighting vs. proportional abundance weighting

####Better separation?

Preliminary results show that InformationUniFrac has a higher and more closely reproducible effect size in a data set with well characterized separation (saliva vs. stool from the Human Microbiome Project). This effect persists even when the data is filtered to have different Shannon diversity, sparsity, or sequencing depth.

##Files

####Files in main folder

#####[InformationUniFrac.r](InformationUniFrac.r)
Generalized UniFrac script with weighting replaced by a Shannon entropy based weighting (details above).

#####[GUniFrac.r](GUniFrac.r)
Generalized UniFrac script ripped straight from the [GUniFrac R package][1]

##Notes

A previous attempt at making a weighted UniFrac that was also a proper distance measure used weighting by the centered log ratio transform, which appeared to separate data with high sequencing depth better in some cases. However, there was a confounding correlation with read count when this method was used on data with lower sequencing depth. That code and tests for it can be found at {this GitHub repository)[https://github.com/ruthgrace/CLRUniFrac].

[1]: http://cran.r-project.org/web/packages/GUniFrac/index.html
