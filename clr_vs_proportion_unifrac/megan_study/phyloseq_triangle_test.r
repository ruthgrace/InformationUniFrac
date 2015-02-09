################################################################################
# Perform UniFrac on esophagus data
################################################################################
data("esophagus")
(y <- UniFrac(esophagus, TRUE))
UniFrac(esophagus, TRUE, FALSE)
UniFrac(esophagus, FALSE)
################################################################################
# Now try a parallel implementation using doParallel, which leverages the
# new parallel core package in R 2.14.0+
# Note that simply loading the doParallel package is not enough, you must
# call a function that registers the backend. In general, this is pretty easy
# with the doParallel package (or one of the alternative do* packages)
#
# Also note that the esophagus example has only 3 samples, and a relatively small
# tree. This is fast to calculate even sequentially and does not warrant
# parallelized computation, but provides a good quick example for using UniFrac()
# in a parallel fashion. The number of cores you should specify during the
# backend registration, using registerDoParallel(), depends on your system and
# needs. 3 is chosen here for convenience. If your system has only 2 cores, this
# will probably fault or run slower than necessary.
################################################################################
library(doParallel)
data(esophagus)
# For SNOW-like functionality (works on Windows):
cl <- makeCluster(3)
registerDoParallel(cl)
UniFrac(esophagus, TRUE)
# Force to sequential backed:
registerDoSEQ()
# For multicore-like functionality (will probably not work on windows),
# register the backend like this:
registerDoParallel(cores=3)
UniFrac(esophagus, TRUE)
################################################################################



UniFrac() accesses the abundance (otu_table-class) and a phylogenetic tree (phylo-class)
data within an experiment-level (phyloseq-class) object. If the tree and contingency table are
separate objects, suggested solution is to combine them into an experiment-level class using the
phyloseq function. For example, the following code
phyloseq(myotu_table, myTree)
returns a phyloseq-class object that has been pruned and comprises the minimum arguments necessary
for UniFrac().
Parallelization is possible for UniFrac calculated with the phyloseq-package, and is encouraged
in the instances of large trees, many samples, or both. Parallelization has been implemented via the
foreach-package. This means that parallel calls need to be preceded by 2 or more commands that120 UniFrac
register the parallel “backend”. This is acheived via your choice of helper packages. One of the
simplest seems to be the doParallel package.



#getting the otu counts table in

data(esophagus)
x1 = phyloseq(otu_table(esophagus), phy_tree(esophagus))
identical(x1, esophagus)
data(GlobalPatterns)
GP <- GlobalPatterns
phyloseq(sample_data(GP), otu_table(GP))
phyloseq(otu_table(GP), phy_tree(GP))
phyloseq(tax_table(GP), otu_table(GP))
phyloseq(phy_tree(GP), otu_table(GP), sample_data(GP))
phyloseq(otu_table(GP), tax_table(GP), sample_data(GP))
phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))




otu_table(object, taxa_are_rows, errorIfNULL=TRUE)





# get the tree in as normal (phyloseq uses ape phylo class, can use ape read.tree method)

# creating the phyloseq object