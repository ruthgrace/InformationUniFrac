#!/usr/bin/env Rscript

options(error=recover)

library(phangorn)
plotParameters <- par()

triangleTest <- function(triangles,distMat,dataFrame) {
	colNames <- colnames(dataFrame)
	for(i in 1:ncol(triangles)) {

		triangle <- triangles[,i]
		# print(paste("testing sample indices",triangle[1],triangle[2],triangle[3],"distances",distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]]))
		if (distMat[triangle[1],triangle[2]] > (distMat[triangle[2],triangle[3]] + distMat[triangle[1],triangle[3]])) {
			difference <- distMat[triangle[1],triangle[2]] - (distMat[triangle[2],triangle[3]] + distMat[triangle[1],triangle[3]])
			print(paste("samples ",paste(triangle,collapse=" "),"do not fit triangle inequality -- 12 > 23+13: ",distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]]))
			dataFrame <- rbind(dataFrame,c(triangle[1],triangle[2],triangle[3],distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]], difference))
			colnames(dataFrame) <- colNames
		}
		else if (distMat[triangle[2],triangle[3]] > (distMat[triangle[1],triangle[2]] + distMat[triangle[1],triangle[3]])) {
			difference <- distMat[triangle[2],triangle[3]] - (distMat[triangle[1],triangle[2]] + distMat[triangle[1],triangle[3]])
			print(paste("samples ",paste(triangle,collapse=" "),"do not fit triangle inequality -- 23 > 12+13: ",distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]]))
			dataFrame <- rbind(dataFrame,c(triangle[1],triangle[2],triangle[3],distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]], difference))
			colnames(dataFrame) <- colNames
		}
		else if (distMat[triangle[1],triangle[3]] > (distMat[triangle[1],triangle[2]] + distMat[triangle[2],triangle[3]])) {
			difference <- distMat[triangle[1],triangle[3]] - (distMat[triangle[1],triangle[2]] + distMat[triangle[2],triangle[3]])
			print(paste("samples ",paste(triangle,collapse=" "),"do not fit triangle inequality -- 13 > 23+12: ",distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]]))
			dataFrame <- rbind(dataFrame,c(triangle[1],triangle[2],triangle[3],distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]], difference))
			colnames(dataFrame) <- colNames
		}
	}
	return(dataFrame)
}


data.otu.tab <- read.table("./master_table_knight.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
taxa <- rownames(data.otu.tab)
data.otu.tab <- apply(data.otu.tab,2,as.numeric)
rownames(data.otu.tab) <- taxa

#remove taxonomy column to make otu count matrix numeric
taxonomy <- rownames(data.otu.tab)
#data.otu.tab <- data.otu.tab[-length(colnames(data.otu.tab))]
data.otu.tab <- t(as.matrix(data.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(data.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
data.otu.tab <- data.otu.tab[,taxaOrder]



unweightedUnifrac <- read.table("unweightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
weightedUnifrac <- read.table("weightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
eUnifrac <- read.table("eUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")

metadata<- read.table("./meta_analysis_table_oct30_2014 - meta_analysis_table_feb19_fixed.tsv", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

sampleIndices <- match(rownames(metadata),rownames(unweightedUnifrac))
sampleIndices <- sampleIndices[!is.na(sampleIndices)]
metadataOrdered <- metadata[sampleIndices,]

sampleIndices <- match(rownames(data.otu.tab),rownames(unweightedUnifrac))
sampleIndices <- sampleIndices[!is.na(sampleIndices)]
data.otu.tab <- data.otu.tab[sampleIndices,]







samples <- rownames(unweightedUnifrac)

fail <- list()
unweightedText <- "unweighted"
weightedText <- "weighted"
informationText <- "information"

columnNames <- c(paste("index",c(1:3),sep=""),"dist12","dist23","dist13")

fail$unweightedText <- data.frame(matrix(ncol=7,nrow=0))
colnames(fail$unweightedText) <- columnNames
fail$weightedText <- data.frame(matrix(ncol=7,nrow=0))
colnames(fail$weightedText) <- columnNames
fail$informationText <- data.frame(matrix(ncol=7,nrow=0))
colnames(fail$informationText) <- columnNames


for (index in 1:20) {
	triangles <- combn(sample(c(1:length(samples)),100),3)
	#write.table(triangles,file="triangle_combinations.txt",sep="\t")
	#triangles <- read.table("triangle_combinations.txt",sep="\t",check.names=FALSE,header=TRUE,row.names=1)

	print(paste("ROUND",index))

	print("UNWEIGHTED UNIFRAC TRIANGLE TEST")
	fail$unweightedText <- triangleTest(triangles,unweightedUnifrac,fail$unweightedText)

	print("WEIGHTED UNIFRAC TRIANGLE TEST")
	fail$weightedText <- triangleTest(triangles,weightedUnifrac,fail$weightedText)

	print("INFORMATION UNIFRAC TRIANGLE TEST")
	fail$informationText <- triangleTest(triangles,eUnifrac,fail$informationText)
}

save(fail,file="failed_triangle_inequality_with_difference.dat")

































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