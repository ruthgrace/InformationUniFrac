#!/usr/bin/env Rscript

options(error=recover)
library(ape)
library(phangorn)
library(phyloseq)

#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
source("../../GUniFrac.R")
source("../../InformationUniFrac.R")

#### USEFUL FUNCTIONS ####

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

writeDistMat <- function(otuFile,treeFile,metaFile,weightedOutFile,informationOutFile) {

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
	eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable.dat")
	save(brazil.tree,file="processedBrazilTree.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata.dat")
	write.table(gUnifrac,file=weightedOutFile,sep="\t")
	write.table(eUnifrac,file=informationOutFile,sep="\t")

}




writeDistMatNoDiv <- function(otuFile,treeFile,metaFile,weightedOutFile,informationOutFile) {
	source("../../GUniFrac_no_prop_div.R")
	source("../../InformationUniFrac_no_prop_div.r")
	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
	eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable.dat")
	save(brazil.tree,file="processedBrazilTree.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata.dat")
	write.table(gUnifrac,file=weightedOutFile,sep="\t")
	write.table(eUnifrac,file=informationOutFile,sep="\t")

}




writeRuthifracDistMat <- function(otuFile,treeFile,metaFile,unweightedOutFile,weightedOutFile,informationOutFile) {
	source("../../RuthiFrac.r")

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	uwUnifrac <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="unweighted")
	wUnifrac <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="weighted")
	iUnifrac <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="information")


	save(brazil.otu.tab,file="processedRuthBrazilTable.dat")
	save(brazil.tree,file="processedRuthBrazilTree.dat")
	save(MyMetaOrdered,file="processedRuthBrazilMetadata.dat")
	write.table(uwUnifrac,file=unweightedOutFile,sep="\t")
	write.table(wUnifrac,file=weightedOutFile,sep="\t")
	write.table(iUnifrac,file=informationOutFile,sep="\t")

}



writeDistMatNoPruningNoDiv <- function(otuFile,treeFile,metaFile,unweightedOutFile,weightedOutFile,informationOutFile) {
	source("../../GUniFrac_no_tree_pruning_no_prop_div.R")
	source("../../InformationUniFrac_no_tree_pruning_no_prop_div.R")

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	calculatedUnifrac <- GUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))
	uwUnifrac <- calculatedUnifrac$unifrac[,,2]
	wUnifrac <- calculatedUnifrac$unifrac[,,1]
	iUnifrac <- InformationUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable_no_pruning.dat")
	save(brazil.tree,file="processedBrazilTree_no_pruning.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata_no_pruning.dat")
	write.table(uwUnifrac,file=unweightedOutFile,sep="\t")
	write.table(wUnifrac,file=weightedOutFile,sep="\t")
	write.table(iUnifrac,file=informationOutFile,sep="\t")

}



writeDistMatNoPruningProperDiv <- function(otuFile,treeFile,metaFile,unweightedOutFile,weightedOutFile,informationOutFile) {
	source("../../GUniFrac_no_tree_pruning_proper_prop_div.R")
	source("../../InformationUniFrac_no_tree_pruning_proper_prop_div.R")

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	calculatedUnifrac <- GUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))
	uwUnifrac <- calculatedUnifrac$unifrac[,,2]
	wUnifrac <- calculatedUnifrac$unifrac[,,1]
	iUnifrac <- InformationUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable_no_pruning.dat")
	save(brazil.tree,file="processedBrazilTree_no_pruning.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata_no_pruning.dat")
	write.table(uwUnifrac,file=unweightedOutFile,sep="\t")
	write.table(wUnifrac,file=weightedOutFile,sep="\t")
	write.table(iUnifrac,file=informationOutFile,sep="\t")

}



writeDistMatNoPruningSumDiv <- function(otuFile,treeFile,metaFile,unweightedOutFile,weightedOutFile,informationOutFile) {
	source("../../GUniFrac_no_tree_pruning_sum_prop_div.R")
	source("../../InformationUniFrac_no_tree_pruning_sum_prop_div.R")

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	calculatedUnifrac <- GUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))
	uwUnifrac <- calculatedUnifrac$unifrac[,,2]
	wUnifrac <- calculatedUnifrac$unifrac[,,1]
	iUnifrac <- InformationUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable_no_pruning.dat")
	save(brazil.tree,file="processedBrazilTree_no_pruning.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata_no_pruning.dat")
	write.table(uwUnifrac,file=unweightedOutFile,sep="\t")
	write.table(wUnifrac,file=weightedOutFile,sep="\t")
	write.table(iUnifrac,file=informationOutFile,sep="\t")

}




writeDistMatNoPruning <- function(otuFile,treeFile,metaFile,weightedOutFile,informationOutFile) {
	source("../../GUniFrac_no_tree_pruning.R")
	source("../../InformationUniFrac_no_tree_pruning.R")

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- GUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
	eUnifrac <- InformationUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	save(brazil.otu.tab,file="processedBrazilTable_no_pruning_no_prop_div.dat")
	save(brazil.tree,file="processedBrazilTree_no_pruning_no_prop_div.dat")
	save(MyMetaOrdered,file="processedBrazilMetadata_no_pruning_no_prop_div.dat")
	write.table(gUnifrac,file=weightedOutFile,sep="\t")
	write.table(eUnifrac,file=informationOutFile,sep="\t")

}


getProcessedPhyloWeightedDistMat <- function(otu,tree,metadata) {

	brazil.otu.tab <- otu_table(otu, taxa_are_rows=FALSE)
	brazil.tree <- tree
	MyMetaOrdered <- metadata

	phylo.tree <- as.phylo(tree)

	# phyloseqObject <- phyloseq(brazil.otu.tab,brazil.tree)
	phyloseqObject <- phyloseq(brazil.otu.tab,phylo.tree)

	weightedUnifrac <- UniFrac(phyloseqObject, TRUE, TRUE)

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- as.matrix(weightedUnifrac)

	# returnList <- list()
	# returnList$gUnifrac <- gUnifrac
	# return(returnList)
	return(gUnifrac)
}




writePhyloProcessedFiles <- function(otuFile,treeFile,metaFile,weightedOutFile,informationOutFile) {

	# read OTU table and format appropriately for input into UniFrac methods
	brazil.otu.tab <- read.table(otuFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	#remove taxonomy column to make otu count matrix numeric
	taxonomy <- brazil.otu.tab$taxonomy
	brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
	brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

	#sort taxa from most to least abundant
	taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
	taxonomy <- taxonomy[taxaOrder]
	brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

	# read and root tree (rooted tree is required)
	brazil.tree <- read.tree(treeFile)
	brazil.tree <- midpoint(brazil.tree)
	brazil.tree <- reorder(brazil.tree,order="cladewise")

	# read metadata
	MyMeta<- read.table(metaFile, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

	# filter OTU table and metadata so that only samples which appear in both are retained
	otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
	otu_indicies <- otu_indicies[!is.na(otu_indicies)]
	brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
	MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

	save(brazil.otu.tab,file="processedBrazilPhyloTable.dat")
	save(brazil.tree,file="processedBrazilPhyloTree.dat")
	save(MyMetaOrdered,file="processedBrazilPhyloMetadata.dat")

	print("about to run phylo for weighted unifrac")

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- getProcessedPhyloWeightedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)

	print("done phylo for weighted unifrac")

	eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]

	print("done entropy unifrac")

	write.table(gUnifrac,file=weightedOutFile,sep="\t")
	write.table(eUnifrac,file=informationOutFile,sep="\t")
}

writeProcessedDistMat <- function(otu,tree,metadata,weightedOutFile,informationOutFile) {

	#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
	source("../../GUniFrac.R")

	source("../../InformationUniFrac.R")

	brazil.otu.tab <- otu
	brazil.tree <- tree
	MyMetaOrdered <- metadata

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
	eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	write.table(gUnifrac,file=weightedOutFile,sep="\t")
	write.table(eUnifrac,file=informationOutFile,sep="\t")

}



getProcessedDistMat <- function(otu,tree,metadata) {

	brazil.otu.tab <- otu
	brazil.tree <- tree
	MyMetaOrdered <- metadata

	#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
	gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
	eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]


	returnList <- list()
	returnList$gUnifrac <- gUnifrac
	returnList$eUnifrac <- eUnifrac
	return(returnList)
}


writeAllInvalidTriangles <- function(uwUnifrac,wUnifrac,iUnifrac,fileName) {
	samples <- rownames(uwUnifrac)

	fail <- list()
	unweightedText <- "unweighted"
	weightedText <- "weighted"
	informationText <- "information"

	columnNames <- c(paste("index",c(1:3),sep=""),"dist12","dist23","dist13","difference")

	fail$unweightedText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$unweightedText) <- columnNames
	fail$weightedText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$weightedText) <- columnNames
	fail$informationText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$informationText) <- columnNames


	triangles <- combn(c(1:length(samples)),3)
	fail$weightedText <- triangleTest(triangles,uwUnifrac,fail$weightedText)

	fail$weightedText <- triangleTest(triangles,wUnifrac,fail$weightedText)

	fail$informationText <- triangleTest(triangles,iUnifrac,fail$informationText)

	save(fail,file=fileName)

}




writeInvalidTriangles <- function(gUnifrac,eUnifrac,fileName) {
	samples <- rownames(gUnifrac)

	fail <- list()
	unweightedText <- "unweighted"
	weightedText <- "weighted"
	informationText <- "information"

	columnNames <- c(paste("index",c(1:3),sep=""),"dist12","dist23","dist13","difference")

	fail$weightedText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$weightedText) <- columnNames
	fail$informationText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$informationText) <- columnNames


	triangles <- combn(c(1:length(samples)),3)

	fail$weightedText <- triangleTest(triangles,gUnifrac,fail$weightedText)

	fail$informationText <- triangleTest(triangles,eUnifrac,fail$informationText)

	save(fail,file=fileName)

}

getInvalidTriangles <- function(gUnifrac,eUnifrac) {
	samples <- rownames(gUnifrac)

	fail <- list()
	unweightedText <- "unweighted"
	weightedText <- "weighted"
	informationText <- "information"

	columnNames <- c(paste("index",c(1:3),sep=""),"dist12","dist23","dist13","difference")

	fail$weightedText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$weightedText) <- columnNames
	fail$informationText <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail$informationText) <- columnNames


	triangles <- combn(c(1:length(samples)),3)

	fail$weightedText <- triangleTest(triangles,gUnifrac,fail$weightedText)

	fail$informationText <- triangleTest(triangles,eUnifrac,fail$informationText)

	return(fail)

}


getInvalidTrianglesOneTable <- function(gUnifrac) {
	samples <- rownames(gUnifrac)
	
	columnNames <- c(paste("index",c(1:3),sep=""),"dist12","dist23","dist13","difference")

	fail <- data.frame(matrix(ncol=7,nrow=0))
	colnames(fail) <- columnNames

	triangles <- combn(c(1:length(samples)),3)

	fail <- triangleTest(triangles,gUnifrac,fail)

	return(fail)
}

# #### MAKE DISTANCE MATRIX ####

# writeDistMat("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","weightedUnifracDistMat.txt","informationUnifracDistMat.txt")

writeDistMatNoDiv("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","weightedUnifracDistMat_no_prop_div.txt","informationUnifracDistMat_no_prop_div.txt")


# writeDistMatNoPruning("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","weightedUnifracDistMat_no_pruning.txt","informationUnifracDistMat_no_pruning.txt")

# writeRuthifracDistMat("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","unweightedRuthifracDistMat.txt","weightedRuthifracDistMat.txt","informationRuthifracDistMat.txt")

# writeDistMatNoPruningNoDiv("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","unweightedUnifracDistMat_no_pruning_no_prop_div.txt","weightedUnifracDistMat_no_pruning_no_prop_div.txt","informationUnifracDistMat_no_pruning_no_prop_div.txt")

writeDistMatNoPruningProperDiv("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","unweightedUnifracDistMat_no_pruning_proper_prop_div.txt","weightedUnifracDistMat_no_pruning_proper_prop_div.txt","informationUnifracDistMat_no_pruning_proper_prop_div.txt")

writeDistMatNoPruningSumDiv("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","unweightedUnifracDistMat_no_pruning_sum_prop_div.txt","weightedUnifracDistMat_no_pruning_sum_prop_div.txt","informationUnifracDistMat_no_pruning_sum_prop_div.txt")

#### READ IN DISTANCE MATRIX ####

gUnifrac <- read.table("weightedUnifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
eUnifrac <- read.table("informationUnifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

gUnifrac <- read.table("weightedUnifracDistMat_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
eUnifrac <- read.table("informationUnifracDistMat_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)


# gUnifrac <- read.table("weightedUnifracDistMat_no_pruning.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# eUnifrac <- read.table("informationUnifracDistMat_no_pruning.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

# uwUnifrac <- read.table("unweightedRuthifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# wUnifrac <- read.table("weightedRuthifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# iUnifrac <- read.table("informationRuthifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

# uwUnifrac <- read.table("unweightedUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# wUnifrac <- read.table("weightedUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# iUnifrac <- read.table("informationUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

# uwUnifrac <- read.table("unweightedUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# wUnifrac <- read.table("weightedUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
# iUnifrac <- read.table("informationUnifracDistMat_no_pruning_no_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

uwUnifrac <- read.table("unweightedUnifracDistMat_no_pruning_proper_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
wUnifrac <- read.table("weightedUnifracDistMat_no_pruning_proper_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
iUnifrac <- read.table("informationUnifracDistMat_no_pruning_proper_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

uwUnifrac <- read.table("unweightedUnifracDistMat_no_pruning_sum_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
wUnifrac <- read.table("weightedUnifracDistMat_no_pruning_sum_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
iUnifrac <- read.table("informationUnifracDistMat_no_pruning_sum_prop_div.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)


#### CALCULATE ALL TRIANGLES ####

writeInvalidTriangles(gUnifrac,eUnifrac,"invalid_triangles.dat")
writeInvalidTriangles(gUnifrac,eUnifrac,"invalid_triangles_no_prop_div.dat")

# writeAllInvalidTriangles(uwUnifrac,wUnifrac,iUnifrac,"invalid_Ruthifrac_triangles.dat")
# writeAllInvalidTriangles(uwUnifrac,wUnifrac,iUnifrac,"invalid_GUnifrac_no_pruning_no_prop_div_triangles.dat")
# writeAllInvalidTriangles(uwUnifrac,wUnifrac,iUnifrac,"invalid_GUnifrac_no_pruning_proper_prop_div_triangles.dat")
# writeAllInvalidTriangles(uwUnifrac,wUnifrac,iUnifrac,"invalid_GUnifrac_no_pruning_proper_sum_div_triangles.dat")

#### EXPLORE PROBLEMATIC SAMPLES ####

save(invalidTriangles,file="invalid_triangles.rds")
numSamples <- length(rownames(gUnifrac))

# which sample appears the most frequently
sampleNums <- c(fail$weightedText[,1],fail$weightedText[,2],fail$weightedText[,3])
histogram <- hist(sampleNums,breaks=c(1:numSamples))
# max count is 184
#sample with max counts is at index 65
problemSample <- brazil.otu.tab[which(histogram$counts == max(histogram$counts)),]

#does this sample have OTUs that other samples dont have?
maxCount <- apply(brazil.otu.tab,2,max)
which(brazil.otu.tab[65,]==maxCount)
#otu at 17th index ("16") is the one for which sample 65 has maxcount
#	[1] Bacteria;Actinobacteria;Actinobacteria;Coriobacteriales;Coriobacteriaceae;Atopobium;rimae

#try dropping this OTU to see if the number of invalid triangles decreases
treeWithout16 <- drop.tip(tree,"16")
writeProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered,"weightedUnifracDistMatWithoutOTU16.txt","informationUnifracDistMatWithoutOTU16.txt")
gUnifracWithout16 <- read.table("weightedUnifracDistMatWithoutOTU16.txt",sep="\t",header=TRUE,row.names=1,check.names=FALSE)
eUnifracWithout16 <- read.table("informationUnifracDistMatWithoutOTU16.txt",sep="\t",header=TRUE,row.names=1,check.names=FALSE)

writeInvalidTriangles(gUnifracWithout16,eUnifracWithout16,"invalid_triangles_without_OTU16.rds")
load("invalid_triangles_without_OTU16.rds")

# ^ failure. found nothing.

#### TRY REMOVING MOST PROBLEMATIC SAMPLE 66 ####

distMat <- getProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)

desiredSamples <- c(1:65,67:length(rownames(brazil.otu.tab)))
gUnifracWithout65 <- distMat$gUnifrac[desiredSamples,desiredSamples]
eUnifracWithout65 <- distMat$eUnifrac[desiredSamples,desiredSamples]

trianglesWithout65 <- getInvalidTriangles(gUnifracWithout65,eUnifracWithout65)

save(trianglesWithout65,file="invalid_triangles_without_problem_sample_65.dat")	



## remove next most problematic sample ##
numSamples <- length(rownames(gUnifracWithout65))
sampleNums <- c(trianglesWithout65$weightedText[,1],trianglesWithout65$weightedText[,2],trianglesWithout65$weightedText[,3])
histogram <- hist(sampleNums,breaks=c(1:numSamples))
#next sample is index 56
desiredSamples <- c(1:56,58:numSamples)
gUnifracWithout56 <- gUnifracWithout65[desiredSamples,desiredSamples]
eUnifracWithout56 <- eUnifracWithout65[desiredSamples,desiredSamples]

trianglesWithout56 <- getInvalidTriangles(gUnifracWithout56,eUnifracWithout56)

save(trianglesWithout56,file="invalid_triangles_without_problem_sample_56.dat")	



# next problematic sample is not that problematic (only in about 20 of the 80 invalid triangles)
numSamples <- length(rownames(gUnifracWithout56))
sampleNums <- c(trianglesWithout56$weightedText[,1],trianglesWithout56$weightedText[,2],trianglesWithout56$weightedText[,3])
histogram <- hist(sampleNums,breaks=c(1:numSamples))


#### MAKE MINIMAL INVALID TRIANGLE CASE ####
load("invalid_triangles.dat")
load("processedBrazilTable.dat")
load("processedBrazilTree.dat")

biggestFail <- fail$weighted[which(fail$weighted$difference == max(fail$weighted$difference)),]

minCase <- list()
minCase$samples <- brazil.otu.tab[c(biggestFail$index1,biggestFail$index2,biggestFail$index3),]
minCase$tree <- brazil.tree
minCase$distMat <- GUniFrac(minCase$samples,minCase$tree,alpha=c(1))$unifracs[,,"d_1"]
otuSums <- apply(minCase$samples,2,sum)
totalReads <- sum(minCase$samples)
minCase$OnePercentSparsitySamples <- minCase$samples[,which(otuSums >= (0.01*totalReads))]

absentOTUs <- minCase$tree$tip.label[!(minCase$tree$tip.label %in% colnames(minCase$OnePercentSparsitySamples))]
minCase$OnePercentSparsityTree <- drop.tip(minCase$tree, absentOTUs)

save(minCase,file = "min_invalid_triangle_case.dat")

# #### TRY DROPPING TIPS FROM TREE TO SEE IF NUMBER OF INVALID TRIANGLES DECREASE ####

# load("processedBrazilTable.dat")
# load("processedBrazilTree.dat")
# load("processedBrazilMetadata.dat")

# #drop each OTU from the tree and re-run triangle calculation
# invalidTriangles <- list()
# for (index in 1:ncol(brazil.otu.tab)) {
# 	treeMinusOne <- drop.tip(brazil.tree,colnames(brazil.otu.tab)[index])
# 	distMat <- getProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)
# 	invalidTriangles[[index]] <- getInvalidTriangles(distMat$gUnifrac,distMat$eUnifrac)
# }

# save(invalidTriangles,file="invalid_triangles_with_dropped_OTUs.dat")

# # ^fail. same number of invalid triangles

# #### EXPLORE EFFECT OF TREE CHANGES ON INVALID TRIANGLES FOR INFORMATION UNIFRAC ####

# nReplicates <- 10

# polytomyReplicates <- list()

# for(index in 1:nReplicates) {
# 	treeWithPolytomies <- brazil.tree

# 	#introduce polytomy (make 3 random edge lengths zero)
# 	treeWithPolytomies$edge.length[sample(c(1:length(treeWithPolytomies$edge.length)),20)] <- 0

# 	#calculate distance matrix
# 	polytomyDistMat <- getProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)

# 	#run triangle analysis
# 	polytomyReplicates[[index]] <- getInvalidTriangles(polytomyDistMat$gUnifrac,polytomyDistMat$eUnifrac)

# }

# save(polytomyReplicates,file="polytomyInvalidTrianglesReplicates.dat")	

# # ^ fail. the exact same broken triangles still occur (and none for information unifrac)


#### TRY USING PHYLOSEQ UNIFRAC METHOD ####

load("processedBrazilTable.dat")
load("processedBrazilTreeForPhylo.dat")
load("processedBrazilMetadata.dat")

print("test")

writePhyloProcessedFiles("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","weightedUnifracPhyloDistMat.txt","informationUnifracPhyloDistMat.txt")
writeInvalidTriangles(gUnifrac,eUnifrac,"invalid_triangles.dat")


#^ same result as gunifrac

#### TRY PROBABLY BROKEN SELF-WRITTEN UNIFRAC ####
# (unweighted unifrac looks completely different from other methods)

gUnweighted <- read.table("unweighted_general_unifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)
gWeighted <- read.table("weighted_general_unifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)
gInformation <- read.table("information_general_unifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)
rUnweighted <- read.table("unweighted_Ruthifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)
rWeighted <- read.table("weighted_Ruthifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)
rInformation <- read.table("information_Ruthifrac_distance_matrix.txt",sep="\t",check.names=FALSE,row.names=1)

fail <- list()
fail$rUnweighted <- getInvalidTrianglesOneTable(rUnweighted)
fail$rWeighted <- getInvalidTrianglesOneTable(rWeighted)
fail$rInformation <- getInvalidTrianglesOneTable(rInformation)
fail$gUnweighted <- getInvalidTrianglesOneTable(gUnweighted)
fail$gWeighted <- getInvalidTrianglesOneTable(gWeighted)
fail$gInformation <- getInvalidTrianglesOneTable(gInformation)



