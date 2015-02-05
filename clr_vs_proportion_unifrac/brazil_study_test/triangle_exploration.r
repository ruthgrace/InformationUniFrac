#!/usr/bin/env Rscript

options(error=recover)
library(ape)
library(phangorn)

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

	#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
	source("../../GUniFrac.R")

	source("../../InformationUniFrac.R")


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

	#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
	source("../../GUniFrac.R")

	source("../../InformationUniFrac.R")

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



#### MAKE DISTANCE MATRIX ####

writeDistMat("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt","./brazil_study_data/fasttree_all_seed_OTUs.tre","./brazil_study_data/metadata_BVsamplesonly.txt","weightedUnifracDistMat.txt","informationUnifracDistMat.txt")



#### READ IN DISTANCE MATRIX ####

gUnifrac <- read.table("weightedUnifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)
eUnifrac <- read.table("informationUnifracDistMat.txt",sep="\t",quote="",check.names=FALSE,header=TRUE,row.names=1)

#### CALCULATE ALL TRIANGLES ####

writeInvalidTriangles(gUnifrac,eUnifrac,"invalid_triangles.dat")


# #### EXPLORE PROBLEMATIC SAMPLES ####

# save(invalidTriangles,file="invalid_triangles.rds")
# numSamples <- length(rownames(gUnifrac))

# # which sample appears the most frequently
# sampleNums <- c(fail$weightedText[,1],fail$weightedText[,2],fail$weightedText[,3])
# histogram <- hist(sampleNums,breaks=c(1:numSamples))
# # max count is 184
# #sample with max counts is at index 65
# problemSample <- brazil.otu.tab[which(histogram$counts == max(histogram$counts)),]

# #does this sample have OTUs that other samples dont have?
# maxCount <- apply(brazil.otu.tab,2,max)
# which(brazil.otu.tab[65,]==maxCount)
# #otu at 17th index ("16") is the one for which sample 65 has maxcount
# #	[1] Bacteria;Actinobacteria;Actinobacteria;Coriobacteriales;Coriobacteriaceae;Atopobium;rimae

# #try dropping this OTU to see if the number of invalid triangles decreases
# treeWithout16 <- drop.tip(tree,"16")
# writeProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered,"weightedUnifracDistMatWithoutOTU16.txt","informationUnifracDistMatWithoutOTU16.txt")
# gUnifracWithout16 <- read.table("weightedUnifracDistMatWithoutOTU16.txt",sep="\t",header=TRUE,row.names=1,check.names=FALSE)
# eUnifracWithout16 <- read.table("informationUnifracDistMatWithoutOTU16.txt",sep="\t",header=TRUE,row.names=1,check.names=FALSE)

# writeInvalidTriangles(gUnifracWithout16,eUnifracWithout16,"invalid_triangles_without_OTU16.rds")
# load("invalid_triangles_without_OTU16.rds")

# # ^ failure. found nothing.

#### TRY DROPPING TIPS FROM TREE TO SEE IF NUMBER OF INVALID TRIANGLES DECREASE ####
load("processedBrazilTable.dat")
load("processedBrazilTree.dat")
load("processedBrazilMetadata.dat")

#drop each OTU from the tree and re-run triangle calculation
invalidTriangles <- list()
for (index in 1:ncol(brazil.otu.tab)) {
	treeMinusOne <- drop.tip(brazil.tree,colnames(brazil.otu.tab)[index])
	distMat <- getProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)
	invalidTriangles[[index]] <- getInvalidTriangles(distMat$gUnifrac,distMat$eUnifrac)
}

save(invalidTriangles,file="invalid_triangles_with_dropped_OTUs.dats")


#### EXPLORE EFFECT OF TREE CHANGES ON INVALID TRIANGLES FOR INFORMATION UNIFRAC ####

treeWithPolytomies <- brazil.tree

#introduce polytomy (make 3 random edge lengths zero)
treeWithPolytomies$edge.lengths[sample(c(1:length(treeWithPolytomies$edge.lengths),3))] <- 0

#calculate distance matrix
polytomyDistMat <- getProcessedDistMat(brazil.otu.tab,brazil.tree,MyMetaOrdered)

#run triangle analysis
polytomyInvalidTriangles <- getInvalidTriangles(polytomyDistMat$gUnifrac,polytomyDistMat$eUnifrac)

save(polytomyInvalidTriangles,file="polytomyInvalidTriangles.dat")