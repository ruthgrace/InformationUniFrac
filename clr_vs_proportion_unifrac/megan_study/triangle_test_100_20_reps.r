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



# unweightedUnifrac <- read.table("unweightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
# weightedUnifrac <- read.table("weightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
# eUnifrac <- read.table("eUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")


unweightedUnifrac <- read.table("phyloseq_unweightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
weightedUnifrac <- read.table("phyloseq_weightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
eUnifrac <- read.table("phyloseq_eUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")



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

