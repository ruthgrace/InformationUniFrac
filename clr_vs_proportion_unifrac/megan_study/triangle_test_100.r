#!/usr/bin/env Rscript

options(error=recover)

library(phangorn)
plotParameters <- par()

triangleTest <- function(triangles,distMat) {
	for(i in 1:ncol(triangles)) {

		triangle <- triangles[,i]
		# print(paste("testing sample indices",triangle[1],triangle[2],triangle[3],"distances",distMat[triangle[1],triangle[2]],distMat[triangle[2],triangle[3]],distMat[triangle[1],triangle[3]]))
		if (distMat[triangle[1],triangle[2]] > (distMat[triangle[2],triangle[3]] + distMat[triangle[1],triangle[3]])) {
			print(paste("samples ",triangle,"do not fit triangle inequality"))
		}
		else if (distMat[triangle[2],triangle[3]] > (distMat[triangle[1],triangle[2]] + distMat[triangle[1],triangle[3]])) {
			print(paste("samples ",triangle,"do not fit triangle inequality"))
		}
		else if (distMat[triangle[1],triangle[3]] > (distMat[triangle[1],triangle[2]] + distMat[triangle[2],triangle[3]])) {
			print(paste("samples ",triangle,"do not fit triangle inequality"))
		}
	}
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

MyMeta<- read.table("./meta_analysis_table_oct30_2014 - meta_analysis_table_feb19_fixed.tsv", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

sampleIndices <- match(rownames(MyMeta),rownames(unweightedUnifrac))
sampleIndices <- sampleIndices[!is.na(sampleIndices)]
MyMetaOrdered <- MyMeta[sampleIndices,]

sampleIndices <- match(rownames(data.otu.tab),rownames(unweightedUnifrac))
sampleIndices <- sampleIndices[!is.na(sampleIndices)]
data.otu.tab <- data.otu.tab[sampleIndices,]

samples <- rownames(unweightedUnifrac)
triangles <- combn(c(1:100),3)
write.table(triangles,file="triangle_combinations.txt",sep="\t")

#triangles <- read.table("triangle_combinations.txt",sep="\t",check.names=FALSE,header=TRUE,row.names=1)

print("UNWEIGHTED UNIFRAC TRIANGLE TEST")
triangleTest(triangles,unweightedUnifrac)

print("WEIGHTED UNIFRAC TRIANGLE TEST")
triangleTest(triangles,weightedUnifrac)

print("INFORMATION UNIFRAC TRIANGLE TEST")
triangleTest(triangles,eUnifrac)



hasOneOrLessOTU <- function(data) {
	if(length(which(data!=0))<=1) {
		return(TRUE)
	}
	return(FALSE)
}

oneOrLessOTUSamples <- apply(data.otu.tab,1,hasOneOrLessOTU)
moreThanOneOTUSamples <- which(oneOrLessOTUSamples==FALSE)
oneOrLessOTUSamples <- which(oneOrLessOTUSamples==TRUE)
if (length(oneOrLessOTUSamples)>0) {
	data.otu.tab <- data.otu.tab[moreThanOneOTUSamples,]
	MyMetaOrdered <- MyMetaOrdered[moreThanOneOTUSamples,]
	unweightedUnifrac <- unweightedUnifrac[moreThanOneOTUSamples,moreThanOneOTUSamples]
	weightedUnifrac <- weightedUnifrac[moreThanOneOTUSamples,moreThanOneOTUSamples]
	eUnifrac <- eUnifrac[moreThanOneOTUSamples,moreThanOneOTUSamples]
}





