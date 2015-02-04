#!/usr/bin/env Rscript

#figuring out why eunifrac produces NA

library(phangorn)

source("../../InformationUniFracTest.r")

#read in everything
unweightedUnifrac <- read.table("unweightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
weightedUnifrac <- read.table("weightedUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")
eUnifrac <- read.table("eUniFracDistanceMatrix.txt",sep="\t",check.names=FALSE,quote="")

MyMeta<- read.table("./meta_analysis_table_oct30_2014 - meta_analysis_table_feb19_fixed.tsv", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

otu_indicies <- match(rownames(MyMeta),rownames(unweightedUnifrac))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
MyMetaOrdered <- MyMeta[otu_indicies,]

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
data.otu.tab <- data.otu.tab[otu_indicies,]


# groups <- MyMetaOrdered$HM_GROUP
# condition <- "heavy_metal"

groups <- MyMetaOrdered$GROUP
condition <- "nutritional_status"

originalUnweighteUnifrac <- unweightedUnifrac
originalWeightedUnifrac <- weightedUnifrac
originalEUnifrac <- eUnifrac

data.tree <- read.tree("./gg_13_5_unannotated.tre")
data.tree <- midpoint(data.tree)

e <- InformationUniFrac(data.otu.tab,data.tree,alpha=c(1))


#figure out which sample copmarisons are problematic


# hasNA <- function(data) {
# 	if (length(which(is.na(data))) > 0) {
# 		return(TRUE)
# 	}
# 	return(FALSE)
# }

# problemSampleIndices <- apply(eUnifrac,1,hasNA)
# problemSampleIndices <- which(problemSampleIndices == TRUE)
# problemSamples <- data.otu.tab[problemSampleIndices,]

# hasOneOTU <- function(data) {
# 	if(length(which(data!=0))==1) {
# 		return(TRUE)
# 	}
# 	return(FALSE)
# }

# hasTwoOTU <- function(data) {
# 	if(length(which(data!=0))==2) {
# 		return(TRUE)
# 	}
# 	return(FALSE)
# }

# hasNoOTU <- function(data) {
# 	if(length(which(data!=0))==0) {
# 		return(TRUE)
# 	}
# 	return(FALSE)
# }

