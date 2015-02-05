#!/usr/bin/env Rscript

#TEST THIS SCRIPT TO MAKE SURE IT WORKS particularly doing mixed otus

options(error=recover)

library(ape)
library(phangorn)
library(vegan)

originalPar <- par()


source("../metrics.r")

replicates <- 5
extraReplicates <- 10

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.3)

aquamarine <- col2rgb("aquamarine4")
transparentaquamarine <- rgb(aquamarine[1]/255,aquamarine[2]/255,aquamarine[3]/255,0.3)


printPcoaData <- function(otu,groups,tree,fileName) {
	unifrac <- GUniFrac(otu, tree, alpha = c(1))
	uwUnifrac <- unifrac$unifrac[,,1]
	wUnifrac <- unifrac$unifrac[,,3]
	eUnifrac <- InformationUniFrac(otu, tree, alpha = c(1))$unifrac[,,1]

	uwUnifrac.pcoa <- pcoa(uwUnifrac)$vectors
	wUnifrac.pcoa <- pcoa(wUnifrac)$vectors
	eUnifrac.pcoa <- pcoa(eUnifrac)$vectors

	write.table(uwUnifrac.pcoa,file=paste(fileName,"unweighted_Unifrac",sep="_"),sep="\t")
	write.table(wUnifrac.pcoa,file=paste(fileName,"weighted_Unifrac",sep="_"),sep="\t")
	write.table(eUnifrac.pcoa,file=paste(fileName,"information_Unifrac",sep="_"),sep="\t")
}



rootTree <- function(tree) {
	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
	}
	return(tree)
}

removeTreeTipLabelSingleQuotes <- function(tree) {
	tree$tip.label <- gsub("'","",tree$tip.label)
	return(tree)
}




low.data <- read.table("low_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
med.data <- read.table("med_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
high.data <- read.table("high_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)

#get metadata (stool vs saliva)
low.groups <- as.factor(gsub("_[^_]*$", "", colnames(low.data)))
med.groups <- as.factor(gsub("_[^_]*$", "", colnames(med.data)))
high.groups <- as.factor(gsub("_[^_]*$", "", colnames(high.data)))

low.tree <- read.tree("./low_sequencing_depth_subtree.tre")
med.tree <- read.tree("./med_sequencing_depth_subtree.tre")
high.tree <- read.tree("./high_sequencing_depth_subtree.tre")

#get rid of extra quotes on OTU labels
low.tree <- removeTreeTipLabelSingleQuotes(low.tree)
med.tree <- removeTreeTipLabelSingleQuotes(med.tree)
high.tree <- removeTreeTipLabelSingleQuotes(high.tree)

#root tree by midpoint if not rooted
low.tree <- rootTree(low.tree)
med.tree <- rootTree(med.tree)
high.tree <- rootTree(high.tree)

#source("../../CLRUniFrac.R")
source("../../GUniFrac.R")
source("../../InformationUniFrac_original.R")
#source("../../CLRDirichletUniFrac.R")

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
low.data.t <- t(low.data)
med.data.t <- t(med.data)
high.data.t <- t(high.data)

#get rid of any OTUs that aren't in the tree (only a couple reads discarded in total)
low.otu.unordered <- low.data.t[,which(colnames(low.data.t) %in% low.tree$tip.label)]
med.otu.unordered <- med.data.t[,which(colnames(med.data.t) %in% med.tree$tip.label)]
high.otu.unordered <- high.data.t[,which(colnames(high.data.t) %in% high.tree$tip.label)]

#order otus by abundance (least to most)
taxaOrder <- rev(order(apply(low.otu.unordered,2,sum)))
low.otu <- low.otu.unordered[,taxaOrder]
taxaOrder <- rev(order(apply(med.otu.unordered,2,sum)))
med.otu <- med.otu.unordered[,taxaOrder]
taxaOrder <- rev(order(apply(high.otu.unordered,2,sum)))
high.otu <- high.otu.unordered[,taxaOrder]

#get rid of samples that aren't analyzeable (low read count or monoculture)
low.samples <- getAnalyzableSamples(low.data.t)
med.samples <- getAnalyzableSamples(med.data.t)
high.samples <- getAnalyzableSamples(high.data.t)

low.data.t <- low.data.t[low.samples,]
med.data.t <- med.data.t[med.samples,]
high.data.t <- high.data.t[high.samples,]

low.groups <- low.groups[low.samples]
med.groups <- med.groups[med.samples]
high.groups <- high.groups[high.samples]


printPcoaData(low.otu,low.groups,low.tree,"low_sequencing_depth_pcoa")
printPcoaData(med.otu,med.groups,med.tree,"med_sequencing_depth_pcoa")
printPcoaData(high.otu,highgroups,high.tree,"high_sequencing_depth_pcoa")




