#!/usr/bin/env Rscript


options(error=recover)

library(ape)
library(phangorn)
library(vegan)

originalPar <- par()


source("../metrics.r")
removeTreeTipLabelSingleQuotes <- function(tree) {
	tree$tip.label <- gsub("'","",tree$tip.label)
	return(tree)
}

rootTree <- function(tree) {
	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
	}
	return(tree)
}

replicates <- 5
extraReplicates <- 10

darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.3)

aquamarine <- col2rgb("aquamarine4")
transparentaquamarine <- rgb(aquamarine[1]/255,aquamarine[2]/255,aquamarine[3]/255,0.3)


low.data <- read.table("low_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
med.data <- read.table("med_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)
high.data <- read.table("high_sequencing_depth_hmp_data.txt",sep="\t",header=TRUE,row.names=1)

#get metadata (stool vs saliva)
low.groups <- as.factor(gsub("_.*", "", colnames(low.data)))
med.groups <- as.factor(gsub("_.*", "", colnames(med.data)))
high.groups <- as.factor(gsub("_.*", "", colnames(high.data)))

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
source("../../InformationUniFrac.R")
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


low.uwUnifrac.pcoa <- read.table(paste("low_sequencing_depth_pcoa","unweighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
low.wUnifrac.pcoa <- read.table(paste("low_sequencing_depth_pcoa","weighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
low.eUnifrac.pcoa <- read.table(paste("low_sequencing_depth_pcoa","information_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)

med.uwUnifrac.pcoa <- read.table(paste("med_sequencing_depth_pcoa","unweighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
med.wUnifrac.pcoa <- read.table(paste("med_sequencing_depth_pcoa","weighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
med.eUnifrac.pcoa <- read.table(paste("med_sequencing_depth_pcoa","information_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)

high.uwUnifrac.pcoa <- read.table(paste("high_sequencing_depth_pcoa","unweighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
high.wUnifrac.pcoa <- read.table(paste("high_sequencing_depth_pcoa","weighted_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)
high.eUnifrac.pcoa <- read.table(paste("high_sequencing_depth_pcoa","information_Unifrac",sep="_"),sep="\t",header=TRUE,row.names=1)


low.kmeans <- list()
med.kmeans <- list()
high.kmeans <- list()

for (i in c(1:10)){
	# low.kmeans[[i]] <- kmeansClustering(low.otu,low.groups,low.tree,"low_sequencing_depth",uwUnifrac.pcoa,wUnifrac.pcoa,eUnifrac.pcoa)
	# med.kmeans[[i]] <- kmeansClustering(med.otu,med.groups,med.tree,"med_sequencing_depth",uwUnifrac.pcoa,wUnifrac.pcoa,eUnifrac.pcoa)
	# high.kmeans[[i]] <- kmeansClustering(high.otu,high.groups,high.tree,"high_sequencing_depth",uwUnifrac.pcoa,wUnifrac.pcoa,eUnifrac.pcoa)

	low.kmeans[[i]] <- kmeansClustering(low.groups,"low_sequencing_depth",low.uwUnifrac.pcoa,low.wUnifrac.pcoa,low.eUnifrac.pcoa)
	med.kmeans[[i]] <- kmeansClustering(med.groups,"med_sequencing_depth",med.uwUnifrac.pcoa,med.wUnifrac.pcoa,med.eUnifrac.pcoa)
	high.kmeans[[i]] <- kmeansClustering(high.groups,"high_sequencing_depth",high.uwUnifrac.pcoa,high.wUnifrac.pcoa,high.eUnifrac.pcoa)


}

save(low.kmeans,file="low_kmeans_replicates.dat")
save(med.kmeans,file="med_kmeans_replicates.dat")
save(med.kmeans,file="med_kmeans_replicates.dat")

printKmeansMetric <- function(kmeans,metric) {
	uw2 <- metric(kmeans$uwUnifrac[[1]])
	w2 <- metric(kmeans$uwUnifrac[[1]])
	e2 <- metric(kmeans$uwUnifrac[[1]])

	for (i in 1:length(kmeans$uwUnifrac)) {
		print(paste("kmeans of",(i+1),"clusters divided by cluster of 2"))
		print("unweighted")
		print(metric(kmeans$uwUnifrac[[i]])/uw2)

		print("weighted")
		print(metric(kmeans$wUnifrac[[i]])/w2)
		
		print("entropy")
		print(metric(kmeans$eUnifrac[[i]])/e2)

	}

	return("")
	
}

metric <- function(kmeansObject) {
	return(kmeansObject$tot.withinss/kmeansObject$totss)
}

