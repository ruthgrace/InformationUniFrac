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

pcoaLabels <- c(rep("pcoa1",3),rep("pcoa12",3),rep("pcoa123",3))
unifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),3)
plotDataColNames <- paste(pcoaLabels, unifracLabels, sep = ".")

pcoaText <- c(rep("PCoA first component",3),rep("PCoA first 2 components",3),rep("PCoA first 3 components",3))
unifracText <- rep(c("Unweighted UniFrac","Weighted UniFrac","Information UniFrac"),3)
plotDataTextLabels <- paste(unifracText,pcoaText, sep = "\n")

screeText <- paste("Axis",rep(c(1:5),3))
screeUnifracTest <- rep(c("Unweighted UniFrac","Weighted UniFrac","Information UniFrac"),5)
screeLabels <- paste(screeUnifracTest,screeText)



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

runReplicate <- function(otu,groups,tree,nSamples) {
	
	
	#sample 50 samples from condition 1
	group1.indices <- which(groups==levels(groups)[1])
	group1.rand <- otu[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2
	group2.indices <- which(groups==levels(groups)[2])
	group2.rand <- otu[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups)[1],nSamples),rep(levels(groups)[2],nSamples))
	return(getAllPcoaMetrics(data,newGroups,tree))
}

runMixedReplicate <- function(otu1,otu2,groups1,groups2,tree,nSamples) {	
	#sample 50 samples from condition 1,group1
	group1.indices <- which(groups1==levels(groups1)[1])
	group1.rand <- otu1[as.integer(sample(group1.indices,nSamples,replace=FALSE)),]
	
	#sample 50 samples from condition 2,group2
	group2.indices <- which(groups2==levels(groups2)[2])
	group2.rand <- otu2[as.integer(sample(group2.indices,nSamples,replace=FALSE)),]
	
	if (ncol(group1.rand)!=ncol(group2.rand)) {
		mixedOTUs <- addDisimilarOTUs(group1.rand,group2.rand,tree)
		group1.rand <- mixedOTUs[[1]]
		gropu2.rand <- mixedOTUs[[2]]
	}

	#concatenate
	data <- rbind(group1.rand,group2.rand)

	#make groups
	newGroups <- c(rep(levels(groups1)[1],nSamples),rep(levels(groups2)[2],nSamples))
	return(getAllPcoaMetrics(data,newGroups,tree))
}

addDisimilarOTUs <- function(otu1,otu2,tree) {
	# remove OTUs not in tree
	otu1.tree.otu <- otu1[,which(colnames(otu1) %in% tree$tip.label)]
	otu2.tree.otu <- otu2[,which(colnames(otu2) %in% tree$tip.label)]

	# find indices of otus common and unique to both lists
	otu1.common.otus <- which(colnames(otu1.tree.otu) %in% colnames(otu2.tree.otu))
	otu2.common.otus <- which(colnames(otu2.tree.otu) %in% colnames(otu1.tree.otu))
	otu1.unique.otus <- which(!(colnames(otu1.tree.otu) %in% colnames(otu1.tree.otu)[otu1.common.otus]))
	otu2.unique.otus <- which(!(colnames(otu2.tree.otu) %in% colnames(otu2.tree.otu)[otu2.common.otus]))

	# construct blank otus from otu2 to add to otu1
	otu2.unique.blanks <- otu2.tree.otu[,otu2.unique.otus]
	otu2.unique.blanks[,] <- 0
	if (nrow(otu2)>nrow(otu1)) {
		otu2.unique.blanks <- otu2.unique.blanks[1:nrow(otu1),]
	}
	else {
		otu2.unique.blanks <- otu2.unique.blanks[c(c(1:nrow(otu2)),rep(1,(nrow(otu1)-nrow(otu2)))),]
	}
	newotu1 <- data.frame(otu1.tree.otu,otu2.unique.blanks)

	# construct blank otus from otu1 to add to otu2
	otu1.unique.blanks <- otu1.tree.otu[,otu1.unique.otus]
	otu1.unique.blanks[,] <- 0
	if (nrow(otu1)>nrow(otu2)) {
		otu1.unique.blanks <- otu1.unique.blanks[1:nrow(otu2),]
	}
	else {
		otu1.unique.blanks <- otu1.unique.blanks[c(c(1:nrow(otu1)),rep(1,(nrow(otu2)-nrow(otu1)))),]
	}
	newotu2 <- data.frame(otu2.tree.otu,otu1.unique.blanks[1:nrow(otu2),])

	#return both new otu tables
	returnList <- list()
	returnList[[1]] <- newotu1
	returnList[[2]] <- newotu2
	return(returnList)
}

compare <- function(replicateMethod, otuList, groupList, comparisonList, treeList, comparisonTitleList, fileName, nSamples) {
	pdf(fileName)
	print(paste("length of comparisonList",length(comparisonList)))
	for (i in 1:length(comparisonList)) {
		print(paste("comparison list",i))
		print(comparisonList[[i]])
		compare <- comparisonList[[i]]
		index1 <- compare[1]
		otu1 <- otuList[[index1]]
		group1 <- groupList[[index1]]
		plotTitle <- comparisonTitleList[i]
		tree <- treeList[[i]]
		if (identical(replicateMethod,runReplicate)) {
			getReplicate(replicateMethod,otu1,group1,tree,plotTitle,nSamples)
		}
		else {
			index2 <- compare[2]
			otu2 <- otuList[[index2]]
			group2 <- groupList[[index2]]
			getReplicate(replicateMethod,otu1,group1,tree,plotTitle,nSamples,group2=group2,otu2=otu2)
		}
	}

	dev.off()

}

getReplicate <- function(replicateMethod,otu1,group1,tree,plotTitle,nSamples,group2=NULL,otu2=NULL) {
	#do comparisons
	reps <- list()
	#columns are unifrac, weighted unifrac, info unifrac for each of separation on component 1, 1&2, 1&2&3
	plot.data <- data.frame(matrix(nrow=5,ncol=9))
	colnames(plot.data) <- plotDataTextLabels

	dist <- data.frame(matrix(nrow=5,ncol=9))
	colnames(dist) <- plotDataTextLabels
	sd <- data.frame(matrix(nrow=5,ncol=9))
	colnames(sd) <- plotDataTextLabels

	scree <- data.frame(matrix(nrow=5,ncol=15))
	colnames(scree) <- screeLabels


	for (i in 1:replicates) {
		if (identical(replicateMethod,runReplicate)) {
			reps[[i]] <- runReplicate(otu1,group1,tree,nSamples)
		}
		else {
			reps[[i]] <- runMixedReplicate(otu1,otu2,group1,group2,tree,nSamples)
		}
		plot.data[i,] <- unlist(data.frame(t(reps[[i]]$effect)))
		dist[i,] <- c(reps[[i]]$meanDist.SD[[1]]$meanDist,reps[[i]]$meanDist.SD[[2]]$meanDist,reps[[i]]$meanDist.SD[[3]]$meanDist)
		sd[i,] <- c(reps[[i]]$meanDist.SD[[1]]$error,reps[[i]]$meanDist.SD[[2]]$error,reps[[i]]$meanDist.SD[[3]]$error)
		scree[i,] <- c(reps[[i]]$screeData$uwUnifrac[1:5],reps[[i]]$screeData$wUnifrac[1:5],reps[[i]]$screeData$eUnifrac[1:5])
	}

	#make plots
	originalPar <- par()

	#effect size plot
	par(mar=c(13, 4, 4, 2) + 0.1)
	par(cex.lab=1.3)
	par(cex.main=1.5)
	stripchart(plot.data,vertical=TRUE,main="Sparsity filter at 0.1%",group.names=colnames(plot.data),pch=19,col=transparentdarkorchid,las=2,ylab="Effect size")
	par(originalPar)

	#mean separation with standard deviation plot (averaged over 5 replicates .....)
	par(mar=c(13, 6, 4, 2) + 0.1)

	meanDist <- unlist(lapply(dist,mean))
	meanSd <- unlist(lapply(sd,mean))
	myBarPlot <- barplot(meanDist,col=transparentdarkorchid,las=2,ylim=c(0,1.2),ylab="Difference between mean positions on\nfirst PCoA component between groups",main="Sparsity filter at 0.1%")
	segments(myBarPlot, meanDist - meanSd, myBarPlot, meanDist + meanSd, lwd=2)
	segments(myBarPlot - 0.1, meanDist - meanSd, myBarPlot + 0.1, meanDist - meanSd, lwd=2)
	segments(myBarPlot - 0.1, meanDist + meanSd, myBarPlot + 0.1, meanDist + meanSd, lwd=2)

	# scree plots
	par(originalPar)
	par(mar=c(12, 6, 4, 2) + 0.1)
	screePlotData <- apply(scree,2,mean)
	screePlotError <- apply(scree,2,sd)
	myBarPlot <- barplot(screePlotData,col=transparentdarkorchid,las=2,ylim=c(0,1),ylab="Variation explained by each axis of PCoA",main="Sparsity filter at 0.1%")
	segments(myBarPlot, screePlotData - screePlotError, myBarPlot, screePlotData + screePlotError, lwd=2)
	segments(myBarPlot - 0.1, screePlotData - screePlotError, myBarPlot + 0.1, screePlotData - screePlotError, lwd=2)
	segments(myBarPlot - 0.1, screePlotData + screePlotError, myBarPlot + 0.1, screePlotData + screePlotError, lwd=2)

	par(originalPar)
	# pcoa plots, only plot first replicate in each comparison
	palette(c(transparentdarkorchid,transparentaquamarine,"blue","black"))
	pcoaGroups <- as.factor(c(rep(1,nSamples),rep(2,nSamples)))
	plot(reps[[1]]$pcoa$uwUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$uwUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$uwUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(reps[[1]]$pcoa$wUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$wUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$wUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(reps[[1]]$pcoa$eUnifrac$vectors[,1],reps[[1]]$pcoa$uwUnifrac$vectors[,2], type="p",col=pcoaGroups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(reps[[1]]$screeData$eUnifrac[1],digits=3),"variance explained"),ylab=paste("Second Component", round(reps[[1]]$screeData$eUnifrac[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)


	par(originalPar)
	# return replicate
	return(reps)
}

#all CLR DIRICHLET commented out while the method is being fixed.
#attempt at sparsity filter instead -- 30 counts minimum is what is stable


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



selfComparisonList3 <- list()
selfComparisonList3[[1]] <- c(1)
selfComparisonList3[[2]] <- c(2)
selfComparisonList3[[3]] <- c(3)

mixedComparisonList3 <- list()
mixedComparisonList3[[1]] <- c(1,2)
mixedComparisonList3[[2]] <- c(1,3)
mixedComparisonList3[[3]] <- c(2,3)

selfComparisonList2 <- list()
selfComparisonList2[[1]] <- c(1)
selfComparisonList2[[2]] <- c(2)

mixedComparisonList2 <- list()
mixedComparisonList2[[1]] <- c(1,2)
mixedComparisonList2[[2]] <- c(2,1)

# sparsityLabels <- c(rep("sparsity.001",3),rep("sparsity.0001",3),rep("sparsity.00001",3))
# plotSparsityDataColNames <- paste(sparsityLabels, unifracLabels, sep = ".")

# diversityLabels <- c(rep("low.diversity",3),rep("high.diversity",3))
# shortUnifracLabels <- rep(c("uwUnifrac","wUnifrac","iUnifrac"),2)
# plotDiversityDataColNames <- paste(diversityLabels, shortUnifracLabels, sep = ".")


# SEQUENCING DEPTH TEST

depth <- list()

depth$low <- low.otu
depth$med <- med.otu
depth$high <- high.otu

depth.groups <- list()
depth.groups$low <- low.groups
depth.groups$med <- med.groups
depth.groups$high <- high.groups

depth.tree <- list()
depth.tree$low <- low.tree
depth.tree$med <- med.tree
depth.tree$high <- high.tree

comparisonTitleList <- c("Sequencing depth < 3000 reads/sample","Sequencing depth 3000-6000 reads/sample","Sequencing depth > 6000 reads/sample")
fileName <- "sequencingDepthPlots.pdf"
nSamples <- 50

compare(runReplicate, depth, depth.groups, selfComparisonList3, depth.tree, comparisonTitleList, fileName, nSamples)


# SEQUENCING DEPTH DIFFERENCE TEST

depth.diff.tree <- list()
depth.diff.tree[[1]] <- high.tree
depth.diff.tree[[2]] <- high.tree
depth.diff.tree[[3]] <- high.tree

comparisonTitleList <- c("Sequencing depth < 3000 vs. 3000-6000 reads/sample","Sequencing depth 3000-6000 vs. > 6000 reads/sample","Sequencing depth 3000-6000 vs. > 6000 reads/sample")
fileName <- "sequencingDepthDifferenceTestPlots.pdf"
nSamples <- 50

compare(runMixedReplicate, depth, depth.groups, mixedComparisonList3, depth.diff.tree, comparisonTitleList, fileName, nSamples)



# # SPARSITY TEST

#remove OTUs rarer than a thresh hold throughout all samples
high.otu.sum <- apply(high.otu,2,sum)
high.total.sum <- sum(high.otu)

sparse <- list()

sparse$otu.001 <- high.otu[,(which(high.otu.sum >= (0.001*high.total.sum)))]
sparse$otu.0001 <- high.otu[,(which(high.otu.sum >= (0.0001*high.total.sum)))]
sparse$otu.00001 <- high.otu[,(which(high.otu.sum >= (0.00001*high.total.sum)))]

sparse.groups <- list()
sparse.groups[[1]] <- sparse.groups[[2]] <- sparse.groups[[3]] <- high.groups

sparse.tree <- list()
sparse.tree[[1]] <- sparse.tree[[2]] <- sparse.tree[[3]] <- high.tree

comparisonTitleList <- c("Sparsity filter at 0.1%","Sparsity filter at 0.01%","Sparsity filter at 0.001%")
fileName <- "sparsityTestPlots.pdf"
nSamples <- 50

compare(runReplicate, sparse, sparse.groups, selfComparisonList3, sparse.tree, comparisonTitleList, fileName, nSamples)



# SPARSITY DIFFERENCE TEST

sparse.diff <- list()

sparse.diff$otu.001 <- high.otu
sparse.diff$otu.001[,(which(high.otu.sum < (0.001*high.total.sum)))] <- 0
sparse.diff$otu.0001 <- high.otu
sparse.diff$otu.0001[,(which(high.otu.sum < (0.0001*high.total.sum)))] <- 0
sparse.diff$otu.00001 <- high.otu
sparse.diff$otu.00001[,(which(high.otu.sum < (0.00001*high.total.sum)))] <- 0


sparse.diff.groups <- list()
sparse.diff.groups[[1]] <- sparse.diff.groups[[2]] <- sparse.diff.groups[[3]] <- high.groups

sparse.diff.tree <- list()
sparse.diff.tree[[1]] <- sparse.diff.tree[[2]] <- sparse.diff.tree[[3]] <- high.tree

comparisonTitleList <- c("Sparsity filter at 0.1% vs 0.01%","Sparsity filter at 0.1% vs 0.001%","Sparsity filter at 0.01% vs 0.001%")
fileName <- "sparsityDifferenceTestPlots.pdf"
nSamples <- 50

compare(runMixedReplicate, sparse.diff, sparse.diff.groups, mixedComparisonList3, sparse.diff.tree, comparisonTitleList, fileName, nSamples)



# SHANNON DIVERSITY TEST

#note that average diversity isn't the same in the different conditions
#	medians are 6.216 for saliva, 5.624 for stool
high.otu.diversity <- diversity(high.otu)

#low/high diversity cutoffs set so that there are at least 10 samples in each condition
#	less samples -> extra replicates
high.diversity.indices <- which(high.otu.diversity > 6)
low.diversity.indices <- which(high.otu.diversity < 5.7)
high.diversity <- high.otu[high.diversity.indices,]
low.diversity <- high.otu[low.diversity.indices,]
high.diversity.groups <- high.groups[high.diversity.indices]
low.diversity.groups <- high.groups[low.diversity.indices]


diversity <- list()

diversity$low <- low.diversity
diversity$high <- high.diversity

diversity.groups <- list()
diversity.groups$low <- low.diversity.groups
diversity.groups$high <- high.diversity.groups

diversity.tree <- list()
diversity.tree[[1]] <- diversity.tree[[2]] <- high.tree

comparisonTitleList <- c("Shannon diversity < 5.7","Shannon diversity > 6")
fileName <- "diversityTestPlots.pdf"
nSamples <- 10

compare(runReplicate, diversity, diversity.groups, selfComparisonList2, diversity.tree, comparisonTitleList, fileName, nSamples)


# SHANNON DIVERSITY DIFFERENCE TEST

comparisonTitleList <- c("Saliva diversity < 5.7 vs. stool diversity > 6.1","Stool diversity < 5.7 vs. saliva diversity > 6.1")
fileName <- "diversityDifferenceTestPlots.pdf"
nSamples <- 10

compare(runMixedReplicate, diversity, diversity.groups, mixedComparisonList2, diversity.tree, comparisonTitleList, fileName, nSamples)
