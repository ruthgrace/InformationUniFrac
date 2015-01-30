#make overlap matrix

library(vegan)

getOverlap <- function (otu.tab) {	

	# Convert into CLR
	otu.tab <- as.matrix(otu.tab)
	otu.tab.original <- otu.tab
	otu.tab[otu.tab==0] <- 0.5
	otu.tab.clr <- apply(otu.tab, 1, function(x){log2(x) - mean(log2(x))})
	numCols <- ncol(otu.tab.clr)

	overlap <- matrix(nrow=numCols,ncol=numCols)

	minval <- min(otu.tab.clr)

	otu.tab.clr.positive <- otu.tab.clr - minval

	for(i in 1:numCols) {
		for (j in i:numCols) {
			compare <- data.frame(otu.tab.clr.positive[,i],otu.tab.clr.positive[,j])
			overlap[i,j] <- sum(apply(compare,1,min)) / sum(apply(compare,1,max))
			overlap[j,i] <- overlap[i,j]
		}
	}

	return(overlap)
}


averageReadCount <- function(otu.tab) {
	otu.tab <- as.matrix(otu.tab)
	totalReadCount <- rowSums(otu.tab,na.rm=TRUE)
	numSamples <- nrow(otu.tab)

	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	return(averageGeneric(numSamples,totalReadCount))
}

averageGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- mean(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

maxGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- max(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

minGeneric <- function(numSamples,totalReadCount) {
	avgReadCount <- matrix(nrow=numSamples,ncol=numSamples)

	for(i in 1:numSamples) {
		for (j in i:numSamples) {
			avgReadCount[i,j] <- min(c(totalReadCount[i],totalReadCount[j]))
			avgReadCount[j,i] <- avgReadCount[i,j]
		}
	}
	return(avgReadCount)
}

printSeparation <- function(ruthClrUnifrac.pcoa,gUnifrac.pcoa,eUnifrac.pcoa,condition1,condition2,groups) {

	print(paste("COMPARING",condition1,"and",condition2))

	#calculate separation clrunifrac
	group1.indices <- which(groups==condition1)
	group2.indices <- which(groups==condition2)

	sd.1 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,1]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,1]))
	sd.2 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,2]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,2]))
	sd.3 <- max(sd(ruthClrUnifrac.pcoa$vectors[group1.indices,3]),sd(ruthClrUnifrac.pcoa$vectors[group2.indices,3]))

	group1.1 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,1])/sd.1
	group2.1 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,1])/sd.1

	group1.2 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,2])/sd.2
	group2.2 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,2])/sd.2

	group1.3 <- mean(ruthClrUnifrac.pcoa$vectors[group1.indices,3])/sd.3
	group2.3 <- mean(ruthClrUnifrac.pcoa$vectors[group2.indices,3])/sd.3

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)



	print(paste("CLRUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("CLRUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("CLRUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))


	#calculate separation gunifrac
	group1.1 <- mean(gUnifrac.pcoa$vectors[group1.indices,1])
	group2.1 <- mean(gUnifrac.pcoa$vectors[group2.indices,1])

	group1.2 <- mean(gUnifrac.pcoa$vectors[group1.indices,2])
	group2.2 <- mean(gUnifrac.pcoa$vectors[group2.indices,2])

	group1.3 <- mean(gUnifrac.pcoa$vectors[group1.indices,3])
	group2.3 <- mean(gUnifrac.pcoa$vectors[group2.indices,3])

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)

	print(paste("gUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("gUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("gUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))


	#calculate separation eunifrac
	group1.1 <- mean(eUnifrac.pcoa$vectors[group1.indices,1])
	group2.1 <- mean(eUnifrac.pcoa$vectors[group2.indices,1])

	group1.2 <- mean(eUnifrac.pcoa$vectors[group1.indices,2])
	group2.2 <- mean(eUnifrac.pcoa$vectors[group2.indices,2])

	group1.3 <- mean(eUnifrac.pcoa$vectors[group1.indices,3])
	group2.3 <- mean(eUnifrac.pcoa$vectors[group2.indices,3])

	group1.12 <- sqrt(group1.1^2 + group1.2^2)
	group2.12 <- sqrt(group2.1^2 + group2.2^2)

	group1.123 <- sqrt(group1.12^2 + group1.3^2)
	group2.123 <- sqrt(group2.12^2 + group2.3^2)

	print(paste("eUniFrac 1st component separation: ",abs(group1.1-group2.1)))
	print(paste("eUniFrac 1st and 2nd component separation: ",abs(group1.12-group2.12)))
	print(paste("eUniFrac 1,2,3 component separation: ",abs(group1.123-group2.123)))

}

getShannonDiversityDiffMat <- function(otu) {
	diversityList <- diversity(otu)
	#put diversity into rows
	diversityList <- t(diversityList)
	diversityList <- t(diversityList)
	distMat <- dist(diversityList,method="manhattan",diag=TRUE,upper=TRUE)
	return(distMat)
}

getAvgShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(averageGeneric(nrow(otu),diversityList))
}

getMinShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(minGeneric(nrow(otu),diversityList))
}

getMaxShannonDiversity <- function(otu) {
	diversityList <- diversity(otu)
	return(maxGeneric(nrow(otu),diversityList))
}

getDataSetSep <- function(otu,groups,tree) {
	#take in otu table, rows are samples, cols are OTUs
	# return list of 1) unweighted UniFrac 2) weighted UniFrac 3) iUniFrac
	unifrac <- GUniFrac(otu, tree, alpha = c(1))
	uwUnifrac <- unifrac$unifrac[,,1]
	wUnifrac <- unifrac$unifrac[,,3]
	eUnifrac <- InformationUniFrac(otu, tree, alpha = c(1))$unifrac[,,1]

	uwUnifrac.pcoa <- pcoa(uwUnifrac)
	wUnifrac.pcoa <- pcoa(wUnifrac)
	eUnifrac.pcoa <- pcoa(eUnifrac)

	uwUnifrac.sep <- getPCoASep(uwUnifrac.pcoa,groups)
	wUnifrac.sep <- getPCoASep(wUnifrac.pcoa,groups)
	eUnifrac.sep <- getPCoASep(eUnifrac.pcoa,groups)

	returnList <- data.frame(t(uwUnifrac.sep),t(wUnifrac.sep),t(eUnifrac.sep))
	colnames(returnList) <- c("uwUnifrac","wUnifrac","eUnifrac")
	return(returnList)
}

getAllPcoaMetrics <- function(otu,groups,tree) {
	unifrac <- GUniFrac(otu, tree, alpha = c(1))
	uwUnifrac <- unifrac$unifrac[,,1]
	wUnifrac <- unifrac$unifrac[,,3]
	eUnifrac <- InformationUniFrac(otu, tree, alpha = c(1))$unifrac[,,1]

	uwUnifrac.pcoa <- pcoa(uwUnifrac)
	wUnifrac.pcoa <- pcoa(wUnifrac)
	eUnifrac.pcoa <- pcoa(eUnifrac)

	uwUnifrac.sep <- getPCoASep(uwUnifrac.pcoa,groups)
	wUnifrac.sep <- getPCoASep(wUnifrac.pcoa,groups)
	eUnifrac.sep <- getPCoASep(eUnifrac.pcoa,groups)

	returnList <- list()
	returnList$effect <- data.frame(t(uwUnifrac.sep),t(wUnifrac.sep),t(eUnifrac.sep))
	colnames(returnList$effect) <- c("uwUnifrac","wUnifrac","eUnifrac")

	uwUnifrac.meanDist <- getMeanDistanceWithErrorData(uwUnifrac.pcoa,groups)
	wUnifrac.meanDist <- getMeanDistanceWithErrorData(wUnifrac.pcoa,groups)
	eUnifrac.meanDist <- getMeanDistanceWithErrorData(eUnifrac.pcoa,groups)

	returnList$meanDist.SD <- list()
	returnList$meanDist.SD$uwUnifrac <- uwUnifrac.meanDist
	returnList$meanDist.SD$wUnifrac <- wUnifrac.meanDist
	returnList$meanDist.SD$eUnifrac <- eUnifrac.meanDist

	uwUnifrac.screeData <- 	getScreePlotData(uwUnifrac.pcoa)
	wUnifrac.screeData <- 	getScreePlotData(wUnifrac.pcoa)
	eUnifrac.screeData <- 	getScreePlotData(eUnifrac.pcoa)

	returnList$screeData <- list()
	returnList$screeData$uwUnifrac <- uwUnifrac.screeData
	returnList$screeData$wUnifrac <- wUnifrac.screeData
	returnList$screeData$eUnifrac <- eUnifrac.screeData

	returnList$pcoa <- list()
	returnList$pcoa$uwUnifrac <- uwUnifrac.pcoa
	returnList$pcoa$wUnifrac <- wUnifrac.pcoa
	returnList$pcoa$eUnifrac <- eUnifrac.pcoa

	return(returnList)
}



# kmeansClustering <- function(otu,groups,tree,filename,uwUnifrac.pcoa,wUnifrac.pcoa,eUnifrac.pcoa) {
kmeansClustering <- function(groups,filename,uwUnifrac.pcoa,wUnifrac.pcoa,eUnifrac.pcoa) {
	# pdf(paste(filename,".pdf",sep=""))


	# #Generate UniFrac distance matrices
	# unifrac <- GUniFrac(otu, tree, alpha = c(1))
	# uwUnifrac <- unifrac$unifrac[,,1]
	# wUnifrac <- unifrac$unifrac[,,3]
	# eUnifrac <- InformationUniFrac(otu, tree, alpha = c(1))$unifrac[,,1]

	# write.table(uwUnifrac,file=paste(filename,"unweighted_distance.mat",sep="_"),sep="\t",quote=FALSE)
	# write.table(wUnifrac,file=paste(filename,"weighted_distance.mat",sep="_"),sep="\t",quote=FALSE)
	# write.table(eUnifrac,file=paste(filename,"entropy_distance.mat",sep="_"),sep="\t",quote=FALSE)





	uwUnifrac.totalVarExplained <- sum(apply(uwUnifrac.pcoa,2,function(x) sd(x)*sd(x)))
	uwUnifrac.varEx <- apply(uwUnifrac.pcoa,2,function(x) sd(x)*sd(x)*uwUnifrac.totalVarExplained)

	wUnifrac.totalVarExplained <- sum(apply(wUnifrac.pcoa,2,function(x) sd(x)*sd(x)))
	wUnifrac.varEx <- apply(wUnifrac.pcoa,2,function(x) sd(x)*sd(x)*wUnifrac.totalVarExplained)

	eUnifrac.totalVarExplained <- sum(apply(eUnifrac.pcoa,2,function(x) sd(x)*sd(x)))
	eUnifrac.varEx <- apply(eUnifrac.pcoa,2,function(x) sd(x)*sd(x)*eUnifrac.totalVarExplained)

	uwUnifrac.pcoa <- t(t(uwUnifrac.pcoa)*uwUnifrac.varEx)
	wUnifrac.pcoa <- t(t(wUnifrac.pcoa)*wUnifrac.varEx)
	eUnifrac.pcoa <- t(t(eUnifrac.pcoa)*eUnifrac.varEx)


	# #plot original groups
	# palette(c(transparentdarkorchid,transparentaquamarine,"blue","black"))
	# groups <- as.factor(groups)





	# plot(uwUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(wUnifrac.pcoa[,1],wUnifrac.pcoa[,2], type="p",col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(eUnifrac.pcoa[,1],eUnifrac.pcoa[,2], type="p",col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	#kmeans
	

	#plot 2 clusters
	uwUnifrac.kmeans2 <- kmeans(uwUnifrac.pcoa,centers=2)
	wUnifrac.kmeans2 <- kmeans(wUnifrac.pcoa,centers=2)
	eUnifrac.kmeans2 <- kmeans(eUnifrac.pcoa,centers=2)

	# plot(uwUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(uwUnifrac.kmeans2$cluster),main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(wUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(wUnifrac.kmeans2$cluster),main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(eUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(eUnifrac.kmeans2$cluster),main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	# print("kmeans 2 within sum of squares:")
	# print(uwUnifrac.kmeans2$tot.withinss)
	# print(wUnifrac.kmeans2$tot.withinss)
	# print(eUnifrac.kmeans2$tot.withinss)

	#plot 3 clusters
	uwUnifrac.kmeans3 <- kmeans(uwUnifrac.pcoa,centers=3)
	wUnifrac.kmeans3 <- kmeans(wUnifrac.pcoa,centers=3)
	eUnifrac.kmeans3 <- kmeans(eUnifrac.pcoa,centers=3)

	# plot(uwUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(uwUnifrac.kmeans3$cluster),main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(wUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(wUnifrac.kmeans3$cluster),main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(eUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(eUnifrac.kmeans3$cluster),main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	# print("kmeans 3 within sum of squares:")
	# print(uwUnifrac.kmeans3$tot.withinss)
	# print(wUnifrac.kmeans3$tot.withinss)
	# print(eUnifrac.kmeans3$tot.withinss)

	#plot 4 clusters
	uwUnifrac.kmeans4 <- kmeans(uwUnifrac.pcoa,centers=4)
	wUnifrac.kmeans4 <- kmeans(wUnifrac.pcoa,centers=4)
	eUnifrac.kmeans4 <- kmeans(eUnifrac.pcoa,centers=4)

	# plot(uwUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(uwUnifrac.kmeans4$cluster),main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(wUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(wUnifrac.kmeans4$cluster),main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(eUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(eUnifrac.kmeans4$cluster),main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	# print("kmeans 4 within sum of squares:")
	# print(uwUnifrac.kmeans4$tot.withinss)
	# print(wUnifrac.kmeans4$tot.withinss)
	# print(eUnifrac.kmeans4$tot.withinss)


	#plot 5 clusters
	uwUnifrac.kmeans5 <- kmeans(uwUnifrac.pcoa,centers=5)
	wUnifrac.kmeans5 <- kmeans(wUnifrac.pcoa,centers=5)
	eUnifrac.kmeans5 <- kmeans(eUnifrac.pcoa,centers=5)

	# plot(uwUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(uwUnifrac.kmeans5$cluster),main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(uwUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(uwUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(wUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(wUnifrac.kmeans5$cluster),main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(wUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(wUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	# plot(eUnifrac.pcoa[,1],uwUnifrac.pcoa[,2], type="p",col=as.factor(eUnifrac.kmeans5$cluster),main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	# print("kmeans 5 within sum of squares:")
	# print(uwUnifrac.kmeans5$tot.withinss)
	# print(wUnifrac.kmeans5$tot.withinss)
	# print(eUnifrac.kmeans5$tot.withinss)



	# dev.off()

	#get kmeans metric
	returnList <- list()

	returnList$uwUnifrac <- list()
	returnList$uwUnifrac[[1]] <- uwUnifrac.kmeans2
	returnList$uwUnifrac[[2]] <- uwUnifrac.kmeans3
	returnList$uwUnifrac[[3]] <- uwUnifrac.kmeans4
	returnList$uwUnifrac[[4]] <- uwUnifrac.kmeans5
	
	returnList$wUnifrac <- list()
	returnList$wUnifrac[[1]] <- wUnifrac.kmeans2
	returnList$wUnifrac[[2]] <- wUnifrac.kmeans3
	returnList$wUnifrac[[3]] <- wUnifrac.kmeans4
	returnList$wUnifrac[[4]] <- wUnifrac.kmeans5
	
	returnList$eUnifrac <- list()
	returnList$eUnifrac[[1]] <- eUnifrac.kmeans2
	returnList$eUnifrac[[2]] <- eUnifrac.kmeans3
	returnList$eUnifrac[[3]] <- eUnifrac.kmeans4
	returnList$eUnifrac[[4]] <- eUnifrac.kmeans5
	

	return(returnList)
}

getScreePlotData <- function(pcoa){
	varExplained <- sum(apply(pcoa$vector,2,function(x) sd(x)*sd(x)))
	varExplainedByComponent <- apply(pcoa$vector,2,function(x) sd(x)*sd(x)/varExplained)
	return(varExplainedByComponent)
}

getMeanDistanceWithErrorData <- function(pcoa,groups) {
	groups <- as.factor(groups)
	#given pcoa and metadata
	# returns separations on axis 1, 1&2, 1&2&3
	group1.1 <- pcoa$vectors[which(groups==levels(groups)[1]),1]
	group2.1 <- pcoa$vectors[which(groups==levels(groups)[2]),1]
	diff.1 <- abs(mean(group1.1) - mean(group2.1))
	sd.1 <- sd(pcoa$vector[,1])

	group1.2 <- pcoa$vectors[which(groups==levels(groups)[1]),2]
	group2.2 <- pcoa$vectors[which(groups==levels(groups)[2]),2]
	diff.2 <- abs(mean(group1.2) - mean(group2.2))
	sd.2 <- sd(pcoa$vector[,2])

	group1.3 <- pcoa$vectors[which(groups==levels(groups)[1]),3]
	group2.3 <- pcoa$vectors[which(groups==levels(groups)[2]),3]
	diff.3 <- abs(mean(group1.3) - mean(group2.3))
	sd.3 <- sd(pcoa$vector[,3])

	diff.12 <- sqrt((diff.1^2) + (diff.2^2))
	diff.123 <- sqrt((diff.12^2) + (diff.3^2))

	diff.12.uncondensed <- sqrt(pcoa$vector[,1]^2 + pcoa$vector[,2]^2)
	diff.123.uncondensed <- sqrt(diff.12.uncondensed^2 + pcoa$vector[,3]^2)
	sd.12 <- sd(diff.12.uncondensed)
	sd.123 <- sd(diff.123.uncondensed)

	returnList <- list()

	separation <- data.frame(c(diff.1,diff.12,diff.123))
	rownames(separation) <- c("separationOn1","separationOn12","separationOn123")
	separation <- t(separation)
	returnList$meanDist <- separation

	error <- data.frame(c(sd.1,sd.12,sd.123))
	rownames(error) <- c("standardDeviationOn1","standardDeviationOn12","standardDeviationOn123")
	error <- t(error)
	returnList$error <- error
	return(returnList)
}

getPCoASep <- function(pcoa,groups) {
	groups <- as.factor(groups)
	#given pcoa and metadata
	# returns separations on axis 1, 1&2, 1&2&3
	group1.1 <- pcoa$vectors[which(groups==levels(groups)[1]),1]
	group2.1 <- pcoa$vectors[which(groups==levels(groups)[2]),1]
	diff.1 <- abs(mean(group1.1) - mean(group2.1))/sd(pcoa$vector[,1])

	group1.2 <- pcoa$vectors[which(groups==levels(groups)[1]),2]
	group2.2 <- pcoa$vectors[which(groups==levels(groups)[2]),2]
	diff.2 <- abs(mean(group1.2) - mean(group2.2))/sd(pcoa$vector[,2])

	group1.3 <- pcoa$vectors[which(groups==levels(groups)[1]),3]
	group2.3 <- pcoa$vectors[which(groups==levels(groups)[2]),3]
	diff.3 <- abs(mean(group1.3) - mean(group2.3))/sd(pcoa$vector[,3])

	diff.12 <- sqrt((diff.1^2) + (diff.2^2))
	diff.123 <- sqrt((diff.12^2) + (diff.3^2))

	returnList <- data.frame(c(diff.1,diff.12,diff.123))
	rownames(returnList) <- c("separationOn1","separationOn12","separationOn123")
	returnList <- t(returnList)
	return(returnList)
}

getAnalyzableSamples <- function(otu.tab) {
	hasMoreThanOneOTU <- function(data) {
		if(length(which(data!=0))>1) {
			return(TRUE)
		}
		return(FALSE)
	}
	indices <- c(1:length(rownames(otu.tab)))
	MoreThanOneOTU <- apply(otu.tab,1,hasMoreThanOneOTU)
	if ( (length(rownames(otu.tab))-length(which(MoreThanOneOTU))) >0) {
		indices <- which(MoreThanOneOTU)
	}


	otu.tab <- as.matrix(otu.tab)
	row.sum <- rowSums(otu.tab)

	if (length(which(apply(otu.tab,1,sum) <= 100)) > 0) {
		properReadCountIndices <- which(apply(otu.tab,1,sum) > 100)
		indices <- indicies[which(indicies %in% properReadCountIndices)]
	}

	return(indices)

}