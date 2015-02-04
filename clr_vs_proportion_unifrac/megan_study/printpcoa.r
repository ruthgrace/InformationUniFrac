#!/usr/bin/env Rscript

options(error=recover)

library(phangorn)
plotParameters <- par()

printPCOA <- function (unweightedUnifrac,weightedUnifrac,eUnifrac,groups, condition,data.otu.tab) {

#	groups <- MyMetaOrdered$HM_GROUP
	groups <- as.factor(groups)

	otuSum <- apply(data.otu.tab,1,sum)




	unweightedUnifrac <- unweightedUnifrac[which(!is.na(groups)),which(!is.na(groups))]
	weightedUnifrac <- weightedUnifrac[which(!is.na(groups)),which(!is.na(groups))]
	eUnifrac <- eUnifrac[which(!is.na(groups)),which(!is.na(groups))]

	# caculate pcoa vectors
	unweightedUnifrac.pcoa <- pcoa(unweightedUnifrac)
	weightedUnifrac.pcoa <- pcoa(weightedUnifrac)
	eUnifrac.pcoa <- pcoa(eUnifrac)

	# calculate total variance explained
	unweightedUnifrac.varExplained <- sum(apply(unweightedUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
	weightedUnifrac.varExplained <- sum(apply(weightedUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
	eUnifrac.varExplained <- sum(apply(eUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

	# calculate proportion of variance explained by first component
	unweightedUnifrac.pc1.varEx <- sd(unweightedUnifrac.pcoa$vector[,1])*sd(unweightedUnifrac.pcoa$vector[,1])/unweightedUnifrac.varExplained
	#calculate proportion of variance explained by second component
	unweightedUnifrac.pc2.varEx <- sd(unweightedUnifrac.pcoa$vector[,2])*sd(unweightedUnifrac.pcoa$vector[,2])/unweightedUnifrac.varExplained
	unweightedUnifrac.pc3.varEx <- sd(unweightedUnifrac.pcoa$vector[,3])*sd(unweightedUnifrac.pcoa$vector[,3])/unweightedUnifrac.varExplained

	weightedUnifrac.pc1.varEx <- sd(weightedUnifrac.pcoa$vector[,1])*sd(weightedUnifrac.pcoa$vector[,1])/weightedUnifrac.varExplained
	weightedUnifrac.pc2.varEx <- sd(weightedUnifrac.pcoa$vector[,2])*sd(weightedUnifrac.pcoa$vector[,2])/weightedUnifrac.varExplained
	weightedUnifrac.pc3.varEx <- sd(weightedUnifrac.pcoa$vector[,3])*sd(weightedUnifrac.pcoa$vector[,3])/weightedUnifrac.varExplained

	eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
	eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained
	eUnifrac.pc3.varEx <- sd(eUnifrac.pcoa$vector[,3])*sd(eUnifrac.pcoa$vector[,3])/eUnifrac.varExplained

	#save plots as PDF
	pdf(paste("megan_data",condition,"pcoa_plots.pdf",sep="_"))



	#choose colors for each condition
	#palette(c("red","purple","cyan","blue","orange","black"))


	otuSum <- apply(data.otu.tab,1,sum)
	otuProp <- apply(data.otu.tab,2,function(x) x/otuSum)
	otuEntropy <- apply(otuProp,1:2,function(x) - x*log2(x))

	getTaxonomyLabel <- function(x) {
		taxaInfo <- strsplit(as.character(x),";")
		taxaInfo <- taxaInfo[[1]]
		return(paste(taxaInfo[length(taxaInfo)-1],taxaInfo[length(taxaInfo)]))
		# if (paste(taxaInfo[length(taxaInfo)-1],taxaInfo[length(taxaInfo)]) %in% interestingTaxa) {
		# 	return(taxaInfo[length(taxaInfo)-1])
		# }
		# else {
		# 	return("")
		# }
	}

	getSampleName <- function(x) {
		sampleInfo <- strsplit(as.character(x),".")
		sampleInfo <- sampleInfo[[1]]
		return(sampleInfo[length(sampleInfo)])
	}

	biplotRownames <- lapply(rownames(unweightedUnifrac), getSampleName)
	biplotTaxanames <- taxonomy

	#plot pcoa plots with legend
	plot(unweightedUnifrac.pcoa$vectors[,1],unweightedUnifrac.pcoa$vectors[,2], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(unweightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(unweightedUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend(0,0,levels(groups),col=palette(),pch=19)

	plot(unweightedUnifrac.pcoa$vectors[,2],unweightedUnifrac.pcoa$vectors[,3], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinate analysis",xlab=paste("Second Component", round(unweightedUnifrac.pc2.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(unweightedUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(unweightedUnifrac.pcoa$vectors[,1],unweightedUnifrac.pcoa$vectors[,3], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(unweightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(unweightedUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)


	plot(weightedUnifrac.pcoa$vectors[,1],weightedUnifrac.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(weightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(weightedUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(weightedUnifrac.pcoa$vectors[,2],weightedUnifrac.pcoa$vectors[,3], col=groups,main="Weighted UniFrac\nprincipal coordinate analysis",xlab=paste("Second Component", round(weightedUnifrac.pc2.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(weightedUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(weightedUnifrac.pcoa$vectors[,1],weightedUnifrac.pcoa$vectors[,3], col=groups,main="Weighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(weightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(weightedUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	#legend(-0.1,-0.12,levels(groups),col=palette(),pch=19)

	rownames(weightedUnifrac.pcoa$vectors) <- biplotRownames
	colnames(otuProp) <- biplotTaxanames
	biplot(weightedUnifrac.pcoa,otuProp[which(!is.na(groups)),])


	par(plotParameters)
	
	plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(eUnifrac.pcoa$vectors[,2],eUnifrac.pcoa$vectors[,3], col=groups,main="Information UniFrac\nprincipal coordinate analysis",xlab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(eUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,3], col=groups,main="Information UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Third Component", round(eUnifrac.pc3.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	#legend(-0.15,-0.4,levels(groups),col=palette(),pch=19)


	# #change appearance of points
	rownames(eUnifrac.pcoa$vectors) <- biplotRownames
	colnames(otuEntropy) <- biplotTaxanames
	otuEntropy[which(is.nan(otuEntropy))] <- 0
	biplot(eUnifrac.pcoa,otuEntropy[which(!is.na(groups)),])


	unweightedUnifrac.vector <- unlist(unweightedUnifrac[lower.tri(unweightedUnifrac,diag=TRUE)])
	weightedUnifrac.vector <- unlist(weightedUnifrac[lower.tri(weightedUnifrac,diag=TRUE)])
	eUnifrac.vector <- unlist(eUnifrac[lower.tri(eUnifrac,diag=TRUE)])

	par(plotParameters)
	
	plot(weightedUnifrac.vector,unweightedUnifrac.vector,main="weighted vs. unweighted unifrac")
	plot(weightedUnifrac.vector,eUnifrac.vector,main="weighted unifrac vs eunifrac")


	dev.off()

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






originalUnweightedUnifrac <- unweightedUnifrac
originalWeightedUnifrac <- weightedUnifrac
originalEUnifrac <- eUnifrac

parameterUnweightedUnifrac <- originalUnweightedUnifrac
parameterWeightedUnifrac <- originalWeightedUnifrac
parameterEUnifrac <- originalEUnifrac


groups <- MyMetaOrdered$HM_GROUP
condition <- "heavy_metal"

printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)



parameterUnweightedUnifrac <- originalUnweightedUnifrac
parameterWeightedUnifrac <- originalWeightedUnifrac
parameterEUnifrac <- originalEUnifrac


groups <- MyMetaOrdered$GROUP
condition <- "nutritional_status"

printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)



parameterUnweightedUnifrac <- originalUnweightedUnifrac
parameterWeightedUnifrac <- originalWeightedUnifrac
parameterEUnifrac <- originalEUnifrac

groups <- MyMetaOrdered$TREATMENT
condition <- "probiotic_treatment"
printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)




parameterUnweightedUnifrac <- originalUnweightedUnifrac
parameterWeightedUnifrac <- originalWeightedUnifrac
parameterEUnifrac <- originalEUnifrac


groups <- MyMetaOrdered$SAMPLE_TYPE
condition <- "body_site"
printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)




parameterUnweightedUnifrac <- originalUnweightedUnifrac
parameterWeightedUnifrac <- originalWeightedUnifrac
parameterEUnifrac <- originalEUnifrac


groups <- MyMetaOrdered$ADULT_INFANT
condition <- "adult_infant"
printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)


adultInfantGroups <- MyMetaOrdered$ADULT_INFANT
bodySiteGroups <- MyMetaOrdered$SAMPLE_TYPE

ageSiteGroups <- paste(as.character(adultInfantGroups),as.character(bodySiteGroups))

groups <- ageSiteGroups
condition <- "age_site"
printPCOA(parameterUnweightedUnifrac,parameterWeightedUnifrac,parameterEUnifrac,groups,condition,data.otu.tab)
