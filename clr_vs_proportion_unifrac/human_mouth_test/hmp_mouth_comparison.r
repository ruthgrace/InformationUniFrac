options(error=recover)

library(ape)
library(phangorn)

varianceExplained <- function(pcoa) {
	# calculate total variance explained
	varExplained <- sum(apply(pcoa$vector,2,function(x) sd(x)*sd(x)))

	#calculate proportion of variance explained by each component
	propVarExplained <- apply(pcoa$vector,1,function(x) sd(x)*sd(x)/varExplained)

	return(propVarExplained)
}


# test read count correlations
plotReadCount <- function(unifrac,otu.tab,comparisonTitle,unifracType) {
	source("../metrics.r")

	#get average read count for each comparison in unifrac matrix
	avg <- averageReadCount(otu.tab)

	#put metrics matricies into single dimensional vectors for plotting
	avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])
	unifrac.vector <- unlist(unifrac[lower.tri(unifrac,diag=TRUE)])

	#plot
	plot(unifrac.vector,avg.vector,main=paste("sequencing depth and",comparisonTitle,"comparison by",unifracType),xlab="UniFrac distance",ylab="Average Total Read Count",col="palegreen",cex.lab=1.4,cex.main=2)
	abline(fit <- lm(avg.vector ~ unifrac.vector),col="darkorchid4")
	print(paste(comparisonTitle,"R squared fit:",summary(fit)$r.squared))
}

plotWeightedVsUnweighted <- function(pairPcoa) {
	weightedUnifrac.vector <- unlist(pairPcoa$unifrac$weighted[lower.tri(pairPcoa$unifrac$weighted,diag=TRUE)])
	unweightedUnifrac.vector <- unlist(pairPcoa$unifrac$unweighted[lower.tri(pairPcoa$unifrac$unweighted,diag=TRUE)])
	plot(weightedUnifrac.vector,unweightedUnifrac.vector,main=paste("weighted vs unweighted Unifrac",pairPcoa$comparisonTitle,sep="\n"))
}
plotWeightedVsInformation <- function(pairPcoa) {
	weightedUnifrac.vector <- unlist(pairPcoa$unifrac$weighted[lower.tri(pairPcoa$unifrac$weighted,diag=TRUE)])
	informationUnifrac.vector <- unlist(pairPcoa$unifrac$information[lower.tri(pairPcoa$unifrac$information,diag=TRUE)])
	plot(weightedUnifrac.vector,informationUnifrac.vector,main=paste("proportion weighted vs entropy weighted Unifrac",pairPcoa$comparisonTitle,sep="\n"))	
}


plotPcoa <- function(pairPcoa) {

	unweighted <- pairPcoa$pcoa$unweighted
	weighted <- pairPcoa$pcoa$weighted
	information <- pairPcoa$pcoa$information

	groups <- pairPcoa$groups
	comparisonTitle <- pairPcoa$comparisonTitle

	#plot pcoa plots with legend
	varEx <- varianceExplained(unweighted)
	plot(unweighted$vectors[,1],unweighted$vectors[,2], type="p",col=groups,main=paste("unweighted UniFrac PCoA for\n",comparisonTitle,sep=""),xlab=paste("First Component", round(varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("bottomright", legend=unique(groups), col=unique(groups),pch=19, title="Mouth sites")
	
	varEx <- varianceExplained(weighted)
	plot(weighted$vectors[,1],weighted$vectors[,2], type="p",col=groups,main=paste("weighted UniFrac PCoA for\n",comparisonTitle,sep=""),xlab=paste("First Component", round(varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("bottomright", legend=unique(groups), col=unique(groups), pch=19, title="Mouth sites")

	varEx <- varianceExplained(information)
	plot(information$vectors[,1],information$vectors[,2], type="p",col=groups,main=paste("information UniFrac PCoA for\n",comparisonTitle,sep=""),xlab=paste("First Component", round(varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("bottomright", legend=unique(groups), col=unique(groups), pch=19, title="Mouth sites")

}





outputFileName <- "hmp_mouth_sites_pairwise_comparison.pdf"

mouth.otu <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("buccal mucosa",20),rep("tongue dorsum",20),rep("attached keratinized gingiva",20),rep("hard palate",20),rep("saliva",20)))

# read and root tree (rooted tree is required)
mouth.tree <- read.tree("./hmp_mouth_subtree.tre")
if (!is.rooted(mouth.tree)) {
	mouth.tree <- midpoint(mouth.tree)
}
#get rid of extra quotes on OTU labels
mouth.tree$tip.label <- gsub("'","",mouth.tree$tip.label)

#get rid of extra OTUs in tree
# absent <- mouth.tree$tip.label[!(mouth.tree$tip.label %in% rownames(mouth.otu))]
# if (length(absent) != 0) {
# 		mouth.tree <- drop.tip(mouth.tree, absent)
# 		write.tree(mouth.tree,file="rep_set_v35_subtree.tre")
# }

source("../../GUniFrac.R")
source("../../InformationUniFrac.R")

#format otu table for input into unifrac methods (rownames are samples, colnames are OTUs)
mouth.otu <- t(mouth.otu)

#get rid of the 8 OTUs that aren't in the tree
mouth.otu <- mouth.otu[,which(colnames(mouth.otu) %in% mouth.tree$tip.label)]

#RAREFY
mouth.original <- mouth.otu
# rarefiedData <- Rarefy(mouth.otu,depth=2000)
# mouth.otu <- rarefiedData[[1]]
# mouth.original <- mouth.original[match(rownames(mouth.otu),rownames(mouth.original)),]


#order otus by abundance (least to most)
taxaOrder <- rev(order(apply(mouth.otu,2,sum)))
mouth.otu <- mouth.otu[,taxaOrder]

#sparsity filter
#remove all OTUs for which the minimum count is < 30
# mouth.otu.min <- apply(mouth.otu,2,min)
# mouth.otu <- mouth.otu[,which(mouth.otu.min) >= 30)]
# mouth.original <- mouth.original[,which(mouth.otu.min) >= 30)]

calculatedUnifrac <- GUniFrac(mouth.otu, mouth.tree, alpha = c(1))

unweightedUnifrac <- calculatedUnifrac$unifrac[,,2]
#weightedUnifrac
gUnifrac <- calculatedUnifrac$unifrac[,,1]
eUnifrac <- InformationUniFrac(mouth.otu, mouth.tree, alpha = c(1))$unifrac[,,1]

pairwisePcoa <- list()
index <- 1

for (i in 1:length(levels(groups))) {
	iplus1 <- i+1
	for (j in iplus1:length(levels(groups))) {
		if (i >= j || j > length(levels(groups))) {
			break
		}
		condition1 <- levels(groups)[i]
		condition2 <- levels(groups)[j]
		samplesOfInterest <- which(groups == condition1 | groups==condition2)
		pairwisePcoa[[index]] <- list()

		pairwisePcoa[[index]]$groups <- groups[samplesOfInterest]

		pairwisePcoa[[index]]$unifrac <- list()
		pairwisePcoa[[index]]$unifrac$unweighted <- unweightedUnifrac[samplesOfInterest,samplesOfInterest]
		pairwisePcoa[[index]]$unifrac$weighted <- gUnifrac[samplesOfInterest,samplesOfInterest]
		pairwisePcoa[[index]]$unifrac$information <- eUnifrac[samplesOfInterest,samplesOfInterest]

		pairwisePcoa[[index]]$comparisonTitle <- paste(condition1,"vs.",condition2)

		pairwisePcoa[[index]]$pcoa <- list()
		pairwisePcoa[[index]]$pcoa$unweighted <- pcoa(pairwisePcoa[[index]]$unifrac$unweighted)
		pairwisePcoa[[index]]$pcoa$weighted <- pcoa(pairwisePcoa[[index]]$unifrac$weighted)
		pairwisePcoa[[index]]$pcoa$information <- pcoa(pairwisePcoa[[index]]$unifrac$information)

		pairwisePcoa[[index]]$otu.tab <- mouth.otu[samplesOfInterest,]
		
		index <- index + 1
	}
}

### PRINT OUT BAR PLOT (MAKES PDF LOADING VERY SLOW)

# #get otu proportions for barplot
# mouth.prop <- t(apply(mouth.otu,1,function(x) x/sum(x)))

# #get otu total read counts
# mouth.sum <- apply(mouth.original,1,sum)

# #convert to dist structure
# unweightedUnifrac.dist <- as.dist(unweightedUnifrac)
# gUnifrac.dist <- as.dist(gUnifrac)
# eUnifrac.dist <- as.dist(eUnifrac)

# #"average" is most similar to UPGMA, apparently
# unweightedUnifrac.dendo <- hclust(unweightedUnifrac.dist, method="average")
# gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")
# eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")



# get default par
plotParameters <- par()

originalPalette <- palette()



#save to pdf
pdf("hmp_mouth_comparison_pcoa_high_read_count_no_bar_plots.pdf")


for (i in pairwisePcoa) {
	plotWeightedVsUnweighted(i)
	plotWeightedVsInformation(i)

	colors <- c("pink","red","purple","blue","orange","black")
	palette(colors)
	par(xpd=TRUE)

	plotPcoa(i)

	par(plotParameters)
	palette(originalPalette)
}









#plot dendogram with bar plots
# colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
# palette(colors)

# #GG legacy code. Fix size, margins, position
# par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1)
# plot(unweightedUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[unweightedUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(gUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[gUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(mouth.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)

#plot(clrDirichletUnifrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
#barplot(t(mouth.prop[clrDirichletUnifrac.dendo$order,]), space=0,col=colors, las=2)





# #plot pcoa first component vs. read count
# par(plotParameters)

# plot(unweightedUnifrac.pcoa$vectors[,1],mouth.sum,main="clr combination weights vs first pcoa")
# lines(lowess(unweightedUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)

# plot(gUnifrac.pcoa$vectors[,1],mouth.sum,main="gunifrac vs first pcoa")
# lines(lowess(gUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)

# plot(eUnifrac.pcoa$vectors[,1],mouth.sum,main="eunifrac vs first pcoa")
# lines(lowess(eUnifrac.pcoa$vectors[,1],mouth.sum), col="yellow") # lowess line (x,y)



# ### PLOT AGAINST SHANNON DIVERSITY METRICS

# #get shannon diversity average matrices
# diversity <- getAvgShannonDiversity(mouth.original)
# diversity.vector <- unlist(diversity[lower.tri(diversity,diag=TRUE)])

# #get shannon diversity difference matrices
# diversity.diff <- getShannonDiversityDiffMat(mouth.original)
# #put into single dimensional vector for plotting
# diversity.diff.vector <- unlist(diversity.diff[lower.tri(diversity.diff,diag=TRUE)])


# darkorchid <- col2rgb("darkorchid4")
# transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.1)
# #plot unifrac vs. shannon diversity distance matrix
# plot(unweightedUnifrac.vector,diversity.diff.vector,main="clrunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
# plot(gUnifrac.vector,diversity.diff.vector,main="gunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)
# plot(eUnifrac.vector,diversity.diff.vector,main="eunifrac vs shannon diversity difference",col=transparentdarkorchid, pch=19)

# plot(unweightedUnifrac.vector,diversity.vector,main="clrunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
# plot(gUnifrac.vector,diversity.vector,main="gunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)
# plot(eUnifrac.vector,diversity.vector,main="eunifrac vs shannon diversity",col=transparentdarkorchid, pch=19)

par(plotParameters)
palette(originalPalette)


dev.off()


#factor analysis of entropy weights


otuSum <- apply(mouth.otu,1,sum)
otuProp <- apply(mouth.otu,2,function(x) x/otuSum)
otuEntropy <- apply(otuProp,1:2,function(x) - x*log2(x))




