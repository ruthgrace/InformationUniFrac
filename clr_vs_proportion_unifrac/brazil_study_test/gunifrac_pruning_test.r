#!/usr/bin/env Rscript

options(error=recover)
#pcoa test script

library(ape)
library(phangorn)

# get default par
plotParameters <- par()


#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
source("../../GUniFrac.R")
source("../../InformationUniFrac.R")
source("../../GUniFrac_no_tree_pruning.R")
source("../../InformationUniFrac_no_tree_pruning.r")

# read OTU table and format appropriately for input into UniFrac methods
brazil.otu.tab <- read.table("./brazil_study_data/td_OTU_tag_mapped_RDPlineage_blastcorrected_vvcfilter_tempgenera.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- brazil.otu.tab$taxonomy
brazil.otu.tab <- brazil.otu.tab[-length(colnames(brazil.otu.tab))]
brazil.otu.tab <- t(as.matrix(brazil.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(brazil.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
brazil.otu.tab <- brazil.otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
brazil.tree <- read.tree("./brazil_study_data/fasttree_all_seed_OTUs.tre")
brazil.tree <- midpoint(brazil.tree)

# read metadata
MyMeta<- read.table("./brazil_study_data/metadata_BVsamplesonly.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(brazil.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
brazil.otu.tab <- brazil.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(brazil.otu.tab),rownames(MyMeta)),]

#run IUniFrac and GUniFrac for comparison, puts distance matrix in eUnifrac and gUnifrac
calculatedUnifrac_prune <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))
gUnweighted_prune <- calculatedUnifrac_prune$unifrac[,,2]
gWeighted_prune <- calculatedUnifrac_prune$unifrac[,,1]
gInformation_prune <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]

calculatedUnifrac <- GUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))
gUnweighted <- calculatedUnifrac$unifrac[,,2]
gWeighted <- calculatedUnifrac$unifrac[,,1]
gInformation <- InformationUniFracNoPrune(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
# rUnweighted <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="unweighted",verbose=TRUE)
# rWeighted <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="weighted",verbose=TRUE)
# rInformation <- getDistanceMatrix(brazil.otu.tab,brazil.tree,method="information",verbose=TRUE)

write.table(gUnweighted_prune,file="unweighted_general_unifrac_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(gWeighted_prune,file="weighted_general_unifrac_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(gInformation_prune,file="information_general_unifrac_distance_matrix.txt",sep="\t",quote=FALSE)

write.table(gUnweighted,file="unweighted_general_unifrac_no_prune_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(gWeighted,file="weighted_general_unifrac_no_prune_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(gInformation,file="information_general_unifrac_no_prune_distance_matrix.txt",sep="\t",quote=FALSE)


# write.table(rUnweighted,file="unweighted_Ruthifrac_distance_matrix.txt",sep="\t",quote=FALSE)
# write.table(rWeighted,file="weighted_Ruthifrac_distance_matrix.txt",sep="\t",quote=FALSE)
# write.table(rInformation,file="information_Ruthifrac_distance_matrix.txt",sep="\t",quote=FALSE)


#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
groups <- MyMetaOrdered$n_status #levels bv, i, n
originalgroups <- groups

# change conditions so that samples which are more than 50% one taxa are colored by that taxa
otuSum <- apply(brazil.otu.tab,1,sum)
otuMax <- apply(brazil.otu.tab,1,max)
otuWhichMax <- apply(brazil.otu.tab,1,which.max)
otuDominated <- which(otuMax > otuSum/2)


otuMaxTax <- taxonomy[otuWhichMax]
#otuDominated <- c(otuDominated[which(as.numeric(otuMaxTax[otuDominated])==32)],otuDominated[which(as.numeric(otuMaxTax[otuDominated])==33)])

taxonomyGroups <- as.character(groups)
otuDominatedNoLacto <- otuDominated[which(!(c(1:length(otuDominated) %in% grep("Lactobacillaceae",otuMaxTax[otuDominated]))))]
otuDominatedLacto <- otuDominated[grep("Lactobacillaceae",otuMaxTax[otuDominated])]
taxonomyGroups[otuDominatedNoLacto] <- as.character(otuMaxTax[otuDominatedNoLacto])

taxonomyGroups <- as.factor(taxonomyGroups)

groups <- taxonomyGroups


# assign appropriate names to single taxa dominated groups
newLevels <- levels(taxonomyGroups)
splittaxa <- strsplit(levels(taxonomyGroups),split=";")

for (i in 1:length(splittaxa)) {
	if (length(splittaxa[[i]])>1) {
		#newLevels[i] <- paste(splittaxa[[i]][length(splittaxa[[i]])-2],splittaxa[[i]][length(splittaxa[[i]])-1],splittaxa[[i]][length(splittaxa[[i]])])
		newLevels[i] <- splittaxa[[i]][length(splittaxa[[i]])-2]
	}
	else {
		newLevels[i] <- splittaxa[[i]][1]
	}
}

levels(taxonomyGroups) <- newLevels


# caculate pcoa vectors
gUnweighted.pcoa <- pcoa(gUnweighted)
gWeighted.pcoa <- pcoa(gWeighted)
gInformation.pcoa <- pcoa(gInformation)

gUnweighted_prune.pcoa <- pcoa(gUnweighted_prune)
gWeighted_prune.pcoa <- pcoa(gWeighted_prune)
gInformation_prune.pcoa <- pcoa(gInformation_prune)
# rUnweighted.pcoa <- pcoa(rUnweighted)
# rWeighted.pcoa <- pcoa(rWeighted)
# rInformation.pcoa <- pcoa(rInformation)

getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}

gUnweighted.varEx <- getVarExplained(gUnweighted.pcoa$vectors)
gWeighted.varEx <- getVarExplained(gWeighted.pcoa$vectors)
gInformation.varEx <- getVarExplained(gInformation.pcoa$vectors)

gUnweighted_prune.varEx <- getVarExplained(gUnweighted_prune.pcoa$vectors)
gWeighted_prune.varEx <- getVarExplained(gWeighted_prune.pcoa$vectors)
gInformation_prune.varEx <- getVarExplained(gInformation_prune.pcoa$vectors)
# rUnweighted.varEx <- getVarExplained(rUnweighted.pcoa$vectors)
# rWeighted.varEx <- getVarExplained(rWeighted.pcoa$vectors)
# rInformation.varEx <- getVarExplained(rInformation.pcoa$vectors)


gUnweighted.vector <- unlist(gUnweighted[lower.tri(gUnweighted,diag=TRUE)])
gWeighted.vector <- unlist(gWeighted[lower.tri(gWeighted,diag=TRUE)])
gInformation.vector <- unlist(gInformation[lower.tri(gInformation,diag=TRUE)])

gUnweighted_prune.vector <- unlist(gUnweighted_prune[lower.tri(gUnweighted_prune,diag=TRUE)])
gWeighted_prune.vector <- unlist(gWeighted_prune[lower.tri(gWeighted_prune,diag=TRUE)])
gInformation_prune.vector <- unlist(gInformation_prune[lower.tri(gInformation_prune,diag=TRUE)])
# rUnweighted.vector <- unlist(rUnweighted[lower.tri(rUnweighted,diag=TRUE)])
# rWeighted.vector <- unlist(rWeighted[lower.tri(rWeighted,diag=TRUE)])
# rInformation.vector <- unlist(rInformation[lower.tri(rInformation,diag=TRUE)])


pdf("GUnifrac_tree_pruning_test_plots.pdf")

#plot pcoa plots
plot(gUnweighted_prune.pcoa$vectors[,1],gUnweighted_prune.pcoa$vectors[,2], col=groups,main="Unweighted GUniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(gUnweighted_prune.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gUnweighted_prune.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.06,-0.11,levels(taxonomyGroups),col=palette(),pch=19,xpd=TRUE)
plot(gWeighted_prune.pcoa$vectors[,1],gWeighted_prune.pcoa$vectors[,2], col=groups,main="Weighted GUniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(gWeighted_prune.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gWeighted_prune.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(gInformation_prune.pcoa$vectors[,1],gInformation_prune.pcoa$vectors[,2], col=groups,main="Information GUniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(gInformation_prune.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gInformation_prune.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

plot(gUnweighted.pcoa$vectors[,1],gUnweighted.pcoa$vectors[,2], col=groups,main="Unweighted GUniFrac No Prune\nprincipal coordinates analysis",xlab=paste("First Component", round(gUnweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gUnweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.06,-0.1,levels(taxonomyGroups),col=palette(),pch=19,xpd=TRUE)
plot(gWeighted.pcoa$vectors[,1],gWeighted.pcoa$vectors[,2], col=groups,main="Weighted GUniFrac No Prune\nprincipal coordinates analysis",xlab=paste("First Component", round(gWeighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gWeighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(gInformation.pcoa$vectors[,1],gInformation.pcoa$vectors[,2], col=groups,main="Information GUniFrac No Prune\nprincipal coordinates analysis",xlab=paste("First Component", round(gInformation.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(gInformation.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
# plot(rUnweighted.pcoa$vectors[,1],rUnweighted.pcoa$vectors[,2], col=groups,main="Unweighted RuthiFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(rUnweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(rUnweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
# plot(rWeighted.pcoa$vectors[,1],rWeighted.pcoa$vectors[,2], col=groups,main="Weighted RuthiFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(rWeighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(rWeighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
# plot(rInformation.pcoa$vectors[,1],rInformation.pcoa$vectors[,2], col=groups,main="Information RuthiFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(rInformation.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(rInformation.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)



plot(gUnweighted.vector,gInformation.vector,main="unweighted vs. information GUnifrac")
plot(gWeighted.vector,gInformation.vector,main="weighted vs. information GUnifrac")
plot(gUnweighted.vector,gWeighted.vector,main="unweighted vs. weighted GUnifrac")

plot(gUnweighted_prune.vector,gInformation_prune.vector,main="unweighted vs. information GUnifrac No Prune")
plot(gWeighted_prune.vector,gInformation_prune.vector,main="weighted vs. information GUnifrac No Prune")
plot(gUnweighted_prune.vector,gWeighted_prune.vector,main="unweighted vs. weighted GUnifrac No Prune")

plot(gUnweighted.vector,gUnweighted_prune.vector,main="unweighted GUniFrac vs. unweighted GUnifrac No Prune")
plot(gWeighted.vector,gWeighted_prune.vector,main="weighted GUniFrac vs. weighted GUnifrac No Prune")
plot(gInformation.vector,gInformation_prune.vector,main="information GUniFrac vs. information GUnifrac No Prune")


dev.off()








