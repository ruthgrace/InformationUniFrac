options(error=recover)
#pcoa test script

library(ape)
library(phangorn)

# get default par
plotParameters <- par()


#this GUniFrac script was ripped straight from the GUniFrac package, with no changes.
source("../../GUniFrac.R")

source("../../InformationUniFrac.R")


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
gUnifrac <- GUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]
eUnifrac <- InformationUniFrac(brazil.otu.tab, brazil.tree, alpha = c(1))$unifrac[,,1]

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
taxonomyGroups[otuDominated] <- as.character(otuMaxTax[otuDominated])
taxonomyGroups <- as.factor(taxonomyGroups)

groups <- taxonomyGroups

# assign appropriate names to single taxa dominated groups
newLevels <- levels(taxonomyGroups)
splittaxa <- strsplit(levels(taxonomyGroups),split=";")

for (i in 1:length(splittaxa)) {
	if (length(splittaxa[[i]])>1) {
		newLevels[i] <- paste(splittaxa[[i]][length(splittaxa[[i]])-2],splittaxa[[i]][length(splittaxa[[i]])-1],splittaxa[[i]][length(splittaxa[[i]])])
	}
	else {
		newLevels[i] <- splittaxa[[i]][1]
	}
}

levels(taxonomyGroups) <- newLevels


# caculate pcoa vectors
gUnifrac.pcoa <- pcoa(gUnifrac)
eUnifrac.pcoa <- pcoa(eUnifrac)

# calculate total variance explained
gUnifrac.varExplained <- sum(apply(gUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.varExplained <- sum(apply(eUnifrac.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
gUnifrac.pc1.varEx <- sd(gUnifrac.pcoa$vector[,1])*sd(gUnifrac.pcoa$vector[,1])/gUnifrac.varExplained
#calculate proportion of variance explained by second component
gUnifrac.pc2.varEx <- sd(gUnifrac.pcoa$vector[,2])*sd(gUnifrac.pcoa$vector[,2])/gUnifrac.varExplained

eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained



# test overlap & read count correlations
source("../metrics.r")

overlap <- getOverlap(brazil.otu.tab)
avg <- averageReadCount(brazil.otu.tab)


#put metrics matricies into single dimensional vectors for plotting
overlap.vector <- unlist(overlap[lower.tri(overlap,diag=TRUE)])
avg.vector <- unlist(avg[lower.tri(avg,diag=TRUE)])

#put distance matrices into single dimensional vectors for plotting
gUnifrac.vector <- unlist(gUnifrac[lower.tri(gUnifrac,diag=TRUE)])
eUnifrac.vector <- unlist(eUnifrac[lower.tri(eUnifrac,diag=TRUE)])

#get shannon diversity average matrices
diversity <- getAvgShannonDiversity(brazil.otu.tab)
diversity.vector <- unlist(diversity[lower.tri(diversity,diag=TRUE)])

#get shannon diversity difference matrices
diversity.diff <- getShannonDiversityDiffMat(brazil.otu.tab)
#put into single dimensional vector for plotting
diversity.diff.vector <- unlist(diversity.diff[lower.tri(diversity.diff,diag=TRUE)])

diversity.max <- getMaxShannonDiversity(brazil.otu.tab)
diversity.max.vector <- unlist(diversity.max[lower.tri(diversity.max,diag=TRUE)])

diversity.min <- getMinShannonDiversity(brazil.otu.tab)
diversity.min.vector <- unlist(diversity.min[lower.tri(diversity.min,diag=TRUE)])


#convert to dist structure
gUnifrac.dist <- as.dist(gUnifrac)
eUnifrac.dist <- as.dist(eUnifrac)

#"average" is most similar to UPGMA, apparently
gUnifrac.dendo <- hclust(gUnifrac.dist, method="average")
eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")

#get otu proportions for barplot
brazil.prop <- t(apply(brazil.otu.tab,1,function(x) x/sum(x)))


#save plots as PDF
pdf("test_plots_with_brazil_study_data_no_bar_plots.pdf")

#plot dendogram with bar plots

#GG legacy code. Fix size, margins, position
#par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1,cex=0.3)

# #order the barplot 
# colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")

# plot(gUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[gUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)


#par(plotParameters)


#for this data set, the colors represent
#  chocolate4: gardnerella vaginalis
#  darkolivegreen: prevotella bivia
#  cyan: lactobacillus crispatus
#  dodgerblue: lactobacillus iners
#  navy: lactobacillus gasseri(johnsonii)
#  magenta: streptococcus (unclassified)
#samples not dominated 50% or more by a single species:
#  red: bacterial vaginosis
#  orange: intermediate
#  blue: normal/healthy
palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","aquamarine"))


#plot overlap vs gunifrac distance
plot(gUnifrac.vector,overlap.vector,main="Proportional abundance\nweighted UniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ gUnifrac.vector),col="darkorchid4")
print("weighted unifrac vs overlap")
print(summary(fit)$r.squared)
#lines(lowess(gUnifrac.vector,overlap.vector), col="darkorchid4") # lowess line (x,y)

plot(eUnifrac.vector,overlap.vector,main="Entropy weighted\nUniFrac vs. overlap",col="palegreen",xlab="UniFrac distance",ylab="Overlap",cex.lab=1.4,cex.main=2)
abline(fit <- lm(overlap.vector ~ eUnifrac.vector),col="darkorchid4")
print("entropy weighted unifrac vs overlap")
print(summary(fit)$r.squared)


palette(c("cyan","dodgerblue","red","orange","blue","black"))


#plot number of reads vs unifrac distances (checking for read count bias)
plot(gUnifrac.vector,avg.vector,main="Proportional abundance weighted\nUniFrac vs. sequencing depth",col="palegreen",xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ gUnifrac.vector),col="darkorchid4")
print("weighted unifrac vs sequencing depth")
print(summary(fit)$r.squared)

plot(eUnifrac.vector,avg.vector,main="Entropy weighted\nUniFrac vs. sequencing depth",col=rgb(.1,1,.1,0.1), pch=19,xlab="UniFrac distance",ylab="Average Total Read Count",cex.lab=1.4,cex.main=2)
#lines(lowess(gUnifrac.vector,avg.vector), col="darkorchid4") # lowess line (x,y)
abline(fit <- lm(avg.vector ~ eUnifrac.vector),col="darkorchid4")
print("entropy weighted vs sequencing depth")
print(summary(fit)$r.squared)


#plot pcoa plots with legend
plot(gUnifrac.pcoa$vectors[,1],gUnifrac.pcoa$vectors[,2], col=groups,main="Proportional abundance weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(gUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(gUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(taxonomyGroups),col=palette(),pch=19)

plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Entropy weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.3,levels(taxonomyGroups),col=palette(),pch=19)


plot(gUnifrac.vector,eUnifrac.vector,main="gunifrac vs eunifrac")


#change color palette for qiime output data (no otu count information for dominant taxa)
# red for bacterial vaginosis, orange for intermediate, blue for normal/healthy
palette(c("red","orange","blue","black"))

# read in unifrac distances from qiime
unifracWeights <- read.table("./brazil_study_data/weighted_unifrac_dm_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
#match up with metadata
indicies <- match(rownames(MyMeta),rownames(unifracWeights))
indicies <- indicies[!is.na(indicies)]
unifracWeightsFiltered <- unifracWeights[indicies,indicies]
meta_indicies <- match(colnames(unifracWeightsFiltered),rownames(MyMeta))
meta_indicies <- meta_indicies[!is.na(meta_indicies)]
#set condition (bv/i/n)
groups <- MyMeta[meta_indicies,]$n_status

#calculate pcoa
unifracWeightsFiltered.ape.pcoa <- pcoa(unifracWeightsFiltered)

#calculate proportion of variance explained by 1st and 2nd components
unifracWeightsFiltered.varExplained <- sum(apply(unifracWeightsFiltered.ape.pcoa$vector,2,function(x) sd(x)*sd(x)))
unifracWeightsFiltered.pc1.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,1])*sd(unifracWeightsFiltered.ape.pcoa$vector[,1])/unifracWeightsFiltered.varExplained
unifracWeightsFiltered.pc2.varEx <- sd(unifracWeightsFiltered.ape.pcoa$vector[,2])*sd(unifracWeightsFiltered.ape.pcoa$vector[,2])/unifracWeightsFiltered.varExplained

#plot qiime unifrac distances pcoa with legend
plot(unifracWeightsFiltered.ape.pcoa$vectors[,1],unifracWeightsFiltered.ape.pcoa$vectors[,2], col=groups,main="pcoa from qiime unifrac distances",xlab=paste("First Component", unifracWeightsFiltered.pc1.varEx,"variance explained"),ylab=paste("Second Component", unifracWeightsFiltered.pc2.varEx,"variance explained"))
legend(-0.6,0.4,levels(groups),col=palette(),pch=1)

#match qiime pcoa vectors with metadata
qiimePCOA <- read.table("./brazil_study_data/weighted_unifrac_pc_from_qiime.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
qiimePCOA <- qiimePCOA[-nrow(qiimePCOA)+1:-nrow(qiimePCOA),]
indicies <- match(rownames(qiimePCOA),rownames(MyMeta))
indicies <- indicies[!is.na(indicies)]
groups <- MyMeta[indicies,]$n_status
qiimePCOA <- qiimePCOA[match(rownames(MyMeta[indicies,]),rownames(qiimePCOA)),]

#calculate proportion of variance explained
qiimePCOA.varExplained <- sum(apply(qiimePCOA,2,function(x) sd(x)*sd(x)))
qiimePCOA.pc1.varEx <- sd(qiimePCOA[,1])*sd(qiimePCOA[,1])/qiimePCOA.varExplained
qiimePCOA.pc2.varEx <- sd(qiimePCOA[,2])*sd(qiimePCOA[,2])/qiimePCOA.varExplained

#plot qiime pecoa vectors with legend
plot(qiimePCOA[,1],qiimePCOA[,2], col=groups,main="qiime pcoa",xlab=paste("First Component", qiimePCOA.pc1.varEx,"variance explained"),ylab=paste("Second Component", qiimePCOA.pc2.varEx,"variance explained"))
legend(0.2,0.3,levels(groups),col=palette(),pch=1)


#plot pcoa first component vs. read count

plot(gUnifrac.pcoa$vectors[,1],otuSum,main="weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(gUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)

plot(eUnifrac.pcoa$vectors[,1],otuSum,main="Entropy weighted UniFrac vs. PCoA first component",col="palegreen",xlab="UniFrac distance",ylab="Average total read count")
lines(lowess(eUnifrac.pcoa$vectors[,1],otuSum), col="darkorchid4") # lowess line (x,y)



palette(c("chocolate4","darkolivegreen","cyan","dodgerblue","navy","magenta","red","orange","blue","aquamarine","darkorchid4"))



darkorchid <- col2rgb("darkorchid4")
transparentdarkorchid <- rgb(darkorchid[1]/255,darkorchid[2]/255,darkorchid[3]/255,0.1)

#plot unifrac vs. shannon diversity distance matrix
plot(gUnifrac.vector,diversity.max.vector,main="gunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.max.vector,main="eunifrac vs max shannon diversity",col=transparentdarkorchid, pch=19)

plot(gUnifrac.vector,diversity.min.vector,main="gunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)
plot(eUnifrac.vector,diversity.min.vector,main="eunifrac vs min shannon diversity",col=transparentdarkorchid, pch=19)


par(plotParameters)


#choose colors for each condition
palette(c("blue4","cornflowerblue","darkcyan","deepskyblue","cyan","dodgerblue","red","orange","blue","black"))

#create and plot iUniFrac PCoA with just intermediate + BV

gUnifrac.ibv <- gUnifrac[which(!(originalgroups=="n")),which(!(originalgroups=="n"))]
eUnifrac.ibv <- eUnifrac[which(!(originalgroups=="n")),which(!(originalgroups=="n"))]

gUnifrac.ibv.pcoa <- pcoa(gUnifrac.ibv)
eUnifrac.ibv.pcoa <- pcoa(eUnifrac.ibv)

ibvGroups <- taxonomyGroups[which(!(originalgroups=="n"))]

# calculate total variance explained
gUnifrac.ibv.varExplained <- sum(apply(gUnifrac.ibv.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.ibv.varExplained <- sum(apply(eUnifrac.ibv.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
gUnifrac.ibv1.varEx <- sd(gUnifrac.ibv.pcoa$vector[,1])*sd(gUnifrac.ibv.pcoa$vector[,1])/gUnifrac.ibv.varExplained
#calculate proportion of variance explained by second component
gUnifrac.ibv2.varEx <- sd(gUnifrac.ibv.pcoa$vector[,2])*sd(gUnifrac.ibv.pcoa$vector[,2])/gUnifrac.ibv.varExplained

eUnifrac.ibv1.varEx <- sd(eUnifrac.ibv.pcoa$vector[,1])*sd(eUnifrac.ibv.pcoa$vector[,1])/eUnifrac.ibv.varExplained
eUnifrac.ibv2.varEx <- sd(eUnifrac.ibv.pcoa$vector[,2])*sd(eUnifrac.ibv.pcoa$vector[,2])/eUnifrac.ibv.varExplained



#repeat with biplot



brazil.ibv <- brazil.otu.tab[which(!(originalgroups=="n")),]


otuSum <- apply(brazil.ibv,1,sum)
otuProp <- apply(brazil.ibv,2,function(x) x/otuSum)
otuEntropy <- apply(otuProp,1:2,function(x) - x*log2(x))

getTaxonomyLabel <- function(x) {
	taxaInfo <- strsplit(as.character(x),";")
	taxaInfo <- taxaInfo[[1]]
	return(paste(taxaInfo[length(taxaInfo)-1],taxaInfo[length(taxaInfo)]))
}

colnames(otuProp) <- lapply(taxonomy,getTaxonomyLabel)
colnames(otuEntropy) <- lapply(taxonomy,getTaxonomyLabel)


plot(gUnifrac.ibv.pcoa$vector[,1],gUnifrac.ibv.pcoa$vector[,2], col=ibvGroups,main="weighted UniFrac PCoA\nBV and Intermediate samples",xlab=paste("First Component", round(gUnifrac.ibv1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(gUnifrac.ibv2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(-0.55,0.3,levels(bvGroups),col=palette(),pch=19)

biplot(gUnifrac.ibv.pcoa,otuProp)

par(plotParameters)

plot(eUnifrac.ibv.pcoa$vector[,1],eUnifrac.ibv.pcoa$vector[,2], col=ibvGroups,main="information UniFrac PCoA\nBV and Intermediate samples",xlab=paste("First Component", round(eUnifrac.ibv1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.ibv2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

biplot(eUnifrac.ibv.pcoa,otuEntropy)



par(plotParameters)


#create and plot iUniFrac PCoA with just BV

gUnifrac.bv <- gUnifrac[which(originalgroups=="bv"),which(originalgroups=="bv")]
eUnifrac.bv <- eUnifrac[which(originalgroups=="bv"),which(originalgroups=="bv")]

gUnifrac.bv.pcoa <- pcoa(gUnifrac.bv)
eUnifrac.bv.pcoa <- pcoa(eUnifrac.bv)

bvGroups <- taxonomyGroups[which(originalgroups=="bv")]

# calculate total variance explained
gUnifrac.bv.varExplained <- sum(apply(gUnifrac.bv.pcoa$vector,2,function(x) sd(x)*sd(x)))
eUnifrac.bv.varExplained <- sum(apply(gUnifrac.bv.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by first component
gUnifrac.bv1.varEx <- sd(gUnifrac.bv.pcoa$vector[,1])*sd(gUnifrac.bv.pcoa$vector[,1])/gUnifrac.bv.varExplained
#calculate proportion of variance explained by second component
gUnifrac.bv2.varEx <- sd(gUnifrac.bv.pcoa$vector[,2])*sd(gUnifrac.bv.pcoa$vector[,2])/gUnifrac.bv.varExplained

eUnifrac.bv1.varEx <- sd(eUnifrac.bv.pcoa$vector[,1])*sd(eUnifrac.bv.pcoa$vector[,1])/eUnifrac.bv.varExplained
eUnifrac.bv2.varEx <- sd(eUnifrac.bv.pcoa$vector[,2])*sd(eUnifrac.bv.pcoa$vector[,2])/eUnifrac.bv.varExplained


#repeat with biplot


brazil.bv <- brazil.otu.tab[which(originalgroups=="bv"),]


otuSum <- apply(brazil.bv,1,sum)
otuProp <- apply(brazil.bv,2,function(x) x/otuSum)
otuEntropy <- apply(otuProp,1:2,function(x) - x*log2(x))

getTaxonomyLabel <- function(x) {
	taxaInfo <- strsplit(as.character(x),";")
	taxaInfo <- taxaInfo[[1]]
	return(paste(taxaInfo[length(taxaInfo)-1],taxaInfo[length(taxaInfo)]))
}

colnames(otuProp) <- lapply(taxonomy,getTaxonomyLabel)
colnames(otuEntropy) <- lapply(taxonomy,getTaxonomyLabel)



plot(gUnifrac.bv.pcoa$vector[,1],gUnifrac.bv.pcoa$vector[,2], col=bvGroups,main="weighted UniFrac PCoA\nBV samples",xlab=paste("First Component", round(gUnifrac.bv1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(gUnifrac.bv2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(-0.55,0.3,levels(bvGroups),col=palette(),pch=19)

biplot(gUnifrac.bv.pcoa,otuProp)

par(plotParameters)

plot(eUnifrac.bv.pcoa$vector[,1],eUnifrac.bv.pcoa$vector[,2], col=bvGroups,main="information UniFrac PCoA\nBV samples",xlab=paste("First Component", round(eUnifrac.bv1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.bv2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

biplot(eUnifrac.bv.pcoa,otuEntropy)


par(plotParameters)

dev.off()
