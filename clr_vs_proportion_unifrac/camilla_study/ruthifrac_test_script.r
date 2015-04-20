options(error=recover)
#pcoa test script

library(ape)
library(phangorn)

# commenting out incorrect clr dirichlet

# get default par
plotParameters <- par()


#this weightedUnifrac script was ripped straight from the weightedUnifrac package, with no changes.
source("../../Ruthifrac.R")


# read OTU table and format appropriately for input into UniFrac methods
breastmilk.otu.tab <- read.table("./camilla_data/td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- breastmilk.otu.tab$taxonomy
breastmilk.otu.tab <- breastmilk.otu.tab[-length(colnames(breastmilk.otu.tab))]
breastmilk.otu.tab <- t(as.matrix(breastmilk.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(breastmilk.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
breastmilk.otu.tab <- breastmilk.otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
breastmilk.tree <- read.tree("./camilla_data/fasttree_all_seed_OTUs.tre")
breastmilk.tree <- midpoint(breastmilk.tree)

# read metadata
MyMeta<- read.table("./camilla_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove infected sample S38I
#MyMeta <- MyMeta[(which(rownames(MyMeta)!="S38I")),]

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(breastmilk.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
breastmilk.otu.tab <- breastmilk.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(breastmilk.otu.tab),rownames(MyMeta)),]


unweightedUnifrac <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="unweighted",verbose=TRUE,pruneTree=TRUE,normalize=TRUE)
weightedUnifrac <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="weighted",verbose=TRUE,pruneTree=TRUE,normalize=TRUE)
eUnifrac <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="information",verbose=TRUE,pruneTree=TRUE,normalize=TRUE)

#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
# groups <- MyMetaOrdered$Gestation #levels ""    "G"   "P_a" "P_b" "T" 
# originalgroups <- groups
# levels(groups) <- c(levels(groups),"Infected")
# groups[which(rownames(MyMetaOrdered)=="S38I")] <- levels(groups)[6]

groups <- rep("Not Infected",length(MyMetaOrdered$Gestation))
groups[which(rownames(MyMetaOrdered)=="S38I")] <- "Infected"
groups <- as.factor(groups)

otuSum <- apply(breastmilk.otu.tab,1,sum)


# # change conditions so that samples which are more than 50% one taxa are colored by that taxa
# otuMax <- apply(breastmilk.otu.tab,1,max)
# otuWhichMax <- apply(breastmilk.otu.tab,1,which.max)
# otuDominated <- which(otuMax > otuSum/2)


# otuMaxTax <- taxonomy[otuWhichMax]
# otuDominated <- c(otuDominated[which(as.numeric(otuMaxTax[otuDominated])==32)],otuDominated[which(as.numeric(otuMaxTax[otuDominated])==33)])

# taxonomyGroups <- as.character(groups)
# taxonomyGroups[otuDominated] <- as.character(otuMaxTax[otuDominated])
# taxonomyGroups <- as.factor(taxonomyGroups)

# groups <- taxonomyGroups

# # assign appropriate names to single taxa dominated groups
# newLevels <- levels(taxonomyGroups)
# splittaxa <- strsplit(levels(taxonomyGroups),split=";")

# for (i in 1:length(splittaxa)) {
# 	if (length(splittaxa[[i]])>1) {
# 		newLevels[i] <- paste("L.",splittaxa[[i]][length(splittaxa[[i]])])
# 	}
# 	else {
# 		newLevels[i] <- splittaxa[[i]][1]
# 	}
# }

# levels(taxonomyGroups) <- newLevels


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

weightedUnifrac.pc1.varEx <- sd(weightedUnifrac.pcoa$vector[,1])*sd(weightedUnifrac.pcoa$vector[,1])/weightedUnifrac.varExplained
weightedUnifrac.pc2.varEx <- sd(weightedUnifrac.pcoa$vector[,2])*sd(weightedUnifrac.pcoa$vector[,2])/weightedUnifrac.varExplained

eUnifrac.pc1.varEx <- sd(eUnifrac.pcoa$vector[,1])*sd(eUnifrac.pcoa$vector[,1])/eUnifrac.varExplained
eUnifrac.pc2.varEx <- sd(eUnifrac.pcoa$vector[,2])*sd(eUnifrac.pcoa$vector[,2])/eUnifrac.varExplained

#save plots as PDF
pdf("camilla_breastmilk_ruthifrac_plots.pdf")


# MAKE BAR PLOTS

# #convert to dist structure
# unweightedUnifrac.dist <- as.dist(unweightedUnifrac)
# weightedUnifrac.dist <- as.dist(weightedUnifrac)
# eUnifrac.dist <- as.dist(eUnifrac)

# #"average" is most similar to UPGMA, apparently
# unweightedUnifrac.dendo <- hclust(unweightedUnifrac.dist, method="average")
# weightedUnifrac.dendo <- hclust(weightedUnifrac.dist, method="average")
# eUnifrac.dendo <- hclust(eUnifrac.dist, method="average")

# #get otu proportions for barplot
# brazil.prop <- t(apply(breastmilk.otu.tab,1,function(x) x/sum(x)))

#plot dendogram with bar plots

# #GG legacy code. Fix size, margins, position
# par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1,cex=0.3)

# plot(unweightedUnifrac.dendo, axes=F, ylab=NULL, ann=F,hang=-1)

# #order the barplot 
# colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
# barplot(t(brazil.prop[unweightedUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(weightedUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[weightedUnifrac.dendo$order,]), space=0,col=colors, las=2)

# plot(eUnifrac.dendo, axes=F, ylab=NULL, ann=F)
# #order the barplot 
# barplot(t(brazil.prop[eUnifrac.dendo$order,]), space=0,col=colors, las=2)

#plot(clrDirichletUniFrac.dendo, axes=F, ylab=NULL, ann=F)
#order the barplot 
#barplot(t(brazil.prop[clrDirichletUniFrac.dendo$order,]), space=0,col=colors, las=2)

#par(plotParameters)



#choose colors for each condition
palette(c("red","black","cyan","dodgerblue","blue","orange"))



#plot pcoa plots with legend
plot(unweightedUnifrac.pcoa$vectors[,1],unweightedUnifrac.pcoa$vectors[,2], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(unweightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(unweightedUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#placement with S38I included
#legend(-0.1,-0.055,levels(groups),col=palette(),pch=19)
# #placement with S38I excluded
# legend(0.055,0.15,levels(groups),col=palette(),pch=19)

plot(weightedUnifrac.pcoa$vectors[,1],weightedUnifrac.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(weightedUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(weightedUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#legend(-0.1,-0.12,levels(groups),col=palette(),pch=19)

plot(eUnifrac.pcoa$vectors[,1],eUnifrac.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinate analysis",xlab=paste("First Component", round(eUnifrac.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Component", round(eUnifrac.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#placement with S38I included
#legend(-0.15,-0.4,levels(groups),col=palette(),pch=19)
# #placement with S38I excluded
# legend(0.4,-0.15,levels(groups),col=palette(),pch=19)


# plot(weightedUnifrac.vector,unweightedUnifrac.vector,main="weighted vs. unweighted unifrac")
# plot(weightedUnifrac.vector,eUnifrac.vector,main="weighted unifrac vs eunifrac")


otuSum <- apply(breastmilk.otu.tab,1,sum)
otuProp <- apply(breastmilk.otu.tab,2,function(x) x/otuSum)
otuEntropy <- apply(otuProp,1:2,function(x) - x*log2(x))

interestingTaxa <- c("Pseudomonas |92","Staphylococcus |77","Pasteurella |93","Escherichia-Shigella |98")

getTaxonomyLabel <- function(x) {
	return("")
	# taxaInfo <- strsplit(as.character(x),";")
	# taxaInfo <- taxaInfo[[1]]
	# if (paste(taxaInfo[length(taxaInfo)-1],taxaInfo[length(taxaInfo)]) %in% interestingTaxa) {
	# 	return(taxaInfo[length(taxaInfo)-1])
	# }
	# else {
	# 	return("")
	# }
}

#the ones i am interested in are:
#Pseudomonas |92
#Staphylococcus |77
#Pasteurella |93
#Escherichia−Shigella |98

#change appearance of points
rownames(eUnifrac.pcoa$vectors) <- rep("o",length(rownames(eUnifrac.pcoa$vectors)))
rownames(eUnifrac.pcoa$vectors)[44] <- "i"

colnames(otuEntropy) <- lapply(taxonomy,function(x) { return ("")})
biplot(eUnifrac.pcoa,otuEntropy)


dev.off()

