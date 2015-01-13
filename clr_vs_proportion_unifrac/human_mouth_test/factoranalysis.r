#!/usr/bin/env Rscript

library(phangorn)

mouth.otu <- read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)

groups <- as.factor(c(rep("buccal mucosa",20),rep("tongue dorsum",20),rep("attached keratinized gingiva",20),rep("hard palate",20),rep("saliva",20)))

# read and root tree (rooted tree is required)
mouth.tree <- read.tree("./hmp_mouth_subtree.tre")
if (!is.rooted(mouth.tree)) {
	mouth.tree <- midpoint(mouth.tree)
}
#get rid of extra quotes on OTU labels
mouth.tree$tip.label <- gsub("'","",mouth.tree$tip.label)

mouth.otu <- t(mouth.otu)

#get rid of the 8 OTUs that aren't in the tree
mouth.otu <- mouth.otu[,which(colnames(mouth.otu) %in% mouth.tree$tip.label)]

#RAREFY
mouth.original <- mouth.otu

taxaOrder <- rev(order(apply(mouth.otu,2,sum)))
mouth.otu <- mouth.otu[,taxaOrder]

fa <- factanal(mouth.otu,factors=1)