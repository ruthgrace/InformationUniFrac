# ruth's unifrac method

#in short, this method calculates the sum of (branch length of node to parent node * |proportional abundance of OTUs in node subtree for sample A - proportional abundance of OTUs in node subtree for sample B| )


library(phangorn)
library(ape)

#valid methods are unweighted, weighted, information. Any other method will result in a warning and the unweighted analysis

getDistanceMatrix <- function(otuTable,tree,method="weighted",verbose=FALSE)  {

	# make sure tree is in correct format (rooted, postorder)

	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
		if(verbose) { print("Rooting tree by midpoint") }
	}

	if (attributes(tree)$order!="postorder") {
		tree <- reorder(tree,order="postorder")
		if (verbose) { print("Reordering tree as postorder for distance calculation algorithm") }
	}

	# get proportions
	readsPerSample <- apply(otuTable,1,sum)
	otu.prop <- otuTable/readsPerSample
	rownames(otu.prop) <- rownames(otuTable)
	colnames(otu.prop) <- colnames(otuTable)

	if(verbose) {	print("calculated proportional abundance")	}

	##get cumulative proportional abundance for the nodes (nodes are ordered same as in the phylo tree representation)

	#cumulative proportional abundance stored in weights
	#weights <- data.frame(matrix(0,ncol=(length(tree$edge.length) + 1),nrow=nrow(otuTable)))
	weights <- data.frame(matrix(NA,ncol=(length(tree$edge.length) + 1),nrow=nrow(otuTable)))
	#each row is a sample
	rownames(weights) <- rownames(otuTable)
	#each column is the abundance weighting for a node in the phylogenetic tree
	colnames(weights) <- c(1:(length(tree$edge.length)+1))

	treeLeaves <- c(1:length(tree$tip.label))

	#loop through edges
	#if child node of edge (of which there is only one) has not been seen before, it is a leaf
	#	(this is a property of a postorder tree -- children are listed before parents)
	if(verbose) {	print("calculating weights...") }
	for (i in c(1:nrow(tree$edge))) {
		parentNode <- tree$edge[i,1]
		childNode <- tree$edge[i,2]

		# if (childNode %in% treeLeaves) {
		# 	otuName <- tree$tip.label[childNode]
		# 	otuIndex <- which(colnames(otu.prop) == otuName)[1]
		# 	weights[,childNode] <- otu.prop[,otuIndex]
		# }

		if (length(which(is.na(weights[,childNode]))) > 0 ) { #if node is all NA, node has not been seen before, and node is a leaf (ie. an OTU)
			#put OTU abundance in weights
			otuName <- tree$tip.label[childNode]
			otuIndex <- which(colnames(otu.prop) == otuName)[1]
			weights[,childNode] <- otu.prop[,otuIndex]
		}

		if (length(which(is.na(weights[,parentNode]))) > 0 ) {
			# initialize parentNode with counts of zero
			weights[,parentNode] <- 0
		}


		# print("new child weights")
		# print(str(weights[,childNode]))
		# print("old parent weights")
		# print(str(weights[,parentNode]))

		#add child node abundance to parent node abundance
		weights[,parentNode] <- weights[,parentNode] + weights[,childNode]
		# print("new parent weights")
		# print(str(weights[,parentNode]))
		
	}


	if(verbose) {	print("done calculating weights")	}

	if (method=="information") {
		if(verbose) {	print("information entropy transform")	}
		#information entropy transform
		weights[] <- (-1) * weights[] * log2(weights[])
		weights <- as.matrix(weights)
		weights[which(is.na(weights))] <- 0
		weights <- data.frame(weights,check.names=FALSE)
	}


	nSamples <- length(rownames(otuTable))
	distanceMatrix <- data.frame(matrix(ncol=nSamples,nrow=nSamples))
	rownames(distanceMatrix) <- rownames(otuTable)
	colnames(distanceMatrix) <- rownames(otuTable)
	branchLengths <- tree$edge.length

	branchLengths <- branchLengths[order(tree$edge[,2])]
	weights <- weights[,which(!is.na(match(colnames(weights),tree$edge[,2])))]


	if(verbose) {	print("calculating pairwise distances...")	}

	for (i in 1:nSamples) {
		for (j in i:nSamples) {

			#remove branch lengths that aren't in either sample
			comparisonBranchLengths <- branchLengths
			excludeBranchLengths <- which( (weights[i,] <= 0) & (weights[j,] <= 0) )
			comparisonBranchLengths[excludeBranchLengths] <- 0

			#turn branch lengths into proportions
			totalBranchLength <- sum(comparisonBranchLengths)
			comparisonBranchLengths <- comparisonBranchLengths/totalBranchLength

			if (method=="weighted" || method=="information") {
				# the formula is sum of (branch lengths * | proportional abundance for sample A - proportional abundance for sample B| )
				distance <- sum( comparisonBranchLengths * abs(weights[i,] - weights[j,]) )
			}
			else {
				if (method!="unweighted") {
					warning(paste("Invalid method",method,", using unweighted Unifrac"))
				}
				# the formula is sum of (branch lengths * (1 if one sample has counts and not the other, 0 otherwise) )
				#	i call the (1 if one sample has counts and not the other, 0 otherwise) xorBranchLength
				xorBranchLength <- as.numeric(xor( weights[i,] > 0, weights[j,] > 0))
				distance <- sum( comparisonBranchLengths *  xorBranchLength)
			}
			distanceMatrix[i,j] <- distance
			distanceMatrix[j,i] <- distance

		}
	}

	#distanceMatrix <- as.matrix(distanceMatrix)

	if(verbose) {	print("done")	}

	return(distanceMatrix)

}