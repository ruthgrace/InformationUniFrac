
#manually calculating minimum case triangle inequality violation distance matrix

load("minimum_case_for_weighted_triangle_violation.dat")

tree <- minCase$FivePercentSparsityTree
samples <- minCase$FivePercentSparsitySamples

source("../../GUniFrac.R")
source("../../RuthiFrac.r")

library(phangorn)

tree <- reorder(tree,order="postorder")

tipIndices <- c(1:length(tree$tip.label))

weights <- data.frame(matrix(ncol=(length(tree$node.label) + length(tree$tip.label)),nrow=length(rownames(samples))))

