# manual mincase test

sampleSums <- apply(minCase$FivePercentSparsitySamples,1,sum)

#set up data structure with branch lenghts and proportional abundances of all samples
weights <- data.frame(matrix(ncol=10,nrow=3))
colnames(weights) <- minCase$FivePercentSparsityTree$edge[,2]
rownames(weights) <- c("branchLengths",rownames(minCase$FivePercentSparsitySamples))

weights[1,] <- minCase$FivePercentSparsityTree$edge.length
weights <- weights[,order(as.numeric(colnames(weights)))]

#assign raw abundance to tip nodes
weights[2:4,1] <- minCase$FivePercentSparsitySamples[,1]
weights[2:4,2] <- minCase$FivePercentSparsitySamples[,3]
weights[2:4,3] <- minCase$FivePercentSparsitySamples[,6]
weights[2:4,4] <- minCase$FivePercentSparsitySamples[,2]
weights[2:4,5] <- minCase$FivePercentSparsitySamples[,4]
weights[2:4,6] <- minCase$FivePercentSparsitySamples[,5]

#add up raw abundance for interior nodes
weights[2:4,which(colnames(weights)==11)[1]] <- weights[2:4,which(colnames(weights)==1)[1]] + weights[2:4,which(colnames(weights)==2)[1]]
weights[2:4,which(colnames(weights)==10)[1]] <- weights[2:4,which(colnames(weights)==11)[1]] + weights[2:4,which(colnames(weights)==3)[1]]
weights[2:4,which(colnames(weights)==9)[1]] <- weights[2:4,which(colnames(weights)==10)[1]] + weights[2:4,which(colnames(weights)==6)[1]]
weights[2:4,which(colnames(weights)==8)[1]] <- weights[2:4,which(colnames(weights)==4)[1]] + weights[2:4,which(colnames(weights)==5)[1]]

#sanity check -- should equal to sampleSums
weights[2:4,which(colnames(weights)==8)[1]] + weights[2:4,which(colnames(weights)==9)[1]]

#convert to proportional abundance
weights[2,] <- weights[2,]/sampleSums[1]
weights[3,] <- weights[3,]/sampleSums[2]
weights[4,] <- weights[4,]/sampleSums[3]

d12 <- sum(weights[1,]*(abs(weights[2,]-weights[3,])))
d23 <- sum(weights[1,]*(abs(weights[3,]-weights[4,])))
d13 <- sum(weights[1,]*(abs(weights[2,]-weights[4,])))

