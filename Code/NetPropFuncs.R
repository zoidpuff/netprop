
library(dplyr)
library(igraph)

# Function that takes in a seed list that is a list of 3 vectors, 
#   1st vector is the seed genes, 2nd vector is the seed scores, 
#   3rd vector is a boolean vector that indicates if the gene is a seed gene or not.
# Return the seed vector that is used for network propagation, the indices in the network of the seed genes, and the scores of the seed genes

createSeedVector <- function(network, seedList, binarize = TRUE) {

    # Get the seed genes indices in the network from the seed list
    seedGenes <- seedList[[1]][which(seedList[[3]])]

    seedGeneInds <- match(seedGenes,V(network)$name)

    seedGeneScores <- seedList[[2]][which(seedList[[3]])]

    # Check if any of the seed genes are not in the network
    if(any(is.na(seedGeneInds))) {
        print(paste0(sum(is.na(seedGeneInds))," Seed gene not found in network"))
        print(seedGenes[is.na(seedGeneInds)])

        seedGeneScores <- seedGeneScores[!is.na(seedGeneInds)]
        seedGeneInds <- seedGeneInds[!is.na(seedGeneInds)]
        
    }

    # Create the seed vector 
    seedVecMaster <- rep(0, length(V(network)))

    if(binarize) {
        seedGeneScores <- rep(1,length(seedGeneScores))
    }

    seedVecMaster[seedGeneInds] <- seedGeneScores


return(list(seedVecMaster,seedGeneInds,seedGeneScores))

}


avgAUROC <- function(network, seedList, nRep, recoverSizeVec, binarize = TRUE,NormFunc = NULL,settingsForNormFunc) {

    # Create the seed vector using the seed list
    seeds <- createSeedVector(network, seedList, binarize = binarize)
    seedVecMaster <- seeds[[1]] # Vector of of length V(network) were only seed genes are non zero
    seedGeneInds <- seeds[[2]] # Vector of indices of the seed genes in the network 

    # Get the number of genes that are in the seed list but not in the network (matching with ENSEMBL ids)
    missingGenes <- length(seedList[[1]]) - length(seedGeneInds)

    # Create a vector that will hold the results
    resVec <- c("nseeds"= length(seedGeneInds), "missing" = missingGenes, "failedReplicates" = 0)
    nameTempl <- c("mean","sd","max","min")

    NAcount <- 0
    # For each recover size (proportion of seed genes that are hidden)
    for(recoverSize in recoverSizeVec) {

        # Initialize the vector that will hold the AUROC values
        procTemp <- rep(NA,nRep)

        recoverN <- max(floor(recoverSize*length(seedGeneInds)),1)

        # For each replicate
        for(i in 1:nRep){

            # Copy the master seed vector
            seedVecTemp <- seedVecMaster

            # Pick a random subset of the seed genes indices
            RecoverInds <- sample(seedGeneInds,recoverN)

            # Set the recover genes to 0 in the seed vector (remove them from the seed vector)
            seedVecTemp[RecoverInds] <- 0

            # Unormalized network propagation scores
            netPropRaw <- igraph::page_rank(network, directed = FALSE, damping = 0.85, personalized = seedVecTemp)$vector

            # If a normalization function is supplied, use it
            if(!is.null(NormFunc)) {
                netPropNorm <- NormFunc(network, seedVecTemp, netPropRaw, settingsForNormFunc)
            } else {
                netPropNorm <- netPropRaw
            }

            # Create a vector that indicates which genes are genes were omitted from the seed vector 
            recoverGenesVec <- rep(FALSE, length(V(network)))
            recoverGenesVec[RecoverInds] <- TRUE

            # Remove the seed genes that were not omitted from the seed vector, as they will always score high
            usedSeedGenes <- setdiff(seedGeneInds,RecoverInds)
            recoverGenesVec <- recoverGenesVec[-usedSeedGenes]
            netPropNorm <- netPropNorm[-usedSeedGenes]

            # Check for NA in the score vector (can happen if the normalization function fails)
            if(any(is.na(netPropNorm))) {
                print("NA in score vector")
                NAcount <- NAcount + 1
                next
            }

            # Compute AUROC with the auroc package (should maybe specify the direction of the ROC curve)
            rocTmp <- pROC::roc(response = recoverGenesVec, predictor = netPropNorm, quiet = TRUE)
            procTemp[i] <- pROC::auc(rocTmp)

        }
        
        # Compute the mean, sd, max, and min of the AUROC values for the replicates
        resVecTemp <- c(mean(procTemp),sd(procTemp),max(procTemp),min(procTemp))

        # Add the results to the result vector
        names(resVecTemp) <- paste0(nameTempl,"_",recoverSize)
        resVec <- c(resVec, resVecTemp)
    }

    resVec["failedReplicates"] <- NAcount

    return(resVec) # names vector with mean, sd, max, min for each recover size

}


ECnormalize <- function(network,seedVector,rawNetPropScores,settings) {
        # Like they do here https://www.frontiersin.org/articles/10.3389/fgene.2019.00004/full#h5
        convergeSuccess <- "notTrue"
        counter <- 0
        while(!is.logical(convergeSuccess) & counter < 100) {
            ecScore <- igraph::page_rank(network, directed = FALSE, damping = 1, personalized = seedVector,algo="arpack")
            convergeSuccess <- all.equal(ecScore$value,1)
            counter <- counter + 1
        }

        if(convergeSuccess != TRUE) {
            print("Failed to converge")
            return(rep(NA,length(rawNetPropScores)))
        } 

        # Choose an epsilon value
        epsilon <- 1e-6
        
        # Add epsilon to avoid log(0) and apply logarithmic transformation for numerical stability
        logTransformedScores <- log(rawNetPropScores + epsilon)
        logTransformedECScores <- log(ecScore$vector + epsilon)

        return(logTransformedScores - logTransformedECScores)
    
}



# Given a seet list of genes, run net prop several times on different with different random seed vectors
# settings is a list with the following elements: nSamples, perserveDegree, degreeSampleSmoothing, minBucketSize

permuteTestNormalize <- function(network, seedVector, netPropTRUE, settings) {
    nSamples <- settings[["nSamples"]]
    perserveDegree <- settings[["perserveDegree"]]
    degreeSampleSmoothing <- settings[["degreeSampleSmoothing"]]
    minBucketSize <- settings[["minBucketSize"]]

    ### Internal Helper functions ###

        # Similar to sample, but if input is numeric of length one, it does returns that number
        sampleIndsFix <- function(x) {
            return(x[sample(length(x),1)])
        }

        # Function which finds the index of the minimum value in a vector and if there are multiple minimum values,
        #   it randomly selects one of them instead of always selecting the first one like the which.min function does
        whichMinRand <- function(deltas) {
            inds <- which(deltas == min(deltas))
            return(inds[sample(length(inds),1)])
        }

        # Pick a random index from the n smallest values in a vector
        pickSimilar <- function(deltas,n) {

            deltasOrd <- deltas[order(deltas)]

            probVec <- rep(0,length(deltasOrd))

            probVec[1] <- 1

            for(i in 2:length(deltasOrd)) {
                if(deltasOrd[i] == deltasOrd[i-1]) {
                    probVec[i] <- probVec[i-1]
                } else {
                    probVec[i] <- probVec[i-1] + 1
                }
            }

            return(sample(order(deltas),1,prob = (1/probVec)^n))
            #return(sample(order(deltas),1,prob = ((1/(log(deltas+1)+1))^n)))
            #return(sample(order(deltas)[1:n],1))
        }


    # Init vectors that will hold the number of times a node was sampled and the number of times it was included in the sampling
    countVec <- rep(1, length(netPropTRUE))
    includedVec <- rep(1, length(netPropTRUE))

    # If we want to permute the seed vector while maintaing the degree distribution of the seed genes
    if(perserveDegree) {

        # Keep track of bucketing attempts
        retryCount <- 0	
        
        seedInds <- which(seedVector != 0)

        # Get the degrees of the seed genes
        seedDegrees <- degree(network, seedInds)
        seedScores <- seedVector[seedInds]
        
        # Get degrees of all nodes in network
        allDegrees <- degree(network)
        
        bucketSizes <- c(-1,0,0)

        nodesInBucketsList <- list(1)

        # While any bucket is smaller then minBucketSize or the number of buckets is not equal to the number of seed genes
        while(any(bucketSizes < minBucketSize) | (length(nodesInBucketsList) != length(seedDegrees))) {
       
            if(degreeSampleSmoothing > 0) {
                # Find similar bucket
                nodesInBuckets <- sapply(allDegrees, function(x) pickSimilar(abs(x - seedDegrees),degreeSampleSmoothing))
            } else {
                # Find the bucket that each node belongs to and store in a vector
                nodesInBuckets <- sapply(allDegrees, function(x) whichMinRand(abs(x - seedDegrees)))
            }
        
            # Create a list of vectors that contain the node indices, sorted into n buckets based on the degrees of the n seed nodes
            nodesInBucketsList <- split(1:length(nodesInBuckets), nodesInBuckets)
            bucketSizes <- sapply(nodesInBucketsList,length)

            if(any(bucketSizes < minBucketSize) | (length(nodesInBucketsList) != length(seedDegrees))) {
                print("Atleast one bucket is smaller then minBucketSize, trying again")
                retryCount <- retryCount + 1
            } 

            if(retryCount > 20) {
                print("Failed to create buckets with minBucketSize, exiting")
                # Quit the function
                return(rep(NA,length(netPropTRUE)))
            }
        }

    
         #print(paste0("Standard deviation of bucket size: ", sd(bucketSizes)))
         #print(summary(bucketSizes))

        #diffVec <- rep(1,length(nodesInBucketsList))
        #for(i in 1:length(nodesInBucketsList)) {
            # Print the degree of the seed node corresponding to the bucket
           # print(paste0("Bucket ", i,   " nNodes: ", length(nodesInBucketsList[[i]]),
            #                            " seed degree: ", seedDegrees[i],
             #                           " mean degree: ", signif(mean(allDegrees[nodesInBucketsList[[i]]])),
              #                          " SD: ", signif(sd(allDegrees[nodesInBucketsList[[i]]]))))

         #                               diffVec[i] <- abs(seedDegrees[i] - mean(allDegrees[nodesInBucketsList[[i]]]))

        
        #}
        #print(paste0("SD of difference between seed degree and mean degree: ", sd(diffVec)))
        #print(summary(diffVec))

    }

    
    # For each permutation
    for(i in 1:nSamples) {

        if(perserveDegree) {
            # Sampling procedure that samples nodes with the same degree as the seed genes 
            sampledNodesInds <- sapply(nodesInBucketsList, function(x) sampleIndsFix(x))

            #sampledNodesInds <- match(sampledNodes,V(network)$name)

            permutedSeedVector <- rep(0, length(V(network)$name))

            # Give sample nodes scores based on original seed vector
            permutedSeedVector[sampledNodesInds] <- seedScores

        } else {
            # Create a permuted seed vector with sample()
            permutedSeedVector <- sample(seedVector)
        }

        # Run network propagation with the permuted seed vector
        netPropFALSE <- igraph::page_rank(network, directed = FALSE, damping = 0.85, personalized = permutedSeedVector)

        # Add to the count vector if the permuted seed vector has a higher score then the true seed vector and the node was not a seed gene
        countVec <- countVec + as.numeric((netPropFALSE$vector >= netPropTRUE) & (permutedSeedVector == 0) ) 

        # Add to the included vector if the node was not a seed gene
        includedVec <- includedVec + as.numeric(permutedSeedVector == 0)

    }

    if(any(includedVec == 1)) {
        print("Atleast one node was always sampled, likely because it was alone in a bucket (no other node with similar degree)")
    }

    permutationScores <- countVec/includedVec

    #return(list(permutationScores, countVec, includedVec))
    return(permutationScores)
}

