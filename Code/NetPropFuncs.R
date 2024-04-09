
library(dplyr)
library(igraph)
library(ggplot2)
#library(doParallel)
#library(foreach)

# Function that takes in a seed list that is a list of 3 vectors, 
#   1st vector is the seed genes, 2nd vector is the seed scores, 
#   3rd vector is a boolean vector that indicates if the gene is a seed gene or not.
# Return the seed vector that is used for network propagation, the indices in the network of the seed genes, and the scores of the seed genes

createSeedVector <- function(network, seedList, binarize = TRUE) {

    # Check seed genes are all unique
    if(length(unique(seedList[[1]])) != length(seedList[[1]])) {
        print("Seed genes are not unique")
        return(NULL)
    }

    # Get the seed genes indices in the network from the seed list
    seedGenes <- seedList[[1]][which(seedList[[3]])]

    seedGeneInds <- match(seedGenes,igraph::V(network)$name)

    seedGeneScores <- seedList[[2]][which(seedList[[3]])]

    # Check if any of the seed genes are not in the network
    if(any(is.na(seedGeneInds))) {
        print(paste0(sum(is.na(seedGeneInds))," Seed gene not found in network"))
        print(seedGenes[is.na(seedGeneInds)])

        seedGeneScores <- seedGeneScores[!is.na(seedGeneInds)]
        seedGeneInds <- seedGeneInds[!is.na(seedGeneInds)]
        
    }

    # Create the seed vector 
    seedVecMaster <- rep(0, length(igraph::V(network)))

    if(binarize) {
        seedGeneScores <- rep(1,length(seedGeneScores))
    }

    seedVecMaster[seedGeneInds] <- seedGeneScores


return(list(seedVecMaster,seedGeneInds,seedGeneScores))

}


avgAUROC <- function(network, seedList, nRep, recoverSizeVec, binarize = TRUE,NormFunc = NULL,settingsForNormFunc) {

    #no_cores <- min(detectCores(),16)

    #cl <- makeCluster(no_cores)
    #registerDoParallel(cl)

    # Create the seed vector using the seed list
    seeds <- createSeedVector(network, seedList, binarize = binarize)
    seedVecMaster <- seeds[[1]] # Vector of of length V(network) were only seed genes are non zero
    seedGeneInds <- seeds[[2]] # Vecor of indices of the seed genes in the network 

    # Get the number of genes that are in the seed list but not in the network (matching with ENSEMBL ids)
    missingGenes <- length(seedList[[1]]) - length(seedGeneInds)

    if(length(seedGeneInds) < 3) {
        print("2 or less seed genes found in network, returning NA")
        return(rep(NA,(3+4*length(recoverSizeVec))))
    }

    # Create a vector that will hold the results
    resVec <- c("nseeds"= length(seedGeneInds), "missing" = missingGenes, "failedReplicates" = 0)
    nameTempl <- c("mean","sd","max","min")

    NAcount <- 0

    # For each recover size (proportion of seed genes that are hidden)
    for(recoverSize in recoverSizeVec) {

        # Initialize the vector that will hold the AUROC values
        procTemp <- rep(NA,nRep)

        recoverN <- max(floor(recoverSize*length(seedGeneInds)),1)
        
        #procTemp <- foreach(i = 1:nRep, .combine = 'c') %dopar% {
        for(i in 1:nRep) {

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

            if (sum(is.na(netPropNorm)) > length(netPropNorm)/2) {
                warning("More than half of the network propagation scores are NA, skipping replicate")
                next
            }

            # Create a vector that indicates which genes are genes were omitted from the seed vector 
            recoverGenesVec <- rep(FALSE, length(igraph::V(network)))
            recoverGenesVec[RecoverInds] <- TRUE

            # Remove the seed genes that were not omitted from the seed vector, as they will always score high
            usedSeedGenes <- setdiff(seedGeneInds,RecoverInds)
            recoverGenesVec <- recoverGenesVec[-usedSeedGenes]
            netPropNorm <- netPropNorm[-usedSeedGenes]

            # Compute AUROC with the auroc package (should maybe specify the direction of the ROC curve)
            procTemp[i] <- pROC::roc(response = recoverGenesVec, predictor = netPropNorm, quiet = TRUE)$auc

   
        }
        

        NAcount <- sum(is.na(procTemp)) + NAcount

        if(sum(is.na(procTemp)) == nRep) {
            resVecTemp <- rep(NA,4)
         } else {
            # Compute the mean, sd, max, and min of the AUROC values for the replicates
            resVecTemp <- c(mean(procTemp),sd(procTemp),max(procTemp),min(procTemp))
         }


        # Add the results to the result vector
        names(resVecTemp) <- paste0(nameTempl,"_",recoverSize)
        resVec <- c(resVec, resVecTemp)
    }

    resVec["failedReplicates"] <- NAcount

    # Stop the cluster
    #stopCluster(cl)

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

        if(settings$logtransform) {

            # Choose an epsilon value
            epsilon <- 1e-9
            
            # Add epsilon to avoid log(0) and apply logarithmic transformation for numerical stability
            logTransformedScores <- log(rawNetPropScores + epsilon)
            logTransformedECScores <- log(ecScore$vector + epsilon)

            return(logTransformedScores - logTransformedECScores)

        } else {
            return(rawNetPropScores/ecScore$vector)
        }
    
}

scaleNormalize <- function(network,seedVector,rawNetPropScores,settings) {
    return(as.vector(scale(rawNetPropScores)))
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
    seedInds <- which(seedVector != 0)
    seedScores <- seedVector[seedInds]

    # If we want to permute the seed vector while maintaing the degree distribution of the seed genes
    if(perserveDegree) {

        # Keep track of bucketing attempts
        retryCount <- 0	
        
        # Get the degrees of the seed genes
        seedDegrees <- igraph::degree(network, seedInds)
        
        # Get degrees of all nodes in network
        allDegrees <- igraph::degree(network)
        
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
                warning("Failed to create buckets with minBucketSize, exiting")
                # Print the number of seeds the number of buckets and the sizes of the buckets
                errString <- paste0("Number of seeds: ", length(seedDegrees), " Number of buckets: ", length(nodesInBucketsList), " Sizes of buckets: ", paste0(bucketSizes,collapse = ", "))
                warning(errString)
                return(rep(NA,length(netPropTRUE)))
            }
            # Stop 
            #stopifnot(!any(is.na(bucketSizes)))

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

    numberOfNodes <- length(seedVector)
    numberOfSeeds <- length(seedInds)

    # For each permutation
    for(i in 1:nSamples) {

        permutedSeedVector <- rep(0, numberOfNodes)

        if(perserveDegree) {
            # Sampling procedure that samples nodes with the same degree as the seed genes 
            sampledNodesInds <- sapply(nodesInBucketsList, function(x) sampleIndsFix(x))

            # Give sample nodes scores based on original seed vector
            permutedSeedVector[sampledNodesInds] <- seedScores

        } else {
            # Create a permuted seed vector with sample()
            # I changed this so that I save inds so I dont have to do x != 0 many times 
            sampledNodesInds <- sample(numberOfNodes, numberOfSeeds)
            permutedSeedVector[sampledNodesInds] <- seedScores
        }

        # Run network propagation with the permuted seed vector
        netPropFALSE <- igraph::page_rank(network, directed = FALSE, damping = 0.85, personalized = permutedSeedVector)

        # Add to the count vector if the permuted seed vector has a higher score then the true seed vector and the node was not a seed gene
        temp <- as.numeric(netPropFALSE$vector >= netPropTRUE)
        temp[sampledNodesInds] <- 0

        countVec <- countVec + temp

        # Add to the included vector if the node was not a seed gene
        temp <- rep(1,numberOfNodes)
        temp[sampledNodesInds] <- 0

        includedVec <- includedVec + temp

    }

    if(any(includedVec == 1)) {
        warning("Atleast one node was always sampled, likely because it was alone in a bucket (no other node with similar degree)")
    }

    permutationScores <- countVec/includedVec

    #return(list(permutationScores, countVec, includedVec))
    return(permutationScores)
}

# A function that takes in a network, a dataframe of gene to disease associations, a cuttoff value value and a normalzation function
# Returns a list of lists, each list contains network propagation scores for a 

runNetProp <- function(network, assocData, cutoff = c("value" = 0.5, "number" = 10),
                         binarize = TRUE, damping = 0.85, NormFunc = NULL, settingsForNormFunc = NULL,returnSeedVec = FALSE) {
    
    # Register the parallel backend
   # no_cores <- 4
    #cl <- makeCluster(no_cores)
    #registerDoParallel(cl)

    # Filter the association data for diseases that have at least n seed genes with a score >= cutoff
    assocDataFilt <- assocData %>% 
                            filter(score >= cutoff[["value"]]) %>% 
                            group_by(diseaseId) %>% 
                            filter(n() >= cutoff[["number"]]) %>%
                            ungroup()

    # Ensure that the association data is unique
    if( nrow(unique(assocDataFilt[,c("diseaseId","targetId")])) != nrow(assocDataFilt)) {
        print("Association contains duplicates diseaseId-targetId pairs")
        print("Merging duplicates (mean of scores)")
        assocDataFilt <- assocDataFilt %>% 
                        group_by(diseaseId,targetId) %>% 
                        summarise(score = mean(score)) %>% 
                        ungroup()

    }

    # For each disease in the association data
    resList <- list()
      for(trait in unique(assocDataFilt$diseaseId)) {

        # Filter the association data for the current disease
        assocDataFiltTemp <- assocDataFilt %>% filter(diseaseId == trait)

        # Create a seed list from the association data
        seedList <- list(assocDataFiltTemp$targetId, assocDataFiltTemp$score, rep(TRUE,nrow(assocDataFiltTemp)))

        # Create seed vector using the seed list
        seeds <- createSeedVector(network, seedList, binarize = binarize)

        if(returnSeedVec) {
            resList[[trait]] <- seeds[[1]]
            next
        }

        if(all(seeds[[1]] == 0)) {
            print(paste0("No seed genes found for ",trait))
            next
        }


        # Run network propagation with the seed vector
        netProp <- igraph::page_rank(network, directed = FALSE, damping = damping, personalized = seeds[[1]])$vector

        if(!is.null(NormFunc)) {
            netProp <- NormFunc(network, seeds[[1]], netProp, settingsForNormFunc)
        }

        resList[[trait]] <- netProp


    }


    return(do.call(rbind,resList))

}





# Function that takes as input a df of network propagation scores, distance function, and a list of related disease name pairs
# It computes the distance between the network propagation vectors for each disease pair and for all possible pairs of diseases df


compareDistanceMetric <- function(netPropScores, distFunc, distFuncSettings, diseasePairs, randomSetSize,  compareALL = FALSE, diseasesDataFrame, returnDist = FALSE, diseaseShortestPaths = NULL) {
    # Check if diseasePairs is a list, if not convert into a list
    if(is.data.frame(diseasePairs) | is.matrix(diseasePairs)) {
        diseasePairs <- list("_" = diseasePairs)
    }

    # Create a list to hold the results
    resList <- list()

    # Make sure the netPropScores is a matrix
    netPropScores <- as.matrix(netPropScores)

    # Compute the distance matrix for all diseases
    allMetric <- distFunc(netPropScores,distFuncSettings)

    # Return the distance matrix if requested
    if(returnDist) {
        resList[["Rawdist"]] <- allMetric
    }

    # Compute the summary statistics of the distances for all diseases
    allMetricVec <- as.vector(allMetric)
    resList[["summaryAll"]] <- c(summary(allMetricVec),"SD." = sd(allMetricVec))
    resList[["DensityAll"]] <- density(allMetricVec,na.rm = TRUE)
    resList[["HistAll"]] <- hist(allMetricVec,plot=FALSE)
    rm(allMetricVec)

    # Function for going from matrix indices to 1d dist object indices
    distdex <-function(i,j,n) {
        return(n*(i-1) - i*(i-1)/2 + j-i)
    }


    excludeFromRandomInds <- c()
    relatedVecComb <- c()

    #print(str(diseasePairs))

    for(pair in names(diseasePairs)) {

        diseasePairsTemp <- diseasePairs[[pair]]

        # Create a vector that maps the rownames of netPropScores to the row indices of the distance matrix
        indsToRow <- setNames(1:nrow(netPropScores),rownames(netPropScores))
        
        # Compute the indices of the related diseases in the distance matrix
        relatedDistInds <- rep(NA,nrow(diseasePairsTemp))
        
        for (i in 1:nrow(diseasePairsTemp)) {
            relatedDistInds[i] <- distdex(indsToRow[diseasePairsTemp[i,1]],indsToRow[diseasePairsTemp[i,2]],length(indsToRow))
        }

        excludeFromRandomInds <- c(excludeFromRandomInds,relatedDistInds)
        
        # Save results for the related pairs
        relatedVec <- as.vector(allMetric[relatedDistInds])
        relatedVecComb <- c(relatedVecComb,relatedVec)

        resList[[paste0("summaryRelated_",pair)]] <- c(summary(relatedVec),"SD." = sd(relatedVec), "n" = length(relatedVec))

        resList[[paste0("related_",pair)]] <- relatedVec

    }

    ## Sample random distances from allMetric 
    excludeFromRandomInds <- unique(excludeFromRandomInds)
    randomVec <- sample(as.vector(allMetric)[-excludeFromRandomInds],length(excludeFromRandomInds)*randomSetSize)

    # Save results for the random pairs
    resList[["summaryRandom"]] <- c(summary(randomVec),"SD." = sd(randomVec), "n" = length(randomVec))

    # Estimate the density of the all related pairs and random pairs together
    densityEst <- density(c(relatedVecComb,randomVec),na.rm = TRUE)
    densityEstRand <- density(randomVec,na.rm = TRUE,from = min(densityEst$x),to = max(densityEst$x),bw = densityEst$bw)
    resList[["densityEstRand"]] <- densityEstRand

    for(pair in names(diseasePairs)) {
        relatedVec <- resList[[paste0("related_",pair)]]
        resList[[paste0("ratios_",pair)]] <- c("mean" = mean(relatedVec)/mean(randomVec), "sd" = sd(relatedVec)/sd(randomVec), "median" = median(relatedVec)/median(randomVec))
            # Compute the auroc for the related and random pairs
        resList[[paste0("AUROC_",pair)]] <- as.numeric(pROC::roc(response = c(rep(1,length(relatedVec)),
                                                    rep(0,length(randomVec))),
                                                    predictor = c(relatedVec,randomVec), quiet = TRUE)$auc)

        # Estimate the density of the related and random pairs separately using the same x values and bandwidths 
        densityEstRel <- density(relatedVec,na.rm = TRUE,from = min(densityEst$x),to = max(densityEst$x),bw = densityEst$bw)

        # Compute the Jensen Shannon divergence and store the results
        resList[[paste0("JSD_",pair)]] <- philentropy::JSD(rbind(densityEstRel$y,densityEstRand$y))

        resList[[paste0("densityEstRel_",pair)]] <- densityEstRel

    }


    ##########################################################
    ######## CLUSTERING AND DIMENSIONALITY REDUCTION #########
    ##########################################################

    # Set the preferred number of clusters to the square root of the number of diseases
    preferredClusterNumber <- round(sqrt(nrow(netPropScores))*1.5) 
    kparam <- round(sqrt(nrow(netPropScores)))

    # Check if the distFunc gave back disance matrix and convert it if not (for downstream tasks) 
    if(!("returnDist" %in% names(distFuncSettings))) {
        # Convert the similarity matrix to a distance matrix
        print("Converting similarity matrix to distance matrix")
        print("Stats before conversion:")
        print(summary(allMetric))
        allMetric <- proxy::as.dist(allMetric) # assumes that the input is a similarity matrix
        print("Stats after conversion:")
        print(summary(allMetric))

    }

    if(any(allMetric <  1e-10)) {
        print("Atleast one distance is 0, setting to the smallest non zero distance")
        allMetric[allMetric < 1e-10] <- 1e-10 # I dont know what value will make isoMDS not fail
        print(summary(allMetric))

    }


    # Compute the correlation between the disease shortest paths and the distance matrix
    if(!is.null(diseaseShortestPaths)) {

        sharedDiseases <- base::intersect(rownames(netPropScores), attr(diseaseShortestPaths,"Labels"))

        # Convert to matrices and reconvert to dists (this is dumb)
        distForCorr <- as.dist(as.matrix(allMetric)[sharedDiseases,sharedDiseases])
        shortPathForCorr <- as.dist(as.matrix(allMetric)[sharedDiseases,sharedDiseases])

        resList[["DistGraphPCC"]] <- cor(distForCorr,shortPathForCorr, method = "pearson")
        reresList[["DistGraphSCC"]] <- cor(distForCorr,shortPathForCorr, method = "spearman")

    }

    #str(as.dist(as.matrix(dist(t(iris[,1:4])))[colnames(iris)[1:2],colnames(iris)[1:2]]))

    # Compute low dim embeddings
    resList[["cMDS"]] <- cmdscale(allMetric,k = 2)
    resList[["isoMDS"]] <- MASS::isoMDS(allMetric,k = 2)$points
    resList[["UMAP"]] <- umap::umap(as.matrix(allMetric),input = "dist")$layout	

    resList[["HClust"]] <- hclust(allMetric)
    resList[["HClustKcut"]] <- cutree(resList[["HClust"]],k = preferredClusterNumber)

    # Run dynamic tree cut with various minClusterSize and take the one with the number of clusters closest to preferredClusterNumber
    #resList[["HClustDTC"]] <- dynamicTreeCut::cutreeDynamic(resList[["HClust"]],minClusterSize = 5)

    DTClist <- list()
    nclust <- c()
    index <- 1
    for(minClusterSize in seq(4,20,1)) {
        DTClist[[as.character(minClusterSize)]] <- dynamicTreeCut::cutreeDynamic(resList[["HClust"]],minClusterSize = minClusterSize)
        nclust[index] <- length(unique(DTClist[[as.character(minClusterSize)]]))
        index <- index + 1
    }
    bestMinClusterSize <- as.character(seq(4,20,1)[which.min(abs(nclust - preferredClusterNumber))])
    resList[["HClustDTC"]] <- DTClist[[bestMinClusterSize]]


    # Run leiden with different resolutions and take the one with the number of clusters closest to preferredClusterNumber
    leidList <- list()
    nclust <- c()
    index <- 1
    snnGraph <- dbscan::sNN(allMetric,k = kparam,kt=kparam/3) %>%
                                    dbscan::adjacencylist() %>%
                                    igraph::graph_from_adj_list() %>%
                                    igraph::as.undirected() %>% # Add edge weights to the graph from similarity matrix
                                    igraph::simplify()
    

    # Add edge weights to the graph from similarity matrix
    #E(snnGraph)$weight <- as.matrix(proxy::as.simil(allMetric))[igraph::get.edgelist(snnGraph)]


    for(res in seq(0.01,2.01,0.05)) {
        leidList[[as.character(res)]] <- snnGraph %>%
                                igraph::cluster_leiden(resolution_parameter = res,n_iterations = 13) %>%
                                igraph::membership() 
                                
        nclust[index] <- length(unique(leidList[[as.character(res)]]))
        index <- index + 1
        }
    
    bestResLeiden <- as.character(seq(0.01,1.01,0.05)[which.min(abs(nclust - preferredClusterNumber))])
    resList[["Leiden"]] <- leidList[[bestResLeiden]]


    # Create plots of the 

    
    coords <- c("cMDS","isoMDS","UMAP")
    clustAlgos <- c("HClustKcut","Leiden","HClustDTC")

    for(clustAlgo in clustAlgos){

        # Get the frequent ancestors of the diseases in each cluster
        resList[[paste0("commonAncestors_",clustAlgo)]] <- getCommonClusterAncestors(diseasesDataFrame,
                                                                    resList[[clustAlgo]],
                                                                    rownames(netPropScores))

        # Create a legend plot that shows the colors of the clusters and the frequent ancestors of the diseases in each cluster
        ## Create dataframe with cluster number and concatenated ancestors
        legendDF <- data.frame("Clusters" = min(resList[[clustAlgo]]):max(resList[[clustAlgo]]),
                                "Ancestors" = sapply(resList[[paste0("commonAncestors_",clustAlgo)]],
                                                                function(x) paste(names(x[1:min(7,length(x))]),collapse = ", ")),
                                "Counts" = as.vector(table(resList[[clustAlgo]])),
                                "Entropy" = sapply(resList[[paste0("commonAncestors_",clustAlgo)]],
                                                    function(x) DescTools::Entropy(x)/log2(length(x)))
                                            )


        legendDF <- legendDF[order(legendDF$Counts,decreasing = TRUE),]
        legendDF$yCoord <- 1:nrow(legendDF)
        legendDF$LineYs <- legendDF$yCoord + 0.5


        # Add a \n to the ancestors strings if the nchar exceeds 60. Add only one \n to the string at the first comma after 60 characters
        legendDF$Ancestors <- gsub("(.{100},)","\\1\n",legendDF$Ancestors)

        sortClustLabs <- names(sort(table(resList[[clustAlgo]]),decreasing = TRUE))

        resList[[paste0("legendDF_",clustAlgo)]] <- legendDF
    

        for(coord in coords) {
            coordDF <- as.data.frame(resList[[coord]]) 
            colnames(coordDF) <- c("Dim1","Dim2")
            rownames(coordDF) <- rownames(netPropScores)
            pairDF <- createPairDF(coordDF,diseasePairs)
            coordDF$Cluster <- factor(resList[[clustAlgo]],levels = sortClustLabs)

            resList[[paste0("coorDF_",clustAlgo,"_",coord)]] <- coordDF
            resList[[paste0("pairDF_",clustAlgo,"_",coord)]] <- pairDF

        }
    }

    # Remove any occurences of "__" in the names of the items in resList
    names(resList) <- gsub("__","",names(resList))


    return(resList)	

}

# Function that takes in a matrix coordinates and a list of pairs of rows and returns a dataframe where each row contains the coordinates of the two rows in the pair

createPairDF <- function(coords,pairs) {
    # Combine pairs list into a single dataframe and remove any duplicated rows and rows that contain ids that are not in the coords dataframe
    pairs <- do.call(rbind,pairs) %>%
                distinct() %>%
                filter(.[[1]] %in% rownames(coords) & .[[2]] %in% rownames(coords))

    resDF <- data.frame(matrix(NA,nrow = nrow(pairs),ncol = ncol(coords)*2))
    for(i in 1:nrow(pairs)) {
        resDF[i,] <- unlist(c(coords[pairs[i,1],],coords[pairs[i,2],]))
    }
    colnames(resDF) <- c("x1","y1","x2","y2")
    return(resDF)
}

diffOne <- function(x) {
    return(abs(1-x))
}

    
computeDistance <- function(matrix, distFuncSettings) {
    method <- distFuncSettings$method
    # Corrlation functions sans kendall
    if (method %in% c("pearson", "spearman")) {
        if ("returnDist" %in% names(distFuncSettings)){
            return(as.dist(abs(1 - cor(t(matrix), method = method))))
        } else {
            return(proxy::as.simil(cor(t(matrix), method = method)))
        }
    } else if (method == "cosine") {
        temp <- proxy::simil(matrix, method = method)
        # Sharpen the cosine similarity
        if ("p" %in% names(distFuncSettings)) {
            temp <- temp^distFuncSettings$p
        }
        # Check if the distance
        if ("returnDist" %in% names(distFuncSettings)){
            return(proxy::as.dist(temp,FUN = diffOne))
        } else {
            return(temp)
        }
    # This is a faster implementation of the kendall
    } else if (method == "kendall") {
        if ("returnDist" %in% names(distFuncSettings)){
            return(as.dist(abs(1 - pcaPP::cor.fk(t(matrix)))))
        } else {
            return(proxy::as.simil(pcaPP::cor.fk(t(matrix))))
        }
    } else if (method == "jsd") {
        # Jenson shannon distance (not divergence!)
        if(nrow(matrix) == 2) {
            val <- unname(philentropy::JSD(matrix))
            return(sqrt(as.dist(matrix(c(0,val,val,0),nrow = 2))))
        } else {
            return(sqrt(as.dist(philentropy::JSD(matrix))))
        }
    # Just the default distance function in R with p serving as the minkowski dist parameter
    } else {
        if ("p" %in% names(distFuncSettings)) {
            return(dist(matrix, method = method, p = distFuncSettings$p))
        } else {
            return(dist(matrix, method = method))
        }
    }
}

safeRound <- function(x, digits = 0) {
    if (is.null(x)) {
        return(NA)
    } else {
        return(round(x, digits))
    }
}


generatePlotsFromDistCompareResults <- function(resList,diseaseMapping = NULL, densityOnly = FALSE) {


    plotList <- list()

    plotList[["DensityPlot"]] <- ggplot(data.frame("x" = resList[["densityEstRand"]]$x, "y" = resList[["densityEstRand"]]$y), aes(x = x, y = y)) +
            geom_line() +
            #geom_line(data = data.frame("x" = resList[["densityEstRand"]]$x, "y" = resList[["densityEstRand"]]$y), aes(x = x, y = y), col = "blue") +
            geom_line(data = data.frame("x" = resList[["densityEstRel"]]$x, "y" = resList[["densityEstRel"]]$y), aes(x = x, y = y), col = "red") +
            theme_classic() +
            labs(title = "Density plot of related and random diseases",
                x = "Distance",
                y = "Density") + 
            # Add AUROC and JSD to the plot
            annotate("text", x = median(resList[["densityEstRand"]]$x), y =  mean(range(resList[["densityEstRand"]]$y)),
                    label = paste("AUROC:",safeRound(resList[["AUROC"]],2),"\n",
                                    "JSD:",safeRound(resList[["JSD"]],2)), size = 5, color = "black") 
                                 #   xlim(min(resList[["densityEstRand"]]$x),max(resList[["densityEstRand"]]$x))
    if(densityOnly) {
        return(plotList)
    }
    
    coords <- c("cMDS","isoMDS","UMAP")
    clustAlgos <- c("HClustKcut","Leiden","HClustDTC")


    if(is.null(diseaseMapping)) {
        diseaseMapping <- data.frame("id" = rownames(resList[[paste0("coorDF_",clustAlgos[1],"_",coords[1])]]), "name" = rownames(resList[[paste0("coorDF_",clustAlgos[1],"_",coords[1])]]))
        idToName <- setNames(diseaseMapping$name, diseaseMapping$id)
    } else {
        idToName <- setNames(diseaseMapping$name, diseaseMapping$id)
    }


    for(clustAlgo in clustAlgos) {

        # Generat a kables table of the legendDF
        legendDF <- resList[[paste0("legendDF_",clustAlgo)]]

        #[[paste0("legendTable_",clustAlgo)]] <- kable(legendDF[,c("Clusters","Counts","Entropy","Ancestors")],caption = paste0("Legend for ",clustAlgo," Clustering"))

        for(coord in coords) {
            coordDF <- resList[[paste0("coorDF_",clustAlgo,"_",coord)]]
            pairDF <- resList[[paste0("pairDF_",clustAlgo,"_",coord)]]
            clusterCenters <- coordDF %>% group_by(Cluster) %>% summarise(Dim1 = median(Dim1),Dim2 = median(Dim2)) %>% ungroup()

            coordDF$Cluster <- as.factor(as.character(coordDF$Cluster))

            coordDF$labels <- paste0(rownames(coordDF), " (", idToName[rownames(coordDF)], ")")

            plotList[[paste0("Plot_",clustAlgo,"_",coord)]]  <- ggplot() +
                    geom_point(data = coordDF,aes_string(x = "Dim1",y = "Dim2",color = "Cluster",text="labels")) +
                    geom_segment(data = pairDF,aes_string(x = "x1",y = "y1",xend = "x2",yend = "y2"),alpha = 0.03) +
                    geom_text(data = clusterCenters,aes_string(x = "Dim1",y = "Dim2",label = "Cluster"),size = 3) +
                    theme_classic() +
                    labs(title = paste0(coord, " with ", clustAlgo, " Clustering"))  
                    #coord_cartesian(xlim = quantile(coordDF$Dim1,c(0.01,0.99)),ylim = quantile(coordDF$Dim2,c(0.01,0.99))) 

             }
    }




    return(plotList)
}


# Function that takes in the disease R object that contains ancestry and descendent info for diseases, and disease clustering info.
# It returns the top n most common ancestry terms of all the diseases in each cluster

getCommonClusterAncestors <- function(diseases,clusterMembership,diseaseNames) {

    if(any(stringr::str_detect(diseaseNames,"-"))) {
        diseaseNames <- sapply(diseaseNames,function(x) strsplit(x,"-")[[1]][1])
    }

    # Create a list that will hold the results
    resList <- list()

    idToName <- setNames(diseases$name, diseases$id)

    # For each cluster
    for(i in min(clusterMembership):max(clusterMembership)) {
        # Get the diseases in the cluster
        clusterDiseases <- diseaseNames[clusterMembership == i]

        # Get the ancestors of the diseases in the cluster
        clusterAncestors <- lapply(clusterDiseases,function(x) getAncestors(diseases,x))

        # If no disease in cluster has ancestors
        if(all(sapply(clusterAncestors,length) == 0)) {
            resList[[as.character(i)]] <- table(idToName[clusterDiseases])
            next
        }

        # Get the most common ancestors of the diseases in the cluster
        sortedAncestors <- sort(table(unlist(clusterAncestors)),decreasing = TRUE)
    	
        # Change from ids to names
        names(sortedAncestors) <- idToName[names(sortedAncestors)]

        names(sortedAncestors) <- paste0(names(sortedAncestors)," (",round(sortedAncestors/length(clusterDiseases),2),")")

        resList[[as.character(i)]] <- sortedAncestors 

        }

    return(resList)


}

getAncestors <- function(diseases,disease) {
    return(unlist(diseases$ancestors[[which(diseases$id == disease)]]))
}   
    


preprocessNetPropDF <- function(netPropDF,varianceBottomQuantileCuttoff,center,scale) {
    # Scale if needed
    if(center | scale) {
        netPropDF <- scale(netPropDF,center = center,scale = scale)
    }
    # Remove genes with low variance based on quantile
    if(varianceBottomQuantileCuttoff > 0) {
        featureVar <- apply(netPropDF,2,var)
        varCutOff <- quantile(featureVar,varianceBottomQuantileCuttoff)
        netPropDF <- netPropDF[,featureVar > varCutOff]
    }
    return(as.data.frame(netPropDF))
}





permuteTestParalell <- function(network, seedVector, netPropTRUE, settings,damping = 0.85) {
    nSamples <- settings[["nSamples"]]
    ncore <- settings[["ncore"]]

    seedInds <- which(seedVector != 0)
    seedScores <- seedVector[seedInds]

    numberOfNodes <- length(seedVector)
    numberOfSeeds <- length(seedInds)

    # For each permutation
    resMat <- foreach(i = 1:nSamples, .combine = rbind, .options.nws= list(chunkSize=nSamples/ncore)) %dopar% {

        permutedSeedVector <- rep(0, numberOfNodes)
 
        sampledNodesInds <- sample(numberOfNodes, numberOfSeeds)
        permutedSeedVector[sampledNodesInds] <- seedScores

        # Run network propagation with the permuted seed vector
        netPropFALSE <- igraph::page_rank(network, directed = FALSE, damping = damping, personalized = permutedSeedVector)

        # Add to the count vector if the permuted seed vector has a higher score then the true seed vector and the node was not a seed gene
        countVec <- as.numeric(netPropFALSE$vector >= netPropTRUE)
        countVec[sampledNodesInds] <- 0

        countVec

    }

    permutationScores <- as.vector(colSums(resMat)/nSamples)

    return(permutationScores)
}



plotDensityGrid <- function(distFuncs, varCuttOffs, normMethod, dataset, masterRes) {

    listOfPlots <- list()
    for(cuttoff in varCuttOffs) {
        if(cuttoff == "0Center") {
            cuttoff <- "0"
            center <- "TRUE"
        } else {
            center <- "FALSE"
        }
        for(dist in distFuncs) {
            #print(paste(cuttoff,center,dist))
            plots <- generatePlotsFromDistCompareResults( masterRes[[dataset]][[normMethod]][[cuttoff]][[center]][["FALSE"]][[dist]] )
            if(center == "TRUE") {
                cuttoffLab <- paste0(cuttoff," Centered")
            } else {
                cuttoffLab <- cuttoff
            }
            listOfPlots[[paste(cuttoffLab,dist)]] <- plots[[1]] + 
            theme(axis.title.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) 
            if(cuttoffLab == varCuttOffs[1]) {
                listOfPlots[[paste(cuttoffLab,dist)]] <- listOfPlots[[paste(cuttoffLab,dist)]] + ylab(dist)
            } else {
                listOfPlots[[paste(cuttoffLab,dist)]] <- listOfPlots[[paste(cuttoffLab,dist)]] + ylab("")
            }
            if(dist == distFuncs[1]) {
                listOfPlots[[paste(cuttoffLab,dist)]] <- listOfPlots[[paste(cuttoffLab,dist)]] + ggtitle(cuttoffLab)
            } else {
                listOfPlots[[paste(cuttoffLab,dist)]] <- listOfPlots[[paste(cuttoffLab,dist)]] + ggtitle("") }
        }
    }

        plotT <- cowplot::plot_grid(plotlist = listOfPlots, ncol = length(varCuttOffs),byrow = FALSE)
        #title <- cowplot::ggdraw() + cowplot::draw_label(paste0(dataset, " ", normMethod), fontface='bold')
       

        #return(cowplot::plot_grid(title, plotT, ncol=1, rel_heights=c(0.05, 1))) # rel_heights values control title margins)
        return(plotT)


}



# Standalone function that analyzes a given distance matrix 

analyzeDistance <- function(allMetric,diseasePairs,randomSetSize,preferredClusterNumber,kparam,traitIDs,diseasesDataFrame) {

    distdex <-function(i,j,n) {
        return(n*(i-1) - i*(i-1)/2 + j-i)
    }


    indsToRow <- setNames(1:nrow(allMetric),traitIDs)


    relatedDistInds <- rep(NA,nrow(diseasePairs))

    for (i in 1:nrow(diseasePairs)) {
        relatedDistInds[i] <- distdex(indsToRow[diseasePairs[i,1]],indsToRow[diseasePairs[i,2]],length(indsToRow))
    }

    # For each disease pair
    relatedVec <- as.vector(allMetric[relatedDistInds])

    #(head(as.vector(allMetric)[-relatedDistInds]))
    #print(head(relatedDistInds))

    # Sample random distances from allMetric 
    randomVec <- sample(as.vector(allMetric)[-relatedDistInds],length(relatedVec)*randomSetSize)	


    resList <- list()
     # Compute the mean and sd of the distances for the disease pairs and the random pairs
    resList[["summaryRelated"]] <- c(summary(relatedVec),"SD." = sd(relatedVec), "n" = length(relatedVec))
    resList[["summaryRandom"]] <- c(summary(randomVec),"SD." = sd(randomVec), "n" = length(randomVec))
    resList[["ratios"]] <- c("mean" = mean(relatedVec)/mean(randomVec), "sd" = sd(relatedVec)/sd(randomVec), "median" = median(relatedVec)/median(randomVec))


    # Compute the auroc for the related and random pairs
    resList[["AUROC"]] <- as.numeric(pROC::roc(response = c(rep(1,length(relatedVec)),
                                                 rep(0,length(randomVec))),
                                                 predictor = c(relatedVec,randomVec), quiet = TRUE)$auc)
    
    # Compute the Jensen Shanonn divergence between the related and random pairs

    ## Estimate the probability density function of the related and random pairs
    ### First estimate the density of the related and random pairs together
    densityEst <- density(c(relatedVec,randomVec),na.rm = TRUE)

    ### Then estimate the density of the related and random pairs separately using the same x values and bandwidths 
    densityEstRel <- density(relatedVec,na.rm = TRUE,from = min(densityEst$x),to = max(densityEst$x),bw = densityEst$bw)
    densityEstRand <- density(randomVec,na.rm = TRUE,from = min(densityEst$x),to = max(densityEst$x),bw = densityEst$bw)

    ### Compute the Jensen Shannon divergence and store the results
    resList[["JSD"]] <- philentropy::JSD(rbind(densityEstRel$y,densityEstRand$y))

    ### Store the density estimates and histograms of the related and random pairs
    resList[["densityEstRel"]] <- densityEstRel
    resList[["densityEstRand"]] <- densityEstRand
    resList[["histRel"]] <- hist(relatedVec,plot=FALSE)
    resList[["histRand"]] <- hist(randomVec,plot=FALSE)

    allMetricVec <- as.vector(allMetric)
    resList[["summaryAll"]] <- c(summary(allMetricVec),"SD." = sd(allMetricVec))
    resList[["DensityAll"]] <- density(allMetricVec,na.rm = TRUE)
    resList[["HistAll"]] <- hist(allMetricVec,plot=FALSE)
    rm(allMetricVec)

    if(any(allMetric <  1e-10)) {
        print("Atleast one distance is 0, setting to the smallest non zero distance")
        allMetric[allMetric < 1e-10] <- 1e-10 # I dont know what value will make isoMDS not fail
        print(summary(allMetric))

    }

    resList[["cMDS"]] <- cmdscale(allMetric,k = 2)
    resList[["isoMDS"]] <- MASS::isoMDS(allMetric,k = 2)$points
    resList[["UMAP"]] <- umap::umap(as.matrix(allMetric),input = "dist")$layout	

    resList[["HClust"]] <- hclust(allMetric)
    resList[["HClustKcut"]] <- cutree(resList[["HClust"]],k = preferredClusterNumber)

    # Run dynamic tree cut with various minClusterSize and take the one with the number of clusters closest to preferredClusterNumber
    #resList[["HClustDTC"]] <- dynamicTreeCut::cutreeDynamic(resList[["HClust"]],minClusterSize = 5)

    DTClist <- list()
    nclust <- c()
    index <- 1
    for(minClusterSize in seq(4,20,1)) {
        DTClist[[as.character(minClusterSize)]] <- dynamicTreeCut::cutreeDynamic(resList[["HClust"]],minClusterSize = minClusterSize)
        nclust[index] <- length(unique(DTClist[[as.character(minClusterSize)]]))
        index <- index + 1
    }
    bestMinClusterSize <- as.character(seq(4,20,1)[which.min(abs(nclust - preferredClusterNumber))])
    resList[["HClustDTC"]] <- DTClist[[bestMinClusterSize]]


    # Run leiden with different resolutions and take the one with the number of clusters closest to preferredClusterNumber
    leidList <- list()
    nclust <- c()
    index <- 1
    snnGraph <- dbscan::sNN(allMetric,k = kparam,kt=kparam/3) %>%
                                    dbscan::adjacencylist() %>%
                                    igraph::graph_from_adj_list() %>%
                                    igraph::as.undirected() %>% # Add edge weights to the graph from similarity matrix
                                    igraph::simplify()
    

    # Add edge weights to the graph from similarity matrix
    #E(snnGraph)$weight <- as.matrix(proxy::as.simil(allMetric))[igraph::get.edgelist(snnGraph)]


    for(res in seq(0.01,2.01,0.05)) {
        leidList[[as.character(res)]] <- snnGraph %>%
                                igraph::cluster_leiden(resolution_parameter = res,n_iterations = 13) %>%
                                igraph::membership() 
                                
        nclust[index] <- length(unique(leidList[[as.character(res)]]))
        index <- index + 1
        }
    
    bestResLeiden <- as.character(seq(0.01,1.01,0.05)[which.min(abs(nclust - preferredClusterNumber))])
    resList[["Leiden"]] <- leidList[[bestResLeiden]]


    # Create plots of the 

    
    coords <- c("cMDS","isoMDS","UMAP")
    clustAlgos <- c("HClustKcut","Leiden","HClustDTC")

    for(clustAlgo in clustAlgos){

        # Get the frequent ancestors of the diseases in each cluster
        resList[[paste0("commonAncestors_",clustAlgo)]] <- getCommonClusterAncestors(diseasesDataFrame,
                                                                    resList[[clustAlgo]],
                                                                    traitIDs)

        # Create a legend plot that shows the colors of the clusters and the frequent ancestors of the diseases in each cluster
        ## Create dataframe with cluster number and concatenated ancestors
        legendDF <- data.frame("Clusters" = min(resList[[clustAlgo]]):max(resList[[clustAlgo]]),
                                "Ancestors" = sapply(resList[[paste0("commonAncestors_",clustAlgo)]],
                                                                function(x) paste(names(x[1:min(7,length(x))]),collapse = ", ")),
                                "Counts" = as.vector(table(resList[[clustAlgo]])),
                                "Entropy" = sapply(resList[[paste0("commonAncestors_",clustAlgo)]],
                                                    function(x) DescTools::Entropy(x)/log2(length(x)))
                                            )


        legendDF <- legendDF[order(legendDF$Counts,decreasing = TRUE),]
        legendDF$yCoord <- 1:nrow(legendDF)
        legendDF$LineYs <- legendDF$yCoord + 0.5


        # Add a \n to the ancestors strings if the nchar exceeds 60. Add only one \n to the string at the first comma after 60 characters
        legendDF$Ancestors <- gsub("(.{100},)","\\1\n",legendDF$Ancestors)

        sortClustLabs <- names(sort(table(resList[[clustAlgo]]),decreasing = TRUE))

        resList[[paste0("legendDF_",clustAlgo)]] <- legendDF
    

        for(coord in coords) {
            coordDF <- as.data.frame(resList[[coord]]) 
            colnames(coordDF) <- c("Dim1","Dim2")
            rownames(coordDF) <- traitIDs
            pairDF <- createPairDF(coordDF,diseasePairs)
            coordDF$Cluster <- factor(resList[[clustAlgo]],levels = sortClustLabs)

            resList[[paste0("coorDF_",clustAlgo,"_",coord)]] <- coordDF
            resList[[paste0("pairDF_",clustAlgo,"_",coord)]] <- pairDF

        }
    }

    return(resList)


}