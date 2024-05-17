library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)


consensus_clustering_walktrap <- function(g, n_iter) {
  # Initialize partition matrix
  partitions <- list()
  
  # Make a new graph with same nodes and no edges
  consensus_graph <- make_empty_graph(vcount(g))
  
  for (i in 1:n_iter) {
    # Perform walktrap clustering
    communities <- igraph::cluster_walktrap(g)
    partitions[[i]] <- communities$membership

   }
   partitions <- do.call(rbind,partitions)

    # Update the consensus graph with the weighted edges only if the edge is present in more than 50% of the partitions
    for (v1 in 1:vcount(g)) {
        for (v2 in v1:vcount(g)) {
            if (v1 != v2) {
                val <- sum(partitions[,v1] == partitions[,v2])/n_iter
                if (val > 0.5) {
                    consensus_graph <- add_edges(consensus_graph, c(v1, v2), weight = val)
                } 
            } 
        }
    }
  consensus_clustering <- igraph::cluster_walktrap(consensus_graph)
  return(consensus_clustering$membership)
}

clusterLargeClusters <- function(graph, clustering, n, consensus = FALSE) {
    largeClusters <- names(table(clustering))[table(clustering) > n]
    modifiedClustering <- clustering
    
    for (cluster in largeClusters) {
        subGraph <- igraph::induced_subgraph(graph, which(clustering == cluster))
        if(consensus) {
            subMembership <- consensus_clustering_walktrap(subGraph, 50)
        } else {
            subMembership <- igraph::cluster_walktrap(subGraph)$membership
        }
        
        # Update the modified clustering with the subcluster membership
        modifiedClustering[clustering == cluster] <- paste0(cluster, ".", subMembership)
    }
    
    return(modifiedClustering)
}

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(11, detectCores())
cl <- makeCluster(no_cores)

registerDoParallel(cl)





# Load the data

load(file = paste0(netPropPath,"/data/intGraphSub.rdata"))
ensgToEntrez <- data.table::fread(paste0(netPropPath,"/data/geneToEnsmbl.txt"))
ensgToEntrez <- setNames( ensgToEntrez$`NCBI gene (formerly Entrezgene) ID`, ensgToEntrez$`Gene stable ID`)

communities <- igraph::cluster_leiden(intGraphSub,n_iterations = 200,resolution = 0.056)
consensusClustersReclustered <- clusterLargeClusters(intGraphSub, communities$membership, 300, consensus = TRUE)	
consensusClustersReclustered <- clusterLargeClusters(intGraphSub, consensusClustersReclustered, 300, consensus = TRUE)
communities$membership <- as.character(consensusClustersReclustered)

communityDF <- data.frame(community = unique(communities$membership),
                        genes = paste0(sapply(unique(communities$membership), function(x) V(intGraphSub)$name[communities$membership == x]), collapse = ","),
)

#communityDF <- read.csv(paste0(netPropPath,"/data/communityDF.csv"))
#load(file = paste0(netPropPath,"/data/communities.rdata"))

# Load the gene set data

enrichRes <- foreach(i = communityDF$community, .packages = c("org.Hs.eg.db"),.combine = list,.errorhandling = "pass") %dopar% {
    genes <- igraph::V(intGraphSub)$name[communities$membership == i]
    returnList <- list()
    if(length(genes) < 3) {
        return(NULL)
    }
    ego1 <- clusterProfiler::enrichGO(genes,
                                    keyType = "ENSEMBL",
                                    OrgDb = org.Hs.eg.db,
                                    ont = "BP",
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH",
                                    minGSSize = 5)

    # Filter out p.adjusted values that are greater than 0.05
    returnDFbp <- ego1@result[ego1@result$p.adjust < 0.05,]
    returnDFbp$cluster <- i
    returnDFbp$db <- "GOBP"

    ego2 <- clusterProfiler::enrichGO(genes,
                                keyType = "ENSEMBL",
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize = 5)

    # Filter out p.adjusted values that are greater than 0.05
    returnDFcc <- ego2@result[ego2@result$p.adjust < 0.05,]
    returnDFcc$cluster <- i
    returnDFcc$db <- "GOCC"

    ego3 <- clusterProfiler::enrichGO(genes,
                                keyType = "ENSEMBL",
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize = 5)
    
    # Filter out p.adjusted values that are greater than 0.05
    returnDFmf <- ego3@result[ego3@result$p.adjust < 0.05,]
    returnDFmf$cluster <- i
    returnDFmf$db <- "GOMF"

    returnDF <- rbind(returnDFbp,returnDFcc,returnDFmf)

    genes <- ensgToEntrez[genes]
    genes <- genes[!is.na(genes)]
    if(length(genes) < 3) {
        return(returnDF)
    }
    epa <-  ReactomePA::enrichPathway(genes,
                  minGSSize = 5,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")

    # Filter out p.adjusted values that are greater than 0.05
    tempDF <- epa@result[epa@result$p.adjust < 0.05,]
    tempDF$cluster <- i
    tempDF$db <- "Reactome"

    returnDF <- rbind(returnDF,tempDF)
    return(returnDF)
}

save(enrichRes, file = paste0(netPropPath,"/data/enrichRes.rdata"))