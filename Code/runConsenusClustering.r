
library(igraph)
load(paste0(netPropPath,"/data/intGraphSub.rdata"))


clusterLargeClusters <- function(graph, clustering, n, consensus = FALSE) {
    largeClusters <- names(table(clustering))[table(clustering) > n]
    modifiedClustering <- clustering
    
    for (cluster in largeClusters) {
        subGraph <- igraph::induced_subgraph(graph, which(clustering == cluster))
        if(consensus) {
            subMembership <- consensus_clustering_walktrap(subGraph, 30)
        } else {
            subMembership <- igraph::cluster_walktrap(subGraph)$membership
        }
        
        # Update the modified clustering with the subcluster membership
        modifiedClustering[clustering == cluster] <- paste0(cluster, ".", subMembership)
    }
    
    return(modifiedClustering)
}

consensus_clustering <- function(g, n_iter, subsample_size, resolution) {
  # Initialize partition matrix
  partitions <- list()
  
  # Make a new grapgh with same nodes and no edges
  consensus_graph <- make_empty_graph(vcount(g))
  
  for (i in 1:n_iter) {
    # Randomly subsample edges
    edges <- sample(E(g), round(subsample_size*ecount(g)))
    subgraph <- subgraph.edges(g, eids = edges,delete.vertices = FALSE)

    #Ensure that the subgraph is has same number of nodes as the original graph
    if(vcount(subgraph) != vcount(g)){stop("Subgraph has different number of nodes than the original graph") }
    # Randomly sample resolution parameter for leiden clustering
    #sampleRes <- -1
    #while(sampleRes < 0) {
    #    sampleRes <- rnorm(1,mean = resolution,sd = 0.005)
    #}

    res <- sample(resolution,1)
    res <- res + rnorm(1,mean = 0,sd = res/10)
    if(res < 0) {res <- 0.1}
    
    # Perform Leiden clustering
    communities <- igraph::cluster_leiden(subgraph, n_iterations =30, resolution = resolution)

    # Break up large clusters with walktrap clustering (several times)
    communitiesAdj <- clusterLargeClusters(subgraph,communities$membership,300)
    communitiesAdj <- clusterLargeClusters(subgraph,communitiesAdj,300)
    communitiesAdj <- clusterLargeClusters(subgraph,communitiesAdj,200)
    partitions[[i]] <- communitiesAdj

   }
   partitions <- do.call(rbind,partitions)

    # Update the consensus graph with the edge weights (this takes a loong time)
    for (v1 in 1:vcount(g)) {
        for (v2 in v1:vcount(g)) {
            if (v1 != v2) {
                val <- sum(partitions[,v1] == partitions[,v2])/n_iter
                if (val > 0.3) {
                    consensus_graph <- add_edges(consensus_graph, c(v1, v2), weight = val)
                } 
            } 
        }
    }
  #consensus_clustering <- igraph::cluster_leiden(consensus_graph, n_iterations = 100, resolution = resolution)
  #return(list(consensus_clustering,consensus_graph))
  return(consensus_graph)
}



resolutionList <- c(seq(0.001,0.1,0.005),seq(0.1,0.31,0.1))
conMat <- consensus_clustering(intGraphSub, n_iter = 1000, subsample_size = 0.8, resolution = resolutionList)

save(conMat,file = paste0(netPropPath,"/data/conMatBigRes.rdata"))