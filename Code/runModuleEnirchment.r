library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(34, detectCores()) 
cl <- makeCluster(no_cores)

registerDoParallel(cl)

netPropPath <- "/cluster/home/gmagnusson/netprop"

load(paste0(netPropPath,"/results/fgseaData_ecdGwasMouse.rdata"))
genesets <- inputForFGSEA$genesets
netPropRes <- as.data.frame(inputForFGSEA$netPropRes)
netPropRes$traitID <- rownames(inputForFGSEA$netPropRes)
rm(inputForFGSEA)
#test = iter(netPropRes, by = "row")
# Add the rownames if netPropRes as a column


# For each row in the netPropDataFrame run the fsgea test against the genesets
#fgsea <- list()

#for(trait in rownames(netPropRes)) {
#     temp <- fgsea::fgsea(pathways = genesets,
#                              stats = netPropRes[trait,],
#                              minSize = 3, maxSize = Inf,eps = 0)
#                              fgsea[[trait]] <- temp[,c("pathway","padj","pval")]
#                              
#}

# Do parallel processing
 fgsea <- foreach(trait = iter(netPropRes, by = "row"), .combine = rbind) %dopar% {
     statsVec <- setNames(as.numeric(trait[1:(length(trait)-1)]),colnames(trait)[1:(length(trait)-1)])
     temp <- fgsea::fgsea(pathways = genesets,
                  stats = statsVec,
                  minSize = 3, maxSize = Inf,eps = 0)
     temp <- temp[temp$padj < 0.005,c("pathway","padj","pval")]
     temp$traitID <- trait$traitID
        temp
 }


save(fgsea, file = paste0(netPropPath,"/results/fgseaRes_ecdGwasMouse.rdata"))
