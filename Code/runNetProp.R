netPropPath <- '/home/gummi/netprop'
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(50, detectCores())
cl <- makeCluster(no_cores)

registerDoParallel(cl)


########### LOAD GRAPH DATA ###########


idToName <- setNames(diseaseDF$name, diseaseDF$id)

rm(diseaseDF)

convergeSuccess <- "notTrue"
counter <- 0
while(!is.logical(convergeSuccess) & counter < 100) {
    res <- igraph::page_rank(intGraph, directed = FALSE, damping = 1, personalized = rep(1,length(V(intGraph))),algo="arpack")
    convergeSuccess <- all.equal(res$value,1)
    counter <- counter + 1
}

if(convergeSuccess != TRUE) {
    stop("Failed to converge")
} 


normList <- list(#list(NULL,NULL,"noNorm"),
                #list(ECnormalize,list("logtransform" = TRUE,"refVec" =res$vector ),"ECnormLog"),
                #list(ECnormalize,list("logtransform" = TRUE,"refVec" =referenceVec$avgVec ),"AverageVecLOR"),
               # list(permuteTestNormalize,list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm"),
                list(permuteTestParalell,list("nSamples" = 100,"ncore" = 0),"permuteNormDegree")
)

for(BINARIZE in c(TRUE)) {
    for(NORMFUNC in normList) {
        print(paste0("Started binarize: ", BINARIZE, " NormFunc: ", NORMFUNC[[3]] ))
        # Replace the outer for loop with a foreach loop
        aurocs <- foreach(trait = unique(assocDataBySourceDirIndiMergedFiltered$diseaseId), .combine = 'rbind', .packages = 'dplyr',.errorhandling = "remove") %dopar% {

            assocDataFiltTemp <- assocDataBySourceDirIndiMergedFiltered %>% filter(diseaseId == trait)

            seedsInd <- assocDataFiltTemp$score >= 0.1

            seedList <- list(assocDataFiltTemp$targetId, assocDataFiltTemp$score, seedsInd)

            diseaseName <- unname(idToName[trait])
            if(length(diseaseName) == 0) {diseaseName <- NA}

            c("trait" = trait,"name" = diseaseName, 
              avgAUROC(network = intGraph,
                       seedList = seedList,
                       nRep = 25,
                       recoverSizeVec = c(0.25, 0.75, 0.9),
                       binarize = BINARIZE,
                       NormFunc = NORMFUNC[[1]],
                       settingsForNormFunc = NORMFUNC[[2]]))
        }

        print(paste0("Finished binarize: ", BINARIZE, " NormFunc: ", NORMFUNC[[3]] ))
        binned <- ifelse(BINARIZE, "_binarized_", "_weighted_")

        write.csv(aurocs,paste0(netPropPath,"/results/overallPhenos_aurocResults",binned,NORMFUNC[[3]],".csv"))
    }
}

# Stop the cluster
stopCluster(cl)
