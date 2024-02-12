
library(igraph)
library(dplyr)

netpropPath <- "/cluster/home/gmagnusson/netprop"

# Load the netprop functions
source(paste0(netpropPath,"/Code/NetPropFuncs.R"))

# Load the network data

pathToCSV <- paste0(netpropPath,"/data/interactionAll.csv")

intData <- read.csv(pathToCSV, stringsAsFactors = FALSE)

intDataFilt  <- intData %>% filter(!(sourceDatabase == "string" & scoring <= 0.4)) # 0.4 only for stringdb

# Create an unweighted graph from the data
intGraph <- graph.data.frame(intDataFilt[,c("targetA","targetB" )], directed = FALSE)
rm(intData, intDataFilt)

# Remove any self loops and redundant edges
intGraph <- simplify(intGraph, remove.multiple = TRUE, remove.loops = TRUE)

# Load the association data and filter for EFO diseases that have at least 10 seed genes with score >= 0.5

pathToAssocCSV <- paste0(netpropPath,"/data/associationByOverallDirect.csv")
assocData <- read.csv(pathToAssocCSV, stringsAsFactors = FALSE)
assocDataFilt10 <- assocData %>% 
                            filter(stringr::str_split(diseaseId, pattern = "_", simplify = TRUE)[,1] == "EFO") %>%
                            filter(score >= 0.5) %>% 
                            group_by(diseaseId) %>% 
                            filter(n() >= 10) %>%
                            ungroup()

rm(assocData)

# Load the EFO to description mapping
diseaseMappings <- read.csv(paste0(netpropPath,"/data/phenotypeMapping.csv"), stringsAsFactors = FALSE) %>%
                    filter(id %in% unique(assocDataFilt10$diseaseId))

print(paste0("Number of traits:", length(unique(assocDataFilt10$diseaseId))))

print("Starting to run netprop")

settingsPerm1 <- list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1)
settingsPerm2 <- list("nSamples" = 100, "perserveDegree" = TRUE, "degreeSampleSmoothing" = 3, "minBucketSize" = 1)


normList <- list(list(NULL,NULL,"noNorm"),
                list(ECnormalize,list("logtransform" = FALSE),"ECnorm"),
                list(ECnormalize,list("logtransform" = TRUE),"ECnormLog"),
                list(permuteTestNormalize,settingsPerm1,"permuteNorm"),
                list(permuteTestNormalize,settingsPerm2,"permuteNormDegree"))


library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(30, detectCores())
cl <- makeCluster(no_cores)
registerDoParallel(cl)

for(BINARIZE in c(TRUE)) {
    for(NORMFUNC in normList) {
        print(paste0("Started binarize: ", BINARIZE, " NormFunc: ", NORMFUNC[[3]] ))

        # Replace the outer for loop with a foreach loop
        aurocs <- foreach(trait = unique(assocDataFilt10$diseaseId), .combine = 'list', .packages = 'dplyr') %dopar% {

            assocDataFiltTemp <- assocDataFilt10 %>% filter(diseaseId == trait)

            seedsInd <- assocDataFiltTemp$score >= 0.5

            seedList <- list(assocDataFiltTemp$targetId, assocDataFiltTemp$score, seedsInd)

            diseaseName <- diseaseMappings[which(diseaseMappings$id == trait),4] 

            c(trait,diseaseName, 
              avgAUROC(network = intGraph,
                       seedList = seedList,
                       nRep = 25,
                       recoverSizeVec = c(0.25, 0.5, 0.75, 0.9),
                       binarize = BINARIZE,
                       NormFunc = NORMFUNC[[1]],
                       settingsForNormFunc = NORMFUNC[[2]]))
        }

        print(paste0("Finished binarize: ", BINARIZE, " NormFunc: ", NORMFUNC[[3]] ))
        binned <- ifelse(BINARIZE, "_binarized_", "_weighted_")

        do.call(rbind,aurocs) %>% as.data.frame() %>% write.csv(paste0(netpropPath,"/results/overallPhenos_aurocResults",binned,NORMFUNC[[3]],".csv"))
    }
}

# Stop the cluster
stopCluster(cl)
