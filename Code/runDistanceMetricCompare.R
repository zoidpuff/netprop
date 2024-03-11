netPropPath <- '/home/gummi/netprop'
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))


########### LOAD GRAPH DATA ###########
intData <- read.csv(paste0(netPropPath,"/data/interactionAll.csv"), stringsAsFactors = FALSE)

# Filter out any row that have scoring less then 0.75 (perhaps other values should be tested)

intDataFilt  <- intData %>% filter(!(sourceDatabase == "string" & scoring <= 0.4)) # 0.4 only for stringdb

rm(intData)

intGraph <- graph.data.frame(intDataFilt[,c("targetA","targetB" )], directed = FALSE)

# Remove any self loops and redundant edges

intGraph <- simplify(intGraph, remove.multiple = TRUE, remove.loops = TRUE)

# remove orphan nodes

intGraph <- delete.vertices(intGraph, which(degree(intGraph) == 0))


########### LOAD ASSOCIATION DATA ###########

assocDataBySource <- read.csv(paste0(netPropPath,"/data/associationByDatasourceDirect.csv"), stringsAsFactors = FALSE)
assocDataOverall <-  read.csv(paste0(netPropPath,"/data/associationByOverallDirect.csv"), stringsAsFactors = FALSE)

# loads diseaseDF object
load("/home/gummi/netprop/data/diseases.rdata")

relationshipsAll <- read.csv(paste0(netPropPath,"/relationshipsWithNames.csv"), stringsAsFactors = FALSE)


########### Filter the association data ###########

assocDataOverallEFO <-  assocDataOverall %>%
                         filter(grepl("EFO",diseaseId))


assocDataRareDiseases <- assocDataBySource %>% filter(datasourceId %in% c("orphanet","clingen","eva","eva_somatic")) 

assocDataRareDiseases <- assocDataRareDiseases[,c("score","diseaseId","targetId")] %>% 
  group_by(targetId,diseaseId) %>% 
  summarise(score = max(score)) %>%
    ungroup()


# Create a list of the association data

assocDataList <- list("AllTraitsOverallScores" = assocDataOverall,
                        "EFOTraitsOverallScores" = assocDataOverallEFO,
                        "RareDiseasesOrphaEvaClingenScores" = assocDataRareDiseases)
# Create a list of distance metrics

distanceMetricList <- list(
  "kendallDist" = list("method" = "kendall","returnDist" = NA),
  "kendall" = list("method" = "kendall"),
  "euclidean" = list("method" = "euclidean","returnDist" = NA),
  "manhattan" = list("method" = "manhattan","returnDist" = NA),
  "cosine" = list("method" = "cosine"),
  "cosineSharp2" = list("method" = "cosine","p" = 2),
  "cosineSharp3" = list("method" = "cosine","p" = 3),
  "cosineDist" = list("method" = "cosine","returnDist" = NA),
  "cosineSharp2Dist" = list("method" = "cosine","p" = 2,"returnDist" = NA),
  "cosineSharp3Dist" = list("method" = "cosine","p" = 3,"returnDist" = NA),
  "pearsonDist" = list("method" = "pearson","returnDist" = NA),	
  "spearmanDist" = list("method" = "spearman","returnDist" = NA)
)


normList <- list(list(NULL,NULL,"noNorm"),
                list(ECnormalize,list("logtransform" = FALSE),"ECnorm"),
                list(ECnormalize,list("logtransform" = TRUE),"ECnormLog"),
                list(permuteTestNormalize,list("nSamples" = 300, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm")
)



preprocessList <- list(c("0","FALSE","FALSE"),
                        c("0.25","FALSE","FALSE"),
                        c("0.5","FALSE","FALSE"),
                        c("0.75","FALSE","FALSE"),
                        c("0","TRUE","FALSE"))




# Run the experiments in foreach loop

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(2, detectCores())
cl <- makeCluster(no_cores)

registerDoParallel(cl)

for(dataset in names(assocDataList)){
    for(NORMFUNC in normList) {
        for(PREPROCESS in preprocessList) {
            print(paste0("Started dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_")))
            # Run the netprop algorithm with the association data and the network
            netPropDataFrame <- runNetProp(network = intGraph,
                                assocData = assocDataList[[dataset]],
                                cutoff = c("value" = 0.6, "number" = 12),
                                binarize = TRUE,
                                damping = 0.85,
                                NormFunc = NORMFUNC[[1]],
                                settingsForNormFunc = NORMFUNC[[2]])

            print(paste0("Length of netpropDataFrame: ", nrow(netPropDataFrame)))
            netPropDataFrame[min(100,nrow(netPropDataFrame)),]
            
            relationships <- relationshipsAll %>% filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
            relationships <- as.matrix(relationships[,c("term1","term2")])

            # Preprocess the netprop data
            netPropDataFrame <- preprocessNetPropDF(netPropDataFrame,   as.numeric(PREPROCESS[1]), 
                                                                        as.logical(PREPROCESS[2]), 
                                                                        as.logical(PREPROCESS[3]))

#,.errorhandling = "remove"
            listRes <- foreach(distanceMetric = names(distanceMetricList), .combine = list, .packages = c('dplyr',"ggplot2")) %dopar% {
                   #for(distanceMetric in names(distanceMetricList)){
                    res <- compareDistanceMetric(as.matrix(netPropDataFrame),
                            computeDistance,
                            distanceMetricList[[distanceMetric]],	
                            relationships,
                            6,
                            TRUE,
                            diseaseDF)
                            save(res, file = paste0(netPropPath,
                                "/results/compareDist/netpropDistanceMetricCompare_",
                                dataset,"_",
                                NORMFUNC[[3]],"_",
                                paste0(PREPROCESS,collapse = "_"),
                                "_",distanceMetric,".rdata"))

            }
        }
    }
}


# Stop the parallel backend
stopCluster(cl)
