netPropPath <- '/home/gummi/netprop'
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(18, detectCores())
cl <- makeCluster(no_cores)

registerDoParallel(cl)


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
#load("/home/gummi/netprop/data/diseases.rdata")
load(paste0(netPropPath,"/data/diseases.rdata"))

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
               # list(permuteTestNormalize,list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm"),
                list(permuteTestParalell,list("nSamples" = 200,"ncore" = no_cores),"permuteNormDegree")
)



preprocessList <- list(c("0","FALSE","FALSE"),
                        c("0.25","FALSE","FALSE"),
                        c("0.5","FALSE","FALSE"),
                        c("0.75","FALSE","FALSE"),
                        c("0","TRUE","FALSE"))




# Run the experiments in foreach loop



for(dataset in names(assocDataList)){
    for(NORMFUNC in normList) {
            # Run the netprop algorithm with the association data and the network
            print(paste0("Started dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]]))

            netPropDataFrame <- runNetProp(network = intGraph,
                                assocData = assocDataList[[dataset]],
                                cutoff = c("value" = 0.5, "number" = 2),
                                binarize = TRUE,
                                damping = 0.85,
                                NormFunc = NORMFUNC[[1]],
                                settingsForNormFunc = NORMFUNC[[2]])

            print(paste0("Length of netpropDataFrame: ", nrow(netPropDataFrame)))
            #netPropDataFrame[min(100,nrow(netPropDataFrame)),]
            
            relationships <- relationshipsAll %>% filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
            relationships <- as.matrix(relationships[,c("term1","term2")])

            temp <- foreach(PREPROCESS = preprocessList,.combine = list) %:%
                        foreach(distanceMetric = names(distanceMetricList), .combine = list, .packages = c('dplyr',"ggplot2"),.errorhandling = "remove") %dopar% {
                            # Preprocess the netprop data
                            #print(paste0("Started dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_")))
                            netPropDataFramePP <- preprocessNetPropDF(netPropDataFrame, as.numeric(PREPROCESS[1]), 
                                                                                        as.logical(PREPROCESS[2]), 
                                                                                        as.logical(PREPROCESS[3]))

                        #for(distanceMetric in names(distanceMetricList)){
                            cat(file="internalStarted.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_"), " DistanceMetric: ", distanceMetric,"\n"))	
                            res <- compareDistanceMetric(as.matrix(netPropDataFramePP),
                                    computeDistance,
                                    distanceMetricList[[distanceMetric]],	
                                    relationships,
                                    8,
                                    TRUE,
                                    diseaseDF)
                                    save(res, file = paste0(netPropPath,
                                        "/results/compareDist/netpropDistanceMetricCompare_",
                                        dataset,"_",
                                        NORMFUNC[[3]],"_",
                                        paste0(PREPROCESS,collapse = "_"),
                                        "_",distanceMetric,".rdata"))
                            
                            cat(file="internalFinished.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_"), " DistanceMetric: ", distanceMetric,"\n"))	

            
        }
    }
}



# Stop the parallel backend
stopCluster(cl)
