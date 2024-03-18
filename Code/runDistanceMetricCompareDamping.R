netPropPath <- '/home/gummi/netprop'
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(21, detectCores())
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

#assocDataBySource <- read.csv(paste0(netPropPath,"/data/associationByDatasourceDirect.csv"), stringsAsFactors = FALSE)
assocDataOverall <-  read.csv(paste0(netPropPath,"/data/associationByOverallDirect.csv"), stringsAsFactors = FALSE)

# loads diseaseDF object
#load("/home/gummi/netprop/data/diseases.rdata")
load(paste0(netPropPath,"/data/diseases.rdata"))

relationshipsAll <- read.csv(paste0(netPropPath,"/relationshipsWithNames.csv"), stringsAsFactors = FALSE)



# Create a list of the association data

assocDataList <- list("AllTraitsOverallScores" = assocDataOverall)
                        #"EFOTraitsOverallScores" = assocDataOverallEFO,
                        #"RareDiseasesOrphaEvaClingenScores" = assocDataRareDiseases)
# Create a list of distance metrics

distanceMetricList <- list(
  #"kendallDist" = list("method" = "kendall","returnDist" = NA),
  #"kendall" = list("method" = "kendall"),
  #"euclidean" = list("method" = "euclidean","returnDist" = NA),
  "manhattan" = list("method" = "manhattan","returnDist" = NA)
  "cosine" = list("method" = "cosine"),
  #"cosineSharp2" = list("method" = "cosine","p" = 2),
  #"cosineSharp3" = list("method" = "cosine","p" = 3),
  #"cosineDist" = list("method" = "cosine","returnDist" = NA),
  #"cosineSharp2Dist" = list("method" = "cosine","p" = 2,"returnDist" = NA),
  #"cosineSharp3Dist" = list("method" = "cosine","p" = 3,"returnDist" = NA),
  #"pearsonDist" = list("method" = "pearson","returnDist" = NA),	
  "spearmanDist" = list("method" = "spearman","returnDist" = NA)
)


normList <- list(list(NULL,NULL,"noNorm"),
                #list(ECnormalize,list("logtransform" = FALSE),"ECnorm"),
                list(ECnormalize,list("logtransform" = TRUE),"ECnormLog"))
               # list(permuteTestNormalize,list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm"),
                #list(permuteTestParalell,list("nSamples" = 200,"ncore" = no_cores),"permuteNormDegree")




dampingList <- list(c("0.95","FALSE","FALSE"),
                    c("0.9","FALSE","FALSE"),
                    c("0.85","FALSE","FALSE"),
                    c("0.8","FALSE","FALSE"),
                    c("0.75","FALSE","FALSE"),
                    c("0.7","FALSE","FALSE"),
                    c("0.65","FALSE","FALSE"),
                    c("0.6","FALSE","FALSE"),
                    c("0.55","FALSE","FALSE"),
                    c("0.5","FALSE","FALSE"),
                    c("0.45","FALSE","FALSE"),
                    c("0.4","FALSE","FALSE"),
                    c("0.35","FALSE","FALSE"),
                    c("0.3","FALSE","FALSE"),
                    c("0.25","FALSE","FALSE"),
                    c("0.2","FALSE","FALSE"),
                    c("0.15","FALSE","FALSE"),
                    c("0.1","FALSE","FALSE"),
                    c("0.05","FALSE","FALSE"))




# Run the experiments in foreach loop


for(dataset in names(assocDataList)){
    for(NORMFUNC in normList) {
        temp <- foreach(dampingFactor = dampingList, .combine = list, .packages = c('dplyr',"ggplot2"),.errorhandling = "stop") %dopar% {
            netPropDataFrame <- runNetProp(network = intGraph,
                                assocData = assocDataList[[dataset]],
                                cutoff = c("value" = 0.5, "number" = 2),
                                binarize = TRUE,
                                damping = as.numeric(dampingFactor[1]),
                                NormFunc = NORMFUNC[[1]],
                                settingsForNormFunc = NORMFUNC[[2]])
            relationships <- relationshipsAll %>% filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
            relationships <- as.matrix(relationships[,c("term1","term2")])
            for(distanceMetric in names(distanceMetricList)){
                cat(file="internalStarted.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Damping: ",dampingFactor[1], " DistanceMetric: ", distanceMetric,"\n"))	
                res <- compareDistanceMetric(as.matrix(netPropDataFrame),
                        computeDistance,
                        distanceMetricList[[distanceMetric]],	
                        relationships,
                        8,
                        TRUE,
                        diseaseDF)
                        save(res, file = paste0(netPropPath,
                            "/results/compareDistDamping/netpropDistanceMetricCompare_",
                            dataset,"_",
                            NORMFUNC[[3]],"_",
                            dampingFactor[1],"_",
                            "_",distanceMetric,".rdata"))
                
                cat(file="internalFinished.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Damping: ",dampingFactor[1], " DistanceMetric: ", distanceMetric,"\n"))	

            }
        }
    }
}



# Stop the parallel backend
stopCluster(cl)
