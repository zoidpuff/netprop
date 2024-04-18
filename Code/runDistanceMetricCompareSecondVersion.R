netPropPath <- '/home/gummi/netprop'
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))

library(doParallel)
library(foreach)

# Register the parallel backend
no_cores <- min(10, detectCores())
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

intGraph <- delete_vertices(intGraph, which(degree(intGraph) == 0))


########### LOAD ASSOCIATION DATA ###########

load(paste0(netPropPath,"/data/shortestPATHSDiseasesFilt.RData"))


assocDataBySourceDirectFiltered <- read.csv(paste0(netPropPath,"/data/associationByDatasourceDirectFiltered.csv"), stringsAsFactors = FALSE)
assocDataBySourceDirIndiMergedFiltered <- read.csv(paste0(netPropPath,"/data/associationByDatasourceDirIndirMergedFiltered.csv"), stringsAsFactors = FALSE)

#allPossibleDiseases <- unique(c(assocDataBySourceDirectFiltered$diseaseId,assocDataBySourceDirIndiMergedFiltered$diseaseId))
#intersectionOfDisease <- intersect(allPossibleDiseases,attr(shortestPATHS,"Labels"))
#shortestPATHS <- as.dist(as.matrix(shortestPATHS)[intersectionOfDisease,intersectionOfDisease])



# loads diseaseDF object
load(paste0(netPropPath,"/data/diseasesNew.rdata"))

relationshipsAll <- read.csv(paste0(netPropPath,"/relationshipsWithNames.csv"), stringsAsFactors = FALSE)

# Create a list of the association data

assocDataList <- list(#"assocDataBySourceDirectFiltered" = assocDataBySourceDirectFiltered,
                      "assocDataBySourceDirIndiMergedFiltered" = assocDataBySourceDirIndiMergedFiltered)

# Create a list of distance metrics
distanceMetricList <- list(
 # "kendall" = list("method" = "kendall"),
 # "euclidean" = list("method" = "euclidean","returnDist" = NA),
 # "manhattan" = list("method" = "manhattan","returnDist" = NA),
 # "minkowski05" = list("method" = "minkowski","returnDist" = NA,"p" = 0.5),
 # "cosine" = list("method" = "cosine"),
 # "cosineSharp2" = list("method" = "cosine","p" = 2),
 # "pearson" = list("method" = "pearson"),	
 # "spearman" = list("method" = "spearman"),
  "jsd" = list("method" = "jsd")
)


normList <- list(list(NULL,NULL,"noNorm"),
               # list(ECnormalize,list("logtransform" = FALSE),"ECnorm"),
                list(ECnormalize,list("logtransform" = TRUE),"ECnormLog"),
               # list(permuteTestNormalize,list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm"),
                list(permuteTestParalell,list("nSamples" = 200,"ncore" = no_cores),"permuteNormDegree")
)



preprocessList <- list(c("0","FALSE","FALSE"),
                        c("0.25","FALSE","FALSE"),
                        c("0.5","FALSE","FALSE"),
                        c("0.75","FALSE","FALSE"),
                        c("0","TRUE","FALSE"))



# Set the vars for testing
#dataset <- names(assocDataList)[2]
#NORMFUNC <- normList[[1]]
#PREPROCESS <- preprocessList[[1]]
#distanceMetric <- names(distanceMetricList)[9]

gc()
# Run the experiments in foreach loop
for(dataset in names(assocDataList)){
    for(NORMFUNC in normList) {
            # Run the netprop algorithm with the association data and the network
            print(paste0("Started dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]]))

            netPropDataFrame <- runNetProp(network = intGraph,
                    assocData = assocDataList[[dataset]],
                    cutoff = c("value" = 0.5, "number" = 5),
                    binarize = TRUE,
                    damping = 0.85,
                    NormFunc = NORMFUNC[[1]],
                    settingsForNormFunc = NORMFUNC[[2]])

            print(paste0("Length of netpropDataFrame: ", nrow(netPropDataFrame)))
            
            relationships <- relationshipsAll %>% filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
            relationships <- as.matrix(relationships[,c("term1","term2")])
            gc()
            temp <- foreach(PREPROCESS = preprocessList,.combine = list) %:%
                        foreach(distanceMetric = names(distanceMetricList), .combine = list, .packages = c('dplyr',"ggplot2"),.errorhandling = "remove") %dopar% {
                            #for(PREPROCESS in preprocessList){
                             #   for(distanceMetric in names(distanceMetricList)){
                                    # Preprocess the netprop data
                                    netPropDataFramePP <- preprocessNetPropDF(netPropDataFrame, as.numeric(PREPROCESS[1]), 
                                                                                                as.logical(PREPROCESS[2]), 
                                                                                                as.logical(PREPROCESS[3]))
                                    
                                    if(any(0>as.matrix(netPropDataFramePP))){return(NULL)}

                                #for(distanceMetric in names(distanceMetricList)){
                                    cat(file="internalStarted.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_"), " DistanceMetric: ", distanceMetric,"\n"))	
                                    res <- compareDistanceMetric(as.matrix(netPropDataFramePP),
                                            computeDistance,
                                            distanceMetricList[[distanceMetric]],	
                                            relationships,
                                            8,
                                            TRUE,
                                            diseaseDF,
                                            FALSE,
                                            shortestPATHS)

                                            save(res, file = paste0(netPropPath,
                                                "/results/compareDist/netpropDistanceMetricCompare_",
                                                dataset,"_",
                                                NORMFUNC[[3]],"_",
                                                paste0(PREPROCESS,collapse = "_"),
                                                "_",distanceMetric,".rdata"))
                                    
                                    gc()
                                    
                                    cat(file="internalFinished.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " Preprocess: ", paste0(PREPROCESS,collapse = "_"), " DistanceMetric: ", distanceMetric,"\n"))	

                               # }
        }
    }
}




# Stop the parallel backend
stopCluster(cl)
