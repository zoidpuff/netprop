netPropPath <- '/home/gummi/netprop'
netPropPath <- "/cluster/home/gmagnusson/netprop"
library(igraph)
library(dplyr)
source(paste0(netPropPath, '/Code/NetPropFuncs.R'))

# Loads shortestPATHS dist object with 
#load("/home/gummi/netprop/data/shortestPATHSDiseases.RData")


########### LOAD GRAPH DATA ###########

if(file.exists(paste0(netPropPath,"/data/intGraph.rdata"))){
    load(paste0(netPropPath,"/data/intGraph.rdata"))
} else {
    intData <- read.csv(paste0(netPropPath,"/data/interactionAll.csv"), stringsAsFactors = FALSE)

    # Filter out any row that have scoring less then 0.75 (perhaps other values should be tested)

    intDataFilt  <- intData %>% filter(!(sourceDatabase == "string" & scoring <= 0.4)) # 0.4 only for stringdb

    rm(intData)

    intGraph <- graph.data.frame(intDataFilt[,c("targetA","targetB" )], directed = FALSE)

    # Remove any self loops and redundant edgess

    intGraph <- simplify(intGraph, remove.multiple = TRUE, remove.loops = TRUE)

    # remove orphan nodes

    intGraph <- delete.vertices(intGraph, which(degree(intGraph) == 0))

    #save(intGraph, file = paste0(netPropPath,"/data/intGraph.rdata"))
}


########### LOAD ASSOCIATION DATA ###########

assocMouseGWAS <- read.csv(paste0(netPropPath,"/data/impcGWASAssocs0504.csv"), stringsAsFactors = FALSE)
mouseTraits <- unique(assocMouseGWAS$diseaseId)
assocMouseGWAS$diseaseId <- paste0(assocMouseGWAS$diseaseId, "-", assocMouseGWAS$source)

assocRareGWAS <-  read.csv(paste0(netPropPath,"/data/rarevsGWASAssocs0505.csv"), stringsAsFactors = FALSE)
rareTraits <- unique(assocRareGWAS$diseaseId)
assocRareGWAS$diseaseId <- paste0(assocRareGWAS$diseaseId, "-", assocRareGWAS$source)

# loads diseaseDF object
#load("/home/gummi/netprop/data/diseases.rdata")
load(paste0(netPropPath,"/data/diseasesNew.rdata"))

relationshipsAll <- read.csv(paste0(netPropPath,"/relationshipsWithNames.csv"), stringsAsFactors = FALSE)

relationshipsRare <- data.frame(term1 = paste0(rareTraits,"-", unique(assocRareGWAS$source)[1]),
                                term2 = paste0(rareTraits,"-", unique(assocRareGWAS$source)[2]))

relationshipsMouse <- data.frame(term1 = paste0(mouseTraits,"-", unique(assocMouseGWAS$source)[1]),
                                term2 = paste0(mouseTraits,"-", unique(assocMouseGWAS$source)[2]))



# Create a list of the association data

assocDataList <- list( "assocRareGWAS" = assocRareGWAS,
                        "assocMouseGWAS" = assocMouseGWAS
                       )


relationshipList <- list("assocRareGWAS" = relationshipsRare,
                        "assocMouseGWAS" = relationshipsMouse)     

# Create a list of distance metrics
distanceMetricList <- list(
  #"kendallDist" = list("method" = "kendall","returnDist" = NA),
 "pearson" = list("method" = "pearson"),
  "kendall" = list("method" = "kendall")
  #"minkowski05" = list("method" = "minkowski","returnDist" = NA,"p" = 0.5),
  #"euclidean" = list("method" = "euclidean","returnDist" = NA),
  #"manhattan" = list("method" = "manhattan","returnDist" = NA),
  #"cosine" = list("method" = "cosine"),
  #"cosineSharp2" = list("method" = "cosine","p" = 2),
  #"cosineSharp3" = list("method" = "cosine","p" = 3),
  #"cosineDist" = list("method" = "cosine","returnDist" = NA),
  #"cosineSharp2Dist" = list("method" = "cosine","p" = 2,"returnDist" = NA),
  #"cosineSharp3Dist" = list("method" = "cosine","p" = 3,"returnDist" = NA),	
  #"spearman" = list("method" = "spearman")
  #"jsd" = list("method" = "jsd")
)

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



normList <- list(#(list(NULL,NULL,"noNorm"),
                list(ECnormalize,list("logtransform" = TRUE,"refVec" =res$vector ),"ECnormLog")
                #list(ECnormalize,list("logtransform" = TRUE,"refVec" =referenceVec$avgVec ),"AverageVecLOR"),
               # list(permuteTestNormalize,list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 0, "minBucketSize" = 1),"permuteNorm"),
                #list(permuteTestParalell,list("nSamples" = 200,"ncore" = no_cores),"permuteNormDegree")
)

# Run the experiments in foreach loop
for(dataset in names(assocDataList)){
    for(NORMFUNC in normList) {

        print(paste0("Started dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]] ))

        netPropDataFrame <- runNetProp(network = intGraph,
                            assocData = assocDataList[[dataset]],
                            cutoff = c("value" = 0.1, "number" = 1),
                            binarize = TRUE,
                            damping = 0.85,
                            NormFunc = NORMFUNC[[1]],
                            settingsForNormFunc = NORMFUNC[[2]],
                            computeStats = list("otGenes" = unique(rareDiseaseAssocsFromOTgenetics$targetId),
                            "evaGenes" = unique(assocDataRareDiseasesCollapsed$targetId),
                            "mouseGenes" = unique(mouseAssocs$targetId)))
        
        
    #save(netPropDataFrame, file = paste0(netPropPath,"/results/netPropDataFrame_ecdMouse_ecnorm.rdata"))

        
        # Remove suffixes from the rownames of netprop
        cleanedTraits <- sapply(strsplit(rownames(netPropDataFrame), "-"), function(x) x[1]) %>% unique()
        
        # Sample 100 traits from the cleanedTraits

        #cleanedTraits <- sample(cleanedTraits, 100)

        relationships <- relationshipList[[dataset]]

         # Filter out rows from the netPropDataFrame that are not in the cleanedTraits adding back the suffix
        #netPropDataFrame <- as.data.frame(netPropDataFrame) %>% 
         #                   filter(rownames(.) %in% paste0(cleanedTraits, "-", unique(strsplit(relationships$term1[1], "-")[[1]][2])) | 
          #                          rownames(.) %in% paste0(cleanedTraits, "-", unique(strsplit(relationships$term2[1], "-")[[1]][2]))) %>%
           #                 as.matrix()


        print(dim(relationships))
        # Remove rows from the relationships that have terms that are not in the rownames for the netPropDataFrame
        relationships <- relationships %>% 
                    filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
        print(dim(relationships))

        relationshipsFromAll <- relationshipsAll %>% 
                    filter(term1 %in% cleanedTraits & term2 %in% cleanedTraits)

        # Get the suffix that comes afer "-" in the relationships dataframe
        dataset1 <- unique(strsplit(relationships$term1[1], "-")[[1]][2])
        dataset2 <- unique(strsplit(relationships$term2[1], "-")[[1]][2])
        
        # Add suffixes to disease ids in relationshipsFromAll
        relationshipsFromAll1 <- data.frame(term1 = paste0(relationshipsFromAll$term1, "-", dataset1),
                                        term2 = paste0(relationshipsFromAll$term2, "-", dataset1))
        

        relationshipsFromAll2 <- data.frame(term1 = paste0(relationshipsFromAll$term1, "-", dataset2),
                                        term2 = paste0(relationshipsFromAll$term2, "-", dataset2))
        
        relationshipsFromAllAcross <- data.frame(term1 = paste0(relationshipsFromAll$term1, "-", dataset1),
                                        term2 = paste0(relationshipsFromAll$term2, "-", dataset2))

        # Create a list of random trait pairs as a control
        traitsDS1 <- paste0(cleanedTraits, "-", dataset1)
        randomRelationships1 <- data.frame(term1 = sample(traitsDS1, 8*length(traitsDS1), replace = TRUE),
                                        term2 = sample(traitsDS1, 8*length(traitsDS1), replace = TRUE)) %>%
                                        filter(term1 != term2)

        traitsDS2 <- paste0(cleanedTraits, "-", dataset2)
        randomRelationships2 <- data.frame(term1 = sample(traitsDS2, 8*length(traitsDS2), replace = TRUE),
                                        term2 = sample(traitsDS2, 8*length(traitsDS2), replace = TRUE)) %>% 
                                        filter(term1 != term2)

        relationshipsLister <- list()

        relationshipsLister[["across"]] <- relationships
        relationshipsLister[[dataset1]] <- relationshipsFromAll1
        relationshipsLister[[dataset2]] <- relationshipsFromAll2
        relationshipsLister[["Relacross"]] <- relationshipsFromAllAcross
        relationshipsLister[[paste0("random", dataset1)]] <- randomRelationships1
        relationshipsLister[[paste0("random", dataset2)]] <- randomRelationships2

    
        for(distanceMetric in names(distanceMetricList)){
                if(distanceMetric == "jsd" & NORMFUNC[[3]] != "noNorm") {next}       
                cat(file="internalStarted.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " DistanceMetric: ", distanceMetric,"\n"))	
                res <- compareDistanceMetric(as.matrix(netPropDataFrame),
                        computeDistance,
                        distanceMetricList[[distanceMetric]],	
                        relationshipsLister,
                        8,
                        TRUE,
                        diseaseDF,
                        FALSE) 
                        save(res, file = paste0(netPropPath,
                            "/results/netpropDistanceResRAREcommon_",
                            dataset,"RelatedSet_",
                            NORMFUNC[[3]],"_",
                            distanceMetric,".rdata"))
                
                cat(file="internalFinished.txt",append = TRUE,paste0("dataset: ", dataset, " NormFunc: ", NORMFUNC[[3]], " DistanceMetric: ", distanceMetric,"\n"))	
            
        }
    }
}



