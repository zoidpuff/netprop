
library(igraph)
library(dplyr)

netpropPath <- "/home/gummi/netprop"

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
# Run auroc test for each trait for each normalization method
for(BINARIZE in c(TRUE,FALSE)) {
    for(NORMFUNC in list(NULL,ECnormalize)) {
        aurocs <- list()
        print(paste0("Binarize: ", BINARIZE, " NormFunc: ", ifelse(is.null(NORMFUNC),"NULL",deparse(substitute(NORMFUNC))) ))
        for(trait in unique(assocDataFilt10$diseaseId)) {

            assocDataFiltTemp <- assocDataFilt10 %>% filter(diseaseId == trait)

            seedsInd <- assocDataFiltTemp$score >= 0.5

            seedList <- list(assocDataFiltTemp$targetId, assocDataFiltTemp$score, seedsInd)

            diseaseName <- diseaseMappings[which(diseaseMappings$id == trait),4] 

            aurocs[[diseaseName]] <- avgAUROC(intGraph, seedList, 100, c(0.1, 0.25, 0.5, 0.75), binarize = BINARIZE, NormFunc = NORMFUNC, settingsForNormFunc = NULL)

            #print(paste("Finished trait ", diseaseName, " with ", sum(seedsInd), sep = ""))
        }
        print(paste0("Finished binarize: ", BINARIZE, " NormFunc: ", ifelse(is.null(NORMFUNC),"NULL",deparse(substitute(NORMFUNC))) ))
        binned <- ifelse(BINARIZE, "_binarized", "_weighted")
        normed <- ifelse(is.null(NORMFUNC), "_unnormalized", "_ECnormalized")

        do.call(rbind,aurocs) %>% as.data.frame() %>% write.csv(paste0("aurocResults",binned,normed,".csv"))
    }
}

settings1 <- list("nSamples" = 100, "perserveDegree" = TRUE, "degreeSampleSmoothing" = 4, "minBucketSize" = 2)
settings2 <- list("nSamples" = 100, "perserveDegree" = FALSE, "degreeSampleSmoothing" = 4, "minBucketSize" = 2)

for(SETTINGS in list(settings1,settings2)) {
    for(trait in unique(assocDataFilt10$diseaseId)) {

        assocDataFiltTemp <- assocDataFilt10 %>% filter(diseaseId == trait)

        seedsInd <- assocDataFiltTemp$score >= 0.5

        seedList <- list(assocDataFiltTemp$targetId, assocDataFiltTemp$score, seedsInd)

        diseaseName <- diseaseMappings[which(diseaseMappings$id == trait),4] 

        aurocs[[diseaseName]] <- avgAUROC(intGraph, seedList, 20, c(0.1, 0.25, 0.5, 0.75), binarize = TRUE, NormFunc = permuteTestNormalize, settingsForNormFunc = SETTINGS)

        #print(paste("Finished trait ", diseaseName, " with ", sum(seedsInd), sep = ""))
    }
    do.call(rbind,aurocs) %>% as.data.frame() %>% write.csv(paste0("aurocResults_permuteNorm_",ifelse(SETTINGS$perserveDegree,"_degreePreserve",""),".csv"))
}