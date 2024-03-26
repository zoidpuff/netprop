netPropDataFrame <- runNetProp(network = intGraph,
                    assocData = assocDataOverallEFO,
                    cutoff = c("value" = 0.5, "number" = 2),
                    binarize = TRUE,
                    damping = 0.85,
                    NormFunc = scaleNormalize,
                    settingsForNormFunc = list())

print(paste0("Length of netpropDataFrame: ", nrow(netPropDataFrame)))
#netPropDataFrame[min(100,nrow(netPropDataFrame)),]

# check the row wise variance and mean of the netpropDataFrame
print(paste0("Mean: ",mean(rowMeans(netPropDataFrame))))
print(paste0("Variance: ",var(rowMeans(netPropDataFrame))))

relationshipsAll <- read.csv(paste0(netPropPath,"/relationshipsWithNames.csv"), stringsAsFactors = FALSE)
relationships <- relationshipsAll %>% filter(term1 %in% rownames(netPropDataFrame) & term2 %in% rownames(netPropDataFrame))
relationships <- as.matrix(relationships[,c("term1","term2")])


res <- compareDistanceMetric(as.matrix(netPropDataFrame),
        computeDistance,
        list("method" = "manhattan","returnDist" = NA),	
        relationships,
        8,
        TRUE,
        diseaseDF)


densityPlotManhattan <- generatePlotsFromDistCompareResults(res,diseaseMapping = NULL, densityOnly = TRUE) 
densityPlotManhattan[[1]] <- densityPlotManhattan[[1]] + ggtitle("")

res <- compareDistanceMetric(as.matrix(netPropDataFrame),
        computeDistance,
        list("method" = "euclidean","returnDist" = NA),	
        relationships,
        8,
        TRUE,
        diseaseDF)


densityPloteuclidean <- generatePlotsFromDistCompareResults(res,diseaseMapping = NULL, densityOnly = TRUE) 
densityPloteuclidean[[1]] <- densityPloteuclidean[[1]] + ggtitle("")


res <- compareDistanceMetric(as.matrix(netPropDataFrame),
        computeDistance,
        list("method" = "spearman","returnDist" = NA),	
        relationships,
        8,
        TRUE,
        diseaseDF)


densityPlotspearman <- generatePlotsFromDistCompareResults(res,diseaseMapping = NULL, densityOnly = TRUE) 
densityPlotspearman[[1]] <- densityPlotspearman[[1]] + ggtitle("")


cowplot::plot_grid(densityPlotManhattan[[1]],densityPloteuclidean[[1]],densityPlotspearman[[1]],ncol = 3,labels = c("Manhattan","Euclidean","Spearman"))