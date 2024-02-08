# benchmark network propagation while leaving out seed genes

setwd('W:/GROUP/Users/Ellen/NetworkPropagation/')

library(igraph)
library(pROC)

variantCertain = read.csv('20240101_variantsEFO_MONDO_withSource_updatedJBTS.csv')

#load open targets interaction network (IntAct, Reactome, SIGNOR, STRING)
intAll <- read.csv('./Datasets/interaction/interactionAll.csv')

#set threshold for STRING
intString = grep('string', intAll$sourceDatabase)
lowString = intString[intAll$scoring[intString] < 0.4]
intHigh = intAll[-lowString,]

intHigh = intHigh[grep('ENSG0', intHigh$targetA),]
intHigh = intHigh[grep('ENSG0', intHigh$targetB),]
intHigh = intHigh[!duplicated(intHigh[,c("targetA", "targetB")]),]

#create graph
intGraph = graph_from_data_frame(intHigh[,c("targetA", "targetB")], directed = F)
intGraphClean = igraph::simplify(intGraph, remove.loops = T, remove.multiple = T, edge.attr.comb = c(weight = 'max', 'ignore'))

#remove big files to save memory
rm(intAll, intGraph, intString, lowString)
gc()

variantCertain = unique(variantCertain[,c("diseaseId", "targetId")])

variantCertainMin10 = variantCertain[variantCertain$diseaseId %in% names(table(variantCertain$diseaseId))[table(variantCertain$diseaseId) > 9],]

diseases = unique(variantCertainMin10$diseaseId)

### benchmark with ChEMBL dataset

benchmark <-
  function(disease = diseases[1], intGraph = intGraphClean) {
    variantList <-
      unique(variantCertainMin10[variantCertainMin10$diseaseId == paste0(disease), 'targetId'])
    testList <-
      sample(variantList, 0.9*length(variantList))
    
    `%notin%` <- Negate(`%in%`)
    trainList = variantList[variantList %notin% testList]
    
    if (length(variantList) < 10) {
      resAUC = NA
    } else {
      V(intGraph)$targetGene = 0
      V(intGraph)$targetGene[V(intGraph)$name %in% trainList] = 1
      
      pageRank = page_rank(intGraph, personalized = V(intGraph)$targetGene)
      
      pageRankRes = data.frame(
        gene = V(intGraph)$name,
        pageRank = pageRank$vector,
        target = V(intGraph)$targetGene
      )
      pageRankRes$test = 0
      pageRankRes[pageRankRes$gene %in% testList, 'test'] = 1
      
      pageRankRes = pageRankRes[pageRankRes$target != 1, ]
      resROC <- roc(pageRankRes$test, pageRankRes$pageRank)
      resAUC <- resROC$auc
    }
    return(resAUC)
  }

aucAll <- sapply(diseases, FUN = function(x){benchmark(disease = x, intGraph = intGraphClean)})

boxplot(aucAll, xlab = 'test', ylab = 'AUC')
