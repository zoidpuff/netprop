
library(dplyr)
library(sparklyr)
library(sparklyr.nested)


pathToDataFull <- setNames(list.files('/home/gummi/netprop/data/newRawData', full.names = TRUE),
                            list.files('/home/gummi/netprop/data/newRawData'))

config <- spark_config()
config$`sparklyr.shell.driver-memory` <- '6G'
config$`sparklyr.shell.executor-memory` <- '6G'
config$`sparklyr.verbose` <- TRUE


for(folder in names(pathToDataFull)[2]){

    sc <- spark_connect(master = "local",
                    log = "console",
                    config = config )

    pathToData <- pathToDataFull[[folder]]

    assocData <- spark_read_parquet(sc,
                                    path = pathToData)

    assocDF <- collect(assocData)
    write.csv(assocDF, paste0("/home/gummi/netprop/data/newData/",folder,".csv"))

    spark_disconnect(sc)
}

