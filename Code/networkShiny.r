require(shiny)
require(visNetwork)
require(igraph)
require(plotly)
netPropPath <- '/home/gummi/netprop'

load(file = paste0(netPropPath,"/results/communityDF_highres.rdata")) #communityDF
load(file = paste0(netPropPath,"/results/intGraphSubContract_highres.rdata")) #intGraphSubContract
load(file = paste0(netPropPath,"/results/colMat.rdata")) #colMat
load(file = paste0(netPropPath,"/data/communitiesHighRes.rdata")) #communities
load(paste0(netPropPath,"/data/diseasesNew.rdata")) #diseaseDF
load(file = "/home/gummi/netprop/umapPlot.rdata") #umapPlot
idToName <- setNames(diseaseDF$name, diseaseDF$id)
nameToId <- setNames(diseaseDF$id, diseaseDF$name)

rm(diseaseDF)

data <- visNetwork::toVisNetworkData(intGraphSubContract)

# Add a size column the corresponds to the number of nodes in the community
data$nodes$size <- log1p(sapply(data$nodes$id, function(x) sum(communities$membership == x)))*5

server <- function(input, output) {
  output$network_proxy_nodes <- renderVisNetwork({
    visNetwork::visNetwork(nodes = data$nodes, edges = data$edges) %>% 
    visNetwork::visIgraphLayout(layout = "layout_nicely") %>%
    visNetwork::visEdges(color=list("opacity" = 0.05)) %>%
    visNetwork::visOptions(highlightNearest = TRUE) 


  })

    output$umapPlot <- renderPlotly({
        # Retrieve the text labels from the umapPlot ggplot object
        # and add them to the plotly object
        # This is a workaround to get the text labels to show up in the plotly object
        umapPlot <- ggplotly(umapPlot, tooltip = c("Cluster","text"))
        umapPlot$x$data[[1]]$text <- umapPlot$x$data[[1]]$text
        # use ggplot object umapPlot
        ggplotly(umapPlot)%>% 
      event_register("plotly_click")
    })

  observe({
    if(!input$hover) {
      selVar <- "plotly_click"
    } else {
      selVar <- "plotly_hover"
    }
    point <- event_data(event = selVar, priority = "event")
    req(point) # to avoid error if no point is clicked
    
    # Get the trait id from the point clicked
    traitId <- stringr::str_split(point$key, " ")[[1]][1]
    dataUpdate <- data
    dataUpdate$nodes$color <- colMat[traitId,]
    visNetwork::visNetworkProxy("network_proxy_nodes") %>%
      visNetwork::visUpdateNodes(dataUpdate$nodes)
  })
  
  observe({
    # Transform the input into just the trait id in the parentheses
    traitId <- gsub(".*\\((.*)\\)", "\\1", input$trait)
    if(!traitId %in% rownames(colMat)) {
      return(NULL)
    }
    dataUpdate <- data
    dataUpdate$nodes$color <- colMat[traitId,]
    visNetwork::visNetworkProxy("network_proxy_nodes") %>%
      visNetwork::visUpdateNodes(dataUpdate$nodes)
    
  })
  
}

ui <- fluidPage(
    titlePanel("Trait Network Visualization"),
  fluidRow(
    column(
      width = 3,
        wellPanel(
            h4("Select a trait to highlight"),
      selectizeInput(inputId ='trait',
                        label = NULL,
                        choices = paste0(idToName[rownames(colMat)], " (",rownames(colMat),")")),
        # Add checkbox to select turn on feature
        checkboxInput("hover", "Select on hover", value = FALSE)
                        ))),
fluidRow(
    column(
      width = 5,
      plotlyOutput("umapPlot", height = "800px")
    ),
    column(
      width = 7,
      visNetworkOutput("network_proxy_nodes", height = "800px")
    )
  ),
    fluidRow(
        column(
        width = 12,
        tableOutput("click")
        )
    )
)

shinyApp(ui = ui, server = server)



# Do separate for gwas and ecd same trait and look at cluster enrichment
# Add same trait, different modality umap