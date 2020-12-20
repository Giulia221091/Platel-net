###libraries
library(shiny)
library(openxlsx)
library(visNetwork)


#### build regulatory networks
setwd("/Users/home/Desktop/PlatelNet/")
##upload data on correlations
miRNA.mRNA.proteins.PDAC <- openxlsx::read.xlsx("hyp_canonical.miRNA.PDAC.xlsx",sheet=1)
miRNA.premRNA.proteins.PDAC <- openxlsx::read.xlsx("hyp_canonical.miRNA.PDAC.xlsx",sheet=2)
miRNA.mRNA.proteins.Benign <- openxlsx::read.xlsx("hyp_canonical.miRNA.Benign.xlsx",sheet=1)
miRNA.premRNA.proteins.Benign <- openxlsx::read.xlsx("hyp_canonical.miRNA.Benign.xlsx",sheet=2)

isomiRNA.mRNA.proteins.PDAC <- read.csv("hyp_ISO_mRNA.PDAC.txt",sep = " ")
isomiRNA.premRNA.proteins.PDAC <- read.csv("hyp_ISO_IntronmRNA.PDAC.txt",sep = " ")
isomiRNA.mRNA.proteins.Benign <- read.csv("hyp_ISO_mRNA.Benign.txt",sep = " ")
isomiRNA.premRNA.proteins.Benign <- read.csv("hyp_ISO_IntronmRNA.Benign.txt",sep = " ")

### upload data on differential expression
proteins.UP.benign.DDAcounts <- read.table("proteins.UP.benign.DDAcounts.txt")
proteins.UP.PDAC.DDAcounts <- read.table("proteins.UP.PDAC.DDAcounts.txt")
miRNA.UP.Benign.iso <- read.csv("miRNA.UP.Benign.iso.txt", sep=" ")
miRNA.UP.PDAC.iso. <- read.csv("miRNA.UP.PDAC.iso.txt",sep = " ")
miRNA.UP.Benign.canonical <- read.table("miRNA.UP.Benign.canonical.txt")
miRNA.UP.PDAC.canonical <- read.table("miRNA.UP.PDAC.canonical.txt")
mRNA.intron.UP.PDAC_withGeneNames <- read.csv("intron-RNA.UP.PDAC.txt",sep = " ")
mRNA.intron.UP.Benign_withGeneNames <- read.table("intron-RNA.UP.Benign.txt")
mRNA.UP.PDAC <- read.table("mRNA.UP.PDAC.txt")
mRNA.UP.Benign <- read.table("mRNA.UP.Benign.txt")



# adjust proteins data for network --> put ".p" in the names
proteins.UP.benign.DDAcounts$for.netw <- paste(rownames(proteins.UP.benign.DDAcounts),".p",sep = "")
proteins.UP.PDAC.DDAcounts$for.netw <- paste(rownames(proteins.UP.PDAC.DDAcounts),".p",sep = "")

colnames(isomiRNA.mRNA.proteins.PDAC)[2] <- "gene"
colnames(isomiRNA.mRNA.proteins.Benign)[2] <- "gene"
colnames(isomiRNA.premRNA.proteins.PDAC)[2] <- "gene"
colnames(isomiRNA.premRNA.proteins.Benign)[2] <- "gene"

#### try to build first network


#tmp <- isomiRNA.mRNA.proteins.Benign
#tmp.proteins <- data.frame(Gene=unique(tmp$gene),Protein=paste(unique(as.character(tmp$gene)),".p",sep = ""))
#tmp.edges <- data.frame(node1=c(as.character(tmp$miRNA),as.character(tmp.proteins$Gene)),node2=c(as.character(tmp$gene),as.character(tmp.proteins$Protein)))

server <- function(input, output) {
  
  
  output$network <- renderVisNetwork({
    
    if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.premRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- miRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- miRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- miRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- miRNA.premRNA.proteins.Benign
    }
    
    if (input$"text" != ""){
      tmp.v2.gene <- tmp[grep(input$"text" , tmp$gene,ignore.case = T),]
      tmp.v2.mirna <- tmp[grep(input$"text" , tmp$miRNA,ignore.case = T),]
      tmp <- rbind.data.frame(tmp.v2.gene,tmp.v2.mirna)
    }
    
    if (nrow(tmp) == 0){
      visNetwork(tmp, main = "Ooops...! Bad luck..", submain = "Please try with another node")
    } else {
      
      tmp.proteins <- data.frame(Gene=unique(tmp$gene),Protein=paste(unique(as.character(tmp$gene)),".p",sep = ""))
      tmp.edges <- data.frame(node1=c(as.character(tmp$miRNA),as.character(tmp.proteins$Gene)),node2=c(as.character(tmp$gene),as.character(tmp.proteins$Protein)))
      
      
      # minimal example
      nodes <- data.frame(id=c(unique(as.character(tmp$miRNA)),unique(as.character(tmp$gene)),as.character(tmp.proteins$Protein)))
      nodes$id <- as.character(nodes$id)
      nodes$"color" <- "lightgrey"
      nodes$"shape" <- "circle"
      
      
     
      nodes[which(nodes$id %in% tmp$miRNA),]$shape <- "triangle"
      nodes[which(nodes$id %in% tmp$gene),]$shape <- "rectangle"
      
      
      
      try(nodes[which(nodes$id %in% proteins.UP.benign.DDAcounts$for.netw),]$color <- "lightblue",silent = T)
      try(nodes[which(nodes$id %in% proteins.UP.PDAC.DDAcounts$for.netw),]$color <- "red", silent = T)
      if (input$"mirna_type" == "iso"){
        try(nodes[which(nodes$id %in% as.character(miRNA.UP.Benign.iso$mir)),]$color <- "lightblue", silent = T)
        try(nodes[which(nodes$id %in% as.character(miRNA.UP.PDAC.iso.$mir)),]$color <- "red", silent = T)
        
      }else{
        try(nodes[which(nodes$id %in% rownames(miRNA.UP.Benign.canonical)),]$color <- "lightblue", silent = T)
        try(nodes[which(nodes$id %in% rownames(miRNA.UP.PDAC.canonical)),]$color <- "red", silent = T)
        
      }
      if (input$"rna_type" == "premrna"){
        try(nodes[which(nodes$id %in% mRNA.intron.UP.PDAC_withGeneNames$hgnc_symbol),]$color <- "red", silent = T)
        try(nodes[which(nodes$id %in% mRNA.intron.UP.Benign_withGeneNames$hgnc_symbol),]$color <- "lightblue", silent = T)
        
        
      }else{
        try(nodes[which(nodes$id %in% rownames(mRNA.UP.PDAC)),]$color <- "red", silent = T)
        try(nodes[which(nodes$id %in% rownames(mRNA.UP.Benign)),]$color <- "lightblue", silent = T)
        
      }
      
      
      edges <- data.frame(from = as.character(tmp.edges$node1), to = as.character(tmp.edges$node2))
      
      
      
      visNetwork(nodes, edges, width = "100%",hight="400px") %>%
        #visLegend(useGroups = F,addNodes = data.frame(label=c("miRNA","Gene","Protein"),shape=c("triangle","rectangle","circle")),addEdges=data.frame(label="link",color="black")) %>%
        visExport() %>%visInteraction(hover = TRUE) %>%
        visEvents(hoverNode = "function(nodes) {
                  Shiny.onInputChange('current_node_id', nodes);
                  ;}") %>%
        visOptions(highlightNearest = TRUE)
    }})
  
  output$shiny_return <- renderPrint({
    input$current_node_id
  })
  
  output$table_mirna <- renderTable({
    if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.premRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- miRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- miRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- miRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- miRNA.premRNA.proteins.Benign
    }
    
    if (input$"text" != ""){
      tmp.v2.gene <- tmp[grep(input$"text" , tmp$gene,ignore.case = T),]
      tmp.v2.mirna <- tmp[grep(input$"text" , tmp$miRNA,ignore.case = T),]
      tmp <- rbind.data.frame(tmp.v2.gene,tmp.v2.mirna)
    }
    
    if (nrow(tmp) == 0){
      visNetwork(tmp, main = "Ooops...! Bad luck..", submain = "Please try with another node")
    } else {
      
      tmp.proteins <- data.frame(Gene=unique(tmp$gene),Protein=paste(unique(as.character(tmp$gene)),".p",sep = ""))
      tmp.edges <- data.frame(node1=c(as.character(tmp$miRNA),as.character(tmp.proteins$Gene)),node2=c(as.character(tmp$gene),as.character(tmp.proteins$Protein)))
      miRNA.tab <- as.data.frame(table(as.character(tmp$miRNA)))
      miRNA.tab.sort <- miRNA.tab[order(-miRNA.tab$Freq),]
      
      colnames(miRNA.tab.sort) <- c("miRNA","Node Degree")
      head(miRNA.tab.sort,5)}})
  
  output$table_gene <- renderTable({
    if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- isomiRNA.premRNA.proteins.Benign
    } else if (input$"mirna_type" == "iso" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- isomiRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "pdac"){
      tmp <- miRNA.mRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "mrna" & input$"groups" == "benign"){
      tmp <- miRNA.mRNA.proteins.Benign
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "pdac"){
      tmp <- miRNA.premRNA.proteins.PDAC
    } else if (input$"mirna_type" == "mirna" & input$"rna_type" == "premrna" & input$"groups" == "benign"){
      tmp <- miRNA.premRNA.proteins.Benign
    }
    
    if (input$"text" != ""){
      tmp.v2.gene <- tmp[grep(input$"text" , tmp$gene,ignore.case = T),]
      tmp.v2.mirna <- tmp[grep(input$"text" , tmp$miRNA,ignore.case = T),]
      tmp <- rbind.data.frame(tmp.v2.gene,tmp.v2.mirna)
    }
    
    if (nrow(tmp) == 0){
      visNetwork(tmp, main = "Ooops...! Bad luck..", submain = "Please try with another node")
    } else {
      
      tmp.proteins <- data.frame(Gene=unique(tmp$gene),Protein=paste(unique(as.character(tmp$gene)),".p",sep = ""))
      tmp.edges <- data.frame(node1=c(as.character(tmp$miRNA),as.character(tmp.proteins$Gene)),node2=c(as.character(tmp$gene),as.character(tmp.proteins$Protein)))
      
      genes.tab <- as.data.frame(table(as.character(tmp$gene)))
      genes.tab.sort <- genes.tab[order(-genes.tab$Freq),]
      colnames(genes.tab.sort) <- c("Gene","Node Degree")
      head(genes.tab.sort,5)}})
  
  output$hpa_link <- renderUI({
    prot_to_look <- as.character(input$current_node_id)
    if (length(grep("\\.p", prot_to_look)) >0){ prot_to_look <- strsplit(prot_to_look, "\\.p")[[1]] }
    link <- paste0("https://www.proteinatlas.org/search/", prot_to_look)
    url <- a("Click here!", href = link)
    tagList("Human Protein Atlas link:", url)
  })
  
  output$genecards_link <- renderUI({
    prot_to_look <- as.character(input$current_node_id)
    if (length(grep("\\.p", prot_to_look)) >0){ prot_to_look <- strsplit(prot_to_look, "\\.p")[[1]] }
    link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", prot_to_look)
    url <- a("Click here!", href = link)
    tagList("GeneCards link:", url)
    
  })
  
}


ui_prova <- fluidPage(
  # title of the app
  titlePanel(title = "Regulatory Network", windowTitle = "Regulatory Network"),
  
  fluidRow(
      column(2,
             radioButtons("mirna_type", h3("miRNA type:"),
                          choices = list("miRNA" = "mirna", "IsomiRs" = "iso")),
             radioButtons("rna_type", h3("RNA type:"),
                          choices = list("mRNA" = "mrna", "pre-mRNA" = "premrna")),
             radioButtons("groups", h3("Groups:"),
                          choices = list("PDAC" = "pdac", "benign lesions" = "benign"))),
    
    column(8,
           
           visNetworkOutput("network"),
             tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
             textOutput("shiny_return"),
             uiOutput(outputId = "hpa_link"),
           uiOutput(outputId = "genecards_link")),
    
    
    column(2,
           textInput("text", h3("Search"), 
                     value = ""),
           div(tableOutput(outputId = 'table_mirna'), style="font-size:80%"),
           div(tableOutput(outputId = 'table_gene'), style="font-size:80%")
           ))
    
  )


shinyApp(ui = ui_prova, server = server)





