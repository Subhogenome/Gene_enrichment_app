#load packages
library(medulloPackage)
library(BiocManager)
options(repos = BiocManager::repositories())
library(bs4Dash)
library(locfit)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(btools)
library(heatmaply)
library(tidyverse)
library(data.table)
library(plotly)
library(ggfortify)
library("shinyWidgets")
library(dplyr)
library(colourpicker)
library(shinysky)
library(rhandsontable)
library(bslib)
library(fresh)
#custom themes
themes <- list("Light" = theme_light(),
               "Minimal" = theme_minimal(),"Gray" = theme_gray() , "BLack and white" = theme_bw() , "Line Draw" = theme_linedraw() , "Dark" = theme_dark() , "Classic" = theme_classic() , "Void" = theme_void(), "test" = theme_test())
#use bootstrap
library(bslib)
theme <- bs_theme(
  bg = "#0b3d91", fg = "white", primary = "#FCC780",
  base_font = font_google("Space Mono"),
  code_font = font_google("Space Mono")
)

mytheme <- create_theme(
  adminlte_color(
    light_blue = "#E7FF6E",
    black = "#FFFFFF"
  )
)
#ui component
ui <- dashboardPage( 
  dashboardHeader(title = "Visual Studio" ),
  dashboardSidebar( sidebarMenu((menuItem("Upload Data", tabName = "Upload_data",fileInput("FileInput", "Upload Count Table" , multiple = TRUE),
                                          fileInput("MetaInput", "Upload Meta data", multiple = TRUE),
                                          actionButton("go", "Go"))))),
  dashboardBody(fluidPage(use_theme(mytheme),tabItem(tabName = "Upload_data",
                        tabBox(title = "Results",
                               id = "tabset" , height= "1850px" , width = "13",
                               tabPanel("Uploaded Data" , 
                                        fluidPage(h1(textOutput("text1")),
                                                  DT::dataTableOutput("table"),style ="overflow-x: scroll"),
                                        h1(textOutput("text2")),
                                        fluidPage(DT::dataTableOutput("Metatable"))),
                               tabPanel("Interactive Plots",
                                        fluidPage(column(6,h4(textOutput("text3")),br(),
                                                         dropdown(
                                                           
                                                           tags$h3("Customize"),
                                                           
                                                           pickerInput(inputId = 'xcol2',
                                                                       label = 'Features',
                                                                       choices = c("treatment","cell_type","sample_id")),
                                                           pickerInput(inputId = 'theme',
                                                                       label = 'Themes',
                                                                       choices = names(themes)),
                                                           sliderInput(inputId = 'size',
                                                                       label = 'Marker Size',
                                                                       value = 2,
                                                                       min = 1, max = 9),
                                                           options = list(`style` = "btn-info")
                                                           ,
                                                           style = "unite", icon = icon("cog"),
                                                           status = "danger", width = "300px",
                                                           animate = animateOptions(
                                                             enter = animations$fading_entrances$fadeInLeftBig,
                                                             exit = animations$fading_exits$fadeOutRightBig
                                                           )
                                                         ),uiOutput("PCA")),
                                                  column(6,h4(textOutput("text4")),br(),
                                                         dropdown(
                                                           
                                                           tags$h3("Customize"),
                                                           pickerInput(inputId = 'theme2',
                                                                       label = 'Themes',
                                                                       choices = names(themes)),
                                                           sliderInput(inputId = 'percent',
                                                                       label = 'Percentage Coverage',
                                                                       value = 5,
                                                                       min = 1, max = 100),
                                                           options = list(`style` = "btn-info")
                                                           ,
                                                           style = "unite", icon = icon("cog"),
                                                           status = "danger", width = "300px",
                                                           animate = animateOptions(
                                                             enter = animations$fading_entrances$fadeInLeftBig,
                                                             exit = animations$fading_exits$fadeOutRightBig
                                                           )
                                                         ), uiOutput("SPCA"))),
                                        HTML('<hr style="color: purple;">'),
                                        fluidPage(column(6 ,h4(textOutput("text5")),uiOutput("DPCA")),
                                                  column(6,h4(textOutput("text6")),br(),
                                                         dropdown(
                                                           
                                                           tags$h3("Customize"),
                                                           colourInput("col1", "High Correlation", "blue"),
                                                           colourInput("col2", "Low Correlation", "white"),
                                                           sliderInput(inputId = 'top',
                                                                       label = 'Top gene according to pAdj',
                                                                       value = 5,
                                                                       min = 1, max = 25),
                                                           
                                                           options = list(`style` = "btn-info")
                                                           ,
                                                           style = "unite", icon = icon("cog"),
                                                           status = "danger", width = "300px",
                                                           animate = animateOptions(
                                                             enter = animations$fading_entrances$fadeInLeftBig,
                                                             exit = animations$fading_exits$fadeOutRightBig
                                                           )),uiOutput("IH"))),
                                        HTML('<hr style="color: purple;">'),
                                        fluidPage(column(6, h3(textOutput("text11"),br(),
                                                               dropdown(
                                                                 
                                                                 tags$h3("Customize"),
                                                                 pickerInput(inputId = 'theme3',
                                                                             label = 'Themes',
                                                                             choices = names(themes)),
                                                                 style = "unite", icon = icon("cog"),
                                                                 status = "danger", width = "300px",
                                                                 animate = animateOptions(
                                                                   enter = animations$fading_entrances$fadeInLeftBig,
                                                                   exit = animations$fading_exits$fadeOutRightBig
                                                                 )
                                                               ),uiOutput("vop"))))),
                               tabPanel("Customised Plots", fluidPage(column(4,
                                                                             
                                                                             uiOutput("toCol1"),
                                                                             uiOutput("toCol2"),
                                                                             uiOutput("toCol3"),
                                                                             uiOutput("toCol4"),
                                                                             actionButton("addrow", "Add Gene"),
                                                                             actionButton("revrow" , "Remove Gene"),
                                                                             actionButton("plots", "Plots")),
                                                                      column(h3("Select or Unselect Genes"), uiOutput("toCol6"), width = 6)),
                                        HTML('<hr style="color: purple;">'),
                                        fluidPage(column(6 ,h4(textOutput("text8")) ,br(),
                                                         dropdown(
                                                           
                                                           tags$h3("Customize"),
                                                           pickerInput(inputId = 'theme3',
                                                                       label = 'Themes',
                                                                       choices = names(themes)),
                                                           sliderInput(inputId = 'size1',
                                                                       label = 'Marker Size',
                                                                       value = 2,
                                                                       min = 1, max = 9),
                                                           
                                                           options = list(`style` = "btn-info")
                                                           ,
                                                           style = "unite", icon = icon("cog"),
                                                           status = "danger", width = "300px",
                                                           animate = animateOptions(
                                                             enter = animations$fading_entrances$fadeInLeftBig,
                                                             exit = animations$fading_exits$fadeOutRightBig
                                                           )
                                                         ),uiOutput("CIH")),
                                                  column(6,h4(textOutput("text9")),br(),
                                                         dropdown(
                                                           
                                                           tags$h3("Customize"),
                                                           colourInput("col4", "High Correlation", "blue"),
                                                           colourInput("col3", "Low Correlation", "white"),
                                                           
                                                           
                                                           options = list(`style` = "btn-info")
                                                           ,
                                                           style = "unite", icon = icon("cog"),
                                                           status = "danger", width = "300px",
                                                           animate = animateOptions(
                                                             enter = animations$fading_entrances$fadeInLeftBig,
                                                             exit = animations$fading_exits$fadeOutRightBig
                                                           )),uiOutput("HIH"))),
                                        fluidPage( column(6,h3("Customised Volcano") , br(),
                                                          dropdown(
                                                            
                                                            tags$h3("Customize"),
                                                            pickerInput(inputId = 'theme4',
                                                                        label = 'Themes',
                                                                        choices = names(themes)),
                                                            style = "unite", icon = icon("cog"),
                                                            status = "danger", width = "300px",
                                                            animate = animateOptions(
                                                              enter = animations$fading_entrances$fadeInLeftBig,
                                                              exit = animations$fading_exits$fadeOutRightBig
                                                            )
                                                          ), uiOutput("CV")))
                                        
                                        
                                        
                               ),
                               tabPanel("FAQs & Help Box")))
  )
  )
)
#server component
server <- function(input, output, session){
  (datasetInput1 <- reactive({
    infile <- input$FileInput
    if(is.null(infile))
      return(NULL)
    read.table(infile$datapath, header = TRUE)
  }))
  (datasetInput2 <- reactive({
    infile <- input$MetaInput
    if(is.null(infile))
      return(NULL)
    read.csv(infile$datapath, header = TRUE)
  }))
  
  (output$text1 <- renderText({
    req(datasetInput1())
    paste("Counts Table")}))
  
  (output$table = DT::renderDataTable(datasetInput1()))
  
  (output$text2 <- renderText({
    req(datasetInput2())
    paste("Meta Table")}))
  
  (output$Metatable = DT::renderDataTable(datasetInput2()))
  
  
  output$Gene = DT::renderDataTable({
    req(datasetInput1())
    counts = datasetInput1()
    counts = select(counts ,gene_name)
    counts
  })
  
  
  
  values <- reactiveValues()
  values$DT <- data.frame(gene_name = "NA", gene_id = "NA"
  )
  newEntry <- observeEvent(input$addrow, {
    newLine <- c(input$textIn , input$textIn1)
    values$DT <- rbind(values$DT, newLine)
  })
  newEntry <- observeEvent(input$revrow, {
    deleteLine <- values$DT[-nrow(values$DT), ]
    values$DT <- deleteLine
  })
  
  
  
  
  
  
  
  output$gene <- renderTable({
    v= values$DT
    v
  })
  
  
  
  output$toCol6 <- renderUI({
    
    df= values$DT
    df
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    items <- as.character(df[[1]])
    checkboxGroupInput("icons", "Add or drop genes:",
                       items , items
    )
    
  })
  output$toCol8 <- renderUI({
    req(datasetInput1())
    
    df= values$DT
    df
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    items <- as.character(df[[1]])
    selectizeInput("textIn3", "Select Sample IDS", items, selected = "SRR1039508", multiple = FALSE,
                   options = NULL)
    
  })
  
  
  output$toCol1 <- renderUI({
    req(datasetInput1())
    df <- datasetInput1()
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    
    item <- as.character(df[[2]])
    selectizeInput("textIn", "Available Gene", item, selected = NULL, multiple = FALSE,
                   options = NULL)
    
    
    
  })
  
  output$toCol2 <- renderUI({
    req(datasetInput1())
    df <- datasetInput1()
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    
    h6("Use capital letters to specify Gene Names ex- FOX1 ")
    
    
    
  }) 
  output$toCol3 <- renderUI({
    req(datasetInput1())
    df <- datasetInput1()
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    items <- as.character(df[[1]])
    
    selectizeInput("textIn1", "Available Gene IDs", items, selected = "NA", multiple = FALSE,
                   options = NULL)
    
    
    
  }) 
  output$toCol4 <- renderUI({
    req(datasetInput1())
    df <- datasetInput1()
    
    # as.character: to get names of levels and not a numbers as choices in case of factors
    
    h6("Use capital letters to specify Gene id ex- ENSXXXXXXXXX")
    
    
  })
  
  
  
  
  observeEvent(input$plots,{
    output$text8 <- renderText({
      paste("Customized PCA plot")})
    output$CIH <- renderUI({
      req(datasetInput1())
      icons <- paste(input$icons)
      c <- data.frame(gene_name = c(icons))
      v=values$DT
      v
      v = merge(x = v , y = c , by = c("gene_name"))
      v = select(v , -gene_id)
      v
      count = datasetInput1()
      count
      count= merge(x = count , y = v , by = c("gene_name"))
      count
      counts<-select(count , -gene_id , -gene_name)
      counts
      row.names(counts) <- count$gene_name
      row.names(counts) 
      df <- normalize(counts)
      df
      bad <- sapply(df, function(x) all(is.nan(x)))
      bad
      df = df[,!bad]
      df
      pca_res <- prcomp(df , scale = TRUE)
      pca_res
      counts$gene_name <- rownames(counts)
      counts
      plot_theme <- reactive({themes[[input$theme3]]})
      p <- autoplot(pca_res , data = counts ,colour = "gene_name" , size= input$size1)
      p
      
      ggplotly(p + plot_theme())
    })
    
    output$text9 <- renderText({
      paste("Customized HeatMap")})
    
    output$HIH <- renderUI({
      req(datasetInput1(), datasetInput2())
      icons <- paste(input$icons)
      c <- data.frame(gene_name = c(icons))
      v=values$DT
      v
      v = merge(x = v , y = c , by = c("gene_name"))
      v = select(v , -gene_id)
      v
      count = datasetInput1()
      count
      meta = datasetInput2()
      
      
      count = count[order(count$gene_name), ]
      count
      counts<-select(count ,-gene_name , -gene_id )
      counts
      counts <- round(counts)
      counts
      
      row.names(counts) <- count$gene_id
      
      counts$gene_id <- rownames(counts) 
      counts
      counts = counts %>% select(gene_id, everything())
      counts
      dd<-DESeqDataSetFromMatrix(counts,meta,design=~treatment,tidy = TRUE)
      dd
      dd<-DESeq(dd)
      dd
      normc<-counts(dd,normalized=TRUE)
      normc
      row.names(normc)<-rownames(normc, do.NULL = TRUE, prefix = "row")
      row.names(normc)
      re<-results(dd,alpha = 0.05)
      re
      vsd <- rlog(dd)
      vsd
      colData(vsd)
      vsd
      
      mat  <- assay(vsd)
      mat
      mat  <- mat - rowMeans(mat)
      mat
      DF2 = as.data.frame(t(mat))
      DF2
      df_t <- transpose(DF2)
      df_t
      rownames(df_t) <- colnames(DF2)
      colnames(df_t) <- rownames(DF2)
      df_t$gene_id <- rownames(df_t)
      df_t
      ncount <- select(count , gene_id , gene_name)
      ncount
      df_t = merge(x = df_t, y = ncount , by = c("gene_id"))
      df_t
      
      df_t = select(df_t , -gene_id)
      df_t
      df_t= merge(x = df_t , y = v , by = c("gene_name"))
      df_t
      row.names(df_t) <- df_t$gene_name
      df_t
      
      df_t = select(df_t , -gene_name)
      df_t
      df_t <- as.matrix(df_t)
      df_t
      heatmaply(
        df_t, colors = c(input$col3, input$col4)
      )
    })
    output$CV <- renderUI({
      icons <- paste(input$icons)
      icon <- c(icons)
      count = datasetInput1()
      count
      meta = datasetInput2()
      
      count
      count = count[order(count$gene_name), ]
      count = count[order(count$gene_name), ]
      counts<-select(count ,-gene_name , -gene_id )
      counts <- round(counts)
      
      
      row.names(counts) <- count$gene_id
      
      counts$gene_id <- rownames(counts) 
      counts = counts %>% select(gene_id, everything())
      dd<-DESeqDataSetFromMatrix(counts,meta,design=~treatment,tidy = TRUE)
      dd<-DESeq(dd)
      normc<-counts(dd,normalized=TRUE)
      row.names(normc)<-rownames(normc, do.NULL = TRUE, prefix = "row")
      re<-results(dd,alpha = 0.05)
      res_order<-re[order(re$padj),]
      res_order<- as.data.frame(res_order)
      res_order[is.na(res_order)] = 0
      res_order$sig<-ifelse(res_order$padj<=0.05,"yes","no")
      
      
      res_order$gene_id <- rownames(res_order) 
      dft = select(count , gene_id , gene_name)
      
      res_order = merge(x = res_order, y = dft , by = c("gene_id"))
      df = res_order[res_order$gene_name %in% icon, ] 
      df = res_order[res_order$gene_name %in% icon, ] 
      df = select(df , -sig)
      df$sig = c("selected")
      res_new= res_order[ !grepl(paste(icon, collapse="|"), res_order$gene_name),]
      
      res_rr= merge(res_new, df, all = TRUE) 
      
      
      View(res_order)
      p <- ggplot(res_rr,aes(x=log2FoldChange, y = -log10(padj),color=sig , text = gene_name))+geom_point()
      
      plot_theme <- reactive({themes[[input$theme3]]})
      
      ggplotly(p + plot_theme())
      
    })
    output$bar <- renderUI({
      req(datasetInput1() , datasetInput2())
      icons <- paste(input$icons)
      c <- data.frame(gene_name = c(icons))
      v=values$DT
      v
      v = merge(x = v , y = c , by = c("gene_name"))
      v = select(v , -gene_id)
      v
      count = datasetInput1()
      count
      count= merge(x = count , y = v , by = c("gene_name"))
      count
      counts= select(count , -gene_id , -gene_name)
      row.names(counts) <- count$gene_name
      row.names(counts)
      df_t <- transpose(counts)
      df_t
      rownames(df_t) <- colnames(counts)
      colnames(df_t) <- rownames(counts)
      df_t$sample_name<- rownames(df_t)                    
      df_t
      meta = datasetInput2()
      p= ggplot(df_t, aes(x = sample_name, y =input$textIn3 , fill = sample_name))+
        geom_col(position = "dodge") 
      ggplotly(p)
      
    })
  })
  
  
  observeEvent(input$go,{
    output$text3 <- renderText({
      paste("Sample wise PCA plot")})
    
    (output$PCA <- renderUI({
      req(datasetInput1() , datasetInput2())
      
      count = datasetInput1()
      meta = datasetInput2()
      counts<-select(count ,-gene_id , -gene_name)
      row.names(counts) <- count$gene_id
      row.names(counts)
      df_t <- transpose(counts)
      df_t
      rownames(df_t) <- colnames(counts)
      colnames(df_t) <- rownames(counts)
      df_t$gene_name<- rownames(df_t)                    
      df_t
      df = df_t %>% select(-gene_name)
      df_n = normalize(df)
      bad <- sapply(df_n, function(x) all(is.nan(x)))
      df = df_n[,!bad]
      pca_res <- prcomp(df, scale = TRUE)
      colnames(df_t)[colnames(df_t) == "gene_name"] <- "sample_id"
      df_t = merge(x = df_t, y = meta , by = c("sample_id"))
      p <- autoplot(pca_res , data = df_t ,colour = input$xcol2 , geom="points", size= input$size
      ) 
      plot_theme <- reactive({themes[[input$theme]]})
      ggplotly(p + plot_theme()  ) 
      
    })
    )
    output$text4 <- renderText({
      
      paste("Gene Wise PCA plot")})
    (output$SPCA <- renderUI({
      req(datasetInput1() , datasetInput2())
      count = datasetInput1()
      c = ncol(count)
      r=nrow(count)
      sr= round(r*(input$percent/100))
      count<- count[1:sr , 1:c]
      counts<-select(count , -gene_id , -gene_name)
      row.names(counts) <- count$gene_name
      row.names(counts) 
      df <- normalize(counts)
      bad <- sapply(df, function(x) all(is.nan(x)))
      df = df[,!bad]
      pca_res <- prcomp(df , scale = TRUE)
      counts$gene_name <- rownames(counts)
      p <- autoplot(pca_res , data = counts ,colour = "gene_name")
      plot_theme <- reactive({themes[[input$theme2]]})
      ggplotly(p + plot_theme())
    }))
    
    output$text5 <- renderText({
      paste("3D PCA")})
    (output$DPCA <- renderUI({
      req(datasetInput1() , datasetInput2())
      count = datasetInput1()
      meta = datasetInput2()
      counts<-select(count ,-gene_name, -gene_id )
      counts <- round(counts)
      row.names(counts) <- count$gene_id
      counts$gene_id <- rownames(counts) 
      counts = counts %>% select(gene_id, everything())
      dd<-DESeqDataSetFromMatrix(counts,meta,design=~treatment,tidy = TRUE)
      dd<-DESeq(dd)
      normc<-counts(dd,normalized=TRUE)
      row.names(normc)<-rownames(normc, do.NULL = TRUE, prefix = "row")
      re<-results(dd,alpha = 0.05)
      vsd <- vst(dd)
      colData(vsd)
      plotPCA3D(vsd, intgroup = 'treatment', ntop = 500,
                returnData = FALSE)
      
    }))
    
    output$text6 <- renderText({
      paste("HeatMap (pAdj value based)")})
    (output$IH <- renderUI({
      req(datasetInput1() , datasetInput2())
      count = datasetInput1()
      meta = datasetInput2()
      count = count[order(count$gene_name), ]
      counts<-select(count ,-gene_name , -gene_id )
      counts <- round(counts)
      
      
      row.names(counts) <- count$gene_id
      
      counts$gene_id <- rownames(counts) 
      counts = counts %>% select(gene_id, everything())
      dd<-DESeqDataSetFromMatrix(counts,meta,design=~treatment,tidy = TRUE)
      dd<-DESeq(dd)
      normc<-counts(dd,normalized=TRUE)
      row.names(normc)<-rownames(normc, do.NULL = TRUE, prefix = "row")
      re<-results(dd,alpha = 0.05)
      vsd <- vst(dd, blind = FALSE)
      colData(vsd)
      topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE),input$top)
      mat  <- assay(vsd)[ topVarGenes, ]
      mat  <- mat - rowMeans(mat)
      DF2 = as.data.frame(t(mat))
      df_t <- transpose(DF2)
      df_t
      rownames(df_t) <- colnames(DF2)
      colnames(df_t) <- rownames(DF2)
      df_t$gene_id <- rownames(df_t)
      ncount <- select(count , gene_id , gene_name)
      df_t = merge(x = df_t, y = ncount , by = c("gene_id"))
      View(df_t)
      df_t = select(df_t , -gene_id)
      row.names(df_t) <- df_t$gene_name
      df_t = select(df_t , -gene_name)
      df_t <- as.matrix(df_t)
      heatmaply(
        df_t,
        colors = c(input$col2, input$col1)
      )
    }))
    
    output$text11 <- renderText({
      paste("Volcano plot")})
    output$vop <- renderUI({
      req(datasetInput1() , datasetInput2())
      count = datasetInput1()
      meta =datasetInput2()
      count = count[order(count$gene_name), ]
      counts<-select(count ,-gene_name , -gene_id )
      counts <- round(counts)
      
      
      row.names(counts) <- count$gene_id
      
      counts$gene_id <- rownames(counts) 
      counts = counts %>% select(gene_id, everything())
      dd<-DESeqDataSetFromMatrix(counts,meta,design=~treatment,tidy = TRUE)
      dd<-DESeq(dd)
      normc<-counts(dd,normalized=TRUE)
      row.names(normc)<-rownames(normc, do.NULL = TRUE, prefix = "row")
      re<-results(dd,alpha = 0.05)
      res_order<-re[order(re$padj),]
      res_order$sig<-ifelse(res_order$padj<=0.05,"yes","no")
      res_order<- as.data.frame(res_order)
      is.data.frame(res_order)
      View(res_order)
      res_order<-na.omit(res_order)
      res_order$gene_id <- rownames(res_order) 
      dft = select(count , gene_id , gene_name)
      res_order = merge(x = res_order, y = dft , by = c("gene_id"))
      p <- ggplot(res_order,aes(x=log2FoldChange, y = -log10(padj),color=sig , text = gene_name))+geom_point()
      
      plot_theme <- reactive({themes[[input$theme3]]})
      
      ggplotly(p + plot_theme())
      
      
    })
  })
  
}

shinyApp(ui = ui, server = server)
