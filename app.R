if(!require(shiny)){
	install.packages("shiny")
}
library(shiny) 

### Change source file path to working directory

FN <- paste(getwd(),"/R_Functions/GS_Pipeline_Jan_2022_FnsApp.R",sep="")
source(FN)

PN <-  paste(getwd(),"/GSPipeline.png",sep="")


## Set the file upload limit to X MB
options(shiny.maxRequestSize=500*1024^2)

ui <- fluidPage( 
  
  fluidRow(
    column(4),
    column(6, div(
      id = "app-title",
      titlePanel(tags$strong("SOYGEN2 Genomic Selection Demo App")),
      tags$p("A simpler version of SOYGEN2 app to perform genomic predictions given genotypic and phenotypic data for workshop")
    ))),
  
  tabsetPanel(id="inData",
              
    ## Tab for loading data
             
       tabPanel("Load Data",
                       #sidebarLayout(
                         
                         #sidebarPanel(
                       fluidRow(
                          column(4),column(width=4,tags$h3(tags$strong("Load Data Files"))),   
                          column(10),column(width=2,actionButton("Data_Home", "prev"),
                          actionButton("Data_Trait", "next"))), 
                      
                     #  fluidRow(
                     #     column(2),column(width=8,tags$h4("Load data files containing phenotypic data, genotypic data and information on target lines."))),
                       tags$br(),
                       fluidRow(
                         column(1),column(width=3,tags$h4(tags$strong("Phenotypic Data"))),
                         column(5,offset=2,tags$h4(tags$strong("Genotypic Data"))) 
                       ),
                       
                      ####  
                     
                     tags$br(),
                      fluidRow(
                         column(1),column(width=5,fileInput("infileBLUEs", "Choose Phenotype File (CSV)", accept = ".csv")),
                         column(6),column(width=5,fileInput("infileVCF", "Choose Genotype File (VCF)", accept = ".vcf"))
                      ),
                      fluidRow(
                        column(1),column(width=5,checkboxInput("header", "Header", TRUE)),
                        column(6),column(width=5,checkboxInput("header", "Header", TRUE))
                      ),
                     fluidRow(
                       column(1),column(width=5,tags$ul(tags$li("A '.csv' file containing phenotypic information. Phenotype file should have the line ID in the first column and phenotypic trait values after that with the trait ID as column name."))),
                       column(width=5,tags$ul(tags$li("A '.vcf' containing the genotypic information of both the candidate lines used to train the GP model and the target lines whose values need to be predicted.")))
                     ),
                     tags$br(),
                      fluidRow(
                           column(1),column(width=5,tags$h4(tags$strong(textOutput("PhenoHeader")))),
                           column(6),column(width=5,tags$h4(tags$strong(textOutput("GenoHeader"))))
                      ),
                      fluidRow(
                          column(1),column(width=5,tags$h6(tableOutput("PhenoTable"))),
                          column(6),column(width=5,tags$h6(tableOutput("GenoTable")))),
                      tags$br()
               ),  
       
     
## Tab for Trait Selection 
              
              tabPanel("Set Target & Trait",
                       ####  
                       tags$br(),
                       fluidRow(
                         column(2),column(width=4,tags$h4(tags$strong("Set Target Population"))),
                         column(6),column(width=4,tags$h4(tags$strong("Select Trait"))),
                         actionButton("Trait_Data", "prev"),
                         actionButton("Trait_TS", "next")
                       ),
                       tags$br(),
                       fluidRow(
                       column(1),column(width=5,fileInput("infileTargetTable", "Choose Target Data Information File (CSV)", accept = ".csv")),
                       column(5),column(width=5,selectInput(inputId="trait","Choose One or More Traits",choices=NULL,multiple=TRUE))),
                       fluidRow( column(1),column(width=5,checkboxInput("header", "Header", TRUE))),
                       tags$br(),
                       fluidRow(
                         column(1),column(width=5,tags$ul(tags$li(" A '.csv' file containing information on the target data, whose genetic values are to be predicted using the genomic prediction model. 
                       The first column in the target data information file should have the line IDs."))), 
                         column(6),column(width=5,tags$ul(tags$li("Select one trait or multiple traits to train the genomic prediction model.")))),
                       tags$br(), 
                       fluidRow(
                        column(1),column(width=5,div(
                                    tags$h4(tags$strong(textOutput("TargetHeader"))),
                                    tags$h6(tableOutput("TargetTable")))),
                         column(5),column(width=5,div( 
                           tags$h4(tags$strong(textOutput("summaryHeader"))),
                           tags$h6(tableOutput("Summary"))))
                       ),
                       tags$br()
                ),
              

              ### GP Tab 
              tabPanel("Genomic Prediction",
                       sidebarLayout(
                         sidebarPanel(
                           selectInput(inputId="TrainSet","Choose Training Set",c("Complete Input Genotype Set","Optimal Train Set from Step 3","Random Train Set from Step 3")),
                           tags$br(),
                           tags$br(),
                           selectInput(inputId="GPModelST","Choose Prediction Model for Single Trait",c("rrBLUP (rrBLUP)","rrBLUP (bWGR)","BayesB (bWGR)","BayesLASSO (bWGR)")),
                           actionButton("RunPredictionsST", "Predict Single Trait!"),
                           tags$br(),
						               tags$br(),
						               tags$br(),
						               tags$br(),
                           selectInput(inputId="GPModelMT","Choose Prediction Model for Multiple Traits",c("BRR (BGLR)","RKHS (BGLR)","Spike-Slab(BGLR)","Mmer (Sommer)")),
                           actionButton("RunPredictionsMT", "Predict Multiple Traits!"),
                           tags$br(),
						               tags$br()
						   
                         ), mainPanel(
                           fluidRow(
                             column(width=8,tags$h3(tags$strong("Train Genomic Prediction Model"))), actionButton("GP_CVR", "prev"),
                             actionButton("GP_VP", "next")
                           ),
                           fluidRow(
                            
                             column(width=10, tags$p(" Select the statistical method to train the genomic prediction model and predict the values of target lines. rrBLUP method is implemented using the
                                           rrBLUP (Endelman 2011) package. Expectation maximization based RR-BLUP, BayesB and BayesLASSO methods are implemented using the bWGR package (Xavier et al. 2020).
                                           Multi-trait predictions are implemented using the BGLR and Sommer packages. The Multitrait function in BGLR implements Bayesian Ridge Regression, RKHS, and Spike-Slab methods. 
                                           The multi-trait GBLUP is implemented using the 'mmer' function in 'sommer' package."))),
                           tags$br(),
                           tags$br(),
                           tags$h4(textOutput("RankedLinesHeader")),
                           tableOutput("Ranked_Lines_for_SelectionST"),
                           tableOutput("Ranked_Lines_for_SelectionMT")
                          
                         )
                       )
              ),
              #   
              #   ## Tab for Results
              tabPanel("Visualize Predictions",
                       sidebarLayout(
                         sidebarPanel(
                          downloadButton("ExportOut", "Export Output Table")
                          ),
                         mainPanel(
                           
                           fluidRow(
                             column(3),column(width=8,tags$h3(tags$strong("Visualize and Explore Predictions"))), actionButton("VP_GP", "prev")
                           ),
                           tags$br(),
                           # tags$h4(textOutput("RankedLinesHeader")),
                           # tableOutput("Ranked_Lines_for_SelectionST"),
                           # tableOutput("Ranked_Lines_for_SelectionMT"),
                           #tags$h4("Predicted vs Observed Values"),
                           #plotOutput("plots")
                           uiOutput("scatter"),
                           uiOutput("scatterMT")
                         )   
                       ) 
              )
  )
)


###

server <- function(input,output){
  
  #Geno
  
  Geno <- reactive({
    
    genoFile <- input$infileVCF
    ext <- tools::file_ext(genoFile$datapath)
    req(genoFile)
    validate(need(ext == "vcf", "Please upload a vcf file"))
    
   # withProgress(message = 'Reading Data', value = 0, {
   #   NUST_Genotypes_VCF <- read.table(genoFile$datapath)})
    
    withProgress(message = 'Converting VCF to Dataframe', value = 0, {
      gt2d <- VCFtoDF(genoFile$datapath) }) 
    
     gt2d
  })
  
  #Pheno
  Pheno <- reactive({
    
    phenoFile <- input$infileBLUEs
    
    ext <- tools::file_ext(phenoFile$datapath)
    req(phenoFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(phenoFile$datapath, header = input$header)
    
  })
  
  observeEvent(input$infileBLUEs, {
    updateSelectInput(inputId = "trait",choices = colnames(Pheno())[2:ncol(Pheno())])})
  
  
  #Target
  TargetTab <- reactive({
    TargetFile <- input$infileTargetTable
    
    ext <- tools::file_ext(TargetFile$datapath)
    req(TargetFile)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    read.csv(TargetFile$datapath, header = input$header)
  })
  
 
 
  #updateNumericInput(inputId = "noCandidates",value=(nrow(Geno())-nrow(TargetTab())),min =1, max=(nrow(Geno())-nrow(TargetTab())))})
  # observeEvent(input$infileTargetTable, {
  #   updateNumericInput(inputId = "noToSelect",value=(nrow(Geno())-nrow(TargetTab())),min =1, max=(nrow(Geno())-nrow(TargetTab())))})
  # 
  
  phenoHead <- eventReactive(input$infileBLUEs,{ paste("Phenotype Table with ",nrow(Pheno())," lines and ",ncol(Pheno())-1," traits",sep="")})
  output$PhenoHeader <- renderText({phenoHead()})
  output$PhenoTable <- renderTable({as.data.frame((Pheno())[1:5,1:5])})
  
  genoHead <- eventReactive(input$infileVCF,{ paste("Genotype Table with ",ncol(Geno())-5," lines and ",nrow((Geno()))," markers",sep="")})
  output$GenoHeader <- renderText({genoHead()})
  output$GenoTable <- renderTable({as.data.frame((Geno())[1:5,1:5])})
  
  
  TargetHead <- eventReactive(input$infileTargetTable,{ paste("Table with information on ",nrow(TargetTab())," Target lines",sep="")})
  output$TargetHeader <- renderText({TargetHead()})
  output$TargetTable <-  renderTable({as.data.frame(TargetTab()[1:5,c(1,2)])})
  
  
  
  # 
  # 
  # eventReactive(list(input$infileVCF,input$infileBLUEs),{
  #   if(is.null(input$infileTargetTable)){
  #     Pheno()[,1] 
  #     
  #   }
  # })
  
  nTraits <- eventReactive(input$infileBLUEs,{(ncol(Pheno())-1)})
  TargetIDs <- eventReactive(input$infileTargetTable,{as.character(TargetTab()[,1])})
  
## Trait 
  
  Trait <- reactive(input$trait)
  
  summaryHead <- eventReactive(input$trait,{paste("Summary of distribution of selected phenotypic values")})
  output$summaryHeader <- renderText({summaryHead()})
  nSelTraits <- reactive(length(Trait()))
  
  SummaryTxt <- eventReactive(input$trait,{ 
    do.call(rbind,lapply(Trait(),function(x) round(summary(Pheno()[,x]),digits=2)))
  })
  
  observeEvent(input$trait,{ 
  output$Summary <- renderTable({
   
     trTable <- cbind(unlist(Trait()),nSelTraits(),SummaryTxt())
    
    print(trTable)
    })
  })
 
  
  
## GP  Models 
  
  GPModelST <- reactive(input$GPModelST)
  outputList <- eventReactive(input$RunPredictionsST,{ 
    if(nSelTraits()==1){
      outputDF <-  withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait()),GPModelST(),optTS())})
      return(outputDF)
    }
    if(nSelTraits()>1){
      outputDFComb <- c()
      for(nSelT in 1:nSelTraits()){
        i <- reactive(nSelT)
       outputDF <-  withProgress(message = 'Running Computations', value = 0, {
        getRankedPredictedValues(predictionData(),nTraits(),unlist(Trait())[i()],GPModelST(),optTS())})
       
       outputDFComb <- cbind(outputDFComb,outputDF[,2])
      } 
      outputDFComb <- cbind(outputDF[,1],outputDFComb,outputDF[,3])
      colnames(outputDFComb) <- c("GermplasmID",unlist(Trait()),"UpperBound of Reliability")
      return(outputDFComb)
    }
    
 })
  
  
### GPModel MT
  
  GPModelMT <- reactive(input$GPModelMT)
# 
  
  
  outputListMT <- eventReactive(input$RunPredictionsMT,{ 
    
    withProgress(message = 'Running Computations', value = 0, {
      getRankedPredictedValuesMT(predictionData(),nTraits(),unlist(Trait()),GPModelMT(),optTS())
    })
  })
  

  
  
 
  
## Render Table for ST  
 
  output$Ranked_Lines_for_SelectionST <- renderTable({
     
     
      as.data.frame(outputList()[1:20,])
     
  })

## Render Table for MT
  
  output$Ranked_Lines_for_SelectionMT <- renderTable({ if(nSelTraits()>1){
    
      (outputListMT()[1:20,])
   }
  }) 
  #%>% bindCache(outputListMT(),processedData(),nTraits(),unlist(Trait()),GPModelMT()) %>% bindEvent()
  
  
  output$scatter <- renderUI({ 
    if(nSelTraits()==1){
      plotOutput("plots")
    }
  }) 
  output$scatterMT <- renderUI({  
    if(nSelTraits()>1){
      plotOut_List <- lapply(1:nSelTraits(),function(x) { 
      i <- reactive(x)
      plotOutput(noquote(paste("'","plot",i(),"'",sep="")))})
      do.call(tagList,plotOut_List)
    }
  })

      
    
 
    output$plots <- renderPlot({
      if(nSelTraits()==1){
            plot.default(outputList()[,2],outputList()[,3],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values",font.lab=2,cex.lab=1.25)
      } 
    })
    
  
    
  observeEvent(input$RunPredictionsST,{
    if(nSelTraits()>1){
      lapply(1:nSelTraits(),function(i){
        loc_i <- reactive(i)
        output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
          plot.default(outputList()[,(nSelTraits()+1)],outputList()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper Bound of Reliability vs Predicted Values",font.lab=2,cex.lab=1.25)
        })
        
      })
    }
  })  
  
  # }
  #   
  #  
  # local({
   
   #printPlots2 <- reactive({  
    observeEvent(input$RunPredictionsMT,{ 
    
        lapply(1:nSelTraits(),function(i){
               loc_i <- reactive(i)
               output[[noquote(paste("'","plot",loc_i(),"'",sep=""))]] <- renderPlot({
                   plot.default(outputListMT()[,(nSelTraits()+1)],outputListMT()[,(nSelTraits()+2)],type="p",xlab="Predicted Value",ylab="Upper Bound of Reliability",main="Upper bound of Reliability vs Predicted Values",font.lab=2,cex.lab=1.25)
                })
         
        })
    
   }) 
       
          
    
  switch_page <- function(i) {
    sel <- reactive(i)
    PageID <- c("Home","Load Data","Set Target & Trait","Optimize Training Population",
                "Cross Validations","Genomic Prediction","Visualize Predictions")
    updateTabsetPanel(inputId = "inData", selected = PageID[sel()])
  }
  
  
  observeEvent(input$Home_Data,{ switch_page(2)})
  observeEvent(input$Data_Home, switch_page(1))
  
  observeEvent(input$Data_Trait, switch_page(3))
  observeEvent(input$Trait_Data, switch_page(2))
  
  observeEvent(input$Trait_GP, switch_page(4))
  observeEvent(input$GP_Trait, switch_page(3))
 
  observeEvent(input$GP_VP, switch_page(5))
  observeEvent(input$VP_GP, switch_page(4))
     
  
  output$ExportOut <- downloadHandler(
      filename = function() {
        "Predictions.csv"
      },
      content = function(file) {
         if(nSelTraits()==1){
          write.csv(as.data.frame(outputList()), file, row.names = FALSE)
         }
         if(nSelTraits()>1){
           write.csv(as.data.frame(outputListMT()), file, row.names = FALSE)
         }
      }
  )

}

shinyApp(server=server,ui=ui)




