#' Applet to start ddsPLS
#'
#' @param ... No parameter is needed explicitly for this app.
#'
#' @return Nothing due to the nature of the applet.
#'
#' @export
#'
#' @importFrom stats rnorm
#' @import shiny
#'
ddsPLS_App <- function(...) {
  vizu <- c("predict","Q2","criterion", "Q2r", "R2r", "R2", "weightX",
            "weightY","loadingX","loadingY")
  pospos <- c("topleft","topright","bottomright","bottomleft",
              "center","top","bottom",
              "left","right")
  lamCh <- c("lambda.min","lambda.1se")
  get_datas <- function(n=200,p2=100,sd2=3){
    phi <- matrix(rnorm(2*n),nrow = n)
    y <- phi[,1,drop=FALSE] + rnorm(n,sd = 0.3)
    p1_1 <- 50
    p1_2 <- 100
    p1_3 <- 50
    p2_1 <- 100
    x1 <- cbind(matrix(rep(phi[,1,drop=FALSE],p1_1),byrow = F,nrow = n) + rnorm(n*p1_1,sd = 0.4),
                matrix(rep(phi[,1,drop=FALSE]+phi[,2,drop=FALSE],p1_2),byrow = F,nrow = n) + rnorm(n*p1_2,sd = 0.4),
                matrix(rep(phi[,2,drop=FALSE],p1_3),byrow = F,nrow = n) + rnorm(n*p1_3,sd = 0.4))
    x2<- cbind(matrix(rep(phi[,1,drop=FALSE],p2_1),byrow = F,nrow = n) + rnorm(n*p2_1,sd = sd2),
               matrix(rnorm(n*p2,sd=sd2),byrow = F,nrow = n))
    list(x1,x2,y)
  }
  get_datas_NoX1 <- function(n=200,p2=100,sd2=3){
    phi <- matrix(rnorm(2*n),nrow = n)
    y <- phi[,1,drop=FALSE] + rnorm(n,sd = 0.3)
    p1_1 <- 50
    p1_2 <- 100
    p1_3 <- 50
    p2_1 <- 100
    x1 <- cbind(matrix(rep(phi[,1,drop=FALSE],p1_1),byrow = F,nrow = n) + rnorm(n*p1_1,sd = 0.4),
                matrix(rep(phi[,1,drop=FALSE]+phi[,2,drop=FALSE],p1_2),byrow = F,nrow = n) + rnorm(n*p1_2,sd = 0.4),
                matrix(rep(phi[,2,drop=FALSE],p1_3),byrow = F,nrow = n) + rnorm(n*p1_3,sd = 0.4))
    x2<- cbind(matrix(rep(phi[,1,drop=FALSE],p2_1),byrow = F,nrow = n) + rnorm(n*p2_1,sd = sd2),
               matrix(rnorm(n*p2,sd=sd2),byrow = F,nrow = n))
    list(x1[,-c(1:p1_1)],x2,y)
  }
  cols_gps <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
  ui <- fluidPage(
    #=======================================
    navbarPage("ddsPLS (data-driven Sparse PLS)",
               tabPanel("Data analysis",
                        titlePanel("Data analysis"),
                        sidebarLayout(
                          sidebarPanel(
                            h2("Data"),
                            fileInput("fileX", "Choose CSV Files for X",
                                      multiple = TRUE,
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                            fileInput("fileY", "Choose CSV File for Y",
                                      multiple = TRUE,
                                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                            actionButton("startimport","Upload",icon=icon("upload"), inline=T),
                            tags$hr(),
                            h2("Analysis"),
                            h3("Classical analysis"),
                            textInput('lamsAnal', 'Enter a vector of lambdas (comma delimited)', "0.1,0.2"),
                            actionButton("runAnal","Run analysis",icon=icon("play")),
                            tags$hr(),
                            h3("Bootstrap analysis"),
                            numericInput('n_B', 'Number of Bootstrap samples',50,min=50,max=1000,step = 50),
                            numericInput('NCORES', "Number of CPU's to be used",1,min=1,max=15,step = 1),
                            actionButton("runB","Run bootstrap analysis",icon=icon("play")),
                            tags$hr(),
                            tableOutput("summaryShort"),
                            tableOutput("summaryShortAnal"),
                            h2("Test data"),
                            fileInput("fileXTest", "Choose CSV Files for X test",multiple = TRUE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                            actionButton("startimportTest","Upload Test",icon=icon("upload"), inline=T)
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel(
                                title = "The data",
                                h2("Settings upload"),
                                checkboxInput("header", "Header", TRUE),
                                radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ";", inline=T),
                                radioButtons("dec", "Decimal",choices = c(Point = ".",Comma = ","),selected = ",", inline=T),
                                radioButtons("quote", "Quote",choices = c(None = "","Double Quote" = '"',"Single Quote" = "'"),selected = '"', inline=T),
                                tags$hr(),
                                h2("General structure"),
                                textOutput("files_n"),
                                tableOutput("files"),
                                tags$hr(),
                                h2("Block summary"),
                                selectInput('inSelect', 'Choose block', "X 1"),
                                h3("Head"),
                                tableOutput('headerBlock'),
                                h3("Correlation structure"),
                                numericInput("sizeplot2", "Size of Correlation plots (pixels)",600,min=200,max=3000,step = 100),
                                plotOutput('plot1R')
                              ),
                              tabPanel(
                                title = "Model summary",
                                verbatimTextOutput("summary")
                              ),
                              tabPanel(
                                title = "Vizualisations",
                                # numericInput('gamma', 'Fusion coefficient (gamma)',0.0,min=0.0,max=100.0,step = 0.1),
                                selectInput('plo', 'Type of vizualization', vizu),
                                selectInput('pos', 'Legend position', pospos),
                                numericInput('sizeplot', 'Size of plot (pixels)',600,min=200,max=3000,step = 100),
                                plotOutput('plot2')
                              ),
                              tabPanel(
                                title = "Test data results",
                                numericInput("sizeplot3", "Size of plot (pixels)",600,min=200,max=3000,step = 100),
                                selectInput('pos3', 'Legend position', pospos),
                                numericInput("cex3", "Symbol size",1,min=0.01,max=5,step = 0.01),
                                numericInput("cex32", "Text size",1,min=0.01,max=5,step = 0.01),
                                plotOutput('plottest')
                              )
                            )
                          )
                        )
               ),
               tabPanel("Credits",
                        titlePanel("Credits"),
                        tags$div(class="header", checked=NA,
                                 "App. based on the package",
                                 tags$a(href="https://github.com/hlorenzo/ddsPLS",
                                        tags$b("hlorenzo/ddsPLS"))
                        ),
                        br(),
                        "hadrien.lorenzo.2015@gmail.com",
                        br(),
                        "2023 May"
               )
    )
    #=======================================
  )
  server <- function(input, output, session) {
    #=======================================
    #=======================================
    output$files <- reactive({matrix(NA,0,3)})
    fileX <- reactive({
      input$fileX
    })
    fileY <- reactive({
      input$fileY
    })
    fileXTest <- reactive({
      input$fileXTest
    })
    sizeplot <- reactive({
      input$sizeplot
    })
    sizeplot2 <- reactive({
      input$sizeplot2
    })
    sizeplot3 <- reactive({
      input$sizeplot3
    })
    datasR <- eventReactive(input$startimport, {
      K <- nrow(fileX())
      tryCatch(
        {
          dfX_list <- lapply(fileX()$datapath,function(fifi){
            read.csv(fifi,header = input$header,sep = input$sep,dec = input$dec,
                     quote = input$quote)
          })
          dfY <- read.csv(fileY()$datapath,header = input$header,
                          sep = input$sep,dec = input$dec,quote = input$quote)
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
      ps <- unlist(lapply(dfX_list,ncol))
      q <- ncol(dfY)
      outputFiles <- matrix(NA,K+1,3); colnames(outputFiles) <- c("Block","File name","Number of variables")
      for(k in 1:K){
        outputFiles[k,] <- c(paste("X",k),fileX()$name[k],ps[k])
      }
      k <- K+1
      outputFiles[k,] <- c("Y",fileY()$name,q)
      output$files_n <- renderText(paste("Number of observations:",nrow(dfY)))
      output$files <- renderTable(outputFiles)
      blockNames <- outputFiles[,1]
      updateSelectInput(session, "inSelect",
                        choices = blockNames
      )
      list(isSimu=F,Xs = dfX_list,Y=dfY,
           ps=ps,
           outputFiles=outputFiles,
           colsReal=unlist(lapply(1:length(ps),function(k){rep(k,ps[k])})))
    })
    X_test <- eventReactive(input$startimportTest, {
      K <- nrow(fileXTest())
      tryCatch(
        {
          dfX_list <- lapply(fileXTest()$datapath,function(fifi){
            read.csv(fifi,header = input$header,sep = input$sep,dec = input$dec,
                     quote = input$quote)
          })
        },
        error = function(e) {
          stop(safeError(e))
        }
      )
      do.call(cbind,dfX_list)
    })
    output$headerBlock <- renderTable({
      dada <- datasR()
      outputFiles <- dada$outputFiles
      K <- length(dada$Xs)
      posHead <- which(outputFiles[,1]==input$inSelect)
      out <- head(dada$Y)
      if(posHead<=K){
        out <- head(dada$Xs[[posHead]])
      }
      out
    })
    output$plot1R <- renderPlot({
      dada <- datasR()
      posHead <- which(dada$outputFiles[,1]==input$inSelect)
      if(posHead<nrow(dada$outputFiles)){
        main <- paste("Correlation Y/X",posHead,sep="")
        dada1 <- dada$Xs[[posHead]]
      }else{
        main <- "Autocorrelation of Y"
        dada1 <- dada$Y
      }
      coco <- cor(dada1,dada$Y)
      layout(matrix(c(rep(1,6),2),1))
      cols <- c("#0000FF","#0303FF","#0606FF","#0909FF","#0C0CFF","#0F0FFF","#1212FF","#1515FF"
                ,"#1818FF","#1B1BFF","#1E1EFF","#2121FF","#2424FF","#2727FF","#2A2AFF","#2D2DFF"
                ,"#3030FF","#3333FF","#3636FF","#3A3AFF","#3D3DFF","#4040FF","#4343FF","#4646FF"
                ,"#4949FF","#4C4CFF","#4F4FFF","#5252FF","#5555FF","#5858FF","#5B5BFF","#5E5EFF"
                ,"#6161FF","#6464FF","#6767FF","#6A6AFF","#6D6DFF","#7070FF","#7474FF","#7777FF"
                ,"#7A7AFF","#7D7DFF","#8080FF","#8383FF","#8686FF","#8989FF","#8C8CFF","#8F8FFF"
                ,"#9292FF","#9595FF","#9898FF","#9B9BFF","#9E9EFF","#A1A1FF","#A4A4FF","#A7A7FF"
                ,"#ABABFF","#AEAEFF","#B1B1FF","#B4B4FF","#B7B7FF","#BABAFF","#BDBDFF","#C0C0FF"
                ,"#C3C3FF","#C6C6FF","#C9C9FF","#CCCCFF","#CFCFFF","#D2D2FF","#D5D5FF","#D8D8FF"
                ,"#DBDBFF","#DEDEFF","#E1E1FF","#E5E5FF","#E8E8FF","#EBEBFF","#EEEEFF","#F1F1FF"
                ,"#F4F4FF","#F7F7FF","#FAFAFF","#FDFDFF","#FFFDFD","#FFFAFA","#FFF7F7","#FFF4F4"
                ,"#FFF1F1","#FFEEEE","#FFEBEB","#FFE8E8","#FFE5E5","#FFE1E1","#FFDEDE","#FFDBDB"
                ,"#FFD8D8","#FFD5D5","#FFD2D2","#FFCFCF","#FFCCCC","#FFC9C9","#FFC6C6","#FFC3C3"
                ,"#FFC0C0","#FFBDBD","#FFBABA","#FFB7B7","#FFB4B4","#FFB1B1","#FFAEAE","#FFABAB"
                ,"#FFA7A7","#FFA4A4","#FFA1A1","#FF9E9E","#FF9B9B","#FF9898","#FF9595","#FF9292"
                ,"#FF8F8F","#FF8C8C","#FF8989","#FF8686","#FF8383","#FF8080","#FF7D7D","#FF7A7A"
                ,"#FF7777","#FF7474","#FF7070","#FF6D6D","#FF6A6A","#FF6767","#FF6464","#FF6161"
                ,"#FF5E5E","#FF5B5B","#FF5858","#FF5555","#FF5252","#FF4F4F","#FF4C4C","#FF4949"
                ,"#FF4646","#FF4343","#FF4040","#FF3D3D","#FF3A3A","#FF3636","#FF3333","#FF3030"
                ,"#FF2D2D","#FF2A2A","#FF2727","#FF2424","#FF2121","#FF1E1E","#FF1B1B","#FF1818"
                ,"#FF1515","#FF1212","#FF0F0F","#FF0C0C","#FF0909","#FF0606","#FF0303","#FF0000")
      image(coco,xaxt="n",yaxt="n",main=main,zlim=c(-1,1),col=cols)
      if(nrow(coco)>1){
        axis(1,at = (0:(nrow(coco)-1) )/(nrow(coco)-1),labels = colnames(dada1))
      }else{
        axis(1,at = 0,labels = colnames(dada1))
      }
      if(ncol(coco)>1){
        axis(2,at = (0:(ncol(coco)-1) )/(ncol(coco)-1),labels = colnames(dada$Y))
      }else{
        axis(2,at = 0,labels = colnames(dada$Y))
      }
      image(t(seq(-1,1,length.out = 24*4)),xaxt="n",las=2,zlim=c(-1,1),col=cols,yaxt="n",main="Legend")
      axis(2,at = (0:10)/10,labels = seq(-1,1,length.out = 11),las=2)
    },width = sizeplot2)
    modelAnal <- eventReactive(input$runAnal, {
      req(fileX(),fileY())
      x <- as.matrix(do.call(cbind,datasR()$Xs))
      y <- as.matrix(datasR()$Y)
      lams <- as.numeric(unlist(strsplit(input$lamsAnal,",")))
      mo <- ddsPLS(x,y,verbose=F,doBoot = F,lambdas = lams)
      return(mo)
    })
    model <- eventReactive(input$runB, {
      req(fileX(),fileY())
      x <- as.matrix(do.call(cbind,datasR()$Xs))
      y <- as.matrix(datasR()$Y)
      mo <- ddsPLS(x,y,
                   verbose=F,doBoot = T,NCORES = input$NCORES,
                   lambdas = NULL,n_B = input$n_B)
      return(mo)
    })
    output$summary <- renderPrint({
      summary(model())
    })
    output$summaryShortAnal <- renderTable({
      req(modelAnal())
      R <- modelAnal()$R
      out <- "No component built"
      if(R>0){
        expl_variance <- round(modelAnal()$varExplained$Cumu[R])
        out <- matrix(c(R,expl_variance),nrow = 1)
        colnames(out) <- c("Components","Explained variance (%)")
        # paste("ddsPLS model built on ",R," component(s).\n Explains ",expl_variance,"% of variance of Y",sep="")
      }
      out
    })
    output$summaryShort <- renderTable({
      req(model())
      R <- model()$R
      if(R>0){
        expl_variance <- round(model()$varExplained$Cumu[R])
        out <- matrix(c(R,expl_variance),nrow = 1)
        colnames(out) <- c("Components","Explained variance (%)")
        # paste("ddsPLS model built on ",R," component(s).\n Explains ",expl_variance,"% of variance of Y",sep="")
      }else{
        out <- "No component built"
      }
      out
    })
    output$plot2 <- renderPlot({
      mo <- model()
      if(is.null(mo)) mo <- modelAnal()
      noModel <- mo$R==0
      if(!noModel){
        colo <- datasR()$colsReal
        plot(mo,type = input$plo,legend.position =input$pos,col=colo)
      }else{
        plot(0,0,col="white",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
        text(x = 0,y=0,"Nothing to be plotted because model is empty.")
      }
    },height = sizeplot)
    output$plottest <- renderPlot({
      x_test <- X_test()
      if(is.data.frame(x_test)){
        x_test <- as.matrix(X_test())
      }
      diagnos <- predict(model(),x_test,legend.position =input$pos3,
                         cex=input$cex3,cex.text=input$cex32)
    },height = sizeplot3,width = sizeplot3)
    #=======================================
    #=======================================
  }
  shinyApp(ui, server, ...)
}
