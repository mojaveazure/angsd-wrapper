library(shiny)
library(genomeIntervals)
library(lattice)
library(Hmisc)
library(ape)  # Used for ABBA BABA graphing
library(data.table)
library(DT)  # Allows R data objects to be displayed as tables on HTML pages
options(shiny.maxRequestSize = - 1)

# Define headers for thetas, Fst and intersect data
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)", "Chr","WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")
sfs.headers <- c("Allele_Frequency")

not.loaded <- TRUE

# Define server logic required to draw plots
shinyServer(

  function(input, output) {
    
    # Thetas input data
    dataInputThetas = reactive({
      data <- input$userThetas
      path <- as.character(data$datapath)
      thetas <- fread(input = path, sep = " ", header = TRUE)
      #setnames(thetas, thetas.headers)
      return(thetas)
    })
    
    # SFS input data
    dataInputSFS = reactive({
      data <- input$userSFS
      path <- as.character(data$datapath)
      Derived <- as.matrix(fread(input = path, header = FALSE))
      sfs <- as.data.frame(t(Derived))
      setnames(sfs, sfs.headers)
      return(sfs) 
    })

    # Fst input data
    dataInputFst = reactive({
      data <- input$userFst
      path <- as.character(data$datapath)
      fst <- fread(file = path, header = TRUE)
      return(fst)
    })

    # Intersect input data
    dataInputIntersect = reactive({
      data <- input$userIntersect
      path <- as.character(data$datapath)
      intersect <- read.table(file=path,
                              sep = "\t",
                              col.names=intersect.headers
      )
      return(intersect)
    })

    # Admixture input data
    dataInputAdmix = reactive({
      data <- input$userAdmix
      path <- as.character(data$datapath)
      admix <- t(as.matrix(read.table(path, header = FALSE, fill = TRUE)))
      return(admix)
    })
    
    # ABBA BABA test input data
    dataInputABBABABA = reactive({
      data <- input$userABBABABA
      path <- as.character(data$datapath)
      ABBABABA <- read.table(path, sep = "\t", header = T)
      return(ABBABABA)
    })

    # PCA plot input data
    dataInputPCA = reactive({
      data <- input$userPCA
      path <- as.character(data$datapath)
      PCA <- read.table(path, header=F)
      return(PCA)
    })

    # GFF file input
    gffInput = reactive({
      data <- input$userAnnotations
      path <- as.character(data$datapath)
      gff <- readGff3(path)
      return(gff)
    })
    
    # Thetas output data
    output$thetaChroms = renderUI({
      if(is.null(input$userThetas)){
        choices <- 10
      }
      else{
      thetas <- dataInputThetas()
      choices <- unique(thetas$Chr)
      }
      selectInput('thetaChrom', 'Chromosome to plot', 
                  choices,
                  multiple = TRUE,
                  selected = choices
                  )
    })
 
    # Create zoomable thetas plots
    ranges <- reactiveValues(x = NULL, y = NULL)
    
    # Create reactive thetaPlot1 to click and drag to select area
    output$thetaPlot1 <- renderPlot({
      #  error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, 
      error = function(err) {
        thetas <- fread(input = "BKN_Thetas.graph.me",
                        sep = " ",
                        header = TRUE)
      }
      )

      thetas <- subset(thetas, Chr == input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type = "gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name == thetas$Chr[1])
        if(length(gff.df.chr$seq_name) == 0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type == "gene")
      }

      # What we are plotting
      thetas.plot <- thetas

      # remove nsites = 0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data
      if(input$rm.nsites) {
        thetas.plot <- subset(thetas.plot, nSites > input$nsites)
      }
      # Divide thetas by the number of sites in each window
      thetas.plot$tW <- thetas.plot$tW/thetas.plot$nSites
      thetas.plot$tP <- thetas.plot$tP/thetas.plot$nSites
      thetas.plot$tF <- thetas.plot$tF/thetas.plot$nSites
      thetas.plot$tH <- thetas.plot$tH/thetas.plot$nSites
      thetas.plot$tL <- thetas.plot$tW/thetas.plot$nSites

      data <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetas.plot$tW,
                     "Pairwise Theta" = thetas.plot$tP,
                     "Fu and Li's Theta" = thetas.plot$tF,
                     "Fay's Theta" = thetas.plot$tH,
                     "Maximum likelihood (L) Theta" = thetas.plot$tL
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t = "p", pch = 19,col = rgb(0,0,0,0.5),
             xlab = "Position (bp)",
             ylab = paste(input$thetaChoice, "Estimator Value"),
             main = paste("Estimators of theta along chromosome", input$thetasText, sep = " ")
        )

        # Different representation of the data to plot
        rug(rect(gff.df.gene$X1, -1e2, 
                 gff.df.gene$X2, 0, 
                 col = rgb(0.18,0.55,0.8,0.75), 
                 border = NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter, 
                                           data, f=0.1), 
                                    col="red")}
      }
      else {
        plot(thetas.plot$WinCenter,
             data, t = "p", pch = 19,col = rgb(0,0,0,0.5),
             xlab = "Position (bp)",
             ylab = paste(input$thetaChoice, "Estimator Value"),
             main = paste("Estimators of theta along chromosome", input$thetasText, sep = " ")
        )
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter, 
                                           data, f = 0.1), 
                                    col="red")}
      }
    })
    
    # Create reactive plot that zooms in below thetaPlot1
    output$thetaPlot2 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
          thetas <- fread(input = "BKN_Thetas.graph.me",
                          sep = " ",
                          header = TRUE)
      }
      )
      
      thetas <- subset(thetas, Chr == input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name == thetas$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type == "gene")
      }
      
      # What we are plotting
      thetas.plot <- thetas
      
      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetas.plot <- subset(thetas.plot, nSites > input$nsites)
      }
      # Divide thetas by the number of sites in each window
      thetas.plot$tW <- thetas.plot$tW/thetas.plot$nSites
      thetas.plot$tP <- thetas.plot$tP/thetas.plot$nSites
      thetas.plot$tF <- thetas.plot$tF/thetas.plot$nSites
      thetas.plot$tH <- thetas.plot$tH/thetas.plot$nSites
      thetas.plot$tL <- thetas.plot$tW/thetas.plot$nSites
      
      # Choose one of the following data sets to plot
      data <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetas.plot$tW,
                     "Pairwise Theta" = thetas.plot$tP,
                     "Fu and Li's Theta" = thetas.plot$tF,
                     "Fay's Theta" = thetas.plot$tH,
                     "Maximum likelihood (L) Theta" = thetas.plot$tL
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t = "p", pch = 19, col = rgb(0,0,0,0.5),
             xlab = "Position (bp)",
             ylab = paste(input$thetaChoice,"Estimator Value"),
             paste("Estimators of theta along chromosome", input$thetasText, sep = " ")
        )
        
        rug(rect(gff.df.gene$X1, -1e2, 
                 gff.df.gene$X2, 0, 
                 col=rgb(0.18,0.55,0.8,0.75), 
                 border=NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter, data, f=0.1), 
                                    col="red")}
      }
      else {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main = paste("Estimators of theta along chromosome", input$thetasText, sep = " "),
             xlim = ranges$x, ylim = ranges$y
        )
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter, 
                                           data, f = 0.1), col = "red")}
      }
    })

    # Creating hover function to output table of x and y values of thetasPlot1 and thetasPlot2
    output$thetaPlot2_hoverinfo <- renderPrint({
      cat("Theta Plot Hover Function")
      str(input$thetaPlot2_hover)
    })
    
    # Creating zoom function in thetasPlot1 and thetasPlot2
    observe({
      brush <- input$thetaPlot1_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    
    # Creating reactive plot for selectionPlot1 and selectionPlot2
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    
    # selectionPlot1
    output$selectionPlot1 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- fread(
            input = "BKN_Thetas.graph.me",
            sep = " ",
            header = TRUE
        )
      }
      )

      thetas <- subset(thetas, Chr == input$thetaChrom)

      # What we are plotting
      thetas.plot <- thetas

      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)

      data <- switch(input$selectionChoice,
                     "Tajima's D" = thetas.plot$Tajima,
                     "Fu and Li's F" = thetas.plot$fuf,
                     "Fu and Li's D" = thetas.plot$fud,
                     "Fay and Wu's H" = thetas.plot$fayh,
                     "Zeng's E" = thetas.plot$zeng
      )
      
      plot(thetas.plot$WinCenter,
           data, t = "p", pch = 19, col = rgb(0,0,0,0.5),
           xlab = "Position (bp)",
           ylab = input$selectionChoice,
           main = paste("Neutrality test statistics along chromosome", input$thetasText, sep = " ")
      )
      if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter, 
                                             data, f=0.1), col="red")}
    })
    
    # selectionPlot2
    output$selectionPlot2 <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
          thetas <- fread(
              input = "BKN_Thetas.graph.me",
              sep = " ",
              header = TRUE
          )
      }
      )
  
      thetas <- subset(thetas, Chr == input$thetaChrom)
      
      # What we are plotting
      thetas.plot <- thetas
      
      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      
      data <- switch(input$selectionChoice,
                     "Tajima's D" = thetas.plot$Tajima,
                     "Fu and Li's F" = thetas.plot$fuf,
                     "Fu and Li's D" = thetas.plot$fud,
                     "Fay and Wu's H" = thetas.plot$fayh,
                     "Zeng's E" = thetas.plot$zeng
      )
      
      # Output plot of thetas
      plot(thetas.plot$WinCenter,
           data, t = "p", pch = 19, col = rgb(0,0,0,0.5),
           xlab = "Position (bp)",
           ylab = input$selectionChoice,
           main = paste("Neutrality test statistics along chromosome", input$thetasText, sep = " "),
           xlim = ranges2$x, 
           ylim = ranges2$y
      )
      if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter, 
                                             data, f = 0.1), col = "red")}
    })
    
    # Creating hover function to output table of x and y values of selectionPlot1 and selectionPlot2
    output$selectionPlot2_hoverinfo <- renderPrint({
      cat("Theta Plot Hover Function")
      str(input$selectionPlot2_hover)
    })
    
    # Creating zoom function in selectionPlot1 and selectionPlot2
    observe({
      brush <- input$selectionPlot1_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })

    # Create zoomable PCA plots
    ranges3 <- reactiveValues(x = NULL, y = NULL)
    

    # SFS Plot output
    output$SFSPlot <- renderPlot({
      sfs <- tryCatch({
        dataInputSFS()

      }, error = function(err) {
        Derived <- as.matrix(fread("output_barley1_PH_DerivedSFS", header = FALSE))
        sfs <- as.data.frame(t(Derived))
        setnames(sfs, sfs.headers)
      })

      # Graph SFS here
      sfs.AFreq <- sfs$Allele_Frequency
      # Throw out 0th class and nth class
      alleles <- sfs.AFreq[seq(2, nrow(sfs)-1)]
      # Proportion to plot
      alleles.prop <- alleles/sum(alleles)
      # Generate barplot
      sfs.bp <- barplot((alleles/sum(alleles)),
              xaxt = "n",
              xlab = "Derived Allele Frequency",
              ylab = "Proportion of SNPs",
              main = paste(input$sfsText, "Site Frequency Spectrum", sep = " "),
              offset = 0,
              xlim = NULL,
              ylim = NULL,
              axes = TRUE,
              names = 1:length(alleles),
              las = 1,
              pch = 18,
              xpd = TRUE,
              col = "#5F9ED1")
      # Label x-axis
      lab <- c(1:length(sfs.bp))
      axis(1, at = sfs.bp, labels = lab)
    })

    # Admixture plot
    output$admixPlot <- renderPlot({
      admix <- tryCatch({
        dataInputAdmix()
      },error = function(err) {
        admix <- t(as.matrix(read.table("ngsadmix_example.txt")))
      })
      barplot(admix, col=c("#006BA4","#FF800E","#A2C8EC",
                           "#898989","#ABABAB","#595959",
                           "#5F9ED1","#CFCFCF","#FFBC79","#C85200"),
              space = 0,
              border = NA,
              xlab = "Individuals",
              ylab = "Admixture proportion")
    })
    
    # ABBA BABA plot output
    output$ABBABABATree <- renderPlot({
      ABBABABA <- tryCatch({
        dataInputABBABABA()
      }, error = function(err){
        ABBABABA <- read.table("abbababa.test", sep="\t", header=T)
      })
      d.current <- subset(ABBABABA, H2 == input$h2 & H3 == input$h3)
      tree <- read.tree(text=paste("(Outgroup,(", input$h3, ",(", input$h2, ",Taxon)));", sep=""))
      plot(tree, type = "cladogram", edge.width = 2, direction='downwards')

    })
    # Graph dotplot
    output$ABBABABAPlot <- renderPlot({
      ABBABABA <- tryCatch({
        dataInputABBABABA()
      }, error = function(err){
        ABBABABA <- read.table("abbababa.test", 
                               sep="\t", 
                               header=T)
      })
      d.c <- subset(ABBABABA, H2 == input$h2 & H3 == input$h3)
      d.current <- as.vector(d.c)
      mypanel.Dotplot <- function(x, y, ...) {
        panel.Dotplot(x, y, ...)
        tips <- attr(x, "other")
        panel.abline(v = 0, lty = 3)
        trellis.par.set(mfrow = c(2, 1))
        panel.arrows(x0 = tips[,1], y0 = y,
                     x1 = tips[,2], y1 = y,
                     length = 0.05, unit = "native",
                     angle = 90, code = 3)
      }
      # Creating plot with error bars 
      Dotplot(factor(d.current$H1) ~ Cbind(d.current$Dstat, 
                                           d.current$Dstat-d.current$SE, 
                                           d.current$Dstat+d.current$SE), 
              col = "blue", pch = 20, panel = mypanel.Dotplot, 
              xlab = "D", ylab = "Taxon", 
              title = paste("D statistic comparison where H2=", input$h2, " and H3=", input$h3, sep = ""))
    })
    # Output table of ABBA BABA values
    output$ABBABABATable <- renderDataTable({
      ABBABABA <- tryCatch({
        dataInputABBABABA()
      }, error = function(err){
        ABBABABA <- read.table("abbababa.test", sep="\t", header=T)
      })
      abTable <- as.data.table(ABBABABA)
      ABBABABATable <- setkey(abTable)
      colnames(ABBABABATable)[10] <- "Pvalue"
      pvalue <- as.data.table(2*pnorm(-abs(ABBABABATable$Z)))
      as.headers.pval <- c("Pvalue")
      setnames(pvalue, as.headers.pval)
      ABBABABATable$Pvalue <- pvalue
      return(ABBABABATable)
    })

    # Create zoomable PCA plots
    ranges4 <- reactiveValues(x = NULL, y = NULL)
    
    # Create PCAPlot1 you can select areas to zoom in on top
    output$PCAPlot1 <- renderPlot({
      PCA <- tryCatch({
        dataInputPCA()
      }, error = function(err) {
        PCA <- read.table("all.pop.covar", header=F)
      })
      eig <- eigen(PCA, symm=TRUE);  # computes eigenvalues of matrices
      # Plot proportion of values
      eig$val <- eig$val/sum(eig$val);
      PC <- as.data.frame(eig$vectors)
      plot(PC$V1, PC$V2, 
           pch=19, col=rgb(0,0,0,0.4), 
           xlab="PC1", ylab="PC2", 
           main="ngsCovar Results",
           asp=1)
    })
    
    # Create PCAPlot2 that will react to PCAPlot1 selected area
    output$PCAPlot2 <- renderPlot({
      PCA <- tryCatch({
        dataInputPCA()
      }, error = function(err) {
        PCA <- read.table("all.pop.covar", header=F)
      })
      eig <- eigen(PCA, symm=TRUE);  # computes eigenvalues of matrices
      eig$val <- eig$val/sum(eig$val);
      PC <- as.data.frame(eig$vectors)
      plot(PC$V1, PC$V2, 
           pch=19, col=rgb(0,0,0,0.4), 
           xlab="PC1", ylab="PC2", 
           main="ngsCovar Results",
           xlim = ranges4$x, ylim = ranges4$y,
           asp=1)
    })
    # Creating hover function in PCAPlot2
    output$PCAPlot2_hoverinfo <- renderPrint({
      cat("PCA Plot Hover Function")
      str(input$PCAPlot2_hover)
    })
    
    # Creating zoom function in PCAPlot1 and PCAPlot2
    observe({
      brush <- input$PCAPlot1_brush
      if(!is.null(brush)) {
        ranges4$x <- c(brush$xmin, brush$xmax)
        ranges4$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges4$x <- NULL
        ranges4$y <- NULL
      }
    })
  })