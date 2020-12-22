geneOnt <-
  function(InputFile = "",
           IDType = "SYMBOL",
           OutputName = "") {
    #Take Input of Gene List and Fold Change
    PracticeInput <- read.csv(file = InputFile)
    
    #Convert the ID type to EntrezID
    geneList = as.list(PracticeInput[1])
    glLength <- lengths(geneList)
    geneList = unlist(geneList)
    geneList = as.character(geneList)
    geneList = mapIds(org.Hs.eg.db, geneList, 'ENTREZID', IDType)
    
    
    #create a dataframe of upregulated and downregulated genes
    logFC <- PracticeInput[2]
    dataframe <- data.frame(geneList, logFC)
    dataframe <- dataframe[!duplicated(dataframe$geneList),]
    dataframe <- na.omit(dataframe)
    # print(dim(dataframe))
    dataframe <- data.frame(dataframe, row.names = "geneList")
    upreg <- filter(dataframe, logFC > 0)
    downreg <- filter(dataframe, logFC < 0)
    
    #create data frame containing total enrichment analysis and filter by p values out
    totalNames <- rownames(dataframe)
    tGO <- goana(totalNames)
    tLengths <- lengths(dataframe)
    z <- as.integer(tLengths / 10)
    colnames(tGO) <- c("a", "b", "c", "d", "e")
    filter <- which(tGO$d > z & tGO$e < 0.05)
    tGO <- tGO[filter,]
    colnames(tGO) <-
      c("Term", "Ont", "Tot_Genes", "DE_genes", "P_DE_Genes")
    
    
    
    #goana and filter upreg
    upnames <- rownames(upreg)
    upGO <- goana(upnames)
    colnames(upGO) <-
      c("Term", "Ont", "Tot_Genes", "DE_genes", "P_DE_Genes")
    
    
    
    #goana and filter dowreg
    dnames <- rownames(downreg)
    dGO <- goana(dnames)
    colnames(dGO) <-
      c("Term", "Ont", "Tot_Genes", "DE_genes", "P_DE_Genes")
    
    #merge data frames and fill in the gaps
    new <- merge.data.frame(upGO, dGO, by = 0, all = TRUE)
    new2 <- new[,-1]
    rownames(new2) <- new[, 1]
    new <- new2[, c(4, 9)]
    colnames(new) <- c("Up", "Down")
    final <-
      merge.data.frame(tGO,
                       new,
                       by = 0,
                       all = FALSE,
                       all.x = TRUE)
    final[is.na(final)] <- 0
    rowNames <- final[, c(1)]
    rownames(final) <- rowNames
    final <- final[,-c(1)]
    finCol <- (final[6] - final[7])
    colnames(finCol) <- "Up_Down_Difference"
    
    
    #create a list for yes and no
    lLength <- lengths(final[1])
    lLength <- as.integer(lLength)
    newList <- list(rep("Y", lLength))
    names(newList) <- "Y_or_N"
    zscore <- (finCol / (sqrt(final[4])))
    names(zscore) <- "Z-score"
    final <- data.frame(final, finCol, zscore, newList)
    final <- final[with(final, order(Ont, P_DE_Genes)), ]
    filter(final, Up_Down_Difference != 0)
    write.csv(final, file = OutputName)
  }

ontBar <-
  function(InputFile = "",
           OntologyCutoff = "10",
           BarColor = "lightgreen",
           LineColor = "darkblue") {
    #Load table output from GeneOnt R
    GO <- read.csv(file = InputFile)
    GO <- filter(GO, Y_or_N == "Y")
    
    #Separate objects into BP, CC, and MF and give top 20 Ontologies
    BP <- filter(GO, Ont == "BP")
    CC <- filter(GO, Ont == "CC")
    MF <- filter(GO, Ont == "MF")
    
    BP <- head(BP, n = OntologyCutoff)
    CC <- head(CC, n = OntologyCutoff)
    MF <- head(MF, n = OntologyCutoff)
    
    names(BP)[1] <- paste("GOID")
    names(CC)[1] <- paste("GOID")
    names(MF)[1] <- paste("GOID")
    
    #scale the p value axis to the DEGene axis to set minimum and maximum values
    maxP_BP <- max(BP[6])
    maxGene_BP <- max(BP[5])
    maxCoef_BP <- (maxGene_BP / maxP_BP)
    
    maxP_CC <- max(CC[6])
    maxGene_CC <- max(CC[5])
    maxCoef_CC <- (maxGene_CC / maxP_CC)
    
    maxP_MF <- max(MF[6])
    maxGene_MF <- max(MF[5])
    maxCoef_MF <- (maxGene_MF / maxP_MF)
    
    #BarGraph Gen
    {
      pdf("Bars.pdf")
      
      maxP_BP <- max(BP[6])
      maxGene_BP <- max(BP[7])
      maxCoef_BP <- (maxGene_BP / maxP_BP)
      maxP_BP <- signif(maxP_BP, 1)
      bpBreaks <-
        c(0, (signif(maxP_BP, 2)) / 3, (2 * signif(maxP_BP, 2) / 3), signif(maxP_BP, 2))
      
      maxP_CC <- max(CC[6])
      maxGene_CC <- max(CC[7])
      maxCoef_CC <- (maxGene_CC / maxP_CC)
      maxP_CC <- signif(maxP_CC, 1)
      ccBreaks <-
        c(0, (signif(maxP_CC, 2)) / 3, (2 * signif(maxP_CC, 2) / 3), signif(maxP_CC, 2))
      
      maxP_MF <- max(MF[6])
      maxGene_MF <- max(MF[7])
      maxCoef_MF <- (maxGene_MF / maxP_MF)
      maxP_MF <- signif(maxP_MF, 1)
      mfBreaks <-
        c(0, (signif(maxP_MF, 2)) / 3, (2 * signif(maxP_MF, 2) / 3), signif(maxP_MF, 2))
      
      BPplot <- ggplot(BP) +
        geom_col(
          aes(
            y = DE_genes,
            x = reorder(Term, -P_DE_Genes),
            fill = "Number of Genes"
          ),
          size = 2,
          color = NA
        ) +
        geom_text(aes(x = Term, y = DE_genes, label = DE_genes),
                  hjust = 2.9,
                  size = 3.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        labs(title = "Biological Processes", y = "Number of Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(
          aes(y = maxCoef_BP * P_DE_Genes, x = reorder(Term, -P_DE_Genes)),
          size = 1.5,
          color = LineColor,
          group = 1
        ) +
        scale_y_continuous(sec.axis = sec_axis(~ . / maxCoef_BP, name = "Adjusted P Value", breaks = bpBreaks)) +
        scale_fill_manual(name = NULL,
                          values = c("Number of Genes" = BarColor)) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = LineColor)) +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(vjust = 2)) +
        theme(aspect.ratio = 2 / 3) +
        coord_flip()
      
      
      MFplot <- ggplot(MF) +
        geom_col(
          aes(
            y = DE_genes,
            x = reorder(Term, -P_DE_Genes),
            fill = "Number of Genes"
          ),
          size = 2,
          color = NA
        ) +
        geom_text(aes(x = Term, y = DE_genes, label = DE_genes),
                  hjust = 2.9,
                  size = 3.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        labs(title = "Molecular Function", y = "Number of Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(
          aes(y = maxCoef_MF * P_DE_Genes, x = reorder(Term, -P_DE_Genes)),
          size = 1.5,
          color = LineColor,
          group = 1
        ) +
        scale_y_continuous(sec.axis = sec_axis(~ . / maxCoef_BP, name = "Adjusted P Value", breaks = bpBreaks)) +
        scale_fill_manual(name = NULL,
                          values = c("Number of Genes" = BarColor)) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = LineColor)) +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(vjust = 2)) +
        theme(aspect.ratio = 2 / 3) +
        coord_flip()
      
      
      CCplot <- ggplot(CC) +
        geom_col(
          aes(
            y = DE_genes,
            x = reorder(Term, -P_DE_Genes),
            fill = "Number of Genes"
          ),
          size = 2,
          color = NA
        ) +
        geom_text(aes(x = Term, y = DE_genes, label = DE_genes),
                  hjust = 2.9,
                  size = 3.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        labs(title = "Cellular Component", y = "Number of Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(
          aes(y = maxCoef_CC * P_DE_Genes, x = reorder(Term, -P_DE_Genes)),
          size = 1.5,
          color = LineColor,
          group = 1
        ) +
        scale_y_continuous(sec.axis = sec_axis(~ . / maxCoef_BP, name = "Adjusted P Value", breaks = bpBreaks)) +
        scale_fill_manual(name = NULL,
                          values = c("Number of Genes" = BarColor)) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = LineColor)) +
        theme(legend.position = "bottom") +
        theme(plot.title = element_text(vjust = 2)) +
        theme(aspect.ratio = 2 / 3) +
        coord_flip()
      
      print(CCplot)
      print(BPplot)
      print(MFplot)
      imagelist <- list(BPplot, CCplot, MFplot)
      dev.off()
    }
  }

ontb2bBar <-
  function(InputFile = "",
           OntologyCutoff = "10",
           UpBarColor = "lightgreen",
           DownBarColor = "lightblue",
           PointColor = "darkblue") {
    GO <- read.csv(file = InputFile)
    GO <- filter(GO, Y_or_N == "Y")
    signifP <- signif(GO$P_DE_Genes, 3)
    GO <- data.frame(GO, signifP)
    colnames(GO) <-
      c(
        "X",
        "Term",
        "Ont",
        "Tot_Genes",
        "DE_Genes",
        "P_DE_Genes",
        "Up",
        "Down",
        "Up_Down_Difference",
        "Z.score",
        "Y_or_N",
        "PValue"
      )
    
    
    
    #Separate objects into BP, CC, and MF and give top 20 Ontologies
    BP <- filter(GO, Ont == "BP")
    CC <- filter(GO, Ont == "CC")
    MF <- filter(GO, Ont == "MF")
    
    BP <- head(BP, n = OntologyCutoff)
    CC <- head(CC, n = OntologyCutoff)
    MF <- head(MF, n = OntologyCutoff)
    
    names(BP)[1] <- paste("GOID")
    names(CC)[1] <- paste("GOID")
    names(MF)[1] <- paste("GOID")
    
    maxP_BP <- max(BP[6])
    maxGene_BP <- max(BP[7])
    maxCoef_BP <- (maxGene_BP / maxP_BP)
    maxP_BP <- signif(maxP_BP, 1)
    bpBreaks <-
      c(0, (signif(maxP_BP, 2)) / 3, (2 * signif(maxP_BP, 2) / 3), signif(maxP_BP, 2))
    
    maxP_CC <- max(CC[6])
    maxGene_CC <- max(CC[7])
    maxCoef_CC <- (maxGene_CC / maxP_CC)
    maxP_CC <- signif(maxP_CC, 1)
    ccBreaks <-
      c(0, (signif(maxP_CC, 2)) / 3, (2 * signif(maxP_CC, 2) / 3), signif(maxP_CC, 2))
    
    maxP_MF <- max(MF[6])
    maxGene_MF <- max(MF[7])
    maxCoef_MF <- (maxGene_MF / maxP_MF)
    maxP_MF <- signif(maxP_MF, 1)
    mfBreaks <-
      c(0, (signif(maxP_MF, 2)) / 3, (2 * signif(maxP_MF, 2) / 3), signif(maxP_MF, 2))
    
    {
      pdf("b2bBars.pdf")
      
      z <- (max(BP$Up) + max(BP$Down)) / 8
      z <- round(z, -1)
      if (z == 0) {
        z <- 10
      }
      
      BPb2bplot <- ggplot(BP, aes(x = reorder(Term, -P_DE_Genes))) +
        geom_col(aes(y = Up, fill = "Upregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(x = Term, y = Up, label = Up),
                  hjust = 1,
                  size = 3.5) +
        geom_col(aes(y = -Down, fill = "Downregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(
          x = Term,
          y = -Down,
          label = Down
        ),
        hjust = -0.1,
        size = 3.5) +
        labs(title = "Biological Processes", y = "Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        scale_y_continuous(breaks = seq(round(-max(BP$Down), -1), round(max(BP$Up), -1), z),
                           labels = abs(seq(round(
                             -max(BP$Down), -1
                           ), round(max(
                             BP$Up
                           ), -1), z))) +
        theme(aspect.ratio = 2 / 3) +
        scale_fill_manual(name = NULL, values = (
          c("Upregulated" = UpBarColor, "Downregulated" = DownBarColor)
        )) +
        theme(legend.position = "bottom") +
        scale_y_continuous(sec.axis = sec_axis(~ . / maxCoef_BP, name = "Adjusted P Value", breaks = bpBreaks)) +
        coord_flip() +
        geom_point(aes(
          y = maxCoef_BP * P_DE_Genes,
          x = Term,
          color = "Adjusted P Value"
        ),
        size = 1) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = PointColor))
      
      
      y <- (max(CC$Up) + max(CC$Down)) / 8
      y <- round(y, -1)
      if (y == 0) {
        y <- 10
      }
      
      CCb2bplot <- ggplot(CC, aes(x = reorder(Term, -P_DE_Genes))) +
        geom_col(aes(y = Up, fill = "Upregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(x = Term, y = Up, label = Up),
                  hjust = 1,
                  size = 3.5) +
        geom_col(aes(y = -Down, fill = "Downregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(
          x = Term,
          y = -Down,
          label = Down
        ),
        hjust = -0.1,
        size = 3.5) +
        labs(title = "Cellular Component", y = "Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        scale_y_continuous(breaks = seq(round(-max(CC$Down), -1), round(max(CC$Up), -1), y),
                           labels = abs(seq(round(
                             -max(CC$Down), -1
                           ), round(max(
                             CC$Up
                           ), -1), y))) +
        theme(aspect.ratio = 2 / 3) +
        scale_fill_manual(name = NULL, values = (
          c("Upregulated" = UpBarColor, "Downregulated" = DownBarColor)
        )) +
        theme(legend.position = "bottom") +
        scale_y_continuous(sec.axis = sec_axis(
          ~ . / maxCoef_CC,
          name = "Adjusted P Value",
          breaks = signif(ccBreaks, 2)
        )) +
        coord_flip() +
        geom_point(aes(
          y = maxCoef_CC * P_DE_Genes,
          x = Term,
          color = "Adjusted P Value"
        ),
        size = 1) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = PointColor))
      
      
      
      x <- (max(MF$Up) + max(MF$Down)) / 8
      x <- round(x, -1)
      if (x == 0) {
        x <- 10
      }
      
      MFb2bplot <- ggplot(MF, aes(x = reorder(Term, -P_DE_Genes))) +
        geom_col(aes(y = Up, fill = "Upregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(x = Term, y = Up, label = Up),
                  hjust = 1,
                  size = 3.5) +
        geom_col(aes(y = -Down, fill = "Downregulated"),
                 size = 2,
                 color = NA) +
        geom_text(aes(
          x = Term,
          y = -Down,
          label = Down
        ),
        hjust = -0.1,
        size = 3.5) +
        labs(title = "Molecular Function", y = "Genes", x = "Gene Ontologies") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal"
        ) +
        scale_y_continuous(breaks = seq(round(-max(MF$Down), -1), round(max(MF$Up), -1), x),
                           labels = abs(seq(round(
                             -max(MF$Down), -1
                           ), round(max(
                             MF$Up
                           ), -1), x))) +
        theme(aspect.ratio = 2 / 3) +
        scale_fill_manual(name = NULL, values = (
          c("Upregulated" = UpBarColor, "Downregulated" = DownBarColor)
        )) +
        theme(legend.position = "bottom") +
        coord_flip() +
        scale_y_continuous(sec.axis = sec_axis(
          ~ . / maxCoef_MF,
          name = "Adjusted P Value",
          breaks = signif(mfBreaks, 2)
        )) +
        geom_point(aes(
          y = maxCoef_MF * P_DE_Genes,
          x = Term,
          color = "Adjusted P Value"
        ),
        size = 1) +
        scale_color_manual(name = NULL,
                           values = c("Adjusted P Value" = PointColor))
      
      print(CCb2bplot)
      print(BPb2bplot)
      print(MFb2bplot)
      dev.off()
      }
  }

ontBubble <-
  function(InputFile = "",
           OntologyCutoff = "10",
           BubbleOutline = "darkblue",
           BubbleFill = "lightblue") {
    GO <- read.csv(file = InputFile)
    GO <- filter(GO, Y_or_N == "Y")
    
    #Separate objects into BP, CC, and MF and give top 20 Ontologies
    BP <- filter(GO, Ont == "BP")
    CC <- filter(GO, Ont == "CC")
    MF <- filter(GO, Ont == "MF")
    
    BP <- head(BP, n = OntologyCutoff)
    CC <- head(CC, n = OntologyCutoff)
    MF <- head(MF, n = OntologyCutoff)
    
    names(BP)[1] <- paste("GOID")
    names(CC)[1] <- paste("GOID")
    names(MF)[1] <- paste("GOID")
    
    pdf("Bubbleplots.pdf")
    
    
    BPbubblePlot <-
      ggplot(BP, aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        size = DE_genes
      )) +
      geom_point(alpha = 0.3,
                 color = BubbleOutline,
                 fill = BubbleFill) +
      scale_size(range = c(1, 24), name = "Number of Genes") +
      geom_text(aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        label = Term
      ), size = 3) +
      labs(title = "Biological Processes", y = "Log Adj. P-value", x = "Z-score") +
      theme(plot.title = element_text(hjust = 0.5))
    
    lBoundBP <- min(BP[10])
    lBoundBP <- (lBoundBP) - 0.75
    uBoundBP <- max(BP[10])
    uBoundBP <- (uBoundBP) + 0.75
    
    BPbubblePlot <-
      BPbubblePlot + scale_x_continuous(limits = c(lBoundBP, uBoundBP))
    
    CCbubblePlot <-
      ggplot(CC, aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        size = DE_genes
      )) +
      geom_point(alpha = 0.3,
                 color = BubbleOutline,
                 fill = BubbleFill) +
      scale_size(range = c(1, 24), name = "Number of Genes") +
      geom_text(aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        label = Term
      ), size = 3) +
      labs(title = "Cellular Component", y = "Log Adj. P-value", x = "Z-score") +
      theme(plot.title = element_text(hjust = 0.5))
    
    lBoundCC <- min(CC[10])
    lBoundCC <- (lBoundCC) - 0.75
    uBoundCC <- max(CC[10])
    uBoundCC <- (uBoundCC) + 0.75
    
    CCbubblePlot <-
      CCbubblePlot + scale_x_continuous(limits = c(lBoundCC, uBoundCC))
    
    MFbubblePlot <-
      ggplot(MF, aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        size = DE_genes
      )) +
      geom_point(alpha = 0.3,
                 color = BubbleOutline,
                 fill = BubbleFill) +
      scale_size(range = c(1, 24), name = "Number of Genes") +
      geom_text(aes(
        x = Z.score,
        y = -log10(P_DE_Genes),
        label = Term
      ), size = 3) +
      labs(title = "Molecular Function", y = "Log Adj. P-value", x = "Z-score") +
      theme(plot.title = element_text(hjust = 0.5))
    
    lBoundMF <- min(MF[10])
    lBoundMF <- (lBoundMF) - 0.75
    uBoundMF <- max(BP[10])
    uBoundMF <- (uBoundMF) + 0.75
    
    MFbubblePlot <-
      MFbubblePlot + scale_x_continuous(limits = c(lBoundMF, uBoundMF))
    
    print(BPbubblePlot)
    print(CCbubblePlot)
    print(MFbubblePlot)
    dev.off()
    
    
  }

ontComp <- function(InputOne = "",
                    InputTwo = "") {
  first <- read.csv(InputOne)
  second <- read.csv(InputTwo)
  combined <- rbind.data.frame(first, second)
  common <- merge(first, second, by = "X")
  first <- merge(x = first,
                 y = common,
                 by = "X",
                 all.x = TRUE)
  noDups <- combined[!duplicated(combined$X),]
  vec <- c()
  com_vec <- as.vector(t(common[1]))
  first_vec <- as.vector(t(first[1]))
  second_vec <- as.vector(t(second[1]))
  for (x in as.vector(t(noDups[1]))) {
    if (x %in% com_vec) {
      vec <- c(vec, "Both")
    }
    else if (x %in% first_vec) {
      vec <- c(vec, "Input1")
    }
    else if (x %in% second_vec) {
      vec <- c(vec, "Input2")
    }
  }
  output <- data.frame(noDups[1], noDups[2], noDups[3], vec)
  write.csv(output, file = "out.csv", row.names = FALSE)
}

ontVenn <-
  function(InputFile,
           FirstInputFileName = "Input 1",
           SecondInputFileName = "Input 2",
           FirstInputColor = 'lightgreen',
           SecondInputColor = 'blue',
           Title = "Gene Ontology Comparison")
  {
    input <- read.csv(InputFile)
    set1 <- input$X[input$vec == "Input1" | input$vec == "Both"]
    set2 <- input$X[input$vec == "Input2" | input$vec == "Both"]
    
    
    venn.diagram(
      x = list(set1, set2),
      category.names = c(FirstInputFileName, SecondInputFileName),
      filename = 'VennDiagram.png',
      output = TRUE,
      
      lwd = 2,
      lty = 'blank',
      fill = c(FirstInputColor, SecondInputColor),
      main = Title
    )
  }

ontVennSub <-
  function(InputFile,
           FirstInputFileName = "Input 1",
           SecondInputFileName = "Input 2",
           FirstInputColor = 'lightgreen',
           SecondInputColor = 'blue') {
    input <- read.csv(InputFile)
    bp <- input[input["Ont"] == "BP", ]
    set1 <- bp$X[bp$vec == "Input1" | bp$vec == "Both"]
    set2 <- bp$X[bp$vec == "Input2" | bp$vec == "Both"]
    
    
    venn.diagram(
      x = list(set1, set2),
      category.names = c(FirstInputFileName, SecondInputFileName),
      filename = 'BP.png',
      main = "Biological Processes Ontology Comparison",
      output = TRUE,
      lwd = 2,
      lty = 'blank',
      fill = c(FirstInputColor, SecondInputColor)
    )
    
    # CC
    cc <- input[input["Ont"] == "CC", ]
    set1 <- cc$X[cc$vec == "Input1" | cc$vec == "Both"]
    set2 <- cc$X[cc$vec == "Input2" | cc$vec == "Both"]
    
    
    venn.diagram(
      x = list(set1, set2),
      category.names = c(FirstInputFileName, SecondInputFileName),
      filename = 'CC.png',
      main = "Cellular Component Ontology Comparison",
      output = TRUE,
      lwd = 2,
      lty = 'blank',
      fill = c(FirstInputColor, SecondInputColor)
    )
    
    # MF
    mf <- input[input["Ont"] == "MF", ]
    set1 <- mf$X[mf$vec == "Input1" | mf$vec == "Both"]
    set2 <- mf$X[mf$vec == "Input2" | mf$vec == "Both"]
    
    
    venn.diagram(
      x = list(set1, set2),
      category.names = c(FirstInputFileName, SecondInputFileName),
      filename = 'MF.png',
      main = "Molecular Function Ontology Comparison",
      output = TRUE,
      lwd = 2,
      lty = 'blank',
      fill = c(FirstInputColor, SecondInputColor)
    )
  }
