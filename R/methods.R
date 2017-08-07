## Creator: Wajid Jawaid
## Author: Wajid Jawaid
## Date: 2 December 2014
## Gottgens lab scRNA-seq tool repo

## Methods definitions file

##' Returns eigen values
##'
##' Returns eigen values
##' @title Eigen values of Reduced Dimension object
##' @param object Dimensionality reduced object
##' @return numeric vestor of ordered eigen values
##' @author Wajid Jawaid
##' @export
setMethod("eigenvals", "reducedDim", function(object) object@eigenvalues)

##' Returns eigen vectors
##'
##' Returns eigen values
##' @title Eigen vectors
##' @param object Dimensionality reduced object
##' @return Annotated matrix of eigen vectors, with eigen vectors in columns.
##' @author Wajid Jawaid
##' @export
setMethod("eigenvecs", "reducedDim", function(object) object@eigenvectors)

##' Generates minmium spanning tree
##'
##' Generates minmium spanning tree and adds diameter path to \code{paths} slot
##' @title Minimum Spanning Tree for Reduced Dimensionality Object
##' @param object Dimensionality reduced object
##' @param weighted Default TRUE. Passed to graph.adjacency() function in igraph.
##' @param useDims Default is NULL.
##' @return Dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("mst", "reducedDim", function(object, weighted = TRUE, useDims=NULL) {
    ev <- eigenvecs(object)
    if (is.null(useDims)) useDims <- 1:ncol(ev)
    d2 <- as.matrix(dist(ev[,useDims]))
    gs <- graph.adjacency(d2, mode = "undirected", weighted = weighted)
    gsMST <- minimum.spanning.tree(gs, weights=E(gs)$weight, algorithm="prim")
    object@graph <- gsMST
    diamPath <- get.diameter(gsMST)
    object <- addPath(object, pathName = "DiameterPath",
                      pathDescription = "Cells found on unsupervized diameter path",
                      path = diamPath)
    return(object)
})

##' Returns the minimum spanning tree
##'
##' Returns the minimum spanning tree
##' @title Retrieves the minimum spanning tree from the \code{graph} slot of a
##' Dimensionality reduced object
##' @param object Dimensionality reduced object
##' @return Minimum spanning tree as an igraph object
##' @author Wajid Jawaid
##' @export
setMethod("getMST", "reducedDim", function(object) object@graph)

##' Retrieve paths from a dimensionality reduced object
##'
##' Retrieve paths from a dimensionality reduced object
##' @title Retrieve paths from a dimensionality reduced object
##' @param object Dimensionality reduced object
##' @return The stored paths for the object
##' @author Wajid Jawaid
##' @export
setMethod("getPaths", "reducedDim", function(object) object@paths)

##' 3D plotting function
##'
##' 3D plotting function
##' @title 3D plotting function
##' @param object A Reduced Dimension object 
##' @param colorBy Character or numeric vector. Character vectors will be converted to a factor
##' and cells coloured accordingly. Cell color can be supplied in the plotCols parameter. If
##' not supplied the ggplots2 coloring scheme will be adopted.
##' @param plotCols Colors to use for the plot
##' @param legendSize Controls the size of the legend
##' @param legendOffset Distance from the edge of the axis to place the legend as a proportion
##' of the length of the axis.
##' @param legPlac Controls the axis along which the legend is placed.
##' @param legPos Controls the depth at which the legend is placed. Default is "front".
##' @param negLeg Moves the legend to the other extreme of the axis
##' @param plotLegend Logical value. Default is "True"
##' @param legVert Controls the vertical orientation of the legend
##' @param legTextPos Controls the position of the legend text
##' @param axLabs Prefix with which to label the axes.
##' @param doPath Logical value controlling whether the diameter path is plotted
##' @param doTree Logical calue controlling whether the minimum spanning tree is plotted
##' @param selectedPath The path to plot. Default is NULL and will plot the diameter path.
##' @param useDims The dimensions of the the dimensionality reduced objects to use
##' @param selectedCells The selected cells to plot
##' @param legendTitle Legend title.
##' @param diamWidth Width of diameter path to be plotted
##' @param diamCol Color of diameter path
##' @param returnID Whether to return the rgl plot id. For use with selecting cells.
##' @param na.col Colour for missing values. Default is Black.
##' @param bringToFront To bring higher colours to fore.
##' @param project2D Default NULL. Projects 3D graph into 2D by the given transformation
##' matrix. This matrix must have been produced using \code{getMVP} on an open rgl device.
##' @param opacity NULL. If given will use the data given in this column to set point opacity.
##' Ensure these are numeric and between 0-1
##' @param size Default 5. Sets point size.
##' @param pch Default 20. Sets point size.
##' @param outline Default is NULL. If require an outline set this to required colour.
##' @param ... Parameters to be passed to \code{plot3d} function.
##' @return 3D plot
##' @author Wajid Jawaid
##' @export
setMethod("plot3", "reducedDim", function(object, colorBy = "", plotCols = NULL,
                                          legendSize = 1,
                                          legendOffset = 0.2, legPlac = 1:3, legPos = "front",
                                          negLeg = FALSE, plotLegend = TRUE, legVert = 1,
                                          legTextPos = 1, axLabs = "Component", doPath = FALSE,
                                          doTree = FALSE, selectedPath = NULL, useDims = 1:3,
                                          selectedCells = c(NA, NULL), legendTitle = "",
                                          diamWidth = 2, diamCol = "black", returnID = FALSE,
                                          na.col="Black", bringToFront = TRUE,
                                          project2D = NULL,
                                          opacity=NULL, size=5, pch=20, outline=NULL, ...) {
    ev <- eigenvecs(object)
    .useDims <- useDims
    if (!is.null(project2D)) {
        ev <- applyMVP(ev[,useDims], project2D)
        useDims <- 1:2
    }
    if (any(is.na(match(useDims, 1:ncol(ev))))) {
        stop("Selected plot dimensions and eigen vector matrix not consistent")
    }
    if (length(useDims) > 3) stop("Cannot plot more than 3 dimensions.")
    pTitle <- attr(colorBy, "title")
    if (length(colorBy) != nrow(ev)) {
        cat("Length of ColorBy incorrect. Colour set to ggPlot colours\n")
        plotCols <- NULL
        # colorBy <- rep("", nrow(ev))
    }
    if (is.null(outline)) {
        outline <- FALSE
    } else {
        outline.col <- outline
        outline <- TRUE
    }
    if (outline) {
        if (!any(pch==21:25)) pch=21    # Check that outline compatible pch is used.
    }
    ev <- ev[,useDims]
    if (!is.numeric(colorBy)) {
        colorBy <- as.factor(colorBy)
        nCols <- length(levels(colorBy))        
        if (!is.null(plotCols)) {
            if (length(plotCols) != nCols) stop("Incorrect number of colours supplied.")
        } else {
            ## ggColors <- hcl(h = seq(15, 375 - 360/nCols, length.out = nCols) %% 360, c = 100,
            ##                 l = 65)
            ggColors <- ggCol(nCols, remix = FALSE)
            plotCols <- ggColors
        }
        ## if (nCols ==1) plotCols <- "black"
        pointColor <- plotCols[as.numeric(colorBy)]
    } else {
        if (is.null(plotCols)) plotCols <- c("grey", "red")
        pointColor <- colorGradient(colorBy, plotCols[1], plotCols[2], na.col=na.col)
        if (length(plotCols) > 2) cat("plotCols > 2. Only first two colours used")
    }
    if (outline)
        if (length(outline.col) == 1) outline.col <- rep(outline.col, length(pointColor))
    if (!is.null(opacity)) {
        pointColor <- rgb(t(col2rgb(pointColor)), alpha=opacity, maxColorValue=255)
        if (outline)
            outline.col <- rgb(t(col2rgb(outline.col)), alpha=opacity, maxColorValue=255)
    }
    if (class(axLabs) == "expression") {
        colnames(ev) <- paste(axLabs, "[", useDims, "]", sep="")
    } else if (as.character(axLabs) != "") {
        colnames(ev) <- paste(axLabs, useDims)
    } else colnames(ev) <- rep("", ncol(ev))
    axLabs <- colnames(ev)
    axLabs <- gsub(" ", "_", axLabs)
    if (is.na(selectedCells)) selectedCells <- 1:nrow(ev)
    mlabel <- function(x) bquote(.(parse(text = x)))  # Function to allow symbols in labels
    if (length(useDims) == 2) {
        if (bringToFront) {
            pC <- order(colorBy, decreasing=FALSE)
            evO <- ev[pC,]
            pointCol <- pointColor[pC]
            if (outline) outline.col <- outline.col[pC]
        }
        if (plotLegend && is.numeric(colorBy)) {
            layout(matrix(1:2, 2), heights = c(1,5))
            legCol <- attr(pointColor, "leg")
            origPar <- par(no.readonly = TRUE)
            par(mar=c(.5,2,4,2), mgp=c(3,.22,0), las=1)
            image(as.numeric(legCol[,1]), 1,
                  matrix(seq_along(legCol[,1]),ncol=1), col=legCol[,2], axes=FALSE,
                  xlab="", ylab="", main = pTitle, cex.axis=.5,
                  sub = expression(paste(log[10], "transformed normalised counts")))
            axis(1, cex.axis=.5, tck=-.1)
            par(mar = c(4.1, 4.1, 1.1, 2.1), mgp = origPar$mgp)
        }
        if (!outline) {
            cellPlot <- plot(evO[selectedCells,], col = pointCol, xlab = mlabel(axLabs[1]),
                             ylab = mlabel(axLabs[2]), pch=pch, ...)
        } else {
            cellPlot <- plot(evO[selectedCells,], col = outline.col, xlab = mlabel(axLabs[1]),
                             ylab = mlabel(axLabs[2]), pch=pch, bg = pointCol, ...)
        }
        if (plotLegend && !(is.numeric(colorBy))) {
            if (legPos == "front") legPos <- "topright"
            legend(legPos, pch = 16, legend = levels(colorBy), col = plotCols)
        }
    } else {
        cellPlot <- plot3d(ev[selectedCells,], col = pointColor, xlab = axLabs[1],
                           ylab = axLabs[2], size=size,
                           zlab = axLabs[3], alpha = opacity/255, ...)
        par3d(ignoreExtent = TRUE)
        ## Legend plotting
        legCoords <- apply(ev, 2, range)[,legPlac]
        if (legPos == "front") legCoords[,2] <- rev(legCoords[,2])
        legCoord <- legCoords[2,]
        legMargSign <- 1
        if (negLeg) {
            legCoord[1] <- legCoords[1,1]
            legMargSign <- -1
        }
        legMargin <- diff(legCoords[,1]) * legendOffset
        legCoord[1] <- legCoord[1] + (legMargin * legMargSign)
        if (legVert == -1) legCoord[3] <- legCoords[1,3]
        if (exists("nCols")) {
            if (nCols > 1) {
                lineSpace = legVert * (diff(legCoords[,3]) * legendSize / (nCols - 1))
            } else lineSpace = 1
            legCoordTab <- data.frame(matrix(rep(legCoord, nCols), ncol = 3, byrow = TRUE),
                                      stringsAsFactors = FALSE)
            legCoordTab[,3] <- seq(legCoordTab[1,3], by = -lineSpace, length.out = nCols)
            legCoordTab[,4] <- levels(colorBy)
            legCoordTab[,5] <- plotCols
            legCoordTab[,6] <- legCoordTab[,1] + (legMargin * legMargSign * legTextPos) / 4
            writeLeg <- function(x) {
                points3d(x=x[1], y=x[2], z=x[3], color=x[5], size = 5)
                x[legPlac[1]] <- x[6]
                text3d(x=x[1], y=x[2], z=x[3], texts=x[4], color="black", adj=0)
            }
            if (plotLegend) {
                ordLeg <- order(legPlac)
                legTitle <- legCoordTab[1,]
                legTitle[,3] <- legTitle[,3] + lineSpace
                apply(legCoordTab[,c(ordLeg,4:6)], 1, writeLeg)
                text3d(legTitle[ordLeg], texts=legendTitle, adj=0)
            }
        }
        ## Legend plotting end
    }
    if(!doTree & !doPath & length(useDims) == 3) invisible(cellPlot)
    if (doTree) plotTree3(object, useDims = .useDims, pointColor = pointColor,
                          project2D = project2D)
    if (doPath) {
        if (is.null(selectedPath)) selectedPath <- "DiameterPath"
        for (i in 1:length(selectedPath)) {
            plotPath3(object, diamWidth = diamWidth, diamCol = diamCol[i], useDims = .useDims,
                      selectedPath = selectedPath[i], project2D = project2D)
        }
    }
    cellPlot <- list(cellPlot, legend = levels(colorBy), col = plotCols)
    invisible(cellPlot)
})

setMethod("plotTree3", "reducedDim", function(object, useDims, pointColor, project2D, ...) {
    ev <- eigenvecs(object)
    edgeList <- get.edgelist(object@graph)
    n <- nrow(edgeList)
    dfo <- order(c(seq(from=1, by=2, length.out=n), seq(from=2, by=2, length.out=n)))
    segments_df <- c(edgeList[,1], edgeList[,2])[dfo]
    inds <- match(segments_df, rownames(ev))
    if (!is.null(project2D)) {
        ev <- applyMVP(ev[,useDims], project2D)
        useDims <- 1:2
    }
    segments_df <- cbind.data.frame(segments_df, ev[inds,useDims], row.names=1:length(segments_df))
    if (length(useDims) == 2) {
        segments_df <- cbind(segments_df, rbind(segments_df[-1,2:3], c(0,0)))
        segments_df <- segments_df[-seq(from=2, by=2, to=nrow(segments_df)),]
        colnames(segments_df) <- c("name", "x0", "y0", "x1", "y1")
        segments(x0=segments_df[,2], y0=segments_df[,3], x1=segments_df[,4], y1=segments_df[,5], col = pointColor[inds])
    }
    else segments3d(segments_df[,2:4], color = pointColor[inds])
})

##' 3D path plot
##'
##' 3D path plot
##' @title Plot path in 3D plot. Only to be called after \code{plot3}
##' @inheritParams plotPath3
##' @return Plots 3D path
##' @author Wajid Jawaid
setMethod("plotPath3", "reducedDim", function(object, diamWidth, diamCol, useDims, selectedPath,
                                              project2D) {
    chosenPath <- ifelse(is.null(selectedPath), "DiameterPath", selectedPath)
    pathNames <- do.call(c, lapply(object@paths, "[[", "name"))
    findDiamPath <- pathNames==chosenPath
    if (!any(findDiamPath)) stop("No Diameter Path / Named path Found")
    if (length(which(findDiamPath)) > 1) warning("Multiple Diameter paths found. First used!")
    findDiamPath <- which(findDiamPath)[1]
    ev <- eigenvecs(object)
    if (!is.null(project2D)) {
        ev <- applyMVP(ev[,useDims], project2D)
        useDims <- 1:2
    }
    dpCoords <- ev[object@paths[[findDiamPath]][["path"]],useDims]
    if (length(useDims) == 2) lines(dpCoords, lwd=diamWidth, col=diamCol)
    else lines3d(dpCoords, lwd=diamWidth, col=diamCol)
})


##' Returns loadings from PCA object
##'
##' Returns loadings from PCA object
##' @title PCA loadings
##' @param x Single Cell Dataset object
##' @return Annotated matrix of loadings
##' @author Wajid Jawaid
##' @export
setMethod("loadings", "PCA", function(x) x@rotation)

##' Returns eigen values for Diffusion Map object
##'
##' In the interest of time the diffusion map algorithm only calculates
##' a given number of largest eigen values and vectors using rARPACK.
##' Therefore not all eigen values will be returned.
##' Returns eigen values for Diffusion Map object
##' @title Diffusion Map eigen values
##' @param object DiffusionMap object
##' @return Numeric vector of eigen values for 
##' @author Wajid Jawaid
##' @export
setMethod("eigenvals", "DiffusionMap", function(object) object@eigenvalues)

##' Returns eigen vectors for Diffusion Map object
##'
##' In the interest of time the diffusion map algorithm only calculates
##' a given number of largest eigen values and vectors using rARPACK.
##' Therefore not all eigen vectors will be returned.
##' Returns eigen vectors for Diffusion Map object
##' @title Diffusion Map eigen vectors
##' @param object DiffusionMap object
##' @return Annotated matrix of eigen vectors, with eigen vectors in columns.
##' @author Wajid Jawaid
##' @export
setMethod("eigenvecs", "DiffusionMap", function(object) object@eigenvectors)

##' Initialise a SCD object
##'
##' Initialise a SCD object
##' @title To initialize a SCD object
##' @param .Object Single Cell Dataset object to be created
##' @param experimentType Either "RNAseq" or "qPCR"
##' @param assayData Annotated matrix of gene expression data
##' @param featureData Annotated data frame of gene information. One column must match rownames
##' in \code{assayData} the gene expression annotated matrix.
##' @param phenoData Annotated data frame of phenotype information. One column must match the colnames in \code{assayData} the gene expression annotated matrix.
##' @param spike Default NULL. Include spike-in data e.g. ERCC-92.
##' @param qc Default NULL. Include quality control counts from HTseq for example.
##' @param ... Parameters passed to other Methods
##' @return Returns a newly initialised SCD object
##' @author Wajid Jawaid
setMethod("initialize", "SCD", function(.Object, experimentType, assayData,
                                        featureData, phenoData, spike = NULL, qc = NULL, ...) {
    .Object@useFilter <- TRUE
    geneIds <- rownames(assayData)
    cellIds <- colnames(assayData)
    data <- prepData(geneIds, featureData, "genodata")
    featureData <- new("AnnotatedDataFrame", data=data,
                          varMetadata=data.frame(colnames(data)))
    data <- prepData(cellIds, phenoData, "phenodata")
    phenoData <- new("AnnotatedDataFrame", data=data,
                     varMetadata=data.frame(colnames(data)))
    if (!is.null(spike)) {
        cellsInData <- colnames(spike)
        missingCells <- setdiff(cellIds, cellsInData)
        if (length(missingCells) != 0) {
            warning("Spike-in counts only available for a subset of cells.")
            missingCells <- matrix(nrow=nrow(spike), ncol=length(missingCells),
                                   dimnames=list(rownames(spike), missingCells))
            spike <- cbind(spike, missingCells)
        }
        data <- spike[,cellIds]
        attr(data, "cellsInData") <- cellsInData
        .Object@spike <- data
    }
    if (!is.null(qc)) {
        cellsInData <- colnames(qc)
        missingCells <- setdiff(colnames(assayData), cellsInData)
        if (length(missingCells != 0)) {
            warning("QC counts only available for a subset of cells.")
            missingCells <- matrix(nrow=nrow(qc), ncol=length(missingCells),
                                   dimnames=list(rownames(qc), missingCells))
            qc <- cbind(qc, missingCells)
        }
        data <- qc[,cellIds]
        attr(data, "cellsInData") <- cellsInData
        .Object@qcCounts <- data
    }
    .Object@annotation = experimentType
    callNextMethod(.Object, phenoData=phenoData, featureData = featureData, ...)
})


setMethod("initialize", "SCESet",
  function(.Object, exprsData = NULL, countData = NULL, tpmData = NULL, 
           fpkmData = NULL, cpmData = NULL, phenoData = NULL,
           featureData = NULL, experimentData = NULL,
           is_exprsData = NULL, cellPairwiseDistances = dist(vector()), 
           featurePairwiseDistances = dist(vector()),
           lowerDetectionLimit = NULL, logExprsOffset = NULL, assayData = NULL,
           ...) {
      if (is.null(assayData)) {
          have.data <- NULL
          for (dataname in c("countData", "tpmData", "cpmData", "fpkmData", 
                             "exprsData")) {
              eData <- get(dataname)
              if (!is.null(eData)) {
                  if (!is.null(have.data)) {
                      warning(sprintf("'%s' provided, '%s' will be ignored", 
                                      have.data, dataname))
                      assign(dataname, NULL)
                  }
                  else {
                      assign(dataname, as.matrix(eData))
                      have.data <- dataname
                  }
              }
          }
          if (is.null(have.data)) {
              ## stop("one set of expression values should be supplied")
              countData <- matrix()
              have.data <- "countData"
          }
          if (!is.null(is_exprsData)) {
              if (have.data != "exprsData") {
                  warning(sprintf("'%s' provided, 'is_exprsData' will be ignored", 
                                  have.data))
                  is_exprsData <- NULL
              }
              else {
                  is_exprsData <- as.matrix(is_exprsData)
              }
          }
          if (is.null(logExprsOffset)) {
              logExprsOffset <- 1
              if (have.data != "countData") {
                  warning("'logExprsOffset' should be set manually for non-count data")
              }
          }
          if (is.null(lowerDetectionLimit)) {
              lowerDetectionLimit <- 0
              if (have.data == "exprsData") {
                  warning("'lowerDetectionLimit' should be set manually for log-expression values")
              }
          }
          if (have.data == "countData") {
              exprsData <-
                  ## scater:::.compute_exprs(countData, size_factors = colSums(countData), 
                  ##                         log = TRUE, sum = FALSE,
                  ##                         logExprsOffset = logExprsOffset)
                  countData
              ## dimnames(exprsData) <- dimnames(countData)
          }
          else if (have.data != "exprsData") {
              exprsData <- log2(get(have.data) + logExprsOffset)
          }
          assaydata <- assayDataNew("environment", exprs = exprsData)
          if (!is.null(is_exprsData)) 
              assaydata[["is_exprs"]] <- is_exprsData
          if (!is.null(tpmData)) 
              assaydata[["tpm"]] <- tpmData
          if (!is.null(fpkmData)) 
              assaydata[["fpkm"]] <- fpkmData
          if (!is.null(countData)) 
              assaydata[["counts"]] <- countData
          if (!is.null(cpmData)) 
              assaydata[["cpm"]] <- cpmData
      } else assaydata <- assayData
      if (is.null(phenoData)) 
          phenoData <- annotatedDataFrameFrom(exprsData, byrow = FALSE)
      if (is.null(featureData)) 
          featureData <- annotatedDataFrameFrom(exprsData, byrow = TRUE)
      expData_null <- new("MIAME", name = "<your name here>", lab = "<your lab here>", 
                          contact = "<email address>", title = "<title for this dataset>", 
                          abstract = "An SCESet", url = "<your website here>", 
                          other = list(notes = "This dataset created from ...", 
                                       coauthors = c("")))
      if (!is.null(experimentData)) {
          if (is(experimentData, "MIAME")) 
              expData <- experimentData
          else {
              expData <- expData_null
              warning("'experimentData' is not an 'MIAME' object, setting to an empty object")
          }
      }
      else {
          expData <- expData_null
      }
      storageMode(assaydata) <- "lockedEnvironment"
      .Object@assayData <- assaydata
      .Object@phenoData <- phenoData
      .Object@featureData <- featureData
      .Object@experimentData <- expData
      .Object@cellPairwiseDistances <- cellPairwiseDistances
      .Object@featurePairwiseDistances <- featurePairwiseDistances
      .Object@lowerDetectionLimit <- lowerDetectionLimit
      .Object@logExprsOffset <- logExprsOffset
      .Object@featureControlInfo <- AnnotatedDataFrame()
      .Object@protocolData <- AnnotatedDataFrame(
          data.frame(labelDescription = rownames(phenoData),
                     row.names = rownames(phenoData)))
      .Object@useForExprs <- "exprs"
      .Object
      ## callNextMethod(.Object, assayData = assayData, annotation=experimentType,
      ##                featureData=featureData, phenoData=phenoData, ...)
  })

##' Retrieve gene expression data from Single Cell Dataset object
##'
##' Retrieve gene expression data from Single Cell Dataset object
##' @title Retrieve gene expression data from Single Cell Dataset object
##' @param object SCD
##' @return matrix of expression values
##' @author Wajid Jawaid
##' @export
setMethod("exprs", "SCD", function(object) {
    fD <- pData(featureData(object))
    if (filterGene(object)) {
        x <- assayDataElement(object, "exprs")[fD[,"included"],]
    } else x <- assayDataElement(object, "exprs")
    x <- subsetData(x, pData(phenoData(object)), saveFilters(object))
})

##' Retrieve size factors
##'
##' Retrieve size factors
##' @title Retrieve size factors
##' @param object SCD
##' @param type Default "bio". For "spike-ins" use "tech"
##' @return Returns matrix of size factors
##' @author wj241
setMethod("sf", "SCD", function(object, type="bio"){
    if (is.null(object@technicalNoise$sizeFactors)) stop("Run QC first.")
    x <- object@technicalNoise$sizeFactors
    if (tolower(type) == "bio") {
        x <- x$sf.data
    } else if (tolower(type) == "tech") {
        x <- x$sf.data.ercc
    } else if (tolower(type) == "both") {
        
    } else stop("Incorrect type argument.")
    x <- matrix(x, nrow=1, dimnames=list(NULL, names(x)))
    x <- subsetData(x, pData(phenoData(object)), saveFilters(object))
    return(x)
})

##' Retrieve gene counts from Single Cell Dataset object
##'
##' Retrieve gene counts from Single Cell Dataset object
##' @title Retrieve counts data from Single Cell Dataset object
##' @param object SCD
##' @param ... Additional parameters
##' @return matrix of counts
##' @author Wajid Jawaid
##' @export
setMethod("counts", "SCD", function(object, ...) {
    x <- assayDataElement(object, "counts")
    fD <- pData(featureData(object))
    if (filterGene(object) && useFilter(object)) x <- x[rownames(fD)[fD[,"included"]],]
    x <- subsetData(x, pData(phenoData(object)), saveFilters(object))
})

##' Retrieve spike in counts
##'
##' Retrieve spike in counts
##' @title Retrieve spike in counts
##' @param object SCD
##' @return Matrix of spike counts
##' @author Wajid Jawaid
setMethod("spikes", "SCD", function(object) {
    x <- object@spike
    x <- subsetData(x, pData(phenoData(object)), saveFilters(object))

    ## cellsInData <- attr(eD, "cellsInData")
    ## pD <-rownames(pData(object))
    ## if (filterQC(object)) eD <- eD[, pD]
    ## if (!useFilter(object)) return(eD)
    ## eD <- eD[, pD]
    ## attr(eD, "cellsInData") <- cellsInData
})

##' Retrieve QC in counts
##'
##' Retrieve QC in counts
##' @title Retrieve QC in counts
##' @param object SCD
##' @return Matrix of QC counts
##' @author Wajid Jawaid
setMethod("qc", "SCD", function(object) {
    x <- object@qcCounts
    x <- subsetData(x, pData(phenoData(object)), saveFilters(object))
})

##' Dimensions of SCD
##'
##' Dimensions of SCD
##' @title Dimensions of SCD
##' @param x Single cell data set
##' @return dimension of object
##' @author Wajid Jawaid
##' @export
setMethod("dim", "SCD", function(x) {
    dimN <- c(nrow(fData(x)), nrow(pData(x)))
    cat("Use filter is:\n")
    cat("\tMain filter: ", useFilter(x), "\n")
    cat("\tGene filter: ", filterGene(x), "\n")
    cat("\tCell filter: ", filterCell(x), "\n")
    cat("Dimensions: ", paste(dimN, collapse = "x"), "\n")
    invisible(dimN)
})

##' Retrieve gene information table from Single Cell Dataset object
##'
##' Retrieve gene information table from Single Cell Dataset object
##' This table should include the ENSEMBL ids matching the gene names in the gene expression
##' matrix along with the corresponding human readable gene names.
##' @title Retrieve gene information table from Single Cell Dataset object
##' @param object SCD
##' @return Data frame of ENSEMBL ids and corresponding human readable gene names
##' @author Wajid Jawaid
##' @export
setMethod("fData", "SCD", function(object) {
    fD <- pData(featureData(object))
    if (useFilter(object) && filterGene(object)) return(orderFrame(fD))
    else return(fD)
})

##' Retrieve phenotype information table from Single Cell Dataset object
##'
##' Retrieve phenotype information table from Single Cell Dataset object
##' This table should include the sample/cell ID, sequencing chip and index tag, embryo stage
##' of the sample, the embryo ID and other data - for now including FACS index data
##' @title Retrieve phenotype data
##' @param object SCD
##' @return Data frame of phenotype data
##' @author Wajid Jawaid
##' @export
setMethod("pData", "SCD", function(object) {
    pD <- pData(phenoData(object))
    if (filterQC(object)) {
        pD <- pD[pD[,"passedQC"],]
        pD <- pD[,-match("passedQC", colnames(pD)), drop=FALSE]
    }   
    if (useFilter(object) && filterCell(object)) {
        return(orderFrame(pD))
    } else return(pD)
})

##' To use filter or not
##'
##' To use filter or not. 
##' @title Predicate indicating whether overall filtering is active
##' @param object Single Cell Dataset
##' @return Logical.
##' @author Wajid Jawaid
##' @export
setMethod("useFilter", "SCD", function(object) {
    object@useFilter
})


##' To use gene filter or not
##'
##' To use filter or not.
##' @title Predicate indicating whether gene filtering is active
##' @param object Single Cell Dataset
##' @return Updated Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setMethod("filterGene", "SCD", function(object) {
    object@filterGene
})

##' To use filter or not
##'
##' To use filter or not. Replaces "noprocessing" argument from previous version
##' @title Predicate indicating whether cell filtering is active.
##' @param object Single Cell Dataset
##' @return Logical
##' @author Wajid Jawaid
##' @export
setMethod("filterCell", "SCD", function(object) {
    object@filterCell
})

##' To filter out cells that failed QC.
##'
##' To filter out cells that failed QC. It is highly recommended that this
##' is not disabled as size factor and technical noise estimation may fail.
##' @title Predicate indicating whether filtering of cells that failed QC is active
##' @param object Single Cell Dataset
##' @return Logical
##' @author Wajid Jawaid
##' @export
setMethod("filterQC", "SCD", function(object) {
    object@filterQC
})

##' Set whether to use filter or not
##'
##' Set whether to use filter or not
##' @title Set whether to use filter or not
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setReplaceMethod("useFilter", signature(object="SCD", value="logical"),
                 function(object, value) {
                     object@useFilter <- value
                     return(object)
                 })

##' Save the size factors
##'
##' Save the size factors
##' @title Save the size factors
##' @param object SCD
##' @param type Can be "bio" or "tech"
##' @param value Matrix of values
##' @return UYpdated SCD object
##' @author wj241
setReplaceMethod("sf", signature(object="SCD", type="character", value="matrix"),
                 function(object, type, value) {
                     if (tolower(type) == "bio") {
                         object@technicalNoise$sizeFactors$sf.data <- value
                     } else if (tolower(type) == "tech") {
                         object@technicalNoise$sizeFactors$sf.data.ercc <- value
                     } else stop("Incorrect type argument.")
                     
                 })

##' Set whether to use gene filter or not
##'
##' Set whether to use gene filter or not
##' @title Set gene filter.
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setReplaceMethod("filterGene", signature(object="SCD", value="logical"),
                 function(object, value) {
                     object@filterGene <- value
                     return(object)
                 })

##' Set whether to use cell filter or not
##'
##' Set whether to use cell filter or not
##' @title Set cell filter.
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setReplaceMethod("filterCell", signature(object="SCD", value="logical"),
                 function(object, value) {
                     object@filterCell <- value
                     return(object)
                 })

##' Set whether to use QC filter or not
##'
##' Set whether to use QC filter or not
##' @title Set QC filter.
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setReplaceMethod("filterQC", signature(object="SCD", value="logical"),
                 function(object, value) {
                     object@filterQC <- value
                     return(object)
                 })

##' Add additional phenotype data
##'
##' Add additional phenotype data
##' @title Add additional phenotype data
##' @inheritParams addPheno
##' @return Updated Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setMethod("addPheno", "SCD", function(object, value, noprocessing=FALSE) {
    origFilter <- saveFilters(object)
    useFilter(object) <- !noprocessing
    if (!is.data.frame(value)) stop("Can only add dataframe object to pData")
    if (nrow(value) != nrow(pData(object)))
        stop("Length of new data does not match number of rows in phenoData.")
    useFilter(object) <- FALSE
    filterQC(object) <- FALSE
    phenoData <- pData(object)
    phenoData[,colnames(value)] <- NA
    idx <- match(rownames(value), rownames(phenoData))
    if (any(is.na(idx))) stop("Rownames do not match row names in phenoData")
    for (i in colnames(value)) {
        phenoData[idx, i] <- value[,i]
        lev <- levels(value[,i])
        if (is.factor(value[,i])) {
            if (any(is.na(phenoData))) {
                phenoData[,i] <- factor(phenoData[,i], levels=1:(length(lev)+1),
                                        labels=c(lev,"NA"))
            } else phenoData[,i] <- factor(phenoData[,i], levels=1:length(lev), labels=lev)
        }
    }
    pData(object) <- phenoData
    object <- restoreFilters(object, origFilter)
    return(object)   
})

##' Exclude selected cells
##'
##' Function to remove selected cells from analyses
##' @title Exclude selected cells
##' @inheritParams excludeCells
##' @return modified Single Cell Data object
##' @author Wajid Jawaid
setMethod("excludeCells", "SCD", function(object, regexpr = NULL, cellNames = NULL,
                                          invert=FALSE, reset=TRUE) {
    origFilter <- saveFilters(object)
    object <- deactivateFilters(object, includeQC = TRUE)
    pD <- pData(object)
    if (reset) pD[,"included"] <- rep(TRUE, nrow(pD))
    if (is.null(cellNames)) {
        xcludedCells <- grepl(regexpr, rownames(pD))
        if (invert) xcludedCells <- !xcludedCells
        pD[xcludedCells, "included"] <- FALSE        
    } else {
        nm <- match(cellNames, rownames(pD))
        if (any(is.na(nm))) stop("All cellNames names are not in the dataset")
        if (invert) nm <- -nm
        if (!is.null(cellNames)) pD[nm, "included"] <- FALSE
    }
    pData(object) <- pD
    object <- restoreFilters(object, origFilter)
    return(object)
})

##' Mark cells that have failed QC
##'
##' Mark cells that have failed QC
##' @title Mark cells that have failed QC
##' @param object Single Cell Dataset
##' @param cellNames Character vector of cell names
##' @param reset Default TRUE. Sets all others as PASSED.
##' @return Updated Single Cell Data object.
##' @author wj241
setMethod("markFailedQC", "SCD", function(object, cellNames, reset = TRUE) {
    svFil <- saveFilters(object)
    object <- deactivateFilters(object, includeQC = TRUE)
    pD <- pData(object)
    if (reset) pD[,"passedQC"] <- TRUE
    pD[cellNames,"passedQC"] <- FALSE
    pData(object) <- pD
    return(object)
})

##' To select only highly variable genes for further analysis
##'
##' To select only highly variable genes for further analysis. Give either genes to
##' be included or genes to be excluded BY ENSMBL id only. Do not give both.
##' @title To select only highly variable genes for further analysis
##' @inheritParams selectVariableGenes
##' @return Returns modififed Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setMethod("selectVariableGenes", "SCD",
          function(object, includeGenes=NULL, excludeGenes=NULL, reset=FALSE) {
              if (is.null(includeGenes) && is.null(excludeGenes)) {
                  if (!reset)
                      stop("Must provide either includeGene or excludeGenes.")
                  
              } else if (!is.null(includeGenes) & !is.null(excludeGenes)) {
                  stop("Must provide either includeGene or excludeGenes. Not both.")
              }
              origFilter <- saveFilters(object)
              object <- deactivateFilters(object)
              genoData <- fData(object)
              if (any(is.na(match(c(includeGenes, excludeGenes), rownames(genoData))))) {
                  stop("Gene names do not match Ensembl ids in genoData.")
              }

              if (is.null(includeGenes) && is.null(excludeGenes)) {
                  genoData[,"included"] <- TRUE
              } else if (is.null(excludeGenes)) {
                  if (reset) genoData[,"included"] <- FALSE
                  genoData[includeGenes,"included"] <- TRUE
              } else if (is.null(includeGenes)) {
                  if (reset) genoData[,"included"] <- TRUE
                  genoData[excludeGenes, "included"] <- FALSE
              }
              fData(object) <- genoData
              object <- restoreFilters(object, origFilter)
              return(object)
          })

##' Calculate PCA for Single Cell Dataset
##'
##' Calculate PCA for Single Cell Dataset
##' @title Principal Component Analysis
##' @param object Single Cell Dataset object
##' @param ... Parameters passed to \code{prcomp}
##' @return Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setMethod("runPCA", "SCD", function(object, ...) {
    x <- prcomp(t(exprs(object)), ...)
    object@pca <- new("PCA",
                      standardDeviation = x$sdev,
                      rotation = x$rotation,
                      center = x$center,
                      eigenvalues = x$sdev / sum(x$sdev),
                      eigenvectors =x$x,
                      scaled = x$scale)
    return(object)
})

##' Retrieve PCA data
##'
##' Retrieve PCA data
##' @title Retrieves PCA object
##' @inheritParams getPCA 
##' @return Returns PCA dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("getPCA", "SCD", function(object, ...) object@pca)

##' Retrieve diffusion map object
##'
##' Retrieve diffusion map object
##' @title Retrieves diffusion map object
##' @inheritParams getDiffMap
##' @return Returns DiffMap dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("getDiffMap", "SCD", function(object, ...) object@diffMap)

##' Retrieve t-SNE object
##'
##' Retrieve t-SNE object
##' @title Retrieves t-SNE object
##' @param object Single Cell Dataset object
##' @param ... Additional parameters
##' @return Returns t-SNE dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("getTSNE", "SCD", function(object, ...) object@tsne)

##' Retrieve Isomap object
##'
##' Retrieve Isomap object
##' @title Retrieves Isomap object
##' @param object Single Cell Dataset object
##' @param ... Additional parameters
##' @return Returns Isomap dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("getIsomap", "SCD", function(object, ...) object@isomap)

##' Retrieve diffusion map object
##'
##' Retrieve diffusion map object
##' @title Perform t-distributed stochastic neighbour embedding
##' @param object Single Cell Data Set object
##' @param ndims Number of dimensions for new embedding
##' @param perplexity Perplecity parameter
##' @param theta Theta parameter
##' @param pca Logical, Perform initial PCA reduction
##' @param seed Set random seed
##' @param verbose To show output
##' @param use_dist Default FALSE. Calculate a distance metric on the
##' data and use this in the Rtsne function.
##' @param dist_fun Default NULL. If NULL will use (1 - spearman_correlation(data)) / 2.
##' Otherwise provide a function to use on the data matrix.
##' @param distmat Deafult NULL. Optionally give the distance matrix to use to speed up analysis.
##' @param ... Additional parameters to pass to Rtsne()
##' @return Returns an upadate Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setMethod("runTSNE", "SCD", function(object, ndims = 3, perplexity = 30, theta = 0.5,
                                     pca = FALSE, seed = NULL, verbose = FALSE,
                                     use_dist = FALSE, dist_fun = NULL, distmat = NULL, ...) {
    data <- exprs(object)
    if (use_dist) {
        if (is.null(distmat)) {
            if (is.null(dist_fun)) {
                data <- (1 - cor(data, method="spearman")) / 2
            } else dist_fun(data)
        } else if (all(dim(distmat) == rep(ncol(data), 2))) {
            data <- distmat
        } else stop("Distance matrix not consistent with data!")
    }
    if (is.null(seed)) {
        seed <- Sys.time()
        seed <- as.numeric(seed)
    }
    set.seed(seed)
    tsne <- Rtsne(t(data), dims = ndims, perplexity = perplexity, theta = theta, pca = pca,
                  verbose = verbose, is_distance = use_dist, ...)
    rownames(tsne$Y) <- rownames(pData(object))
    object@tsne <- new("TSNE",
                       eigenvectors = tsne$Y,
                       theta = tsne$theta,
                       perplexity = tsne$perplexity,
                       N = tsne$N,
                       origD = tsne$origD,
                       seed = seed,
                       nnError = tsne$itercosts,
                       algorithm = "bh-tsne"
                       )
    return(object)
})

##' Diffusion map dimensionality reduction
##'
##' Diffusion map dimensionality reduction as performed by Florian Buettner / Fabian Theis
##' Diffusion map dimensionality reduction
##' @title Diffusion Maps
##' @inheritParams diffuse
##' @return Single Cell Data object with diffusion map
##' @author Wajid Jawaid
##' @importFrom roots diffuseMat
##' @export
setMethod("diffuse", "SCD", function(object, ndims = 4, nn = 0.2, sigma = NULL,
                                     removeFirst = TRUE, useARPACK = TRUE,
                                     distfun = NULL) {
    data <- exprs(object)
    decomp <- diffuseMat(data, ndims=ndims, nsig=nn, sigmas=sigma, removeFirst=removeFirst,
                          useARPACK=useARPACK, distfun = distfun)
    object@diffMap <- new("DiffusionMap",
                          eigenvalues = decomp$values,
                          eigenvectors = decomp$vectors,
                          numberConverged = decomp$nconv,
                          numberIterations = decomp$niter,
                          sigma = sigma,
                          nn = decomp$nn,
                          usedARPACK = decomp$usedARPACK)
    cat("Done.\n")
    return(object)
})

##' Returns the minimum spanning tree
##'
##' Returns the minimum spanning tree
##' @title Retrieves the minimum spanning tree from the \code{graph} slot of a
##' Dimensionality reduced object
##' @param object Dimensionality reduced object
##' @param reducedDim The dimensionality reduction technique from which the MST is required
##' @return Minimum spanning tree as an igraph object
##' @author Wajid Jawaid
##' @export
setMethod("getMST", "SCD", function(object, reducedDim) getMST(slot(object,reducedDim)))

##' Given a vector of cells, this will add a path to the paths slot.
##'
##' Given a vector of cells, this will add a path to the paths slot.
##' @title Store path
##' @param object Dimensionality reduced object
##' @param pathName Name given to the path
##' @param pathDescription A description of the path
##' @param path A numeric vector of cells on the path
##' @return Dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("addPath", "reducedDim", function(object, pathName, pathDescription, path) {
    newPath <- list(name = pathName, decription = pathDescription, path = path)
    usedNames <- do.call(c, lapply(object@paths, "[[", "name"))
    if (any(usedNames == pathName)) cat("Path name already in use. Duplicate entry created!\n")
    paths <- c(object@paths, list(newPath))
    names(paths) <- c(usedNames, pathName)
    object@paths <- paths
    return(object)
})

##' Given a vector of cells, this will add a path to the paths slot.
##'
##' Given a vector of cells, this will add a path to the paths slot.
##' @title Store path
##' @param object Single cell dataset or dimensionality reduced object
##' @param reducedDim A dimensionality reduced object
##' @param pathName Name given to the path
##' @param pathDescription A decription of the path 
##' @param path A numeric vector of cells on the path
##' @param ... Parameters passed to the Method \code{addPath}
##' @return Dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setMethod("addPath", "SCD", function(object, reducedDim, pathName, pathDescription,
                                     path, ...) {
    slot(object, reducedDim) <- addPath(slot(object, reducedDim), pathName = pathName,
                                        pathDescription = pathDescription, path = path)
    return(object)
})

##' Finds and saves path between two cells
##'
##' Finds and saves path between two cells
##' @title Find path between two selected cells
##' @param object A single cell dataset object
##' @param reducedDim Reduced diensionality object
##' @param from Numeric value representing the starting cell of path
##' @param to Numeric value representing the ending cell of path
##' @param pathName Name to give path
##' @param pathDescription Description of path
##' @return Altered object of same class
##' @author Wajid Jawaid
##' @export
setMethod("findPath", "SCD", function(object, reducedDim, from, to, pathName,
                                      pathDescription) {
    sp <- get.shortest.paths(getMST(object, reducedDim), from = from, to = to,
                             weights = NA)$vpath[[1]]
    addPath(object, reducedDim= reducedDim, pathName = pathName,
            pathDescription = pathDescription, path = sp)
})

##' Removes paths from object
##'
##' Removes paths from object
##' @title Removes paths from object
##' @param reducedDim The reduced dimension object to use
##' @inheritParams removePath
##' @return SCD
##' @author Wajid Jawaid
setMethod("removePath", "SCD", function(object, reducedDim, paths) {
    slot(object, reducedDim) <- removePath(slot(object, reducedDim), paths)
    return(object)
})

##' Removes paths from object
##'
##' Removes paths from object
##' @title Removes paths from object
##' @inheritParams removePath
##' @return SCD
##' @author Wajid Jawaid
setMethod("removePath", "reducedDim", function(object, paths) {
    ind <- paths
    paths <- getPaths(object)
    if (is.numeric(ind)) {
        if (!all(ind %in% 1:length(paths))) stop("Path indexes not present.")
        paths[ind] <- NULL
    } else if (!is.numeric(ind) && is.character(ind)) {
        if (ind == "ALL") {
            paths <- list()
        } else {
            names(paths) <- sapply(paths, "[[", "name")
            paths[ind] <- NULL
        }
    } else stop("paths argument must be numeric or character.")
    object@paths <- paths
    return(object)
})

##' Generate minimum spanning true
##'
##' General minimum spanning true
##' @title Generate minimum spanning true
##' @param object Single Cell Dataset object
##' @param wR Name of slot containing Dimensionality Reduced object
##' @param weighted Default TRUE. Passed to igraph.
##' @param useDims Deafult is NULL. Sets dimensions to use.
##' @return Returns single cell data object
##' @author Wajid Jawaid
##' @export
setMethod("mst", "SCD", function(object, wR, weighted=TRUE, useDims=NULL) {
    slot(object, wR) <- mst(slot(object, wR), weighted=TRUE, useDims=useDims)
    return(object)
})

##' General Plot function
##'
##' Plot function
##' @title General plot function
##' @param x Single Cell Dataset object
##' @param reduceMethod Dimensionality reduction method to choose. At the moment can be either "pca"
##' or "diffMap"
##' @param dims Number of dimensions
##' @param colorBy Which column of the annotated phenotype data frame to use for colouring points
##' @param useDims Dimensions to be plotted
##' @param logContinuous Take log for continuous variables
##' @param opacity Default is NULL. Set opacity for each point in pData and pass column
##' name as character vector. Values in the range 0:255.
##' @param ... Parameters passed to \code{plot3d} or \code{plot}
##' @return Returns the requested plot
##' @author Wajid Jawaid
##' @inheritParams plot3
##' @export
setMethod("plot", signature(x="SCD", y = "ANY"),
          function(x, reduceMethod = c("pca", "diffMap","tsne","isomap"), dims = c(2,3),
                   colorBy = "embryoStage", useDims = NULL, logContinuous = TRUE,
                   opacity = NULL, ...) {
              reduceMethod <- match.arg(reduceMethod)
              if (length(dims) > 1) dims <- dims[1]
              if (!is.null(useDims)) dims <- length(useDims)
              ind <- match(rownames(eigenvecs(slot(x, reduceMethod))), rownames(pData(x)))
              pCol <- pData(x)[ind, colorBy]
              attr(pCol, "title") <- colorBy
              if (!is.null(opacity)) opacity <- pData(x)[,opacity]
              if (is.numeric(pCol)) {
                  if (logContinuous) pCol <- log(pCol + 1)
              }
              if (dims == 3) {
                  if (is.null(useDims)) useDims <- 1:3
                  plotID <- plot3(slot(x, reduceMethod), colorBy = pCol,
                                  useDims = useDims, opacity = opacity, ...)
              }
              if (dims == 2) {
                  if (is.null(useDims)) useDims <- 1:2
                  plotID <- plot3(slot(x, reduceMethod), colorBy = pCol,
                                  useDims = useDims, opacity = opacity, ...)
              }
              invisible(plotID)
          })

##' Used for selecting cells
##'
##' Used for selecting cells. Generates plot and then allows user to choose cells. Cells
##' are marked on the plot and are returned.
##' @title Select cells function
##' @param object Single Cell Dataset object
##' @param reduceMethod Dimensionality reduction method to choose. At the moment can be either "pca" or "diffMap"
##' @param useDims Dimensions to be plotted
##' @param mark Logical. Default "TRUE" marks the cells with names on the plot
##' @param returnInd Logical. Default "FALSE" returns a numeric vector of named indices
##' @param label.cex Default 0.3. Set for label size
##' @param pdf Defaut is NULL. Else will generate a pdf.
##' @param ... Further parameters passed to the plot function
##' @return Returns a plot with marked cells and a vector of either named indices or just
##' names if "returnId" is set to "FALSE".
##' @author Wajid Jawaid
##' @export
setGeneric("selectCells", function(object, reduceMethod, useDims = 1:3, mark = TRUE,
                                   returnInd = TRUE, label.cex=.3, pdf=NULL, ...)
    standardGeneric("selectCells"))

##' Used for selecting cells
##'
##' Used for selecting cells. Generates plot and then allows user to choose cells. Cells
##' are marked on the plot and are returned.
##' @title Select cells function
##' @inheritParams selectCells
##' @return Returns a plot with marked cells and a vector of either named indices or just
##' names if "returnId" is set to "FALSE"".
##' @author Wajid Jawaid
##' @export
setMethod("selectCells", signature(object="SCD"),
          function(object, reduceMethod, useDims = 1:3, mark = TRUE, returnInd = FALSE,
                   label.cex=.3, pdf=NULL, ...) {
              if (length(useDims)==3) {
                  open3d()
                  plotID <- plot(object, reduceMethod = reduceMethod, useDims = useDims,
                                 returnID = TRUE, ...)
                  readline("Put plot into desired orientation. Press Enter when ready...")
                  chosen <- selectpoints3d(plotID[1], value = FALSE)[,2]
                  cellName <- colnames(exprs(object))[chosen]
                  cellPos <- eigenvecs(slot(object, reduceMethod))[,useDims]
                  cellPos <- cellPos[match(cellName, rownames(cellPos)),,drop=FALSE]
                  if (mark) text3d(x=cellPos[,1], y = cellPos[,2], z=cellPos[,3],
                                   texts = cellName)
                  names(chosen) <- cellName
                  if (!is.null(pdf)) rgl.postscript(pdf, fmt="pdf")
              } else if (length(useDims)==2) {
                  plot(object, reduceMethod = reduceMethod, useDims = useDims, ...)
                  readline("Put plot into desired orientation. Press Enter when ready...")
                  sel <- do.call(cbind, locator(n=2))
                  cellPos <- eigenvecs(slot(object, reduceMethod))[,useDims]
                  xVals <- sort(sel[,1])
                  yVals <- sort(sel[,2])
                  xIn <- cellPos[,1] > xVals[1] & cellPos[,1] < xVals[2]
                  yIn <- cellPos[,2] > yVals[1] & cellPos[,2] < yVals[2]
                  chosen <- which(xIn & yIn)
                  cellPos <- cellPos[chosen,,drop=FALSE]                  
                  text(cellPos[,1], cellPos[,2], labels=rownames(cellPos), cex=label.cex)
                  cellName <- rownames(cellPos)
                  if (!is.null(pdf)) {
                      pdf(pdf)
                      plot(object, reduceMethod = reduceMethod, useDims = useDims, ...)
                      text(cellPos[,1], cellPos[,2], labels=rownames(cellPos), cex=label.cex)
                      dev.off()
                  }
              }
              if (returnInd) cellName <- chosen
              if (!is.null(pdf)) {system(paste(fileOpen, pdf, sep="")); dev.off()}
              return(cellName)
          })

##' Normalises and performs technical noise analysis for single cell RNAseq(Brennecke et al.)
##'
##' Normalises and performs a technical noise analysis for single cell RNAseq.
##' Give raw counts. Data is normalised using size-factor normalisation identical to DESeq2,
##' then a technical noise analysis is perfromed as described in Brennecke et al. Normalised
##' expression values are log10 normalised.
##' @title Technical noise analysis for single cell RNAseq
##' @param object SCD class
##' @param useERCC Default TRUE. To use ERCC for technical noise estimation.
##' @param cvThresh Defalt 0.3. The CV threshold to use for choosing cells to fit. See quant.
##' @param quant Default 0.8. The quantile of data to fit so that linear model is not fitted to
##' plateau at lower gene expression values. See Breenecke et al. for further details
##' @param minBiolDisp Default 0.25. Assumed biological dispersion. See paper.
##' @param fdr Deafult .1. False Discovery Rate for chi-squared test.
##' @param data.lengths Gene lengths in kilobases for PK "per kilobase normalisation"
##' @param ercc.lengths ERCC lengths in kilobases for PK "per kilobase normalisation"
##' @param meanForFit Provide a minimum mean for fit. (OPTIONAL)
##' @return Returns updated SCD
##' @author Wajid Jawaid
##' @export
setMethod("techVar", signature(object = "SCD"),
          function(object, useERCC = TRUE, cvThresh=.3, quant=.8,
                   minBiolDisp=.5^2, fdr=.1, data.lengths = NULL,
                   ercc.lengths = NULL, meanForFit=NULL) {

              origFilter <- saveFilters(object)
              filterGene(object) <- FALSE
              data.ercc <- NULL
              if (useERCC) data.ercc <- spikes(object)
              retVal <- techVarSub(object, data.ercc, cvThresh, quant,
                                   minBiolDisp, fdr, data.lengths, ercc.lengths,
                                   meanForFit)
              object <- selectVariableGenes(object, includeGenes = retVal[["highVarGenes"]],
                                            reset = TRUE)
              object@technicalNoise <- retVal
              object <- restoreFilters(object, origFilter)
              object
})

##' Perform QC on data
##'
##' Perform QC on data
##' @title Perform QC on data
##' @param object SCD
##' @param selectedCells Default "ALL". If data not available on all cells give subset.
##' @param cutoffs Default 2e5, 0.2, 0.2, 0, 0, 1, 0.
##' Set thresholds for:
##' \itemize{
##'  \item{number of mapped reads (greater than threshold), }
##'  \item{ratio of genes to total number of reads (greater than threshold), }
##'  \item{ratio of mitchondrial genes to mitochondrial + nuclear genes (less than threshold),}
##'  \item{number of nuclear genes (greater than threshold), }
##'  \item{number of genes with >10 reads per million (greater than threshold),}
##'  \item{ercc:mapped (less than threshold) and}
##'  \item{nuclearGenes:mapped (greater than threshold).}
##' }
##' @param metaLaneID ID of column that contains information about flow cell lane.
##' @param mitochondrialIdenitifier Regexp identifying mitochondrial gene in geneTable
##' @param pdf Prefix name of pdfs.
##' @param qcFeatures Features from the QC to plot.
##' @return SCD object
##' @author Wajid Jawaid
setMethod("performQC", signature(object="SCD"),
          function(object, selectedCells = "ALL", cutoffs = c(2e5, .2, .2, 0, 0, 1, 0),
                   metaLaneID = "flowCell", mitochondrialIdenitifier = "^mt|^MT",
                   pdf = NULL, qcFeatures = "ALL") {
              origFilter <- saveFilters(object)
              object <- deactivateFilters(object, includeQC = TRUE)
              counts <- counts(object)
              ercc.data <-  spikes(object)
              ercc.data <- ercc.data[,attr(ercc.data, "cellsInData")]
              htseqQC <- qc(object)
              htseqQC <- htseqQC[,attr(htseqQC, "cellsInData")]
              if (length(qcFeatures) > 1) {
                  htseqQC <- htseqQC[qcFeatures,]   
              } else if (qcFeatures != "ALL") {
                  htseqQC <- htseqQC[qcFeatures,]   
              }
              geneTable <- fData(object)
              meta <- pData(object)
              if (selectedCells != "ALL") {
                  counts <- counts[,selectedCells]
                  ercc.data <- ercc.data[,selectedCells]
                  htseqQC <- htseqQC[,selectedCells]
                  meta <- meta[selectedCells,]
              }
              cellIntersect <- intersect(colnames(counts), colnames(ercc.data))
              cellIntersect <- intersect(cellIntersect, colnames(htseqQC))
              missingCells <- setdiff(colnames(counts), cellIntersect)
              if (length(missingCells) > 0) {
                  cat("Spike-in and HTSeq data not available for all cells.\n")
                  counts <- counts[,cellIntersect]
                  ercc.data <- ercc.data[,cellIntersect]
                  htseqQC <- htseqQC[,cellIntersect]
                  meta <- meta[cellIntersect,]
              }
              cat("Performing QC ... ")
              qc <- qcFunc(counts, ercc.data, htseqQC, geneTable, meta,
                           metaLaneID, mitochondrialIdenitifier, pdf, cutoffs)
              object@qcOutput <- qc
              failedQC <- do.call(c, qc[[1]])
              meta <- pData(object)
              meta[failedQC, "passedQC"] <- FALSE
              pData(object) <- meta
              filterQC(object) <- TRUE
              cat("Done.\n")
              ## cat("Estimating size factors ... ")
              ## counts <- counts(object)
              ## spikes <- spikes(object)
              ## sf.data <- estSizeFact2(counts)
              ## data <- t(t(counts) / sf.data)
              ## sf.data.ercc <- estSizeFact2(spikes)
              ## cat("Done.\n")
              ## cat("Applying size factors ... ")
              ## data.ercc <- t(t(spikes) / sf.data.ercc)
              ## cat("Done.\n")
              ## varGenes <- vector("list", 6)
              ## names(varGenes)<- c("normalisedData", "sizeFactors", "highVarGenes",
              ##                     "adjusted.p", "points", "parameters")
              ## varGenes[["normalisedData"]] <- list(n.data.ercc = data.ercc)
              ## varGenes[["sizeFactors"]] <- list(sf.data = sf.data, sf.data.ercc = sf.data.ercc)
              ## object@technicalNoise <- varGenes
              ## ## data <- cbind(data, assayDataElement(object, "counts")[,failedQC])
              ## data <- cbind(data, matrix(NA, nrow = nrow(data), ncol = length(failedQC),
              ##                            dimnames = list(rownames(data), failedQC)))
              ## data <- data[,colnames(assayDataElement(object, "counts"))]
              ## filterQC(object) <- FALSE
              ## exprs(object) <- log10(1 + data)
              ## filterQC(object) <- TRUE
              object <- restoreFilters(object, origFilter)
              object
          })


##' Normalise counts data
##'
##' Normalise counts data
##' @title Normalise counts data
##' @param object SCD
##' @param normMethod DESeq or scran 
##' @param logFnc log10
##' @param pseudocount Default 1. Added to normalised counts before log transformation
##' @param ... Additional parameters
##' @return SCD object
##' @author Wajid Jawaid
##' @importFrom scran computeSumFactors
setMethod("normalize", signature(object="SCD"),
          function(object, normMethod = c("DESeq", "scran"), logFnc = log10,
                   pseudocount = 1, ...) {
              origFilter <- saveFilters(object)
              object <- deactivateFilters(object, includeQC = FALSE)
              object@normalise <- normMethod
              if (normMethod == "DESeq") {
                  cat("Estimating DESeq size factors ... ")
                  counts <- counts(object)
                  spikes <- spikes(object)
                  sf.data <- estSizeFact2(counts)
                  data <- t(t(counts) / sf.data)
                  sf.data.ercc <- estSizeFact2(spikes)
                  cat("Done.\n")
                  cat("Applying size factors ... ")
                  data.ercc <- t(t(spikes) / sf.data.ercc)
                  cat("Done.\n")
                  varGenes <- vector("list", 6)
                  names(varGenes)<- c("normalisedData", "sizeFactors", "highVarGenes",
                                      "adjusted.p", "points", "parameters")
                  varGenes[["normalisedData"]] <- list(n.data.ercc = data.ercc)
                  varGenes[["sizeFactors"]] <- list(sf.data = sf.data, sf.data.ercc = sf.data.ercc)
                  failedQC <- pData(phenoData(object))[,"failedQC", drop = FALSE]
                  failedQC <- rownames(failedQC)[failedQC[,"failedQC"]]
                  data <- cbind(data, matrix(NA, nrow = nrow(data), ncol = length(failedQC),
                                             dimnames = list(rownames(data), failedQC)))
                  data <- data[,colnames(assayDataElement(object, "counts"))]
                  filterQC(object) <- FALSE
                  exprs(object) <- logFnc(pseudocount + data)
              } else if (normMethod == "scran") {
                  cat("Estimating scran size factors ... ")
                  object <- computeSumFactors(object, ...)
                  sf <- sizeFactors(object)
                  if (any(is.na(sf))) {
                      failScran <- names(sf)[is.na(sf)]
                      object <- excludeCells(object, cellNames = failScran, reset = FALSE)
                  }
                  cat("done.\nApplying scran size factors ... ")                  
                  exprs(object) <- logFnc(t(t(counts(object) / sizeFactors(object))) + pseudocount)
                  varGenes <- vector("list", 6)
                  names(varGenes)<- c("normalisedData", "sizeFactors", "highVarGenes",
                                      "adjusted.p", "points", "parameters")
                  varGenes$sizeFactors$sf.data <- sizeFactors(object)
              }
              object@technicalNoise <- varGenes
              object <- restoreFilters(object, origFilter)
              object
          })

setMethod("coerce", signature(from="SCD", to="SCESet"), function(from, to) {
    new("SCESet", countData = assayDataElement(from, "counts"),
        cellPairwiseDistances = from@cellPairwiseDistances,
        featurePairwiseDistances = from@featurePairwiseDistances,
        phenoData = phenoData(from),
        featureData = featureData(from))
})

##' Set size factors
##'
##' Set size factors
##' @title Set size factors
##' @param object SCD Single Cell Dataset
##' @param type optional character argument providing the type or name of the
#' size factors to be accessed or assigned.
##' @param ... Additional parameters
##' @param value A vector to be stored
##' @return SCD object
##' @author Wajid Jawaid
setReplaceMethod("sizeFactors", signature(object = "SCD", value = "numeric"),
                 function(object, type = NULL, ..., value) {
                     ofield <- .construct_sf_field(object, type)
                     newCol <- data.frame(value,
                                          row.names = rownames(pData(object)))
                     colnames(newCol) <- ofield
                     addPheno(object, newCol)
                 })

.construct_sf_field <- function (object, type) {
    ofield <- "size_factor"
    if (!is.null(type)) {
        fc_available <- object@featureControlInfo$name
        if (length(fc_available) == 0L) {
            stop("no named controls specified in the SCESet object")
        }
        type <- match.arg(type, fc_available)
        ofield <- paste0(ofield, "_", type)
    }
    return(ofield)
}

