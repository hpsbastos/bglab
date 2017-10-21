## Author: Wajid Jawaid
## Date: 2 December 2014
## bglab package: Berthold Gottgens Lab scRNA-seq repo

##' @importFrom rARPACK eigs
##' @importFrom RColorBrewer brewer.pal
##' @import grDevices
##' @import graphics
##' @importFrom stats as.dendrogram as.dist cor density dist fitted.values gaussian hclust median order.dendrogram p.adjust pchisq prcomp predict quantile reorder sd var lm
##' @import utils
##' @import methods
##' @importFrom igraph graph.adjacency minimum.spanning.tree get.diameter get.shortest.paths get.edgelist E
##' @import ggplot2
##' @importFrom rgl plot3d par3d rgl.cur r3dDefaults open3d rgl.postscript text3d points3d segments3d lines3d selectpoints3d
##' @import Rtsne
##' @import Biobase
##' @import scater
##' @importMethodsFrom BiocGenerics counts


suppressMessages(require(Biobase, quietly = TRUE))
suppressMessages(require(igraph, quietly = TRUE))
suppressMessages(require(rgl, quietly = TRUE))
suppressMessages(require(ggplot2, quietly = TRUE))
suppressMessages(require(Rtsne, quietly = TRUE))


if (version$os == "linux-gnu") {
    fileOpen <- "xdg-open "
} else fileOpen <- "open "

##' Make new SCD object
##'
##' Construct a new Single cell dataset object for the bglab package
##' @title Make a new SCD object
##' @param experimentType Experiment Type. RNAseq or qPCR.
##' @param assayData Matrix of data with genes in rows and cells in columns
##' @param featureData Data frame of further details for each gene. Rownames must
##' match row names of assayData.
##' @param phenoData Data frame of further details for each cell. Rownames must
##' match column names of assayData.
##' @param counts Raw counts data.
##' @param spike Default NULL. Include spike-in data e.g. ERCC-92.
##' @param qc Default NULL. Include quality control counts from HTseq for example.
##' @param cellPairwiseDistances CPW
##' @param featurePairwiseDistances FPW
##' @param ... Further parameters.
##' @return New Single Cell Dataset object.
##' @author Wajid Jawaid
##' @export
newSCD <- function(experimentType, assayData = NULL, featureData, phenoData,
                   counts, spike = NULL, qc = NULL,
                   cellPairwiseDistances = dist(vector()),
                   featurePairwiseDistances = dist(vector()), ...) {
    if (is.null(assayData)) assayData = counts
    new("SCD", experimentType = experimentType, assayData = assayData,
        featureData = featureData, phenoData = phenoData, countData = counts,
        spike = spike, qc = qc, cellPairwiseDistances = cellPairwiseDistances,
        featurePairwiseDistances = featurePairwiseDistances, ...)
}

subsetData <- function(ncx, pD, filters) {
    ## cellsInData <- attr(x, "cellsInData")
    filterCell <- filters[[1]] && filters[[3]]
    filterQC <- filters[[4]]
    nrp <- nrow(pD)
    ## ncx <- ncol(x)
    if (nrp == ncx) {
        if (filterCell && filterQC){
            x <- pD[,"included"] & pD[,"passedQC"]
        } else if (!filterCell && filterQC) {
            x <- pD[,"passedQC"]
        } else if (filterCell && !filterQC) {
            x <- pD[,"included"]
        } else {
            x <- rep(TRUE, nrp)
        }
    } else {
        filterInc <- rownames(pD)[pD[,"included"]]
        filterPass <- rownames(pD)[pD[,"passedQC"]]
        filterBoth <- intersect(filterInc, filterPass)
        if (length(filter) != ncx) {
            if (filterCell && filterQC){
                x <- filterBoth
            } else if (!filterCell && filterQC) {
                x <- filterPass
            } else if (filterCell && !filterQC) {
                x <- filterInc
            } else 
                stop("Internal error in subsetData().")   
        }
    }
    ## attr(x, "cellsInData") <- cellsInData
    return(x)
}

##' Prepares annotated data frame
##'
##' Prepares annotated data frame
##' @title Prepares data frame
##' @param eD Annotation from gene expression matrix
##' @param aD Annotated data frame to be prepared.
##' @param nm Name to pass in case of error.
##' @return Returns prepared data frame
##' @author Wajid Jawaid
prepData <- function(eD, aD, nm) {
    annData <- split(t(aD), f = colnames(aD))
    whichCol <- lapply(annData, function(x) !(any(is.na(match(eD, x)))))
    whichCol <- do.call(c, whichCol)
    if (!any(whichCol)) stop(paste(nm, "does not match assayData."))
    firstGeneCol <- names(whichCol)[whichCol][1]
    rownames(aD) <- aD[, firstGeneCol]
    aD$included <- rep(FALSE, nrow(aD))
    aD$order <- rep(NA, nrow(aD))
    aDmatch <- match(eD, rownames(aD))
    aD$included[aDmatch] <- TRUE
    aD$order[aDmatch] <- 1:length(eD)
    aD <- aD[order(aD$order),]
    aD <- aD[aD$included,]
    if (nm == "phenodata") aD$passedQC <- aD$included
    return(aD)
}

##' Reorder data frame
##'
##' Reorder data frame
##' @title Reorder frame
##' @param object Annotated data frame
##' @return Reordered frame
##' @author Wajid Jawaid
orderFrame <- function(object) {
    temp <- object[object[, "included"],]
    temp <- temp[order(temp[, "order"]),]
    temp[,-match(c("included", "order"), colnames(temp)),drop=FALSE]
}

saveFilters <- function(object) {
    return(list(useFilter(object), filterGene(object), filterCell(object), filterQC(object)))
}

restoreFilters <- function(object, savedFil) {
    useFilter(object) <- savedFil[[1]]
    filterGene(object) <- savedFil[[2]]
    filterCell(object) <- savedFil[[3]]
    filterQC(object) <- savedFil[[4]]
    return(object)
}

deactivateFilters <- function(object, includeQC = FALSE) {
    useFilter(object) <- FALSE
    filterGene(object) <- FALSE
    filterCell(object) <- FALSE
    if (includeQC) filterQC(object) <- FALSE
    return(object)
}

##' Reads in file
##' Reads in the qRT-PCR file. A tab delimited file containing \eqn{-\delta}CT qRT-PCR
##' data
##'
##' Reads in the qRT-PCR file. A tab delimited file containing normalized
##' \eqn{-\delta}CT qRT-PCR data, assuming the first line is a header containing gene names
##' and the first coloumn contains cell identifiers.
##' @title Reads in qRT-PCR file
##' @param file Path of file and name containing \eqn{-\delta}CT qRT-PCR data
##' @param exclude Include names of genes that should be excluded from the analysis e.g.
##' housekeeping genes that were used to normalize the data or for those that did not
##' work well
##' @param sep Separator used in file containing qPCR data
##' @param cellIDcol Column containing the cell ID to be used
##' @param removeCol Unnecessary columns to be removed
##' @param cellStageCol  Column containing cell stage data
##' @return A matrix where each row represents a cell and each coloumn a gene
##' @author Wajid Jawaid
##' @export
read.qPCR <- function(file, exclude=NULL, sep="\t", cellIDcol = 1, removeCol = NULL, cellStageCol = NULL) {
    qPCR <- readLines(file)
    temp <- strsplit(qPCR, split=sep)
    genes <- temp[[1]]
    temp[[1]] <- NULL
    temp <- t(data.frame(temp, check.names = FALSE, stringsAsFactors = FALSE))
    attr(temp, "dimnames") <- NULL
    cellID <- temp[,cellIDcol]
    colnames(temp) <- genes
    cellStage <- temp[,cellStageCol]
    temp <- temp[,-c(cellIDcol, removeCol, cellStageCol)]
    temp <- apply(temp, 2, as.numeric)
    rownames(temp) <- cellID
    if (!is.null(exclude)) {
        temp <- temp[,-sapply(exclude, function(x) (grep(x, colnames(temp))))]
    }
    if (!is.null(cellStageCol)) attr(temp, "cellStage") <- factor(cellStage)
    return(temp)
}

##' Helper function for \code{plotGene}
##'
##' Helper function for \code{plotGene} to calculate colour gradient
##' @title Helper function for \code{plotGene}
##' @param data Gene expression data from \code{exprs} of single cell dataset object.
##' @param lcol Colour to be taken for low values. Default is "grey".
##' @param hcol Colour to be taken by high values. Default is "red".
##' @param na.col Colour for missing values, Default is Black.
##' @param opacity Opacity value from 0 to 255
##' @param maxVal Default NULL. Set the maximum and minimum color values manually.
##' @param minVal Default NULL. Set the maximum and minimum color values manually.
##' @param channelOrder Default 1:3 indicating red, green, blue.
##' @return Returns a vector of colours for each point to be plotted.
##' @author Wajid Jawaid
##' @export
colorGradient <- function(data, lcol="grey", hcol="red", na.col="Black", opacity = NULL,
                          maxVal = NULL, minVal = NULL, channelOrder = 1:3) {
    if (is.null(nrow(data))) {
        idx <- is.na(data)
        fdata <- data[!idx]
        colfunc <- colorRampPalette(c(lcol, hcol))
        if (is.null(minVal)) {
            lowest <- floor(min(fdata, na.rm=TRUE) * 10)/10    
        } else lowest <- minVal
        if (is.null(maxVal)) {
            highest <- ceiling(max(fdata, na.rm=TRUE) * 10)/10    
        } else highest = maxVal
        shiftCol <- ifelse(highest == lowest, highest /100,  (highest - lowest) / 100)
        ## lowest <- lowest - shiftCol
        ## highest <- highest + shiftCol
        if (highest == lowest) {
            lowest <- 0
            highest <- ifelse(highest == 0, 1, 1 + lowest)
        }
        pCol <- cut(fdata, seq(lowest-shiftCol, highest+shiftCol, length.out=100))
        pointColorT <- colfunc(100)[pCol]
        pointColor <- rep(NA, length(idx))
        pointColor[idx] <- na.col
        pointColor[!idx] <- pointColorT
        if (!is.null(opacity)) {
            pointColor <- rgb(t(col2rgb(pointColor)), alpha=opacity, maxColorValue=255)
        }
        pCol <- cbind(seq(lowest, highest, length.out = 100), colfunc(100))
        attr(pointColor, "leg") <- pCol
        return(pointColor)
    } else {
        e.c <- attr(data, "e.c")
        pointColor <- matrix(rep(e.c, ncol(data)), 4)
        ## gMax <- attr(data, "max")
        ## ndata <- data * (1 - e.c[1:nrow(data)]) / gMax
        ## ndata <- data * (1 - e.c[1:nrow(data)])
        ndata <- ((data / rowMax(data)) * (1 - e.c[1:nrow(data)]))
        pointColor[1:nrow(ndata),] <- pointColor[1:nrow(ndata),] + ndata
        if (!is.null(opacity)) pointColor[4,] <- opacity/255
        pointColor <- pointColor[c(channelOrder,4),]
        pointColor <- apply(pointColor, 2, function(x) rgb(x[1], x[2], x[3], x[4]))
        return(pointColor)
    }
}

##' Plots dimensionality reduced object with points coloured by level of gene expression.
##'
##' Plots dimensionality reduced object with points coloured by level of gene expression.
##' @title Visual representation of gene expression.
##' @param scd single cell dataset object
##' @param gene Gene to be visualised.
##' @param useDims dimensions to be used
##' @param reduceMethod The dimensionality reduction technique to be used
##' @param use.ggPlot Plot using ggPlots or native R graphics.
##' @param main Plot title. Defaults to gene name.
##' @param pch Point caharcter. Default is 16.
##' @param includeHidden Default FALSE, will only include highly variable genes.
##' @param bringToFront Bring high values to fore.
##' @param project2D Default NULL. Projects 3D graph into 2D by the given transformation
##' matrix. This matrix must have been produced using \code{getMVP} on an open rgl device.
##' @param opacity Default NULL. Set opacity.
##' @param legPos Default "topright". Position of legend
##' @param empty.channel Default (0.2,0.2,0.2). Sets the zero state of cell
##' not expressing anything 
##' @param adjust Default (0,0,0). Adjustment for red, green and blue channel,
##' respectively.
##' @param legend Default FALSE. Plot legend or not.
##' @param multigene Default FALSE. For internal use indicating function
##' is being called from the multiGene() function.
##' @param useMeta Default FALSE. To use information from the meta data.
##' @param bgCol Default black. Background colour.
##' @param leg.text.col Default white. Legend text colour.
##' @param leg.inset Default 0. Legend inset.
##' @param plt.mar Default rep(0,4). Plot margins
##' @param channelOrder Default 1:3 indicating red, green, blue.
##' @param ... Additional parameters to pass to plot
##' @inheritParams colorGradient
##' @return Produced plot
##' @author Wajid Jawaid
##' @export
plotGene <- function(scd, gene, useDims=1:2, reduceMethod = c("pca","diffMap","tsne","isomap"),
                     lcol="grey", hcol="red", use.ggPlot=FALSE, main = gene, pch = 20,
                     includeHidden = FALSE, bringToFront = TRUE, project2D = NULL,
                     opacity = NULL, legPos = "topright", empty.channel = c(.2,.2,.2,.5),
                     adjust = 0, minVal = NULL, maxVal = NULL, legend = FALSE,
                     multigene = FALSE, useMeta = FALSE, bgCol = "black",
                     leg.text.col = "white", leg.inset = 0, plt.mar = rep(0,4),
                     channelOrder = 1:3, ...) {
    reduceMethod <- match.arg(reduceMethod)
    mg <- length(gene) > 1
    redDim <- eigenvecs(slot(scd, reduceMethod))
    if (mg) ng <- length(gene)
    origFilter <- saveFilters(scd)
    useFilter(scd) <- TRUE
    filterCell(scd) <- TRUE
    cells <- rownames(pData(scd))
    pD <- pData(scd) 
    filterGene(scd) <- !includeHidden
    ## data <- exprs(scd)
    ## data <- data[,rownames(pD)]
    geneTable <- fData(scd)
    redDim <- eigenvecs(slot(scd, reduceMethod))[cells,]
    if (!is.null(project2D)) {
        redDim <- applyMVP(redDim[,useDims], project2D)
        useDims <- 1:2
    }
    ## data <- data[,rownames(redDim)]
    ## data <- data[,match(rownames(redDim), colnames(data))]
    geneID <- rownames(geneTable[match(gene, geneTable[,2]), ])
    if (any(is.na(geneID)) || any(grepl("NA", geneID))) {
        scd <- restoreFilters(scd, origFilter)
        cat(geneID, "\n")
        stop("Gene not found. Try includeHidden = \"TRUE\".")   
    }
    pointCol <- assayDataElement(scd, "exprs")[geneID,cells]
    ## pointCol <- data[match(geneID, rownames(data)),cells]
    if (!is.null(opacity)) opacity <- pD[,opacity]
    if (bringToFront) {
        if (!mg) pC <- pointCol
        else pC <- colSums(pointCol)
        pC <- order(pC, decreasing=FALSE)
        redDim <- redDim[pC,]
        if (!mg) pointCol <- pointCol[pC]
         else pointCol <- pointCol[,pC]
        opacity <- opacity[pC]
    }
    attr(pointCol, "max") <- max(pointCol)
    attr(pointCol, "min") <- min(pointCol)
    attr(pointCol, "e.c") <- empty.channel
    attr(pointCol, "adjust") <- adjust
    r.c <- empty.channel; r.c[1] <- 1; r.c <- rgb(r.c[1], r.c[2], r.c[3], r.c[4])
    g.c <- empty.channel; g.c[2] <- 1; g.c <- rgb(g.c[1], g.c[2], g.c[3], g.c[4])
    b.c <- empty.channel; b.c[3] <- 1; b.c <- rgb(b.c[1], b.c[2], b.c[3], b.c[4])
    if (length(useDims)==2) {
        if (!use.ggPlot) {
            if (mg) {
                par(bg = bgCol, mar=plt.mar)
            } else {
                if (legend) {
                    if (!multigene) layout(matrix(2:1,2), heights=c(1,5))
                    par(mar=c(2,2,1,2))   
                }
            }
            pCol <- colorGradient(pointCol, lcol, hcol, opacity=opacity,
                                  minVal = minVal, maxVal = maxVal,
                                  channelOrder = channelOrder)
            plot(redDim[, useDims], col = pCol,
                 main = ifelse(legend, "", main), pch=pch, ...)
            if (mg && legend) {
                legend(legPos, legend = gene, fill = c(r.c,g.c,b.c)[1:ng],
                       text.col = leg.text.col, inset = leg.inset)
            } else {
                if (legend) {
                    legCol <- attr(pCol, "leg")
                    origPar <- par(no.readonly = TRUE)
                    par(mar=c(.5,2,4,2), mgp=c(3,.22,0), las=1)
                    image(as.numeric(legCol[,1]), 1,
                          matrix(seq_along(legCol[,1]),ncol=1), col=legCol[,2], axes=FALSE,
                          xlab="", ylab="", main = main, cex.axis=.5, sub = expression(paste(log[10], "transformed normalised counts")))
                    axis(1, cex.axis=.5, tck=-.1)
                    par(origPar)
                }
            }
        } else {
            geneDF <- data.frame(PC1 = redDim[,1], PC2 = redDim[,2], pointColor = pointCol)
            ggplot(data = geneDF, aes_string(x = "PC1", y = "PC2", color="pointColor")) +
                geom_point(size = 2) +
                ggtitle(gene) +
                scale_color_gradient(low=lcol, high=hcol)
        }
    } else if (length(useDims)==3) {
        dimReduc <- slot(scd, reduceMethod)
        plot3(dimReduc,
              colorBy=assayDataElement(scd, "exprs")[geneName2id(scd,gene),
                                                     rownames(eigenvecs(dimReduc))],
              plotCols=c(lcol,hcol), main = main, useDims = useDims, opacity = opacity, ...)
    }
    scd <- restoreFilters(scd, origFilter)
    invisible(list(pCol=pCol, pointCol=pointCol, minVal=minVal,
                   maxVal=maxVal, opacity=opacity, channelOrder=channelOrder,
                   genes = gene))
}

##' Produce pdf of \code{plotGene} output
##'
##' Produce pdf of \code{plotGene} output
##' @title Produce pdf of \code{plotGene} output
##' @param directory Directory where to place file, without trailing forward slash.
##' @inheritParams plotGene
##' @inheritParams colorGradient
##' @return Produces a file of the name "gene.pdf" in \code{directory}.
##' @author Wajid Jawaid
##' @export
printGene <- function(gene, directory, data) {
    pdf(paste(directory, "/", gene, ".pdf", sep=""))
    plotGene(gene, data)
    dev.off()
}

##' Writes data from single cell dataset object for use in Matlab
##'
##' Writes data from single cell dataset object for use in Matlab
##' @title Writes data from single cell dataset object for use in Matlab
##' @param hscd Single cell dataset object
##' @param file Filename (and path) to save file to
##' @param cols Vector of colors
##' @param colourBy Name of the column in phenotype data by which to colour cells
##' @return Generates a ".mat " file
##' @author Wajid Jawaid
##' @importFrom R.matlab writeMat
##' @export
exportToMatlab <- function(hscd, file, cols, colourBy) {
    exprData <- exprs(hscd)
    colorData <- as.numeric(pData(hscd)[,colourBy])
    colorMatrix <- t(col2rgb(cols)) / 255
    stageOrd <- levels(pData(hscd)[,colourBy])
    R.matlab::writeMat(paste(file,".mat",sep=""), exprData=exprData, colorData=colorData,
                       colorMatrix=colorMatrix, stageOrd=stageOrd)
}

##' Returns human readable gene names
##'
##' Returns human readable gene names
##' @title Returns human readable gene names
##' @param scd Single cell Dataset object
##' @param genes Ensembl gene IDs
##' @return Named vector
##' @author Wajid Jawaid
##' @export
geneID2Name <- function(scd, genes) {
    origFilter <- saveFilters(scd)
    useFilter(scd) <- filterGene(scd) <- FALSE
    g <- fData(scd)
    idx <- match(genes, rownames(g))
    gNames <- g[idx,2]
    names(gNames) <- g[idx,1]
    scd <- restoreFilters(scd, origFilter)
    return(gNames)
}

##' Returns Ensembl id given a gene name
##'
##' Returns Ensembl id given a gene name
##' @title Returns Ensembl id given a gene name
##' @param scd Single Cell Dataset
##' @param genes Gene names
##' @param useAgrep To use approximate match
##' @return Named vector
##' @author Wajid Jawaid
##' @export
geneName2id <- function(scd, genes, useAgrep = FALSE) {
    origFilter <- saveFilters(scd)
    useFilter(scd) <- filterGene(scd) <- FALSE
    g <- fData(scd)
    if (!useAgrep) idx <- match(genes, g[,2])
    else idx <- agrep(genes,g[,2])
    geneId <- rownames(g[idx,])
    names(geneId) <- g[idx,2]
    scd <- restoreFilters(scd, origFilter)
    return(geneId)
}

##' Find gene name
##'
##' Find gene name
##' @title Find an approximate gene name
##' @param scd Single Cell data set object
##' @param g approximate gene name or alias
##' @param includeHidden Search genes that are otherwise excluded.
##' @return returns names from loaded ENSMBL gene table
##' @author Wajid Jawaid
##' @export
findGene <- function(scd, g, includeHidden=TRUE) {
    origFilter <- saveFilters(scd)
    useFilter(scd) <- filterGene(scd) <- !includeHidden
    gTab <- fData(scd)
    scd <- restoreFilters(scd, origFilter)
    gTab[agrep(g, gTab[,2], ignore.case = TRUE),2]
}

## ##' Fast vectorised Euclidean distance calculator
## ##'
## ##' Calculates Euclidean distances between vectors arranged as columns in a matrix.
## ##' @title Fast vectorised Euclidean distance calculator
## ##' @return Returns a matrix of pairwise distances
## ##' @author Wajid Jawaid
## ##' @param x Matrix with vectors in columns.
## ##' @param squared Will not perform the square root, i.e. will return the squared `L2-norm'.
## ##' @export
## fastDist <- function(x, squared = FALSE) {
##     a <- colSums(x^2)
##     a <- a * matrix(1, ncol(x), ncol(x))
##     a <- a + t(a)
##     ab <- t(x) %*% x
##     d <- a - 2 * ab
##     diag(d) <- 0
##     if (!squared) d <- sqrt(d)
##     return(d)
## }

##' Decorate plots with log10 values
##'
##' Decorate plots with log10 values. The plot must have prviously been generated.
##' @title Decorate plots with log10 values
##' @return Values added to plot
##' @author Wajid Jawaid
neatLabels <- function() {
    atx <- axTicks(1, log=TRUE)
    temp <- log(atx,10)
    atx <- atx[temp == round(temp)]
    xtlabs <- sapply(atx, function(i) as.expression(bquote(10^ .(log(i,10)))))
    aty <- axTicks(2, log=TRUE)
    temp <- log(aty,10)
    aty <- aty[temp == round(temp)]
    rm(temp)
    ytlabs <- sapply(aty, function(i) as.expression(bquote(10^ .(log(i,10)))))
    axis(1, at=atx, labels=xtlabs)
    axis(2, at=aty, labels=ytlabs)
    box()
}

##' Plots histogram of gene expression distribution
##'
##' Plots histogram of gene expression distribution
##' @title Histogram of gene distribution
##' @param scd Single Cell Dataset object
##' @param gene Gene to be plotted
##' @param fun Function to use on gene expression values. Defaults to log10 normalisedcounts i.e. \eqn{10^{count}}
##' @param ... Further parameters to pass to \code{hist()}
##' @return Gneerates a plot and returns value from \code{hist()} invisibly.
##' @author Wajid Jawaid
##' @export
geneDistribution <- function(scd, gene, fun = function(x) 10^x, ...) {
    origFilter <- saveFilters(scd)
    filterGene(scd) <- FALSE
    if (is.na(g <- geneName2id(scd, gene))) stop("Gene not found")
    if (fData(scd)[g,"included"]==FALSE)
        cat("Note: Gene is not included in analysis.\n")
    h <- hist(fun(exprs(scd)[g, ]), main = geneID2Name(scd, g), ...)
    scd <- restoreFilters(scd, origFilter)
    invisible(h)
}

##' Generates multiple 3D plots with marked gene expression
##'
##' Generates multiple 3D plots with marked gene expression
##' @title Generates multiple 3D plots with marked gene expression
##' @param scd Single cell dataset
##' @param reduceMethod dimensionality reduction method to be plotted
##' @param geneVect character vector of gene names
##' @param nr Number of rows on scree
##' @param nc Number of columns on screen
##' @param res Default 1439 x 893. Resolution of graphics device. CIMR 3200x1024
##' @param useDims Dimensions to be plotted
##' @param height height of display windows. Deafaults to 7.
##' @param width width of display window. Defaults to 7.
##' @param mar Plot margins. Defaults to c(0,0,3.1,0).
##' @param pdf.name To output to pdf file. Provide name of file. Defaults to NULL
##' @param includeFiltered Defaults to FALSE. To include filtered values.
##' @param vp Default NULL. Pass par() or par3d() options to set desired viewpoint.
##' Particularly useful when generating pdf(s).
##' @param arrange Default TRUE. To arrange the windows
##' @param project2D Default NULL. Projects 3D graph into 2D by the given transformation
##' matrix. This matrix must have been produced using \code{getMVP} on an open rgl device.
##' @param layoutMat Default NULL. Provides greater control of the layout matrix.
##' @param legend Default FALSE. To plot legend or not.
##' @param ... Parameters to pass to plotGene.
##' @return Plots 3D plots.
##' @author Wajid Jawaid
##' @export
multiGenes <- function(scd, reduceMethod, geneVect, nr, nc, res = c(1439, 893),
                       useDims=1:3, height=7, width=7, mar=c(0,0,3.1,0), pdf.name=NULL,
                       includeFiltered = FALSE, vp=NULL, arrange=TRUE, project2D = NULL,
                       layoutMat = NULL, legend = FALSE, ...) {
    origFilter <- saveFilters(scd)
    useFilter(scd) <- filterGene(scd) <- !includeFiltered
    if (length(useDims)==3 && is.null(project2D)) {
        if (is.null(vp)) {
            vp$userMatrix <- r3dDefaults$userMatrix
            vp$zoom <- 1
        }
        if (!is.null(pdf.name)) pdf.name <- gsub(".pdf", "", pdf.name)
        goi <- geneVect
        if (!includeFiltered) {
            goi <- goi[fData(scd)[geneName2id(scd, goi), "included"]]
        } else goi <- goi[goi %in% geneID2Name(scd, rownames(fData(scd)))]
        goi <- goi[!is.na(goi)]
        plotPos <- expand.grid(seq(0,res[1],length.out = nc+1)[1:nc], seq(44,res[2],
                                                                          length.out = nr+1)[1:nr])
        size <- round(c(plotPos[2,1] - plotPos[1,1], plotPos[1+nc,2] - plotPos[1,2]))
        for (g in 1:length(goi)) {
            tl <- round(as.numeric(plotPos[g,]))
            if (arrange) vp$windowRect <- c(tl, tl + size)
            open3d(zoom=vp$zoom, userMatrix=vp$userMatrix, windowRect=vp$windowRect)
            plotGene(scd, gene=goi[g], reduceMethod=reduceMethod, useDims=useDims,
                     includeHidden = TRUE, ...)
            if (!is.null(pdf.name)) {
                rgl.postscript(paste(pdf.name,"-",goi[g], "temp.pdf",sep=""), fmt="pdf")
            }
        }
        if (!is.null(pdf.name)) {
            system(paste("pdfjoin *", pdf.name, "*temp.pdf ", "--outfile ", pdf.name,
                         ".pdf", sep=""))
            system(paste("rm ", pdf.name, "*temp.pdf", sep=""))
            system(paste(fileOpen, pdf.name, ".pdf", sep=""))
        }
    } else if (length(useDims)==2 || (length(useDims)==3 && !is.null(project2D))) {
        goi <- geneVect
        if (is.null(pdf.name)) {
            dev.new(height=height, width=width)
        } else {
            pdf(pdf.name, height=height, width=width)
        }
        par(mar=mar)
        if (is.null(layoutMat)) {
            layout(matrix(1:(nc*nr),nr,nc))
        } else {
            convertMat <- function(x) {
                y <- matrix(nrow=2*nrow(x), ncol=ncol(x))
                y[seq(1,nrow(y), by=2),] <- 2 * x
                y[seq(2,nrow(y), by=2),] <- 2 * x - 1
                y
            }
            if (legend) {
                h <- rep.int(c(1,5),nrow(layoutMat))
            } else {
                h <- rep.int(1, nrow(layoutMat))
            }
            layout(convertMat(layoutMat), heights = h)
        }
        for (g in goi) {
            plotGene(scd, gene=g, reduceMethod=reduceMethod, useDims=useDims,
                     includeHidden = TRUE, project2D = project2D, multigene = TRUE,
                     legend = legend, ...)
            if (!is.null(vp)) par(vp)
        }
        if (!is.null(pdf.name)) dev.off()
    }
    scd <- restoreFilters(scd, origFilter)
}

##' Estimates size factors as in DESeq.
##'
##' Estimates size factors as in DESeq.
##' @title Estimates size factors as in DESeq.
##' @param x Matrix of values. Cells in columns. Genes in rows.
##' @return Numeric vector of size factors.
##' @author Wajid Jawaid
##' @export
estSizeFact <- function(x) {
    ln_gm_mean <- rowMeans(log(x))
    apply(x, 2, function(y) exp(median((log(y) - ln_gm_mean)[is.finite(ln_gm_mean)])))
}


##' Normalises and performs technical noise analysis for single cell RNAseq(Brennecke et al.)
##'
##' Normalises and performs a technical noise analysis for single cell RNAseq.
##' Give raw counts. Data is normalised using size-factor normalisation identical to DESeq2,
##' then a technical noise analysis is perfromed as described in brennecke et al.
##' @title Technical noise analysis for single cell RNAseq
##' @param object SCD object
##' @param data.ercc A matrix of counts with cells in columns and ERCC spike-ins in rows.
##' The columns must mutch those in the gene count table passed in "data".
##' @param cvThresh Defalt 0.3. The CV threshold to use for choosing cells to fit. See quant.
##' @param quant Default 0.8. The quantile of data to fit so that linear model is not fitted to
##' plateau at lower gene expression values. See Breenecke et al. for further details
##' @param minBiolDisp Default 0.25. Assumed biological dispersion. See paper.
##' @param fdr Deafult .1. False Discovery Rate for chi-squared test.
##' @param data.lengths Gene lengths in kilobases for PK "per kilobase normalisation"
##' @param ercc.lengths ERCC lengths in kilobases for PK "per kilobase normalisation"
##' @param meanForFit Provide a minimum mean for fit. (OPTIONAL)
##' @param renormalise Default FALSE. Renormalise or not. Will only work with "DESeq"
##' normalisation method
##' @return Returns a list
##' @importFrom statmod glmgam.fit
##' @author Wajid Jawaid
techVarSub <- function(object, data.ercc = NULL, cvThresh=.3, quant=.8, minBiolDisp=.5^2,
                       fdr=.1, data.lengths = NULL, ercc.lengths = NULL, meanForFit=NULL,
                       renormalise = FALSE) {
    tData <- counts(object)
    if (object@normalise == "DESeq" && renormalise) {
        sf.data <- estSizeFact2(tData)
    } else {
        sf.data <- sf(object)
    }
    tData <- t(t(tData) / sf.data)
    aMean <- rowMeans(tData)
    cv2a <- (apply(tData,1,sd) / aMean)^2
    if (is.null(data.ercc)) {
        useERCC <- FALSE
        sf.data.ercc <- sf.data
        data.ercc <- tData
        sMean <- aMean
        cv2s <- cv2a
    } else {
        useERCC <- TRUE
        sf.data.ercc <- sf(object, type = "tech")
        data.ercc <- object@technicalNoise$normalisedData$n.data.ercc[,names(sf.data.ercc)]
        sMean <- rowMeans(data.ercc)
        cv2s <- (apply(data.ercc,1,sd) / sMean)^2
    }
    if (!is.null(data.lengths) && !is.null(ercc.lengths)) {
        if (length(data.lengths) != nrow(tData) || length(ercc.lengths) != nrow(data.ercc)) {
            stop("Gene/ERCC lengths data does not match data given.")
        }
        aMean <- aMean / data.lengths * 1000
        sMean <- sMean / ercc.lengths * 1000
    }
    if (is.null(meanForFit))
        meanForFit <- unname(quantile(sMean[which(cv2s > cvThresh)], quant))
    useForFit <- sMean >= meanForFit
    ## meanForFit
    ## table(useForFit)
    fit <- statmod::glmgam.fit(cbind(a0=1, a1tilde=1/sMean[useForFit]), cv2s[useForFit])
    residual <- var(log(fitted.values(fit)) - log(cv2s[useForFit]))
    total <- var(log(cv2s[useForFit]))
    1 - residual / total
    a0 <- fit$coefficients["a0"]
    a1t <- fit$coefficients["a1tilde"]
    ## lines(xg, a1t/xg + a0, col="red", lwd=2)
    df <- ncol(tData) - 1
                                        #lines(xg, (a1t/xg + a0) * qchisq(.025,df) / df, col="red", lwd=2, lty=2)
    ## minBiolDisp <- .25^2  #### Parameter changed from .5 to .25 ####
    xi <- mean(1/sf.data.ercc)
    m <- ncol(tData)
    if (!useERCC) {
        psi <- xi + (a1t - xi)
        chi2.values <- df * cv2s / (psi / sMean + a0)
        pA <- pchisq(chi2.values, df, lower.tail = FALSE)
    } else {
        psi <- xi + (a1t - xi) * mean(sf.data.ercc/sf.data)
        cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
        testDenom <- (aMean*psi + cv2th*aMean^2) / (1 + cv2th/m)
        pA <- 1 - pchisq(apply(tData,1,var) * (m-1) / testDenom, m-1)
    }
    padj <- p.adjust(pA, "BH")
    highVarGenes <- names(padj)[which(padj<fdr)]
    ## table(padj<.1)
    ## lines(xg, (psi/xg + a0 + minBiolDisp) * qchisq(.975,df) / df, col="red", lwd=2, lty=2)
    varGenes <- vector("list", 6)
    names(varGenes) <- c("normalisedData", "sizeFactors", "highVarGenes", "adjusted.p", "points", "parameters")
    attr(highVarGenes, "fdr") <- fdr
    varGenes[["normalisedData"]] <-
        list(n.data.ercc = object@technicalNoise$normalisedData$n.data.ercc)
    varGenes[["sizeFactors"]] <-
        list(sf.data = object@technicalNoise$sizeFactors$sf.data,
             sf.data.ercc = object@technicalNoise$sizeFactors$sf.data.ercc)
    varGenes[["highVarGenes"]] <- highVarGenes
    varGenes[["adjusted.p"]] <- padj
    varGenes[["points"]] <- list(genes=data.frame(mean=aMean, cv2=cv2a,
                                                  row.names = rownames(tData),
                                                  highVar = padj<fdr),
                                 spikes=data.frame(mean=sMean, cv2=cv2s,
                                                   row.names = rownames(data.ercc)))
    varGenes[["parameters"]] <- list(a0=unname(a0), a1tilde=unname(a1t), xi=unname(xi),
                                     df=unname(df), meanForFit=unname(meanForFit),
                                     psi=unname(psi), useForFit=unname(table(useForFit)[2]),
                                     minBiolDisp=minBiolDisp, usedERCC=useERCC,
                                     cvThresh = cvThresh, quant = quant)
    return(varGenes)
}

##' Plots a technical noise analysis plot
##'
##' Plots a technical noise analysis plot
##' @title Plots a technical noise analysis plot
##' @param varGenes Object returned from \code{techVar()} function
##' @param pdf Default NULL. If anything other than NULL will produce a pdf of that name.
##' @param pch Default 20. Point type to use in plot, if specific points are not defined.
##' See below. Parameter is plossed to base R \code{plot} function
##' @param pch.var Defaults to pch above if not defined. Point type for highly variable genes.
##' @param pch.notvar Defaults to pch above if not defined. Point type for not highly
##' variable genes.
##' @param pch.spike Defaults to pch above if not defined. Point type for spike-ins.
##' @param ... Further parameters to pass to plot
##' @return Produces a plot.
##' @author Wajid Jawaid
##' @export
plotTechVar <- function(varGenes, pdf=NULL, pch=20, pch.var = pch, pch.notvar = pch,
                        pch.spike = pch, ...) {
    if (!varGenes$parameters$usedERCC) pch.spike = ""
    if (!is.null(pdf)) pdf(pdf)
    minBiolDisp = varGenes$parameters[["minBiolDisp"]]
    par(mar=c(5.1,4.6,4.1,2.1))
    g <- varGenes$points$genes
    g <- g[g[,"mean"]>0,]
    s <- varGenes$points$spikes
    s <- s[s[,"mean"]>0,]
    plotRanges <- rbind(g[,c("mean", "cv2")], s[,c("mean", "cv2")])
    plot(g[!g[,"highVar"],c("mean", "cv2")], col="Gray", pch=pch.notvar, log="xy",
         axes=FALSE, ylab=expression(sigma^{2}/mu^{2}), xlab=expression(mu),
         xlim = range(plotRanges[,1]), ylim = range(plotRanges[,2]), ...)
    points(g[g[,"highVar"],c("mean","cv2")], pch=pch.var, col="Red")
    points(s[,c("mean", "cv2")], pch=pch.spike, col="Blue")
    neatLabels()
    xg <- log(par()$xaxp[1:2],10)
    xg <- 10^seq(xg[1],xg[2],length.out = 1000)
    a1t <- varGenes$parameters[["a1tilde"]]
    a0 <- varGenes$parameters[["a0"]]
    psi <- varGenes$parameters[["psi"]]
    xi <- varGenes$parameters[["xi"]]
    df <- varGenes$parameters[["df"]]
    usedERCC <- varGenes$parameters[["usedERCC"]]
    if (usedERCC) {
        lines(xg, a1t/xg + a0, col="red", lwd=2)
        lines(xg, (psi/xg + a0 + minBiolDisp), col="red", lwd=2, lty=2)
    } else lines(xg, (psi/xg + a0), col="red", lwd=2, lty=2)
                                        # lines(xg, (psi/xg + a0 + minBiolDisp) * qchisq(.975,df) / df, col="red", lwd=2, lty=2)
    if (!is.null(pdf)) dev.off()
}

##' Get rendering data from an open rgl plot
##'
##' Get rendering data from an open rgl plot
##' @title Retrieve rendering data from current rgl plot
##' @return Matrix to achieve transformation from model space to the screen
##' @author Wajid Jawaid
##' @export
getMVP <- function() {
    if (rgl.cur()==0) stop("No rgl window open")
    modelMatrix <- par3d()$modelMatrix
    projMatrix <- par3d()$projMatrix
    projMatrix %*% modelMatrix
}

##' Applies transformation from rgl plot to data
##'
##' Applies transformation from rgl plot to data
##' @title Applies transformation from rgl plot to data
##' @param data 3D matrix of value to be plotted, with each dimension in a column
##' @param MVP Default NULL. The transformation matrix of Projection matrix %*% Viewer matrix %*%
##' Model matrix. This can be retrieved using the \code{getMVP} function or if not supplied is
##' retrieved from the current open rgl device.
##' @return Returns the transformed data ready to be plotted using the R base plot function.
##' @author Wajid Jawaid
##' @export
applyMVP <- function(data, MVP=NULL) {
    if (is.null(MVP)) MVP <- getMVP()
    data <- rbind(t(data),1)
    data <- t(MVP %*% data)
    data <- data[,1:3] / data[,4]
}

##' Plots gene profile on trajectory
##'
##' Plots gene profile on trajectory
##' @title Plots gene profile on trajectory
##' @param scd SCD object
##' @param gene Character vector naming gene of interest
##' @param reduceMethod Default "diffMap". Dimensionality reduction method.
##' @param paths Vector of path names. Path names for which gene profiles are to be plotted
##' @param pathlty Vector of line types for different paths. Give in same order as paths. 
##' @param geneCol Vector of colors. Give in same order as genes.
##' @param winLength Integer. Window length to use for sliding average
##' @param legendPos Character vector of length 1. Giving position of legend or
##' optionally "none" if no legend is to be plotted. 
##' @param tpos Default 0. Position where to place the title.
##' @param useGam Use GAM for smoothing
##' @param ... Additional parameter to pass to plot.
##' @return Returns plot and a list of results.
##' @author Wajid Jawaid
##' @export
plotGeneProfile <- function(scd, gene, reduceMethod="diffMap", paths, pathlty=1:length(paths), geneCol =1:length(gene), winLength=10, legendPos="topright", tpos=0, useGam = TRUE,  ...) {
    gpro <- list(all=exprs(scd)[geneName2id(scd, gene),,drop=FALSE])
    .paths <- getPaths(slot(scd, reduceMethod))
    names(.paths) <- sapply(.paths, "[[", "name")
    paths <- .paths[paths]
    for (i in 1:length(paths)) gpro[[paths[[i]]$name]] <- gpro$all[,paths[[i]]$path,drop=FALSE]
    gpro$all <- NULL
    if (useGam) {
        gpro2 <- lapply(gpro, function(x) apply(x, 1, smoothGn))
    } else gpro2 <- lapply(gpro, function(x) apply(x, 1, slidingAv, winLength=winLength))
    ylim <- range(do.call(c, lapply(gpro2, range)))
    xlim <- c(1,max(sapply(gpro2, nrow)))
    plot(0, type="n", xlim = xlim, ylim = ylim, xlab="pseudotime", ylab="", main="", ...)
    ## plot(gpro2[[1]][,1], lty=pathlty[1], type="l", ylim = ylim, xlab="pseudotime", ylab="",
    ##      main="", col = geneCol[1], ...)
    mtext(paste(gene, collapse=" "), side = 3, line = tpos)
    for (i in 1:length(gpro2)) {
        for (j in 1:length(gene)) {
            points(gpro2[[i]][,j], lty=pathlty[[i]], col=geneCol[[j]], type="l")
        }
    }
    if (legendPos != "none") legend(legendPos, legend = c(gene, names(paths)),
                                    lty = c(rep(1,length(gene)), pathlty), col = c(geneCol, rep("gray", length(paths))))
    invisible(list(raw=gpro, slidingAv=gpro2))
}

##' Calculates sliding average
##'
##' Calculates sliding average
##' @title Calculates sliding average
##' @param data Numeric vector for which sliding average is to be calculated
##' @param winLength Integer defining window length to be used
##' @return Numeric vector of sliding average for the given window length
##' @author Wajid Jawaid
##' @export
slidingAv <- function(data, winLength=10) {
    l <- length(data)
    csum <- cumsum(data)
    csumShift <- csum[-c(1:(winLength-1))]
    csum <- c(0,csum)[1:length(csumShift)]
    (csumShift - csum) / winLength
}

##' Gam smoothing of gene expression
##'
##' Gam smoothing of gene expression
##' @title Gam smoothing of gene expression
##' @param data Vector of data
##' @param bs Default "tp", represents thin plate, may be cubic or other.
##' @param family Deafult gaussian().
##' @param n Not used.
##' @param ... Additional parameters to pass to predict
##' @return Returns vector of smoothed values
##' @author Wajid Jawaid
##' @importFrom mgcv gam
smoothGn <- function(data, bs="tp", family=gaussian(), n=length(data), ...) {
    dat <- cbind.data.frame(x=1:length(data), y=data)
    predict(gam(y ~ s(x, bs="tp"),
                data=dat,
                family=gaussian()), ...)
}


##' Differential expression
##'
##' Identify differentially expressed genes using Wilcoxon Rank sum test
##' @title Differential expression using Wilcoxon Ranks sum test
##' @param scd Single cell data set object
##' @param group1 Names of cells in group 1
##' @param group2 Names of cells in group 2
##' @param thr Default 1e-4. False discovery rate threshold
##' @param sdFilter Default TRUE. Remove genes with zero standard deviation
##' @param colSideCols Which column to use for addtional column side colors for heatmap.
##' @param cscCols Colours to use for the different levels
##' @param minCounts Only include genes with more than this number of counts
##' over all the cells used in DE.
##' @return Returns a list object
##' @author Wajid Jawaid
##' @export
##' @importFrom coin wilcox_test
doDE <- function(scd, group1, group2, thr=1e-4, sdFilter=TRUE, colSideCols = NULL, cscCols = "NULL", minCounts = 1) {
    scd <- deactivateFilters(scd)
    data <- exprs(scd)
    ubexpr <- data[,group1,drop=FALSE]
    lbexpr <- data[,group2,drop=FALSE]
    notExpressed <- which(rowSums(ubexpr) + rowSums(lbexpr) <= minCounts)
    ubexpr <- ubexpr[-notExpressed,,drop=FALSE]
    lbexpr <- lbexpr[-notExpressed,,drop=FALSE]
    rmSdZero <- function(x) x[apply(x, 1, sd) > 0,]
    if (sdFilter) {
        notDynamic <- which(apply(ubexpr,1,sd)==0)
        notDynamic <- c(notDynamic, which(apply(lbexpr,1,sd)==0))
        ubexpr <- ubexpr[-notDynamic,,drop=FALSE]
        lbexpr <- lbexpr[-notDynamic,,drop=FALSE]
    }

    p <- vector("numeric", nrow(ubexpr))
    names(p) <- rownames(ubexpr)
    f <- p
    cbexpr <- cbind(ubexpr, lbexpr)
    id <- factor(c(rep("test", length(group1)), rep("control", length(group2))))
    for (i in 1:nrow(ubexpr)) {
        p[i] <- coin::pvalue(wilcox_test(cbexpr[i,] ~ id))
        ## p[i] <- pvalue(wilcox_test(u ~ p, dframe))
        ## p[i] <- wilcox.test(ubexpr[i,], lbexpr[i,])$p.value
        f[i] <- mean(ubexpr[i,]) / mean(lbexpr[i,])
    }
    
    p.orig <- p
    p <- p.adjust(p, method="fdr")

    sum(p<thr, na.rm=TRUE)
    p1 <- p[which(p<thr)]
    f1 <- f[which(p<thr)]
    de <- data[names(p1),c(colnames(ubexpr), colnames(lbexpr)), drop=FALSE]
    csc <- NULL
    if (!is.null(colSideCols)) {
        csc <- pData(scd)[c(colnames(ubexpr), colnames(lbexpr)),colSideCols]
        levels(csc[,colSideCols]) <- cscCols
        csc <- as.matrix(csc)
    }

    temp <- factor(c(rep("u", length(group1)), rep("l", length(group2))), levels=c("u", "l"),
                   labels=c("red", "blue"))
    names(temp) <- c(group1, group2)
    temp <- temp[match(c(colnames(ubexpr), colnames(lbexpr)), names(temp))]
    if(!is.null(csc)) {
        csc <- cbind(csc, group=temp)
    } else {
        csc <- as.matrix(temp)
        colnames(csc) <- "Group"
    }
    rownames(de) <- geneID2Name(scd, rownames(de))
    upGene <- geneID2Name(scd, names(f1)[which(f1 > 1)])
    downGene <- geneID2Name(scd, names(f1)[which(f1 < 1)])
    
    return(list(heatmap=de, p.vals.all=p, p.vals.thr=p1, fold.change.all=f, fold.change.thr=f1,
                ColSideCols=csc, p.orig=p.orig, upGene = upGene, downGene = downGene))
}

##' List variables ls() on steroids
##'
##' List variables ls() on steroids
##' @title List variables ls() on steroids
##' @return List of variables and memory usage
##' @author wj241
##' @export
memls <- function() {
    z <- sapply(ls(.GlobalEnv), function(x) object.size(get(x)))
    as.matrix(rev(sort(z)))
}

##' Estimate size factor supercharged
##'
##' Estimate size factor supercharged
##' @title Estimate size factor supercharged
##' @param x Matrix of gene expression
##' @return Size factors
##' @author wj241
estSizeFact2 <- function(x) {
    expGenes <- !apply(x, 1, function(y) 0 %in% y)
    estSizeFact(x[expGenes,])
}

##' Identify expressed genes
##'
##' Identify expressed genes
##' @title Identify expressed genes
##' @param x Gene expression matrix.
##' @return Logical vector indicating whether gene is expressed.
##' @author wj241
expressedGene <- function(x) {
    rowSums(x) == 0
    ## !apply(x, 1, function (x) all(x==0))
}

##' Quality control function
##'
##' Quality control function
##' @title Quality control function
##' @param counts raw counts
##' @param ercc.data Spike in data
##' @param htseqQC HTseq output
##' @param geneTable Genetable
##' @param meta Metadata dataframe
##' @param metaLaneID Lane that contains flowcell lane information
##' @param mitochondrialIdenitifier Prefix for mitochondrial genes in geneTable
##' @param pdf Name of output pdf
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
##' @param orientation Default "landscape". Can also be "portrait".
##' @return List 
##' @author Wajid Jawaid
qcFunc <- function(counts, ercc.data, htseqQC, geneTable, meta, metaLaneID = "flowCell",
                   mitochondrialIdenitifier = "^mt|^MT", pdf = NULL,
                   cutoffs = c(2e5, .2, .2, 0, 0, 1, 0),
                   orientation = c("landscape", "portrait")) {
    orientation <- match.arg(orientation)
    mitochondrialGenesIDs <- geneTable[grep(mitochondrialIdenitifier,
                                            geneTable[,"Associated_Gene_Name"]),1]
    nuclearGeneIDs <- geneTable[-match(mitochondrialGenesIDs, geneTable[,1]),1]
    nERCC <- colSums(ercc.data)
    nQC <- colSums(htseqQC)
    nMitochondrialGenes <- colSums(counts[mitochondrialGenesIDs,])
    nNuclearGenes <- colSums(counts[c(nuclearGeneIDs),])
    nGenes <- nMitochondrialGenes + nNuclearGenes
    nTotalReads <- nGenes + nERCC + nQC
    nMappedReads <- nGenes + nERCC
    rpm <- t(t(1000000 * counts[nuclearGeneIDs,]) / nNuclearGenes)
    nHighCoverageReads <- apply(rpm, 2, function(x) sum(x > 10))

    qcPlots <- c(list(nGenes), split(htseqQC, rownames(htseqQC))[rownames(htseqQC)])
    qcPlots <- c(list(nMappedReads, nNuclearGenes), lapply(qcPlots, function(x) x / nTotalReads),
                 list(nHighCoverageReads), list(nMitochondrialGenes / nGenes),
                 list(nERCC / nMappedReads), list(nNuclearGenes / nMappedReads))
    names(qcPlots) <- c("nMappedReads", "nNuclearReads", "fGenes:Total",
                        paste0(rownames(htseqQC), ":Total"),
                        "nHighCoverageReads(>10rpm)", "mit:mappedToGene", "ercc:mapped",
                        "nuclearGenes:mapped")
    cells <- unique(meta[,metaLaneID])
    failQC <- vector("list", length(cells))
    names(failQC) <- cells

    for (cell in cells) {
        n <- meta[,metaLaneID] == cell
        n <- rownames(meta)[n]
        n2 <- rep(FALSE, ncol(counts))
        names(n2) <- colnames(counts)
        n2[n] <- TRUE
        n <- n2
        rm(n2)
        countThreshold <- 0
        ind <- as.character(cell)
        failQC[[ind]] <- names(which(nMappedReads[n] < cutoffs[1]))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(nNuclearGenes[n] < cutoffs[4])))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(qcPlots[["fGenes:Total"]][n] <= cutoffs[2])))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(qcPlots[["mit:mappedToGene"]][n] >= cutoffs[3])))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(qcPlots[["nHighCoverageReads(>10rpm)"]][n] <= cutoffs[5])))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(qcPlots[["ercc:mapped"]][n] >= cutoffs[6])))
        failQC[[ind]] <- c(failQC[[ind]],
                           names(which(qcPlots[["nuclearGenes:mapped"]][n] <= cutoffs[7])))
        failQC[[ind]] <- unique(failQC[[ind]])
        colQC <- names(n)[n]
        names(colQC) <- colQC
        colQC <- gsub(".*", "black", colQC)
        colQC[failQC[[ind]]] <- "red"
        nps <- length(qcPlots)
        npsRow <- floor(sqrt(nps))
        npsCol <- ceiling(nps / npsRow)
        if (orientation == "portrait") {
            temp <- npsRow
            npsRow <- npsCol
            npsCol <- temp
            rm(temp)
        }
        if(is.null(pdf)) {
            dev.new()
        } else {
            pdf(paste0(pdf, cell, ".pdf"), height = 4*npsRow, width = 4*npsCol)   
        }
        layout(matrix(c(rep(1,npsCol), 1 + (1:(npsRow * npsCol))), 1+ npsRow, npsCol, TRUE),
               heights = c(.25, rep(2,npsRow)))
        par(mar=c(0,0,0,0))
        plot(0, axes = FALSE, type = "n", xlab="", ylab="", , xlim=c(-1,1), ylim = c(-1,1))
        text(0,-0.5,paste(cell, ", failed =", length(failQC[[ind]])), cex=2)
        par(mar=c(5.1,4.1,1.1,1.1))
        for (i in 1:length(qcPlots)) {
            plot(nTotalReads[n], qcPlots[[i]][n], xlab="Total Reads", xaxt = "n",
                 ylab=names(qcPlots)[i], pch=20, col=colQC)
            axis(1, at = c(1e6, 5e6, 1e7), las=1)
            qnames <- names(qcPlots)[i]
            if (qnames == "mit:mappedToGene") abline(h=cutoffs[3], col = "red", lty=2)
            if (qnames == "fGenes:Total") abline(h=cutoffs[2], col = "red", lty=2)
            if (qnames == "nMappedReads") abline(h=cutoffs[1], col = "red", lty=2)
            if (qnames == "nNuclearReads") abline(h=cutoffs[4], col = "red", lty=2)
            if (qnames == "nHighCoverageReads(>10rpm)") abline(h=cutoffs[5], col = "red", lty=2)
            if (qnames == "ercc:mapped") abline(h=cutoffs[6], col = "red", lty=2)
            if (qnames == "nuclearGenes:mapped") abline(h=cutoffs[7], col = "red", lty=2)
        }
        if (!is.null(pdf)) dev.off()
    }
    exportMeta <- meta[,metaLaneID]
    names(exportMeta) <- rownames(meta)
    return(list(failQC=failQC, qcPlots=qcPlots, totalReads = nTotalReads, meta=meta[,metaLaneID,drop=FALSE]))
}

##' Generate a FACS like plot from single-cell RNAseq data
##'
##' Generate a FACS like plot from single-cell RNAseq data
##' @title Generate a FACS like plot from single-cell RNAseq data
##' @param scd Single Cell Dataset
##' @param gene1 Gene to plot on x-axis
##' @param gene2 Gene to plot on y-axis
##' @param t1 threshold 1
##' @param t2 threshold 2
##' @param filterCol perform only on cells with a unique measure on this feature
##' @param filterVal Linked to filterVal
##' @param nrpoints Deafult Inf. 
##' @param pch Default 0.3.
##' @param nbin Default 200.
##' @param bandwidth Default 0.1.
##' @param cex Default 0.3.
##' @param colramp See smoothScatter() function.
##' @param t.bold Default TRUE. To draw bold black lines for the thresholds.
##' @param t.lty Default 2. For supplementary threshold line type.
##' @param t.col Default "azure". Colour for supplementary threshold line.
##' @param ... Additional parameters passed to smoothScatter()
##' @return generates plot
##' @author Wajid Jawaid
##' @export
facsPlot <- function(scd, gene1, gene2, t1=NULL, t2=NULL, filterCol=NULL, filterVal=NULL,
                     nrpoints = Inf, pch = 16, nbin = 200, bandwidth = .1, cex = .3,
                     colramp=colorRampPalette(brewer.pal(9,"YlOrRd")), t.bold = TRUE,
                     t.lty = 2, t.col = "azure", ...) {
    if (!is.null(filterCol) && !is.null(filterVal)) {
        ind <- grep(filterVal, pData(scd)[,filterCol])
        selCells <- rownames(pData(scd))[ind]
    } else selCells <- rownames(pData(scd))
    origFilter <- saveFilters(scd)
    filterGene(scd) <- FALSE
    data <- exprs(scd)
    if (gene2 %in% colnames(pData(scd))) {
        data <- rbind(data, pData(scd)[,gene2])
        rownames(data)[length(rownames(data))] <- gene2
        genes <- c(geneName2id(scd, gene1), gene2)
    } else genes <- geneName2id(scd, c(gene1, gene2))
    fdat <- t(data[genes,selCells])
    scd <- restoreFilters(scd, origFilter)
    colnames(fdat) <- c(gene1, gene2)
    smoothScatter(fdat, nrpoints = nrpoints, pch = pch, nbin = nbin,
                  bandwidth = bandwidth, cex = cex, colramp = colramp,
                  ...)
    if (t.bold) abline(v = t1, h = t2)
    abline(v = t1, h = t2, lty = t.lty, col = t.col)
    cellNums <- table(fdat[,2]>t2, fdat[,1]>t1)
    rngs <- apply(fdat, 2, range)
    rngs <- rbind(rngs[1,], c(t1,t2), rngs[2,])
    x <- c(rngs[3,1], t1) - diff(rngs[,1]) / 20
    y1 <- (t2 + diff(rngs[,2]) / 10)[1]
    y2 <- (t2 - diff(rngs[,2]) / 10)[1]
    dts <- cbind(expand.grid(c(y1,y2),x), rev(as.vector(cellNums)))[,c(2,1,3)]
    text(dts[,1:2], labels = dts[,3])
    invisible(fdat)
}

##' Include all cells
##'
##' Include all cells
##' @title Include all cells
##' @param scd SCD
##' @return Updated SCD
##' @author wj241
##' @export
resetCells <- function(scd) {
    pD <- pData(phenoData(scd))
    pD[,"included"] <- TRUE
    pData(scd) <- pD
    scd
}

##' Differential expression using DESeq2
##'
##' Differential expression using DESeq2
##' @title Differential expression using DESeq2
##' @param scd SCD
##' @param group1 First group. Query group.
##' @param group2 Second group. Reference group.
##' @param fdr FDR
##' @param gene Vector of genes (ENSEMBL ids) to use for DE.
##' @param rmPattern Remove these genes
##' @param minCounts Remove genes with raw counts less than this
##' @param nMinCells The minimum number of cells that must
##' have the gene expressed (more than 0 counts) for the gene
##' to be included for differential expression.
##' @return List.
##' @author Wajid Jawaid
##' @importFrom DESeq2 DESeqDataSetFromMatrix
##' @export
doDESeq <- function(scd, group1, group2, fdr, gene= NULL, rmPattern = NULL, minCounts = 1,
                    nMinCells = ceiling(min(length(group1), length(group2))/2)) {
    cells <- c(group1, group2)
    if (is.null(gene)) {
        temp <- assayDataElement(scd, "counts")[, cells]   
    } else temp <- assayDataElement(scd, "counts")[gene, cells]
    temp2 <- pData(scd)[cells,]
    numCellsGeneExpr <- apply(temp, 1, function(x) sum(x>0))
    gene <- names(numCellsGeneExpr)[numCellsGeneExpr > nMinCells]
    temp <- temp[gene,]
    condition <- c(rep("test", length(group1)), rep("control", length(group2)))
    temp2$condition <- factor(condition)
    dds <- DESeq2::DESeqDataSetFromMatrix(temp, temp2, design = ~ condition)
    dds <- dds[rowSums(DESeq2::counts(dds)) > minCounts,]
    if (!is.null(rmPattern)) {
        dds <- dds[-grep(rmPattern, geneID2Name(scd, rownames(dds))), ]
    }
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, alpha=fdr)
    resOrdered <- res[order(res$padj),]
    resOrdered <- resOrdered[!is.na(resOrdered$padj),]
    upGeneID <- rownames(resOrdered)[resOrdered$padj < fdr & resOrdered$log2FoldChange >0]
    upGene <- geneID2Name(scd, upGeneID)
    downGeneID <- rownames(resOrdered)[resOrdered$padj < fdr & resOrdered$log2FoldChange < 0]
    downGene <- geneID2Name(scd, downGeneID)
    return(list(heatmap = assayDataElement(scd, "exprs")[c(upGeneID, downGeneID),
                                                         c(group1, group2)],
                ColSideColors = data.frame(Group = c(rep("red", length(group1)),
                                           rep("blue", length(group2))),
                                           row.names = c(group1, group2)),
                DESeq=resOrdered, upGene=upGene, downGene=downGene,
                backgroundGeneSet = gene))
}


## collapseSameName <- function(scd, name) {
##     return(scd)
## }


##' Finds cells from SCD object matching query term
##'
##' Finds cells from SCD object matching query term
##' @title Finds cells from SCD object matching query term
##' @param scd Single Cell Datast object.
##' @param pattern Grep pattern match.
##' @param column Column of pData to which pattern should be matched.
##' @return Names of identified cells.
##' @author wj241
##' @export
findCells <- function(scd, pattern, column) {
    rownames(pData(scd))[grep(pattern, pData(scd)[,column])]
}

##' Copy R content to clipboard
##'
##' Copy R content to clipboard
##' @title Copy R content to clipboard
##' @param x Vector of characters
##' @param sep Default = "\\n". Separator.
##' @param row.names Not used.
##' @param col.names Not used.
##' @return Copies vector to clipboard.
##' @author Wajid Jawaid
##' @export
clipboard <- function(x, sep="\n", row.names=FALSE, col.names=TRUE){
    if (any(grep("darwin", version$os))) {
        con <- pipe("pbcopy", open="w")   
    } else if (any(grep("linux", version$os))) {
        con <- pipe("xclip -selection clipboard -i", open="w")
    } else stop("Unknow operating system.")
     writeLines(x, con, sep=sep)
     close(con)
}

##' Write a list of gene vectors to file
##'
##' Write a list of gene vectors to file
##' @title Write a list of gene vectors to file
##' @param listofGenes Named list containing vectors of genes. Names
##' will be used as a header
##' @param file Filename for export
##' @return Returns data.frame of list with coloumn post-column binding.
##' @author Wajid Jawaid
##' @export
writeGeneList <- function(listofGenes, file) {
    maxL <- max(sapply(listofGenes, length))
    lG <- lapply(listofGenes, function(x) c(x, rep("", maxL - length(x))))
    lG <- do.call(cbind, lG)
    colnames(lG) <- names(listofGenes)
    ## lG <- rbind(colnames(lG), lG) # Required only if writing GMX file.
    write.table(lG, file=file, quote=FALSE, sep="\t", row.names = FALSE)
    invisible(lG)
}


##' Plots QC from SCD
##'
##' Plots QC from SCD
##' @title Plots QC from SCD
##' @param scd Single cell dataset
##' @param catG Assign colour to each cell. Named vector of colours with
##' names corressponding to cells.
##' @param ylims List object with a vector ylim (lower limit, upper limit)
##' in each entry of the list.
##' @return Generates plot.
##' @author wj241
##' @export
plotQC <- function(scd, catG = NULL, ylims = NULL) {
    qcPlots <- scd@qcOutput$qcPlots
    failQC <- do.call(c, scd@qcOutput$failQC)
    nTotalReads <- scd@qcOutput$totalReads
    colQC <- rep("black", length(nTotalReads))
    names(colQC) <- names(nTotalReads)
    colQC[failQC] <- "red"
    if (!is.null(catG)) colQC[names(catG)] <- catG
    nps <- length(qcPlots)
    npsRow <- floor(sqrt(nps))
    npsCol <- ceiling(nps/npsRow)
    layout(matrix(c(rep(1, npsCol), 1 + (1:(npsRow * npsCol))), 
                  1 + npsRow, npsCol, TRUE), heights = c(0.25, rep(2, 
                                                                   npsRow)))
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, type = "n", xlab = "", ylab = "", 
       , xlim = c(-1, 1), ylim = c(-1, 1))
    text(0, -0.5, paste(paste(names(scd@qcOutput$failQC), collapse = " and "),
                        ", failed =", length(failQC)), 
         cex = 2)
    par(mar = c(5.1, 4.1, 1.1, 1.1))
    for (i in 1:length(qcPlots)) {
        plot(nTotalReads, qcPlots[[i]], xlab = "Total Reads", 
             xaxt = "n", ylab = names(qcPlots)[i], pch = 20, 
             col = colQC, ylim = ylims[[i]])
        axis(1, at = c(1e+06, 5e+06, 1e+07), las = 1)
        qnames <- names(qcPlots)[i]
        ## if (qnames == "mit:mappedToGene") 
        ##     abline(h = cutoffs[3], col = "red", lty = 2)
        ## if (qnames == "fGenes:Total") 
        ##     abline(h = cutoffs[2], col = "red", lty = 2)
        ## if (qnames == "nMappedReads") 
        ##     abline(h = cutoffs[1], col = "red", lty = 2)
        ## if (qnames == "nNuclearReads") 
        ##     abline(h = cutoffs[4], col = "red", lty = 2)
    }   
}


##' Plots legend for plotGene() function
##'
##' Plots legend for plotGene() function when using multiple genes
##' @title Plots legend for plotGene() function
##' @param x output from plotgene
##' @param mpar plot parameters
##' @return generates plot
##' @author Wajid Jawaid
##' @export
plotLegend <- function(x, mpar = list(mar=c(4,4,2,2))) {
    data <- x$pointCol
    channelOrder <- x$channelOrder
    genes <- x$genes
    xgene <- x$genes[1]
    ygene <- x$genes[2]
    e.c <- attr(data, "e.c")
    xprMax <- rowMax(data)
    xprMax[xprMax == 0] <- 0.1
    ## xprMax <- ceiling(xprMax * 10) / 10
    ndata <- t(sapply(xprMax, function(x) seq(0,x,length=10)))
    xprSeq <- ndata
    ndata <- ((ndata / rowMax(ndata)) * (1 - e.c[1:nrow(ndata)])) + e.c[1:nrow(ndata)]
    ndata <- expand.grid(ndata[1,], ndata[2,])
    ndata[,3] <- 0
    ndata <- ndata[,channelOrder]
    ndata <- apply(ndata, 1, function(x) rgb(x[1], x[2], x[3]))
    ndata <- matrix(ndata, 10, 10)
    xStep <- mean(diff(xprSeq[1,]))
    yStep <- mean(diff(xprSeq[2,]))
    origPar <- par(no.readonly = TRUE)
    par(mpar)
    plot(NULL, xlim = c(0, xprMax[1] + xStep/2), ylim = c(0, xprMax[2] + yStep/2),
         axes = FALSE, xlab = xgene, ylab = ygene)
    axis(1); axis(2)
    xprSeq <- xprSeq - c(xStep, yStep)/2
    xprSeq <- cbind(xprSeq, xprSeq[,10] + c(xStep, yStep))
    for (i in 1:10) {
        for (j in 1:10) {
            rect(xprSeq[1,i], xprSeq[2,j], xprSeq[1,i+1], xprSeq[2,j+1], col = ndata[i,j],
                 border = "transparent")
        }
    }
    ## par(origPar)
}


##' Find highly variable genes (log linear model)
##'
##' Find highly variable genes (log linear model).
##' Model is fitted to:
##' \deqn{\frac{\sigma^{2}}{\mu^{2}} = a * \mu^{k}}
##' Using the non-linear squares method, if residualsInLogSpace is FALSE. Otherwise
##' using a linear model +/- a quadartic term.
##' @title Find highly variable genes (log linear model)
##' @param scd Single Cell Dataset
##' @param minMean Fit model to genes with mean expression gretaer than minMean
##' @param residualsInLogSpace Use lm instead of nls to optimise residuals in log-space.
##' @param quadratic Fit an order 2 model to log cv-squared against log mean
##' @param se The number of standard errors to return.
##' @inheritParams winsor
##' @return Returns a list of means, cv2, fit object and variable genes
##' @author Wajid Jawaid
##' @export
##' @importFrom stats nls coefficients qnorm predict.lm

logVarGenes <- function(scd, minMean = 0, fraction = 0.05, lower = FALSE,
                        residualsInLogSpace = TRUE, quadratic = TRUE,
                        se = qnorm(p = 0.975)) {
    if (fraction >= 1) {
        stop("Fraction must be [0,1)")
    } else if (fraction > 0) {
        x <- t(apply(counts(scd), 1, function(xx) winsor(xx, fraction, lower)))
    } else x <- counts(scd)
    gc()
    ind <- rowSums(x) == 0
    x <- x[!ind,]
    gc()
    mu = rowMeans(x)
    sigma2 = apply(x, 1, var)
    rm(x)
    cv2 <- sigma2 / mu^2
    
    wmu <- mu > minMean
    limmu <- mu[wmu]
    limcv2 <- cv2[wmu]
    lmu <- 10^(seq(min(log10(mu)), max(log10(mu)), length.out = 10))
    if (residualsInLogSpace) {
        lglimmu <- log10(limmu)
        lglimmu2 = lglimmu^2
        lglimcv2 <- log10(limcv2)
        lglmu = log10(lmu)
        lgmu <- log10(mu)
        
        gfit <- lm(lglimcv2 ~ lglimmu + lglimmu2)
        rlcv2 <- predict.lm(gfit, data.frame(lglimmu = lglmu, lglimmu2 = lglmu^2),
                            se.fit = TRUE)
        rpcv2 <- predict.lm(gfit, data.frame(lglimmu = lgmu, lglimmu2 = lgmu^2), se.fit = TRUE)
        browser()
        lcv2 <- 10^(rlcv2$fit)
        ucilcv2 <- 10^(rlcv2$fit + se * rlcv2$se)
        pcv2 <- 10^(rpcv2$fit + se * rpcv2$se)
    } else {
        gfit <- nls(limcv2 ~ a * limmu^k, start=c(a=20, k = 1))
        lcv2 <- coefficients(gfit)[["a"]] * lmu^coefficients(gfit)[["k"]]
        pcv2 <- coefficients(gfit)[["a"]] * mu^coefficients(gfit)[["k"]]
    }
    
    plot(mu, cv2, pch=20, log = "xy")
    abline(v = minMean, lty = 2, col = "blue")
    lines(lmu, lcv2, col = "orange", lwd = 2)
    if (residualsInLogSpace) lines(lmu, ucilcv2, col = "orange", lwd = 2, lty = 2)
    highVarGenes <- names(cv2)[cv2>pcv2]
    points(mu[highVarGenes], cv2[highVarGenes], col = "red", pch = 20)
    return(list(mean = mu, cv2 = cv2, fit = gfit, varGenes = highVarGenes))

}

##' Generate ggplot colors
##'
##' Generate ggplot colors
##' @title Generate ggplot colors
##' @param x Vector of group allocation
##' @param alpha Default 1. Measure of opacity.
##' @param returnFull Default TRUE. Return value for each vector entry.
##' @param remix Mix the colours so sequential colours are not identical.
##' @return Returns colors for each entry of the vector or for each group allocation.
##' @author Wajid Jawaid
##' @importFrom grDevices hcl
##' @export
ggCol <- function(x, alpha = 1, returnFull = TRUE, remix = TRUE) {
    if(length(x) > 1) {
        nCols <- length(sort(unique(x)))
    } else nCols <- x
    ggColors <- hcl(h = seq(15, 375 - 360/nCols, length.out = nCols) %% 360, c = 100,
                    l = 65, alpha = alpha)
    if (!is.factor(x) && (length(x) > 1)) x <- factor(x)
    if (remix) {
        if (!nCols%%2) {
            ggColors <- as.vector(t(matrix(ggColors, ncol = 2)))
        } else {
            ggColors <- c(ggColors[nCols], as.vector(t(matrix(ggColors[-nCols], ncol = 2))))
        }
    }
    if (returnFull && (length(x) > 1)) {
        retVal <- ggColors[x]
        names(retVal) <- x
    } else {
        retVal <- ggColors
        names(retVal) <- levels(x)
    }
    return(retVal)
}

##' Winsorization function
##'
##' Winsorization function.
##' Transformation limiting spurious extreme values.
##' Will set all outliers to a specific percentile of the data.
##' @title Winsorization function
##' @param x Numeric vector
##' @param fraction The fraction of cells on
##' @param lower By default this function will only winsorise the upper end.
##' If you also require the bottom end to be winsorised set this to TRUE.
##' @return Return winsorized vector
##' @author Wajid Jawaid
winsor <- function (x, fraction=0.05, lower = FALSE) {
    lim <- quantile(x, probs = c(fraction, 1-fraction)) 
    x[ x > lim[2] ] <- lim[2]
    if(lower) x[ x < lim[1] ] <- lim[1]
    x
    
}


## ##' Louvain clustering on transition matrix
## ##'
## ##' Louvain clustering on transition matrix
## ##' @title Louvain clustering on transition matrix
## ##' @param mkv Transition matrix
## ##' @return Returns a list with graph, dataframe and community object
## ##' @author Wajid Jawaid
## ##' @export
## ##' @importFrom igraph cluster_louvain communities
## findLouvain <- function(mkv) {
##     gph <- graph.adjacency(mkv, weighted = TRUE, mode = "undirected")
##     cll <- cluster_louvain(gph)
##     lvnClust <- communities(cll)
##     lvnClust <- lapply(1:length(lvnClust),
##                        function(x) cbind.data.frame(cell=lvnClust[[x]],
##                                                     class=as.character(x)))
##     lvnClust <- do.call(rbind, lvnClust)
##     rownames(lvnClust) <- lvnClust$cell
##     lvnClust
##     lvnClust[,1] <- NULL                # Remove duplicate name column
##     lvnClust <- lvnClust[rownames(mkv),,drop=FALSE]
##     return(list(gph=gph, clust=lvnClust, cll = cll))
##     ## lvnClust <- lvnClust[rownames(pData(scd)),,drop=FALSE]    
## }

##' Perform Louvain clustering
##'
##' Perform Louvain clustering
##' @title Perform Louvain clustering
##' @param scd Single Cell Dataset object
##' @param pcaDims Number of dimensions to use
##' @param nsig knn for automatic sigma calculation
##' @param d2 Optionally supply distance matrix
##' @param sim Optionally supply similiarity matrix
##' @param plotDims Number of dimensions to use on Fruchter-Reingold layout
##' of graph.
##' @param layoutIter Number of iterstions to run Fruchterman-Reingold graph 
##' @return List
##' @author Wajid Jawaid
##' @importFrom igraph layout_with_fr
##' @importFrom stats dist
##' @importFrom roots applyGaussianKernelwithVariableSigma
##' @importFrom roots diffuseMat sparseMarkov findLouvain
##' @export
clustLouvain <- function(scd, pcaDims = 50, nsig = 20, d2 = NULL, sim = NULL, plotDims = 3,
                         layoutIter = 1000) {
    if ((length(eigenvecs(getPCA(scd))) != 0) &&
        nrow(eigenvecs(getPCA(scd))) == nrow(pData(scd))) {
        cat("Using previoulsy calculated PCA.\n")
    } else {
        cat("Performing PCA ... ")
        scd <- runPCA(scd, scale. = TRUE)
        cat("Done.\n")
    }
    numDims <- pcaDims
    plot(eigenvals(getPCA(scd))*100, type = "l", ylab = "% variance explained",
         main = paste0("First ", numDims, " dims explain ",
                       format(sum(eigenvals(getPCA(scd))[1:pcaDims])*100, digits = 3),
                       "% of the variance."))
    abline(v = numDims, lty = 2, col = "red")
    ## y <- diffuseMat(exprs(scd), distfun = function(x) (1-cor(x))^2, nsig = 10,
    ##                  sqdistmat = cosineSqDistAll)
    if (is.null(d2)) {
        cat("Calculating distances ... ")
        d2 <- as.matrix((1-cor(t(eigenvecs(getPCA(scd))[,1:pcaDims]), method = "spearman"))^2)
    } else {
        cat("Using supplied distance matrix.")
    }
    if (is.null(sim)) {
        cat("Done.\nCalculating similarities ... ")
        sigmas <- apply(d2, 1, function(x) sqrt(sort(x)[nsig])/2)
        mkv <- applyGaussianKernelwithVariableSigma(d2, sigmas)
        rownames(mkv) <- colnames(mkv)
        diag(mkv) <- 0
    } else {
        cat("Using supplied simailarity/mkv matrix.")
        mkv <- sim
    }
    cat("Done.\nMake sparse ... ")
    mkv <- sparseMarkov(mkv, knn = nsig)
    mkv <- (mkv + t(mkv)) / 2
    ## sum(apply(mkv, 1, sd) == 0)
    ## sum(apply(mkv, 2, sd) == 0)
    cat("Done.\nFind communities ... ")
    lvnClust <- findLouvain(mkv)
    cat("Done.\nGenerate layout ... ")
    l <- layout_with_fr(lvnClust$gph, dim = plotDims, niter = layoutIter)
    cat("Done.\nAnnotate and Return.\n")
    rownames(l) <- rownames(pData(scd))
    colnames(l) <- paste("D", 1:plotDims)
    return(list(g = lvnClust, dist = d2, layout = l, plotDims = plotDims,
                layoutIter = layoutIter, mkv = mkv))
}
