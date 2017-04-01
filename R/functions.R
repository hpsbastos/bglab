## Author: Wajid Jawaid
## Date: 2 December 2014
## bglab package: Berthold Gottgens Lab scRNA-seq repo

##' @importFrom rARPACK eigs
##' @importFrom RColorBrewer brewer.pal
##' @import grDevices
##' @import graphics
##' @importFrom stats as.dendrogram as.dist cor density dist filter fitted.values gaussian hclust median order.dendrogram p.adjust pchisq prcomp predict quantile reorder sd var lm
##' @import utils
##' @import methods
##' @importFrom igraph graph.adjacency minimum.spanning.tree get.diameter get.shortest.paths get.edgelist E
##' @import ggplot2
##' @importFrom rgl plot3d par3d rgl.cur r3dDefaults open3d rgl.postscript text3d points3d segments3d lines3d selectpoints3d
##' @import Rtsne
##' @import Biobase
##' @importMethodsFrom BiocGenerics counts


require(Biobase, quietly = TRUE)
require(igraph, quietly = TRUE)
require(rgl, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(Rtsne, quietly = TRUE)


if (version$os == "linux-gnu") {
    fileOpen <- "xdg-open "
} else fileOpen <- "open "

##' Make new SCD object
##'
##' Construct a new Single cell dataset object for the bglab package
##' @title Make a new SCD object
##' @param experimentType Experiment Type. RNAseq or qPCR.
##' @param assayData Matrix of data with genes in rows and cells in columns
##' @param genoData Data frame of gene names with geneIds in rownames. Rownames must match rownames of assayData.
##' @param phenoData Data frame of further details for each cell. Rownames must match column names of assayData.
##' @param counts Raw counts data.
##' @param spike Default NULL. Include spike-in data e.g. ERCC-92.
##' @param qc Default NULL. Include quality control counts from HTseq for example.
##' @param ... Further parameters.
##' @return New Single Cell Dataset object.
##' @author Wajid Jawaid
##' @export
newSCD <- function(experimentType, assayData = NULL, genoData, phenoData, counts, spike = NULL,
                   qc = NULL, ...) {
    if (is.null(assayData)) assayData = counts
    new("SCD", experimentType = experimentType, exprs = assayData, assayData = assayData,
        genoData = genoData, phenoData = phenoData, counts = counts, spike = spike, qc = qc, ...)
}

subsetData <- function(x, pD, filters) {
    cellsInData <- attr(x, "cellsInData")
    filterCell <- filters[[1]] && filters[[3]]
    filterQC <- filters[[4]]
    nrp <- nrow(pD)
    ncx <- ncol(x)
    if (nrp == ncx) {
        if (filterCell && filterQC){
            x <- x[,(pD[,"included"] & pD[,"passedQC"])]
        } else if (!filterCell && filterQC) {
            x <- x[,pD[,"passedQC"]]
        } else if (filterCell && !filterQC) {
            x <- x[,pD[,"included"]]
        }
    } else {
        filterInc <- rownames(pD)[pD[,"included"]]
        filterPass <- rownames(pD)[pD[,"passedQC"]]
        filterBoth <- intersect(filterInc, filterPass)
        if (length(filter) != ncol(x)) {
            if (filterCell && filterQC){
                x <- x[,filterBoth]
            } else if (!filterCell && filterQC) {
                x <- x[,filterPass]
            } else if (filterCell && !filterQC) {
                x <- x[,filterInc]
            } else 
                stop("Internal error in subsetData().")   
        }
    }
    attr(x, "cellsInData") <- cellsInData
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

##' Fast vectorised Euclidean distance calculator
##'
##' Calculates Euclidean distances between vectors arranged as columns in a matrix.
##' @title Fast vectorised Euclidean distance calculator
##' @return Returns a matrix of pairwise distances
##' @author Wajid Jawaid
##' @param x Matrix with vectors in columns.
##' @param squared Will not perform the square root, i.e. will return the squared `L2-norm'.
##' @export
fastDist <- function(x, squared = FALSE) {
    a <- colSums(x^2)
    a <- a * matrix(1, ncol(x), ncol(x))
    a <- a + t(a)
    ab <- t(x) %*% x
    d <- a - 2 * ab
    diag(d) <- 0
    if (!squared) d <- sqrt(d)
    return(d)
}

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

##' Generic diffusion function
##'
##' Generic diffusion function
##' @title Generic diffusion function
##' @param data Matrix of data with genes in rows and cells in columns.
##' @param bluntDensityCorrection Blunt density correction
##' @param normBy Normalise by Column, row or both.
##' @param distfun A different distance function
##' @inheritParams diffuse
##' @return List conatining diffusion map.
##' @author Wajid Jawaid
diffuseMat <- function(data, ndims = 4, nn = 0.2, sigma = 12,
                       removeFirst = TRUE, bluntDensityCorrection = FALSE,
                       useARPACK = TRUE, normBy = c("Column", "Row", "Original"),
                       distfun = NULL) {
    normBy <- tolower(match.arg(normBy))
    nnp <- nn
    nn <- ceiling(ncol(data) * nn)      # Number of nearest neighbours to include
    KT <- sigma^2                       # Diffusion scale parameter
    cat("Calculating distance matrix. Please wait...")
    if (is.null(distfun)) {
        d2 <- as.matrix(fastDist(data, squared = TRUE))       # distance matrix calculation
    } else d2 <- as.matrix(distfun(data))
    cat("Done.\n")
    cat("Applying Gaussian kernel.")
    W <- exp(-d2 / (2*KT))              # Apply Gaussian kernel
    if (nnp < 1 & nnp > 0) {
        cat("Calculating and retaining nearest neighbours.\n")
        R <- apply(d2, 2, function(x) sort(x)[nn])
        R <- matrix(rep(R,ncol(d2)), ncol = ncol(d2)) # Find distance for nn closest neighbour
        W <- (d2<R) * W                     # Only keep nn nearest neighbours
    }
    cat("Calculating and applying local density correction.\n")
    D <- colSums(W)
    q <- D %*% t(D)                     # Calculate local density
    if (bluntDensityCorrection) q <- sqrt(q)
    diag(W) <- 0    
    H <- W / q                          # Correct for local density
    colS <- colSums(H)
    rowS <- rowSums(H)
    if (normBy[1]=="column") {
        eS <- colS == 0                 # Excluded cells with no transitions
        Hp <- t(t(H[!eS,!eS]) / colS[!eS])  # Normalise matrix   
    } else if (normBy[1]=="row") {
        eS <- rowS == 0
        Hp <- H[!eS,!eS] / rowS[!eS]
    } else if (normBy[1]=="original") {
        eS <- colS == 0
        Hp <- H[!eS,!eS] / colS[!eS]
    } else cat("No normalistion performed.\n")
    cat("Calculating eigen vectors. Please wait...")
    n <- nrow(d2)
    chooseDims <- function(E) {     # Sub-function to sort and select largest eigenvecs/vals
        if (ncol(E$vectors) <= ndims) ndims <- ncol(E$vectors) - 1
        eigOrd <- order(Re(E$values), decreasing = TRUE)
        E$values <- Re(E$values[eigOrd][startDim:(ndims+1)])
        E$vectors <- Re(E$vectors[,eigOrd][,startDim:(ndims + 1)]) # Remove first eigen vector
        rownames(E$vectors) <- colnames(data)[!eS]
        withNA <- matrix(NA, ncol(data), ncol(E$vectors),
                         dimnames=list(colnames(data), 1:ncol(E$vectors)))
        withNA[rownames(E$vectors),] <- E$vectors[rownames(E$vectors),]
        E$vectors <- withNA
        return(E)
    }
    if (useARPACK) {
        decomp <- rARPACK::eigs(Hp, which = "LR", k = ndims + 1)
        if (removeFirst) {
            startDim <- 2
        } else startDim <- 1
        decomp <- chooseDims(decomp)
        decomp$usedARPACK <- TRUE
    } else {   
        decomp <- eigen(Hp)                      # Eigen decomposition
        if (removeFirst) {
            startDim <- 2
        } else startDim <- 1
        decomp <- chooseDims(decomp)
        decomp$nconv <- integer(0)
        decomp$niter <- integer(0)
        decomp$usedARPACK <- FALSE
    }
    if (length(which(eS)) != 0) cat(paste("\nWarning:\nCells \"", paste(names(which(eS)), collapse = ", "), "\" have no transitions and have been removed from the analysis\n", sep=""))
    decomp$nn <- nnp
    return(decomp)
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
##' @return Returns a list
##' @importFrom statmod glmgam.fit
##' @author Wajid Jawaid
techVarSub <- function(object, data.ercc = NULL, cvThresh=.3, quant=.8, minBiolDisp=.5^2,
                       fdr=.1, data.lengths = NULL, ercc.lengths = NULL, meanForFit=NULL) {
    tData <- counts(object)
    sf.data <- estSizeFact2(tData)
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

##' Look up available databases on Enrichr
##'
##' Look up available databases on Enrichr
##' @title Look up available databases on Enrichr
##' @return dataframe of available Enrichr databases
##' @author Wajid Jawaid
##' @importFrom httr GET POST
##' @importFrom rjson fromJSON
##' @export
listEnrichrDbs <- function() {
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    dbs <- GET(url="http://amp.pharm.mssm.edu/Enrichr/datasetStatistics")$content
    dbs <- intToUtf8(dbs)
    dbs <- fromJSON(dbs)
    dbs <- lapply(dbs$statistics, function(x) do.call(cbind.data.frame, x))
    dbs <- do.call(rbind.data.frame, dbs)
    options(stringsAsFactors = dfSAF)
    dbs
}

##' Gene enrichment using Enrichr
##'
##' Gene enrichment using Enrichr
##' @title Gene enrichment using Enrichr
##' @param genes Character vector of gene names or dataframe of gene names in
##' in first column and a score between 0 and 1 in the other.
##' @param databases Character vector of databases to search.
##' See http://amp.pharm.mssm.edu/Enrichr/ for available databases.
##' @return Returns a data frame of enrichment terms, p-values, ...
##' @author Wajid Jawaid
##' @importFrom httr GET POST
##' @importFrom rjson fromJSON
##' @export
enrichr <- function(genes, databases = NULL) {
    if (is.null(databases)) {
        dbs <- c("ChEA 2015", "Epigenomics Roadmap HM ChIP-seq",
         "ENCODE and ChEA Consensus TFs from ChIP-X",
         "TF-LOF Expression from GEO", "ENCODE Histone Modifications 2015",
         "Transcription Factor PPIs", "KEGG 2016", "WikiPathways 2016", "CORUM",
         "SILAC Phosphoproteomics", "Humancyc 2016", "NCI-Nature 2016", "Panther 2016",
         "GO Biological Process 2015", "GO Cellular Component 2015",
         "GO Molecular Function 2015",
         "MGI Mammalian Phenotype Level 3", "MGI Mammalian Phenotype Level 4",
         "Human Phenotype Ontology", "OMIM Disease", "OMIM Expanded",
         "Mouse Gene Atlas", "Human Gene Atlas", "Cancer Cell Line Encyclopedia",
         "ESCAPE")
        databases <- gsub(" ", "_", dbs)
    }
    cat("Uploading data to Enrichr... ")
    if (is.vector(genes)) {
        temp <- POST(url="http://amp.pharm.mssm.edu/Enrichr/enrich",
                     body=list(list=paste(genes, collapse="\n")))
    } else if (is.data.frame(genes)) {
        temp <- POST(url="http://amp.pharm.mssm.edu/Enrichr/enrich",
                     body=list(list=paste(paste(genes[,1], genes[,2], sep=","),
                                          collapse="\n")))
    } else {
        warning("genes must be a vector of gene names or a dataframe with genes and score.")
    }
    GET(url="http://amp.pharm.mssm.edu/Enrichr/share")
    cat("Done.\n")
    dbs <- as.list(databases)
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- lapply(dbs, function(x) {
        cat("  Querying ", x, "... ", sep="")
        r <- GET(url="http://amp.pharm.mssm.edu/Enrichr/export",
                 query=list(file="API", backgroundType=x))
        r <- intToUtf8(r$content)
        tc <- textConnection(r)
        r <- read.table(tc, sep = "\t", header = TRUE, quote = "")
        close(tc)
        cat("Done.\n")
        return(r)
    })
    options(stringsAsFactors = dfSAF)
    cat("Parsing results... ")
    names(result) <- dbs
    cat("Done.\n")
    return(result)
}

##' Print Enrichr output.
##'
##' Print Enrichr output to text file.
##' @title Print Enrichr output to text file.
##' @param data Output from Enrichr function.
##' @param file Name of output file.
##' @param sep Default TAB. How to separate fields.
##' @param columns Columns from each entry of data.
##' 1-"Index", 2-"Name", 3-"Adjusted_P-value", 4-"Z-score"         
##' 5-"Combined_Score", 6-"Genes", 7-"Overlap_P-value" 
##' @return Produces file.
##' @author Wajid Jawaid
##' @export
printEnrich <- function (data, file, sep = "\t", columns = c(2,3,6)) {
    enrich <- file(file, "w")
    for (i in 1:length(data)) {
        writeLines(names(data)[i], enrich)
        n <- nrow(data[[i]])
        if (n > 10) n <- 10
        if (n > 0) {
            writeLines(paste(apply(data[[i]][1:n, columns, drop=FALSE], 1,
                                   function(x) paste(x, collapse = sep)),
                             collapse = "\n"), enrich)
        writeLines("\n", enrich)
        }
    }
    close(enrich)
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
##' @return List 
##' @author Wajid Jawaid
qcFunc <- function(counts, ercc.data, htseqQC, geneTable, meta, metaLaneID = "flowCell",
                   mitochondrialIdenitifier = "^mt|^MT", pdf = NULL,
                   cutoffs = c(2e5, .2, .2, 0, 0, 1, 0)) {
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
        temp <- scd@counts[, cells]   
    } else temp <- scd@counts[gene, cells]
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


##' Generic diffusion function using automated
##' individualised sigma calculation
##'
##' Generic diffusion function using automated individualised sigma calculation.
##'
##' A Gaussian kernel is applied to the chosen distance metric producing
##' an \eqn{n \times n} square unnormalised symmetric transition matrix, \eqn{A}.
##' Let \eqn{D} be an \eqn{n \times n} diagonal matrix with row(column) sums of
##' \eqn{A} as entries. The density corrected transition matrix will now
##' be:
##' 
##' \deqn{D^{-1} A D^{-1}}{D^{-1} * A * D^{-1}}
##'
##' and can be normalised:
##' 
##' \deqn{B^{-1} D^{-1} A D^{-1}}{B^{-1} * D^{-1} * A * D{^-1}}
##'
##' where \eqn{B} is an \eqn{n \times n} diagonal matrix with row sums of
##' the density corrected transition matrix as entries. The eigen decomposition of
##' this matrix can be simplified by solving the symmetric system:
##'
##' \deqn{B^{-\frac{1}{2}} D^{-1} A D^{-1} B^{-\frac{1}{2}} R^\prime = %
##'       R^\prime \lambda^\prime}{%
##'       B^{-0.5} * D^{-1} * A * D^{-1} * B^{-0.5} = R^\prime * L^\prime}
##'
##' where \eqn{R^\prime}{R^\prime} is a matrix of the right eigenvectors that solve
##' the system and \eqn{\lambda^\prime}{L'} is the corresponding eigenvalue
##' diagonal matrix. Now the solution of:
##'
##' \deqn{B^{-1} D^{-1} A D^{-1} R = R \lambda}{%
##'       B^{-1} * D^{-1} * A * D^{-1} * R = R * L}
##'
##' in terms of \eqn{R^\prime} and \eqn{B^{-\frac{1}{2}}} is:
##'
##' \deqn{B^{-1} D^{-1} A D^{-1} B^{-\frac{1}{2}} R^\prime = %
##'       B^{-\frac{1}{2}} R^\prime \lambda^\prime}{%
##'       B^{-1} * D^{-1} * A * D^{-1} * B^{-0.5} = B^{-0.5} * R' * L'}
##'
##' and
##'
##' \deqn{R = B^{-\frac{1}{2}} R^\prime}{R = B^{-0.5} * R'}
##'
##' This \eqn{R} without the first eigen vector is returned as the diffusion map.
##' @title Generic diffusion function
##' @param data Matrix of data with genes in rows and cells in columns.
##' @param nsig For automatic sigma calculation
##' @param distfun A different distance function that returns the \strong{squared}
##' distance
##' @param sigmas Manually provide sigma
##' @param sqdistmat \strong{Squared} distance matrix.
##' Give your own squared distance matrix.
##' @inheritParams diffuse
##' @return List output containing:
##' \tabular{rl}{%
##'   \emph{values} \tab Eigenvalues, excluding the first eigenvalue, which should%
##'                       always be 1.\cr
##'    \emph{vectors} \tab Matrix of eigen vectors in columns, first eigen vector%
##'                        removed.\cr
##'    \emph{nconv} \tab Number of eigen vectors/values that converged.\cr
##'    \emph{niter} \tab Iterations taken for Arnoldi algorithm to converge.\cr
##'    \emph{nops} \tab  Number of operations. \cr
##'    \emph{val0} \tab 1st eigen value - should be 1. If not be suspicious!\cr
##'    \emph{vec0} \tab 1st eigen vector - should be \eqn{n^{-\frac{1}{2}}{1/sqrt(n)}},%
##'                      where n is the number of cells/samples.\cr
##'    \emph{usedARPACK} \tab Predicates use of ARPACK for spectral decomposition.\cr
##'    \emph{distfun} \tab Function used to calculate the squared distance.\cr
##'    \emph{nn} \tab Number of nearest neighbours used for calculating \code{sigmas}.\cr
##'    \emph{d2} \tab Matrix of squared distances, returned from \code{distfun}.\cr
##'    \emph{sigmas} \tab Vector of sigmas. Same length as number of cells if individual\cr
##'                  \tab sigmas were calculated, otherwise a scalar if was supplied.\cr
##'    \emph{gaussian} \tab Unnormalised transition matrix after applying Gaussian.\cr
##'    \emph{markov}  \tab Normalised \code{gaussian} matrix.\cr
##'    \emph{densityCorrected} \tab Matrix after applying density correction to %
##'                                 \code{markov}.\cr
##' }
##' @author Wajid Jawaid
##' @references
##' Haghverdi, L., Buettner, F., Theis, F.J., 2015. Diffusion maps for high-dimensional single-cell analysis of differentiation data. Bioinformatics 31, 2989–2998. doi:10.1093/bioinformatics/btv325
##'
##' ##' Haghverdi, L., Büttner, M., Wolf, F.A., Buettner, F., Theis, F.J., 2016. Diffusion pseudotime robustly reconstructs lineage branching. Nat Meth 13, 845–848. doi:10.1038/nmeth.3971
##'
##' Angerer, P., Haghverdi, L., Büttner, M., Theis, F.J., Marr, C., Buettner, F., 2016. destiny : diffusion maps for large-scale single-cell data in R. Bioinformatics 32, 1241–1243. doi:10.1093/bioinformatics/btv715
##' @export
diffuseMat2 <- function(data, ndims = 20, nsig = 5,
                        removeFirst = TRUE, useARPACK = TRUE,
                        distfun = NULL, sigmas = NULL, sqdistmat = NULL) {
    if (!is.null(sqdistmat)) {
        d2 <- sqdistmat
    } else {
        cat("Calculating distance matrix. Please wait... ")
        if (is.null(distfun)) {
            d2 <- as.matrix(fastDist(data, squared = TRUE))       # distance matrix calculation
        } else d2 <- as.matrix(distfun(data))
        cat("Done.\n")
    }
    if (is.null(sigmas)) {
        cat("Calculating sigmas. Please wait... ")
        sigmas <- apply(d2, 1, function(x) sqrt(sort(x)[nsig])/2)
        cat("Done.\n")
        cat("Applying Gaussian kernel.\n")
        W <- applyGaussianKernelwithVariableSigma(d2, sigmas)
    } else {
        W <- exp(-d2 / (2*sigmas^2))
    }
    diag(W) <- 0
    cat("Calculating normalised transition matrix.\n")
    markov <- W / rowSums(W)
    cat("Performing density correction.\n")
    D <- rowSums(W)
    q <- D %*% t(D)                     # Calculate local density
    H <- W / q                          # Correct for local density
    dH <- rowSums(H)
    ## cat("Calculating density corrected normalised transition matrix.\n")
    ## markovDensityCorrected <- H / dH
    cat("Calculating related symmetric system.\n")
    rdH <- 1/sqrt(dH)
    Hp <- H * (rdH %*% t(rdH))
    cat("Calculating eigen vectors for symmetric system. Please wait... ")
    n <- nrow(d2)
    # Sub-function to sort and select largest eigenvecs/vals
    if (useARPACK) {                    # Eigen decomposition
        decomp <- rARPACK::eigs(Hp, which = "LR", k = ndims + 1)
    } else {   
        decomp <- eigen(Hp)
    }
    colnames(decomp$vectors) <- paste0("DC", 0:(ncol(decomp$vectors)-1))
    cat(" Done.\nTransforming to eigen vectors of normalised transition matrix ... ")
    decomp$vectors <- decomp$vectors * rdH
    rownames(decomp$vectors) <- colnames(data)
    ## decomp$vectors <- decomp$vectors / sqrt(colSums(H))
    cat(" Done.\nChoosing dims... ")

    if (removeFirst) {
        startDim <- 2
        decomp$val0 <- Re(decomp$values[1])
        decomp$vec0 <- Re(decomp$vectors[,1])
        decomp$values <- Re(decomp$values[-1])
        decomp$vectors <- Re(decomp$vectors[,-1])
    } else startDim <- 1
    cat(" Done.\n")
    if (!useARPACK) decomp$nconv <- decomp$niter <- integer(0)
    decomp$usedARPACK <- useARPACK
    decomp$distfun <- distfun
    decomp$nn <- nsig
    decomp$d2 <- d2
    decomp$sigmas <- sigmas
    rownames(W) <- colnames(W)
    decomp$gaussian <- W
    rownames(markov) <- colnames(markov)
    decomp$markov <- markov
    rownames(H) <- colnames(H)
    decomp$densityCorrected <- H
    return(decomp)
}

##' Predicts diffusion map projection from new data points
##'
##' Predicts diffusion map projection from new data points
##' @title Predicts diffusion map projection from new data points
##' @param dm Output from diffuseMat2 function
##' @param x Matrix of new data points. Features in rows and cells in
##' columns.
##' @param data Original data used to generate diffusion map
##' @param distfun A distance function that takes new data as first
##' paramter and previous data as second variable returning a squared
##' distance measure, with each sample in the rows and distance to
##' previous data points in columns, e.g. function(x, y) (1 - cor(x, y))^2.
##' @return Returns a matrix with projected diffusion components.
##' @author Wajid Jawaid
##' @export
diffuseProj <- function(dm, x, data, distfun) {
    ## Calculate new distances
    d2 <- distfun(x, data)             # dim(d2) is ncols(x) by ncols(data)
    if (length(dm$sigmas) == 1) {
        W <- exp(-d2 / (2*dm$sigmas^2))    # Gaussian Kernel
    } else {
        #Find distance to nsig cell then / 2.
        sigmas <- apply(d2, 1, function(x) sqrt(sort(x)[dm$nn])/2)
        W <- applyGaussianKernelwithVariableSigma(d2, sigmas, dm$sigmas)
    }
    ## diag(W) <- 0                       # No self-transitions
    W[which(W == 1)] <- 0
    ## Perform density normalisation by dividing by the product of the sums of
    ## incoming and outgoing weights
    H <- W / rowSums(W)
    H <- t(t(H) / rowSums(dm$gaussian))
    ## Make markovian
    Hp <-  H / rowSums(H)
    predMat <- Hp %*% cbind(dm$vec0, dm$vectors) %*% diag(c(1, 1/dm$values))
    colnames(predMat) <- paste0("DC", 0:(ncol(predMat)-1))
    rownames(predMat) <- colnames(x)
    predMat[,-1, drop=FALSE]
}

##' Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##'
##' Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##' @title Apply Gaussian Kernel using Laleh Haghverdi's variable sigma
##' @param d2 Squared distance metric
##' @param rsigmas Sigmas for cells in the rows
##' @param csigmas Sigmas for cells in the columns
##' @return Returns matrix of same size as d2.
##' @author Wajid Jawaid
applyGaussianKernelwithVariableSigma <- function(d2, rsigmas, csigmas = NULL) {
    if (is.null(csigmas)) {
        sigmaMat <- rsigmas %*% t(rsigmas)
        sigmaSqd <- rsigmas^2
        sigmaSum <- expand.grid(1:length(sigmaSqd), 1:length(sigmaSqd))
        sigmaSum <- sigmaSqd[sigmaSum[,1]] + sigmaSqd[sigmaSum[,2]]
        sigmaSum <- matrix(sigmaSum, length(sigmaSqd))
    } else {
        sigmaMat <- rsigmas %*% t(csigmas)
        rsigmaSqd <- rsigmas^2
        csigmaSqd <- csigmas^2
        sigmaSum <- expand.grid(1:length(rsigmaSqd), 1:length(csigmaSqd))
        sigmaSum <- rsigmaSqd[sigmaSum[,1]] + csigmaSqd[sigmaSum[,2]]
        sigmaSum <- matrix(sigmaSum, length(rsigmaSqd))
    }
    W <- sqrt(2*sigmaMat / sigmaSum) * exp(-d2 / sigmaSum)
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
logVarGenes <- function(scd, minMean = 0, fraction = 0.05, lower = FALSE, residualsInLogSpace = TRUE, quadratic = TRUE, se = qnorm(p = 0.975)) {
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
    lmu <- seq(min(mu), max(mu), length.out = 10)
    if (residualsInLogSpace) {
        lglimmu <- log10(limmu)
        lglimcv2 <- log10(limcv2)
        lglmu = log10(lmu)
        lglimmu2 = lglimmu^2
        lgmu <- log10(mu)
        gfit <- lm(lglimcv2 ~ lglimmu + lglimmu2)
        rlcv2 <- 10^(predict.lm(gfit, data.frame(lglimmu = lglmu, lglimmu2 = lglmu^2),
                               se.fit = TRUE))
        rpcv2 <- 10^(predict.lm(gfit, data.frame(lglimmu = lgmu, lglimmu2 = lgmu^2)))
        lcv2 <- rlcv2$fit
        ucilcv2 <- rlcv2$fit + se * rlcv2$se
        pcv2 <- rpcv2$fit + se * rpcv2$se
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
        return(ggColors[x])
    } else return(ggColors)
    
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
