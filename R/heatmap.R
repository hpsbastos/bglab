##' Heatmap function based on heatmap.2
##'
##' Heatmapl function based on heatmap.2
##' @title Heatmap function based on heatmap.2
##' @param x numeric matrix of the values to be plotted.
##' @param Rowv determines if and how the row dendrogram should be reordered.    By default, it is TRUE, which implies dendrogram is computed and reordered based on row means. If NULL or FALSE, then no dendrogram is computed and no reordering is done. If a dendrogram, then it is used "as-is", ie without any reordering. If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
##' @param Colv determines if and how the column dendrogram should be reordered.    Has the options as the Rowv argument above and additionally when x is a square matrix, Colv = \code{Rowv} means that columns should be treated identically to the rows.
##' @param distfun function used to compute the distance (dissimilarity) between both rows and columns. Defaults to \code{\link{dist}}.
##' @param hclustfun function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to \code{\link{hclust}}.
##' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured.
##' @param reorderfun character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured.
##' @param symm logical indicating if x should be treated symmetrically; can only be true when x is a square matrix.
##' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is \code{row} if symm false, and \code{none} otherwise.
##' @param na.rm logical indicating whether NA's should be removed.
##' @param revC logical indicating if the column order should be reversed for plotting, such that e.g., for the symmetric case, the symmetry axis is as usual.
##' @param add.expr expression that will be evaluated after the call to image. Can be used to add components to the plot.
##' @param breaks (optional) Either a numeric vector indicating the splitting points for binning x into colors, or a integer number of break points to be used, in which case the break points will be spaced equally between min(x) and max(x).
##' @param symbreaks Boolean indicating whether breaks should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise.
##' @param col colors used for the image. Defaults to heat colors (heat.colors).
##' 
##' @param colsep (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
##' @param rowsep (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
##' @param sepcolor (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
##' @param sepwidth (optional) Vector of length 2 giving the width (colsep) or height (rowsep) the separator box drawn by colsep and rowsep as a function of the width (colsep) or height (rowsep) of a cell. Defaults to c(0.05, 0.05)
##' @param cellnote     (optional) matrix of character strings which will be placed within each color cell, e.g. p-value symbols.
##' 
##' @param notecex     (optional) numeric scaling factor for cellnote items.
##' 
##' @param notecol (optional) character string specifying the color for cellnote text. Defaults to "green".
##' @param na.color Color to use for missing value (NA). Defaults to the plot background color.
##' @param trace character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'column'.
##' @param tracecol character string giving the color for "trace" line. Defaults to "cyan".
##' @param hline Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if trace is 'row' or 'both'. Vertical lines are only drawn if trace 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
##' @param vline Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if trace is 'row' or 'both'. Vertical lines are only drawn if trace 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
##' @param linecol Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if \code{\link{trace}} is 'row' or 'both'. Vertical lines are only drawn if \code{\link{trace}} is 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
##' @param margins numeric vector of length 2 containing the margins (see \code{\link{par}}(mar= *)) for column and row names, respectively.
##' @param ColSideColors (optional) character vector of length \code{\link{ncol}}(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
##' @param RowSideColors (optional) character vector of length \code{\link{nrow}}(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
##' @param cexRow positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
##' @param cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
##' @param labRow character vectors with row and column labels to use; these default to \code{\link{rownames}}(x) or \code{\link{colnames}}(x), respectively.
##' @param labCol character vectors with row and column labels to use; these default to \code{\link{rownames}}(x) or \code{\link{colnames}}(x), respectively.
##' @param srtRow angle of row/column labels, in degrees from horizontal
##' @param srtCol angle of row/column labels, in degrees from horizontal
##' @param adjRow 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation).
##' @param adjCol 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation).
##' @param offsetRow Number of character-width spaces to place between row/column labels and the edge of the plotting region.
##' @param offsetCol Number of character-width spaces to place between row/column labels and the edge of the plotting region.
##' @param key logical indicating whether a color-key should be shown.
##' @param keysize numeric value indicating the size of the key
##' @param density.info character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key.
##' @param denscol character string giving the color for the density display specified by density.info, defaults to the same value as tracecol.
##' @param symkey Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise.
##' @param densadj Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key. (See the \code{adjust} parameter for the \code{\link{density}} function for details.) Defaults to 0.25.
##' @param key.title character string for the key title
##' @param key.xlab character string for key's x-axis label
##' @param key.ylab character string for key's y-axis label
##' @param key.xtickfun function defining positions of ticks for key's x-axis
##' @param key.ytickfun function defining positions of ticks for key's y-axis
##' @param key.par defines the key's parameters e.g. margins
##' @param main main, x- and y-axis titles; defaults to none.
##' @param xlab main, x- and y-axis titles; defaults to none.
##' @param ylab main, x- and y-axis titles; defaults to none.
##' @param lmat visual layout: position matrix, column height, column width. See below for details
##' @param lhei visual layout: position matrix, column height, column width. See below for details
##' @param lwid visual layout: position matrix, column height, column width. See below for details
##' @param extrafun extra functions
##' @param key.cex size of key text
##' @param key.title.cex size of key title
##' @param ColSideColAxis changes axis on which the group names are printed. Use \code{2} for left and \code{4} for right.
##' @param axis.cex size of text in ColSideColAxis
##' @param ... additional parameters
##' @return Heatmap
##' @author Wajid Jawaid
##' @export
heatmap.minus <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
        w), symm = FALSE, scale = c("none", "row", "column"), 
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
    symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", 
    col = "heat.colors", colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("column", "row", "both", 
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors = NULL, RowSideColors = NULL, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, 
    key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, 
    key.par = list(), main = NULL, xlab = NULL, ylab = NULL, 
    lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, key.cex = NULL,
    key.title.cex = NULL,  ColSideColAxis = 2, axis.cex = NULL, ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none") && !is.null(breaks)) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
            nr)) 
            stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd > 
            nc)) 
            stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$rowName <- rownames(x)[rowInd]
    retval$colName <- colnames(x)[colInd]
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors)) 
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
            nc) 
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                      1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)) 
            stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                           1), 1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
	rsc = RowSideColors[rowInd, ,drop=FALSE]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(colnames(RowSideColors)) > 0) {
            axis(1, 0:(ncol(rsc)-1), colnames(RowSideColors), 
                las = 2, tick = FALSE, cex.axis = axis.cex)
        }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop = FALSE]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
            axis(ColSideColAxis, 0:(ncol(csc)-1),colnames(ColSideColors), 
                las = 2, tick = FALSE, cex.axis = axis.cex)            
            ## axis(ColSideColAxis, (1:(ncol(csc))/(ncol(csc))) - 1,colnames(ColSideColors), 
            ##     las = 2, tick = FALSE, cex.axis = axis.cex)
        }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!gtools::invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol, 
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(2, 2, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab)) 
            mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab)) 
            mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title)) 
            mar[3] <- 1
        if (is.null(key.cex)) pcex <- 1
        else pcex <- key.cex
        par(mar = mar, cex = pcex, mgp = c(2, 1, 0))
        if (length(key.par) > 0) 
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row") 
                key.xlab <- "Row Z-Score"
            else if (scale == "column") 
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, cex = key.cex)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                  0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Density"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = key.cex)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                  labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title)) 
                title(key.title, cex.main = key.title.cex)
            par(cex = key.cex)
            if (is.null(key.ylab)) 
                key.ylab <- "Count"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5, cex = key.cex)
        }
        else title("Color Key", cex.main = key.title.cex)
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun)) 
        extrafun()
    if (!is.null(ColSideColors)) retval$ColSideColors <- ColSideColors
    if (!is.null(RowSideColors)) retval$RowSideColors <- RowSideColors
    invisible(retval)
}

##' Heatmap helper function for leaf re-ordering
##'
##' Heatmap helper function for leaf re-ordering, the Gottgens way.
##' @title Heatmap helper function for leaf re-ordering
##' @param x Matrix
##' @param cor.method method to be used to calculate distance. Default is (1 - spearman rank correlation). This is overridden by \code{distfun}.
##' @param clust.method clustering method to be used. Default is \code{ward.D2}
##' @param distfun function to use for distance calculation.
##' @param col define colour. If defines please ensure that \code{breaks} are passed to \code{\link{heatmap.minus}}.
##' @param breaks colour breaks.
##' @param reorder Logical whether to reorder leaves
##' @param scale Scale by "row", "column" or "none"
##' @param ... Parameters passed to \code{\link{heatmap.minus}}
##' @return Calls \code{\link{heatmap.minus}} and returns list with re-ordering
##' @author Wajid Jawaid
##' @inheritParams heatmap.minus
##' @importFrom cba order.optimal
##' @export
heatmap.gottgens <- function(x, cor.method = "spearman", clust.method = "ward.D2",
                             distfun = NULL, col = NULL, breaks = NULL,
                             reorder = c("both", "row", "column", "none"),
                             Rowv = NULL, Colv = NULL,
                             scale = c("none", "row", "column"), trace = "none",
                             key.cex = .4, ColSideColAxis=4, symbreaks=FALSE, margins=c(3,3),
                             cexRow=0.6, keysize=0.5, key.title="Key",
                             key.par=list(mar=c(2.5,0.5,1.5,.5)), density.info="none",
                             ...) {
    scale <- match.arg(scale)
    reorder <- match.arg(reorder)
    hr <- hc <- list()
    if (reorder != "none") {
        if (is.null(distfun)) {
            if (reorder == "both" || reorder == "row") 
                dr <- as.dist(1 - cor(t(x), method = cor.method))
            if (reorder == "both" || reorder == "column") 
                dc <- as.dist(1 - cor(x, method = cor.method))
        } else {
            if (reorder == "both" || reorder == "row")             
                dr <- distfun(x)
            if (reorder == "both" || reorder == "column")             
                dc <- distfun(t(x))
        }

        if (reorder == "both" || reorder == "row") {
            hr <- hclust(dr, method = clust.method, members = NULL)
            rhr <- cba::order.optimal(dr, hr$merge)
            hr$merge <- rhr$merge
            hr$order <- rhr$order
            Rowv <- as.dendrogram(hr)
        }

        if (reorder == "both" || reorder == "column") {        
            hc <- hclust(dc, method = clust.method, members = NULL)
            rhc <- cba::order.optimal(dc, hc$merge)
            hc$merge <- rhc$merge
            hc$order <- rhc$order
            Colv <- as.dendrogram(hc)
        }
        

    }
    
    if (is.null(col)) {
        ## data(heatmapCols, package = "bglab")
        ## heatmapCols <- get("heatmapCols", pos = globalenv())
        col <- colorRampPalette(rev(heatmapCols))(1000)
        breaks <- seq(min(x),max(x), length.out=1001)
    }

    if (scale == "row" || scale == "column") breaks <- NULL

    retval <- heatmap.minus(x, Rowv = Rowv, Colv = Colv, col = col, breaks = breaks,
                            scale = scale, trace = trace, key.cex = key.cex,
                            ColSideColAxis = ColSideColAxis, symbreaks = symbreaks,
                            margins = margins, cexRow = cexRow, keysize = keysize,
                            key.title = key.title, key.par = key.par,
                            density.info = density.info, ...)
    retval$rowTree <- hr
    retval$colTree <- hc
    invisible(retval)
}

##' Generate heatmaps as subsets of a larger heatmap
##'
##' Generate heatmaps as subsets of a larger heatmap. Call selectionTray
##' before calling this function and select top.
##' left and bottom right corners to chose selection.
##' @title Subset Heatmap.
##' @param scd Single Cell Dataset object.
##' @param retHeat List returned from heatmap original heatmap.
##' @param sTdev The display on which selectionTray is displayed
##' @param pdf Filename to export pdf
##' @param boxLwid Line width of box showing selection
##' @param boxCol Colour of box line
##' @param data If the data to heatmap was altered, enter the new data as a matrix.
##' @return Produces a heatmap to select points and then a submap.
##' @author Wajid Jawaid
##' @export
submap <- function(scd, retHeat, sTdev=dev.cur(), pdf = NULL, boxLwid=3,boxCol="black", data=NULL) {
    if (is.null(data)) {
        hmap2 <- exprs(scd)[rev(retHeat$rowInd), retHeat$colInd]
    } else hmap2 <- data[rev(retHeat$rowInd), retHeat$colInd]
    if (!is.null(retHeat$ColSideColors)) {
        rsc2 <- retHeat$ColSideColors[retHeat$colInd]
    } else rsc2 <- NULL
    if (!is.null(rsc2)) {
        rsc2 <- matrix(rsc2)
        colnames(rsc2) <- "Stage"
    }
    dev.set(sTdev)
    cat("Choose your points.\n")
    sel <- do.call(cbind, locator(n=2))
    sel[sel<0] <- 0
    sel[sel>1] <- 1
    xRa <- round(sel[,1] * ncol(hmap2))
    yRa <- round(sel[,2] * nrow(hmap2))
    selCells <- colnames(hmap2)[xRa[1]:xRa[2]]
    selGenes <- rownames(hmap2)[yRa[1]:yRa[2]]
    ## rownames(hmap2) <- fData(scd)[match(rownames(hmap2), rownames(fData(scd))),
    ##                                "Associated_Gene_Name"]
    if (!is.null(rsc2)) rsc2 <- rsc2[xRa[1]:xRa[2],,drop=FALSE]
    cat("Generating submap.\n")
    if (is.null(pdf)) dev.new()
    else pdf(pdf)
    heatmap.gottgens(hmap2[yRa[1]:yRa[2], xRa[1]:xRa[2]], Rowv=FALSE, Colv=FALSE,
                     dendrogram = "none", trace = "none", cexRow = .5, cexCol = .5,
                     key.title = "Color Key", key.ylab = "", density.info = "none",
                     lmat = rbind(4:3,2:1), lhei = c(1,10), lwid = c(1,15),
                     key.cex = .2, key.title.cex = .2,
                     if(!is.null(rsc2)) ColSideColors = rsc2, margins = c(6,6),
                     ColSideColAxis = 4, axis.cex = 1, reorder = "none", col=retHeat$col)
    selectionTray(scd,retHeat, pdf=!is.null(pdf), data=data)
    rect(sel[1,1], sel[1,2], sel[2,1], sel[2,2], lwd=boxLwid, border=boxCol)
    if (!is.null(pdf)) {dev.off(); system(paste(fileOpen, pdf, sep=""))}
    dev.set(sTdev)
    invisible(x=list(cells=selCells, genes=selGenes, cellInd=xRa, geneInd=yRa))
}

##' Displays heatmap to choose area for submap
##'
##' Displays heatmap to choose area for submap
##' @title Selection Heatmap
##' @inheritParams submap 
##' @return Plots heatmap and returns the display.
##' @author Wajid Jawaid
##' @param pdf Will it be output to a pdf
##' @export
selectionTray <- function(scd, retHeat, pdf=FALSE, data=NULL) {
    if (is.null(data)) {
        hmap2 <- exprs(scd)[rownames(fData(scd)), rownames(pData(scd))]
        hmap2 <- hmap2[rev(retHeat$rowInd), retHeat$colInd]
    } else hmap2 <- data[rev(retHeat$rowInd), retHeat$colInd]
    if (is.null(retHeat$col)) {
        retHeat$col <- NULL
        retHeat$breaks <- NULL}
    if (!pdf) dev.new()
    par(mar=c(1,1,1,1))
    cat("Generating heatmap for selection.\n")
    image(seq(0,1,length.out=ncol(hmap2)), seq(0,1,length.out=nrow(hmap2)), t(hmap2),
          col = retHeat$col, breaks = retHeat$breaks, ylim=c(1,0), axes=FALSE, xlab = "",
          ylab = "")
    sTdev <- dev.cur()
    invisible(dev.cur())
}


##' Alters matrix rownames to only include genes of interest
##'
##' Alters matrix rownames to only include genes of interest
##' @title Alters matrix rownames to only include genes of interest
##' @param scd SCD object to use for translating gene names
##' @param data Matrix with genes in rows
##' @param genes vector of genes of interest
##' @return Matrix with only genes of interest in rownames, other genes contain empty string.
##' @author wj241
##' @export
markGOI <- function(scd, data, genes) {
    rn <- rownames(data)
    rn <- unname(geneID2Name(scd, rn))
    ind <- match(genes, rn)
    notfound <- genes[is.na(ind)]
    warning(paste("Genes:", paste(notfound, collapse=", "), "were not found on the heatmap."))
    ind <- ind[!is.na(ind)]
    rn[-ind] <- ""
    rownames(data) <- rn
    return(data)
}
