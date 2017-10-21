## Creator: Wajid Jawaid
## Author: Wajid Jawaid
## Date: 2 December 2014
## Gottgens lab scRNA-seq tool repo

## Generics definitions file


##' Returns eigen values
##'
##' Returns eigen values
##' @title Eigen values
##' @param object Dimensionality reduced object
##' @return numeric vestor of ordered eigen values
##' @author Wajid Jawaid
##' @export
setGeneric("eigenvals", function(object) standardGeneric("eigenvals"))


##' Returns eigen vectors
##'
##' Returns eigen values
##' @title Eigen vectors
##' @param object Dimensionality reduced object
##' @return Annotated matrix of eigen vectors, with eigen vectors in columns.
##' @author Wajid Jawaid
##' @export
setGeneric("eigenvecs", function(object) standardGeneric("eigenvecs"))



##' Generates minmium spanning tree
##'
##' Generates minmium spanning tree and adds diameter path to \code{paths} slot
##' @title Minimum Spanning Tree for Reduced Dimensionality Object
##' @param object Dimensionality reduced object
##' @param ... Any additional parameters
##' @return Dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("mst", function(object, ...) standardGeneric("mst"))



##' Returns the minimum spanning tree
##'
##' Returns the minimum spanning tree
##' @title Retrieves the minimum spanning tree from the \code{graph} slot of a Single
##' Cell Dataset or a dimensionality reduced object
##' @param object Single Cell Dataset object or dimensionlaity reduced object
##' @param reducedDim The dimensionality reduced object
##' @return Minimum spanning tree as an igraph object
##' @author Wajid Jawaid
##' @export
setGeneric("getMST", function(object, reducedDim) standardGeneric("getMST"))



##' Retrieve paths from a dimensionality reduced object
##'
##' Retrieve paths from a dimensionality reduced object
##' @title Retrieve paths from a dimensionality reduced object
##' @param object Dimensionality reduced object
##' @return The stored paths for the object
##' @author Wajid Jawaid
##' @export
setGeneric("getPaths", function(object) standardGeneric("getPaths"))



##' 3D plotting function
##'
##' 3D plotting function
##' @title 3D plotting function
##' @param object A Reduced Dimension object 
##' @param colorBy Character or numeric vector. Character vectors will be converted to a factor
##' and cells coloured accordingly. Cell color can be supplied in the plotCols parameter. If not
##' supplied the ggplots2 coloring scheme will be adopted.
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
##' @param diamWidth Width of diameter path line
##' @param diamCol Color of diameter path line
##' @param ... Parameters to be passed to \code{plot3d} function.
##' @return 3D plot
##' @author Wajid Jawaid
setGeneric("plot3", function(object, colorBy = "", plotCols = NULL, legendSize = 1,
                             legendOffset = 0.2, legPlac = 1:3, legPos = "front",
                             negLeg = FALSE, plotLegend = TRUE, legVert = 1,
                             legTextPos = 1, axLabs = "Component", doPath = FALSE,
                             doTree = FALSE, selectedPath = NULL, useDims = 1:3,
                             selectedCells = c(NA, NULL), legendTitle = "",
                             diamWidth = 2, diamCol = "black", ...)
           standardGeneric("plot3"))



##' 2D/3D Path plotting function from package \code{bglab}
##'
##' 2D/3D Path plotting function from package \code{bglab}.
##' This function requires \code{plot3} to have been called
##' @title Plots given path
##' @param object Dimensionality reduced object
##' @param useDims Dimensions to be plotted
##' @param pointColor Color for each point
##' @param ... Parameters passed to plot
##' @return A path 
##' @author Wajid Jawaid
setGeneric("plotTree3", function(object, useDims, pointColor, ...) standardGeneric("plotTree3"))



##' 3D path plot
##'
##' 3D path plot
##' @title Plots path in 3D plot. Only to be called after \code{plot3}
##' @param object Dimensionality Reduced object
##' @param diamWidth Single length numeric controlling the diameter path line
##' @param diamCol Single value for color of the diameter path
##' @param useDims Dimensions to be used
##' @param selectedPath Name of path stored in paths slot
##' @param project2D Default NULL. Projects 3D graph into 2D by the given transformation
##' matrix. This matrix must have been produced using \code{getMVP} on an open rgl device.
##' @param ... Additional parameters to pass to plot
##' @return Plots 3D path
##' @author Wajid Jawaid
setGeneric("plotPath3", function(object, diamWidth, diamCol, useDims, selectedPath, project2D,
                                 ...) standardGeneric("plotPath3"))

##' Retrieve size factors
##'
##' Retrieve size factors
##' @title Retrieve size factors
##' @param object SCD
##' @param type Default "bio". For "spike-ins" use "tech"
##' @return Returns matrix of size factors
##' @author wj241
##' @export
setGeneric("sf", function(object, type="bio") standardGeneric("sf"))

##' Retrieve spike in counts
##'
##' Retrieve spike in counts
##' @title Retrieve spike in counts
##' @param object SCD
##' @param ... Additional parameters
##' @return Matrix of spike counts
##' @author Wajid Jawaid
##' @export
setGeneric("spikes", function(object, ...) standardGeneric("spikes"))



##' Retrieve QC counts
##'
##' Retrieve QC counts
##' @title Retrieve QC in counts
##' @param object SCD
##' @param ... Additional parameters
##' @return Matrix of QC counts
##' @author Wajid Jawaid
##' @export
setGeneric("qc", function(object, ...) standardGeneric("qc"))

##' To use filter or not
##'
##' To use filter or not. Replaces "noprocessing" argument from previous version
##' @title Add additional phenotype data
##' @param object Single Cell Dataset
##' @return Logical.
##' @author Wajid Jawaid
##' @export
setGeneric("useFilter", function(object) standardGeneric("useFilter"))

##' To use gene filter or not
##'
##' To use gene filter or not. 
##' @title Predicate indicating whether gene filtering is active
##' @param object Single Cell Dataset
##' @return Logical.
##' @author Wajid Jawaid
##' @export
setGeneric("filterGene", function(object) standardGeneric("filterGene"))

##' To use cell filter or not
##'
##' To use cell filter or not.
##' @title Predicate indicating whether cell filtering is active
##' @param object Single Cell Dataset
##' @return Logical.
##' @author Wajid Jawaid
##' @export
setGeneric("filterCell", function(object) standardGeneric("filterCell"))

##' To filter out cells that failed QC.
##'
##' To filter out cells that failed QC. It is highly recommended that this
##' is not disabled as size factor and technical noise estimation may fail.
##' @title Predicate indicating whether filtering of cells that failed QC is active
##' @param object Single Cell Dataset
##' @return Logical
##' @author Wajid Jawaid
##' @export
setGeneric("filterQC", function(object) standardGeneric("filterQC"))



##' Set whether to use filter or not
##'
##' Set whether to use filter or not
##' @title Set whether to use filter or not
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setGeneric("useFilter<-", function(object,value) standardGeneric("useFilter<-"))

##' Set whether to use filter or not
##'
##' Set whether to use filter or not
##' @title Set whether to use filter or not
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setGeneric("filterGene<-", function(object,value) standardGeneric("filterGene<-"))

##' Set whether to use filter or not
##'
##' Set whether to use filter or not
##' @title Set whether to use filter or not
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setGeneric("filterCell<-", function(object,value) standardGeneric("filterCell<-"))

##' Set whether to use QC or not
##'
##' Set whether to use QC or not
##' @title Set whether to use QC filter or not
##' @param object SCD object
##' @param value Logical value.
##' @return Updated SCD object
##' @author Wajid Jawaid
##' @export
setGeneric("filterQC<-", function(object,value) standardGeneric("filterQC<-"))

##' Save the size factors
##'
##' Save the size factors
##' @title Save the size factors
##' @param object SCD
##' @param type Can be "bio" or "tech"
##' @param value Matrix of values
##' @return UYpdated SCD object
##' @author wj241
##' @export
setGeneric("sf<-", function(object, type, value) standardGeneric("sf<-"))



##' Add additional phenotype data
##'
##' Add additional phenotype data
##' @title Add additional phenotype data
##' @param object Single Cell Dataset
##' @param value Value of pData, must have same length as the number of cells / samples.
##' @param noprocessing Default FALSE. Will return unfiltered/unordered results if set to TRUE
##' @return Updated Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setGeneric("addPheno", function(object, value, noprocessing=FALSE) standardGeneric("addPheno"))



##' Exclude selected cells
##'
##' Function to remove selected cells from analyses
##' @title Exclude selected cells
##' @param object Single Cell Dataset object
##' @param regexpr A regular expression of the cell names to be removed
##' @param cellNames A vector of cell names to be removed
##' @param invert Default FALSE. If TRUE include the given names and exclude other cells.
##' @param reset Default TRUE. If TRUE reset allow choice from all cells in dataset.
##' @return modified Single Cell Data object
##' @author Wajid Jawaid
##' @export
setGeneric("excludeCells", function(object, regexpr = NULL, cellNames = NULL, invert = FALSE,
                                    reset = TRUE)
    standardGeneric("excludeCells"))



##' Mark cells that have failed QC
##'
##' Mark cells that have failed QC
##' @title Mark cells that have failed QC
##' @param object Single Cell Dataset
##' @param cellNames Character vector of cell names
##' @param reset Default TRUE. Sets all others as PASSED.
##' @return Updated Single Cell Data object.
##' @author wj241
##, @export
setGeneric("markFailedQC", function(object, cellNames, reset = TRUE)
    standardGeneric("markFailedQC"))

##' To select only highly variable genes for further analysis
##'
##' To select only highly variable genes for further analysis. Give either genes to
##' be included or genes to be excluded BY ENSMBL id only. Do not give both.
##' @title To select only highly variable genes for further analysis
##' @param object Single Cell dataset object
##' @param includeGenes Genes to be included
##' @param excludeGenes Genes to be excluded
##' @param reset Default is FALSE. If TRUE will set selected genes as provided and all others
##' as the alternative. If FALSE will leave other genes unchanged.
##' @return Returns modififed Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setGeneric("selectVariableGenes", function(object, includeGenes=NULL, excludeGenes=NULL,
                                           reset = FALSE)
    standardGeneric("selectVariableGenes"))



##' Calculate PCA for Single Cell Dataset
##'
##' Calculate PCA for Single Cell Dataset
##' @title Principal Component Analysis
##' @param object Single Cell Dataset object
##' @param ... Parameters passed to \code{prcomp}
##' @return Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setGeneric("runPCA", function(object, ...) standardGeneric("runPCA"))



##' Retrieve PCA data
##'
##' Retrieve PCA data
##' @title Retrieves PCA object
##' @param object Single Cell Dataset object
##' @param ... Additional parameters
##' @return Returns PCA dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("getPCA", function(object, ...) standardGeneric("getPCA"))



##' Retrieve diffusion map object
##'
##' Retrieve diffusion map object
##' @title Retrieves diffusion map object
##' @param object Single Cell Dataset object
##' @param ... Additonal parameters
##' @return Returns DiffMap dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("getDiffMap", function(object, ...) standardGeneric("getDiffMap"))



##' Retrieve t-SNE object
##'
##' Retrieve t-SNE object
##' @title Retrieves t-SNE object
##' @param object Single Cell Dataset object
##' @param ... Additional parameters
##' @return Returns t-SNE dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("getTSNE", function(object, ...) standardGeneric("getTSNE"))



##' Retrieve Isomap object
##'
##' Retrieve Isomap object
##' @title Retrieves Isomap object
##' @param object Single Cell Dataset object
##' @param ... Additional parameters
##' @return Returns Isomap dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("getIsomap", function(object, ...) standardGeneric("getIsomap"))



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
##' @param ... Additional parameters to be passed to Rtsne
##' @return Returns an upadate Single Cell Dataset object
##' @author Wajid Jawaid
##' @export
setGeneric("runTSNE", function(object, ndims = 3, perplexity = 30, theta = 0.5, pca = FALSE,
                               seed = NULL, verbose = FALSE, use_dist = FALSE,
                               dist_fun = NULL, ...) standardGeneric("runTSNE"))



##' Diffusion map dimensionality reduction
##'
##' Diffusion map dimensionality reduction as performed by Florian Buettner / Fabian Theis
##' Diffusion map dimensionality reduction
##' @title Diffusion Maps
##' @param object Single Cell Data object
##' @param ndims Number of dimensions required if using ARPACK otherwise all will be returned
##' @param nn Number of nearest neighbours. Can be set to 1 
##' @param sigma Default NULL. Parameter for Gaussian kernel.
##' @param removeFirst Remove the non informative first eigen vector
##' @param useARPACK By default uses rARPACK but can be set to false to use the eigen
##' function from the base package.
##' @param distfun A distance function.
##' @param ... Other parameters to pass on.
##' @return Single Cell Data object with diffusion map
##' @author Wajid Jawaid
##' @export
setGeneric("diffuse", function(object, ndims = 4, nn = 0.2, sigma = NULL, removeFirst = TRUE,
                               useARPACK = TRUE, distfun = NULL, ...)
    standardGeneric("diffuse"))





##' Given a vector of cells, this will add a path to the paths slot.
##'
##' Given a vector of cells, this will add a path to the paths slot.
##' @title Store path
##' @param object Sngle cell dataset or dimensionality reduced object
##' @param reducedDim A dimensionality reduced object
##' @param pathName Name given to the path
##' @param pathDescription A decription of the path 
##' @param path A numeric vector of cells on the path
##' @param ... Parameters passed to the Method \code{addPath}
##' @return Dimensionality reduced object
##' @author Wajid Jawaid
##' @export
setGeneric("addPath", function(object, reducedDim, pathName, pathDescription, path, ...)
    standardGeneric("addPath"))





##' Finds and saves path between two cells
##'
##' Finds and saves path between two cells
##' @title Find path between two selected cells
##' @param object A single cell dataset object or a reduced dimensionality object
##' @param reducedDim Reduced diensionality object
##' @param from Numeric value representing the starting cell of path
##' @param to Numeric value representing the ending cell of path
##' @param pathName Name to give path
##' @param pathDescription Description of path
##' @return Altered object of same class
##' @author Wajid Jawaid
##' @export
setGeneric("findPath", function(object, reducedDim, from, to, pathName, pathDescription)
    standardGeneric("findPath"))



##' Generic function to remove saved paths
##'
##' Generic function to remove saved paths
##' @title Generic function to remove saved paths
##' @param object A singel cell dataset object
##' @param paths Paths to be removed. Can be numeric or character vector of path names.
##' If "ALL" then all paths will be removed.
##' @param ... Additional parameters
##' @return A single cell dataset object
##' @author Wajid Jawaid
##' @export
setGeneric("removePath", function(object, paths, ...)
    standardGeneric("removePath"))

##' Normalises and performs technical noise analysis for single cell RNAseq(Brennecke et al.)
##'
##' Normalises and performs a technical noise analysis for single cell RNAseq.
##' Give raw counts. Data is normalised using size-factor normalisation identical to DESeq2,
##' then a technical noise analysis is perfromed as described in brennecke et al.
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
setGeneric("techVar", function(object, useERCC = TRUE, cvThresh=.3, quant=.8,
                               minBiolDisp=.5^2, fdr=.1, data.lengths = NULL,
                               ercc.lengths = NULL, meanForFit=NULL)
    standardGeneric("techVar"))

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
##' @param ... Additional parameters
##' @return SCD object
##' @author Wajid Jawaid
##' @export
setGeneric("performQC", function(object, selectedCells = "ALL",
                   cutoffs = c(2e5, .2, .2, 0, 0, 1, 0),
                   metaLaneID = "flowCell", mitochondrialIdenitifier = "^mt|^MT",
                   pdf = NULL, qcFeatures = "ALL", ...) standardGeneric("performQC"))


