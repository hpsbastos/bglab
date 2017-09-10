## Creator: Wajid Jawaid
## Author: Wajid Jawaid
## Date: 2 December 2014
## Gottgens lab scRNA-seq tool repo

## Class definitions file


###############################################################################################
############################## Reduced Dimension (Virtual Class) ##############################
###############################################################################################

##' A virtual S4 class
##'
##' \code{reducedDim} is a virtual S4 class used by \linkS4class{PCA} and
##' \linkS4class{DiffusionMap} classes.
##' @slot eigenvalues A vector of eigen values
##' @slot eigenvectors A matrix with eigen vectors in columns.
##' @slot graph Contains the minimum spanning tree
##' @slot paths List of the various calculated paths
setClass("reducedDim", representation = representation(
                           eigenvalues = "numeric",
                           eigenvectors = "matrix",
                           graph = "ANY",
                           paths = "list",
                           "VIRTUAL"
                           )
         )



#####################################################################
## Other Reduced Dimensionsality (Inherits from Reduced Dimension) ##
#####################################################################

##' Other dimensionality redcution
##' @export
setClass("goggles", representation = representation(
                         "reducedDim",
                         l = "matrix",
                         clust = "list",
                         adj = "matrix",
                         dmat = "matrix",
                         pca = "matrix",
                         sparse = "numeric",
                         nadj = "matrix",
                         nadja = "matrix",
                         seed = "numeric"
                     ))

###############################################################################################
######################### PCA Class (Inherits from Reduced Dimension) #########################
###############################################################################################

##' Princincple Component Analysis on dataset
##'
##' Contains the PCA reduced data. Inheroted from the virtual class \linkS4class{reducedDim}
##' @slot standardDeviation Contains the calculated standard deviations
##' @slot rotation The transformation to rotate the data
##' @slot center Returned from \code{\link{prcomp}}
##' @slot scaled Logical value stores whether the data was scaled prior to PCA
##' @export
setClass("PCA", representation = representation(
                    "reducedDim",
                    standardDeviation = "numeric",
                    rotation = "matrix",
                    center = "ANY", # "numlogic",
                    scaled = "ANY" # "numlogic",
                    )
         )

###############################################################################################
#################### Diffusion Map Class (Inherits from Reduced Dimension) ####################
###############################################################################################

##' Diffusion Map dimensionality reduction on dataset
##'
##' Performed as designed by Florian Buettner
##' @slot numberConverged The number of converged eigenvalues and eigenvectors when using
##' "rARPACK"
##' @slot numberIterations The number of iterations to achieve convergence in "rARPACK".
##' @slot sigmas Value of sigma used for the Gaussian
##' @slot nn Number of nearest neighours used
##' @slot usedARPACK Logical indicating whether rARPACK was used for eigen deconstruction
##' @slot gaussian Post applying Gaussian kernel.
##' @slot markov Gaussian with rows normalised giving stochastic transition matrix.
##' @slot H H matrix.
##' @slot d2 Squared ditance matrix.
##' @slot distfun Function used to calculate distances.
##' @export
setClass("DiffusionMap", representation = representation(
                             "reducedDim",
                             numberConverged = "integer",
                             numberIterations = "integer",
                             sigmas = "numeric",
                             nn = "numeric",
                             usedARPACK = "logical",
                             gaussian = "matrix",
                             markov = "matrix",
                             H = "matrix",
                             d2 = "matrix",
                             distfun = "function"
                             )
         )

###############################################################################################
####################### TSNE Class (Inherits from Reduced Dimension) ##########################
###############################################################################################


##' TSNE dimensionality reduction on dataset
##'
##' Performed in Matlab - Barnes-Hut fast tSNE algorithm used
##' @slot theta Theta parameter
##' @slot perplexity Perplexity parameter
##' @slot N Number of cells/samples
##' @slot origD Original number of cells/samples
##' @slot seed Random seed used, saved for reproducibility
##' @slot nnError The tSNE nearest neighbour error score
##' @slot algorithm The algorithm used e.g. Barnes-Hut
##' @export
setClass("TSNE", representation = representation(
                     "reducedDim",
                     theta = "numeric",
                     perplexity = "numeric",
                     N = "integer",
                     origD = "integer",
                     seed = "numeric",
                     nnError = "numeric",
                     algorithm = "character"
                     )
         )

###############################################################################################
###################### Isomap Class (Inherits from Reduced Dimension) #########################
###############################################################################################


##' Isomap dimensionality reduced object
##'
##' Performed in Matlab
##' @slot nn Number of nearest neighbours used in the algorithm
##' @slot nnMat Adjacency matrix of nearest neighbours
##' @export
setClass("Isomap", representation = representation(
                     "reducedDim",
                     nn = "integer",
                     nnMat = "matrix"
                     )
         )

###############################################################################################
################################## Single Cell Dataset Class ##################################
###############################################################################################

##' Single Cell Dataset
##'
##' This is the main class in the \pkg{bglab} package. It stores all the data. Three datatables
##' are required to generate a new object. Both qPCR and RNAseq data may be used.
##' @slot useFilter Logical. To filter dataset or not
##' @slot filterGene Logical. To filter genes from dataset or not
##' @slot filterCell Logical. To filter cells from dataset or not
##' @slot pca Contains a \linkS4class{PCA} object
##' @slot diffMap Contains a \linkS4class{DiffusionMap} object
##' @slot tsne Contains a \linkS4class{TSNE} object
##' @slot isomap Contains a \linkS4class{Isomap} object
##' @slot counts Contains a matrix of raw counts.
##' @slot spike Spike-in matrix.
##' @slot qcCounts Quality Control matrix.
##' @slot filterQC Default TRUE. Will not return cells that failed QC when set to TRUE.
##' @slot normalise Either "DESeq" or "scran".
##' @slot qcOutput List containing QC output.
##' @slot technicalNoise List containing technical noise output.
##' @slot useForExprs Deprecated - only in for compatibility with scater.
##' Highly recommended.
##' @export
setClass("SCD",
         contains = "SCESet",
         representation = representation(
             useFilter = "logical",
             filterGene = "logical",
             filterCell = "logical",
             pca = "PCA",
             diffMap = "DiffusionMap",
             tsne="TSNE",
             isomap="Isomap",
             ## goggles = "goggles",
             counts="matrix",
             spike="matrix",
             qcCounts="matrix",
             filterQC="logical",
             normalise="character",
             qcOutput="list",
             technicalNoise="list",
             useForExprs="character"
             ),
         prototype = prototype(
             SCESet = new("SCESet", cellPairwiseDistances = dist(vector()),
                          featurePairwiseDistances = dist(vector())),
             useFilter = TRUE,
             filterGene = TRUE,
             filterCell = TRUE,
             pca = new("PCA"),
             diffMap = new("DiffusionMap"),
             tsne=new("TSNE"),
             isomap=new("Isomap"),
             ## goggles = new("goggles"),
             counts=new("matrix"),
             spike=new("matrix"),
             qcCounts=new("matrix"),
             filterQC = TRUE,
             normalise=character(),
             qcOutput = list(),
             technicalNoise = list(),
             useForExprs="exprs"
             ),
         validity = function(object) {
             isValid <- TRUE
             eType <- tolower(object@annotation)
             if (!any(eType == c("rnaseq", "qpcr"))) {
                 isValid <- FALSE
                 cat("Unknown experiment type. Please use \"RNAseq\" or \"qPCR\".\n")
             }
             ## data <- exprs(object)
             geneIds <- rownames(object)
             cellIds <- colnames(object)
             geneData <- split(t(fData(object)), f = colnames(fData(object)))
             whichCol <- lapply(geneData, function(x) !(any(is.na(match(geneIds, x)))))
             whichCol <- do.call(c, whichCol)
             if (!any(whichCol)) {
                 cat("Annotated gene data frame does not contain all gene IDs used in
                      assay data matrix")
                 isValid <- FALSE
             }
             return(isValid)
         }
         )

setValidity("SCD", function(object) {
    msg <- NULL
    valid <- TRUE
    if(valid) TRUE else msg
})


