% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{clustLouvain}
\alias{clustLouvain}
\title{Perform Louvain clustering}
\usage{
clustLouvain(scd, pcaDims = 50, nsig = 20, dist = NULL, plotDims = 3,
  layoutIter = 1000)
}
\arguments{
\item{scd}{Single Cell Dataset object}

\item{pcaDims}{Number of dimensions to use}

\item{nsig}{knn for automatic sigma calculation}

\item{dist}{Provide distance matrix to speed up clustering}

\item{plotDims}{Number of dimensions to use on Fruchter-Reingold layout
of graph.}

\item{layoutIter}{Number of iterstions to run Fruchterman-Reingold graph}
}
\value{
List
}
\description{
Perform Louvain clustering
}
\details{
Perform Louvain clustering
}
\author{
Wajid Jawaid
}
