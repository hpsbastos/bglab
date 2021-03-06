% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{colorGradient}
\alias{colorGradient}
\title{Helper function for \code{plotGene}}
\usage{
colorGradient(data, lcol = "grey", hcol = "red", na.col = "Black",
  opacity = NULL, maxVal = NULL, minVal = NULL, channelOrder = 1:3)
}
\arguments{
\item{data}{Gene expression data from \code{exprs} of single cell dataset object.}

\item{lcol}{Colour to be taken for low values. Default is "grey".}

\item{hcol}{Colour to be taken by high values. Default is "red".}

\item{na.col}{Colour for missing values, Default is Black.}

\item{opacity}{Opacity value from 0 to 255}

\item{maxVal}{Default NULL. Set the maximum and minimum color values manually.}

\item{minVal}{Default NULL. Set the maximum and minimum color values manually.}

\item{channelOrder}{Default 1:3 indicating red, green, blue.}
}
\value{
Returns a vector of colours for each point to be plotted.
}
\description{
Helper function for \code{plotGene}
}
\details{
Helper function for \code{plotGene} to calculate colour gradient
}
\author{
Wajid Jawaid
}
