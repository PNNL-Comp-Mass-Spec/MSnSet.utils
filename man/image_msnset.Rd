% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_msnset.R
\name{image_msnset}
\alias{image_msnset}
\title{Visualize MSnSet}
\usage{
image_msnset(
  m,
  valueName = "value",
  facetBy = NULL,
  sOrderBy = NULL,
  fOrderBy = NULL,
  valRange = NULL,
  maxNRows = 50,
  tile_instead_of_raster = FALSE
)
}
\arguments{
\item{m}{MSnSet object}

\item{valueName}{name of the value to be displayed. Default - "value".}

\item{facetBy}{character of pheno data column containing factor value.
Default is NULL.}

\item{sOrderBy}{character of pheno data telling how to order samples}

\item{fOrderBy}{character of feature data telling how to order features}

\item{valRange}{number for the pseudocolor range from -valRange to +valRange.
Default is NULL. In that case the scale goes from 0.025 to 0.975
quantile.}

\item{maxNRows}{maximum number of rows to display. Default is 50.}

\item{tile_instead_of_raster}{logical: If TRUE, use \code{geom_tile} instead of 
\code{geom_raster}}
}
\description{
Visualize MSnSet
}
