% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_by_feature_property.R
\name{normalize_by_feature_property}
\alias{normalize_by_feature_property}
\title{Normalization of LC-MS/MS Data}
\usage{
normalize_by_feature_property(eset, property, method = c("lowess", "loess"),
  ...)
}
\arguments{
\item{eset}{ExpressionSet/MSnSet object}

\item{property}{character the column in the fData that the
relative intensities regress against.}

\item{method}{character either "lowess" or "loess" at this point.}

\item{...}{passed to \code{lowess} or \code{loess}.}
}
\description{
Fits a nonparametric trend into relative 
abundance data. Any significant (at this point just any)
trend is corrected. 
Converts one ExpressionSet/MSnSet to another ExpressionSet/MSnSet.
}
\note{
So far the only property I have in mind is elution time in 
     label-free LC-MS/MS data.
}
