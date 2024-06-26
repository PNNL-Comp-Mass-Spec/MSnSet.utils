% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAUC.R
\name{plotAUC_gg}
\alias{plotAUC_gg}
\title{Plot AUC (with ggplot2)}
\usage{
plotAUC_gg(
  modelingResult,
  CI = FALSE,
  rectilinear = FALSE,
  no_numeric_policy = c("warning", "plot_blank", "error"),
  seed = 0
)
}
\arguments{
\item{modelingResult}{output of \code{lr_modeling} or \code{rf_modeling}.}

\item{CI}{(logical) whether to plot confidence intervals. Default is
\code{FALSE}.}

\item{rectilinear}{(logical) whether to prevent diagonal lines being formed from jumps
in TPR from CI boundaries. Default is \code{FALSE}.
See technical, mathematical details of this operation under Details.}

\item{no_numeric_policy}{(character) either \code{"warning"},
\code{"plot_blank"}, or \code{"error"}. Defaults to \code{"warning"}. If
\code{modelingResult} does not contain any numeric values,
\code{"plot_blank"} will plot a blank ROC curve, \code{"warning"} will emit
a warning in addition to plotting a blank ROC curve, and \code{"error"}
will throw an error.}

\item{seed}{(numeric) the random seed to use when bootstrapping for ROC confidence
intervals. Passed to \code{\link[base]{set.seed}}.}
}
\description{
Plot AUC (with ggplot2) after LOOCV model evaluation.
}
\details{
If \code{rectilinear = TRUE}, transforms segment \eqn{y = mx + b}
  between \eqn{fpr1} and \eqn{fpr2} to the line \eqn{x = avg(fpr1, fpr2)}
  between \eqn{y = m \times fpr1 + b} and \eqn{y = m \times fpr2 + b}.
  Surrounding horizontal segments are extended to this new vertical segment.
}
\seealso{
\code{\link{plotAUC}}
}
