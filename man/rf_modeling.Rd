% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_modeling.R
\name{rf_modeling}
\alias{rf_modeling}
\title{Random Forest Predictive Models}
\usage{
rf_modeling(
  msnset,
  features,
  response,
  pred.cls,
  K = NULL,
  sel.feat = TRUE,
  sel.alg = c("varSelRF", "Boruta", "top"),
  cores = NULL,
  seed = 0,
  ...
)
}
\arguments{
\item{msnset}{MSnSet object}

\item{features}{character vector features to select from for building
prediction model. The features can be either in featureNames(msnset) or in
pData(msnset).}

\item{response}{factor to classify along. Must be only 2 levels.}

\item{pred.cls}{class to predict}

\item{K}{specifies the cross-validation type. Default NULL means LOOCV.
Another typical value is 10.}

\item{sel.feat}{logical defining if to select features or use the entire set?}

\item{sel.alg}{character.
\itemize{
    \item varSelRF
    \item Boruta
    \item top - just selecting top 3 features
}}

\item{cores}{integer; the number of processes to use for Boruta.
If \code{NULL}, one less than the maximum number of processes will be used.
If this number is greater than 125, a maximum of 125 cores will be used.}

\item{seed}{numeric; the random seed to use during various stages of 
\code{rf_modeling}. Passed to \code{\link[base]{set.seed}}.}

\item{...}{Extra arguments. Currently passed only to Boruta algorithm.}
}
\value{
list
     \describe{
         \item{\code{prob}}{is the probabilities (response) from LOOCV
         that the sample is "case"}
         \item{\code{features}}{list of selected features
         for each iteration of LOOCV}
         \item{\code{top}}{top features over all iterations}
         \item{\code{auc}}{AUC}
         \item{\code{pred}}{prediction perfomance obtained by
                 \code{ROCR::prediction}}
     }
}
\description{
The objective is - given the number of features, select the most
  informative ones and evaluate the predictive random forest model. The
  feature and model selection performed independently for each round of
  LOOCV.
}
\examples{
# Not run
\dontrun{
data(srm_msnset)
head(varLabels(msnset))
head(msnset$subject.type)
# reduce to two classes
msnset <- msnset[,msnset$subject.type != "control.1"]
msnset$subject.type <- as.factor(msnset$subject.type)
plotAUC(out)
# top features consistently recurring in the models during LOOCV
print(out$top)
# the AUC
print(out$auc)
# probabilities of classifying the sample right, if the feature selection
# and model training was performed on other samples
plot(sort(out$prob))
abline(h=0.5, col='red')
}
}
