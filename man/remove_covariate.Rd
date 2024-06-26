% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_covariate.R
\name{remove_covariate}
\alias{remove_covariate}
\alias{correct_batch_effect}
\alias{remove_batch_effect}
\alias{correct_batch_effect_empiricalBayesLM}
\alias{correct_batch_effect_NA}
\title{Removing Covariate Effect form Expression Data}
\usage{
remove_covariate(x, cov_name)

correct_batch_effect(
  m,
  batch_name,
  least_count_threshold = 2,
  BPPARAM = bpparam(),
  ...
)

remove_batch_effect(x, batch_name, ref_level = NULL, subset_by = c(NULL, NULL))

correct_batch_effect_empiricalBayesLM(
  x,
  removed_cov_name,
  retained_cov_name = NULL,
  least_proportion_threshold = 0.5,
  ...
)

correct_batch_effect_NA(
  m,
  batch_name,
  cov_name = NULL,
  par.prior = FALSE,
  prior.plots = FALSE,
  ...
)
}
\arguments{
\item{x}{MSnSet or ExpressionSet object}

\item{cov_name}{covariate name. Must be in pData(x). At this point it can be
only one name.}

\item{m}{an object of class Eset or MSnSet.}

\item{batch_name}{same thing as covariate name. Using "batch" instead of
"covariate" to keep it consistent with `ComBat`. Must be in pData(x).
At this point it can be only one name.}

\item{least_count_threshold}{minimum number of feature observations
required per batch. The default values is 2, the minimum `ComBat` can
handle safely.}

\item{BPPARAM}{BiocParallelParam for parallel operation.
Default is bpparam(). Use bpparam("SerialParam") if you want to restrict
it to only one thread.}

\item{...}{other arguments for \code{\link[sva]{ComBat}}}

\item{ref_level}{In case a certain factor level should be reference
and kept at zero bias. Default is NULL, i.e. none.}

\item{subset_by}{vector of two strings from varLabels(x). First is the
variable name for subsetting the data. Second is the variable value
to retain.}

\item{removed_cov_name}{covariate name to be removed. Must be in pData(x).}

\item{retained_cov_name}{covariate name to be included in the model but
not removed. Must be in pData(x).}

\item{least_proportion_threshold}{minimum proportion of feature observations
required for the whole dataset. The default value is 50\%.}

\item{par.prior}{argument from \code{\link[sva]{ComBat}}. By default we recommend non-parametric.}

\item{prior.plots}{argument from \code{\link[sva]{ComBat}}. No plots by default.}
}
\description{
The main purpose of this function is to remove batch effect from the data.
Batch can be associated with different days of sample processing (as factor)
or with run order (continuous). Can also be used to remove any unwanted
effects from the data.
}
\section{Functions}{
\itemize{
\item \code{correct_batch_effect()}: wrapper around sva::ComBat

\item \code{remove_batch_effect()}: A flexible batch correction function

\item \code{correct_batch_effect_empiricalBayesLM()}: An empirical Bayesian approach to batch correction
with discrete or continuous covariates

\item \code{correct_batch_effect_NA()}: wrapper around ComBat.NA

}}
\note{
The algorithm essentially uses an LM. The reason for re-inventing the
wheel is presense of missing values in proteomics datasets more then usual.
}
\examples{
# example 1
set.seed(1)
means <- rep(c(1,2,3), each=3)
nrows <- 5
e <- matrix(rep(means, nrows), ncol=length(means), byrow=TRUE)
e <- e +
    matrix(rnorm(Reduce(`*`, dim(e)), sd=0.3), ncol=length(means), byrow=TRUE)

# add missing values in increasing frequency
extreme <- 10 # controls how quickly increases propotion of NAs
# 1 - means it will reach 100\% by the end
# 2 - means only 50\% will be missing by the last row
# N - is 1/N-th
freqs <- (1:nrow(e)-1)/(extreme*(nrow(e)-1))
mis <- t(sapply(freqs, rbinom, n=ncol(e), size = 1))

mis[mis == 1] <- NA
e[5,8:9] <- NA
e[4,4:6] <- NA
e <- e + mis
image(e)
library("ggplot2"); library("reshape2")
ggplot(melt(e), aes(x=Var1, y=Var2, fill=value)) + geom_raster()

# generating factors
facs <- gl(length(unique(means)),length(means)/length(unique(means)))
# alternative is correction for continuous variable
cova <- seq_along(means)

library("Biobase")
m <- ExpressionSet(e)
pData(m)$pesky <- facs
pData(m)$runorder <- cova
m2 <- remove_covariate(m, "pesky")
m3 <- remove_covariate(m, "runorder")

image(exprs(m))
image(exprs(m2))
image(exprs(m3))

# Example 2 (real-world)
data(cptac_oca)
# let's test for iTRAQ_Batch effect
res <- eset_lm(oca.set, "~ iTRAQ_Batch", "~ 1")
# not too strong, but there
hist(res$p.value, 50)
image_msnset(oca.set, facetBy="iTRAQ_Batch")
oca.fixed <- remove_covariate(oca.set, "iTRAQ_Batch")
res <- eset_lm(oca.fixed, "~ iTRAQ_Batch", "~ 1")
hist(res$p.value, 50)
image_msnset(oca.fixed, facetBy="iTRAQ_Batch")

# Example for correct_batch_effect
# Not run
\dontrun{
data("cptac_oca") # oca.set object
plot_pca(oca.set, phenotype = 'Batch')
oca.set.2 <- correct_batch_effect(oca.set, batch_name = "Batch")
plot_pca(oca.set.2, phenotype = 'Batch')
}

data("cptac_oca") # oca.set object
plot_pca(oca.set, phenotype = 'Batch')
oca.set.2 <- remove_batch_effect(oca.set,
                                 batch_name = "Batch", ref_level="X14",
                                 subset_by=c("tumor_stage","IIIC"))
plot_pca(oca.set.2, phenotype = 'Batch')

# Example for correct_batch_effect_empiricalBayesLM
data("cptac_oca") # oca.set object
plot_pca(oca.set, phenotype = 'Batch')

# Not run
\dontrun{
oca.set.2 <- correct_batch_effect_empiricalBayesLM(oca.set,
                                 removed_cov_name = "Batch",
                                 retained_cov_name = "tumor_stage")
plot_pca(oca.set.2, phenotype = 'Batch')
}

# Example for correct_batch_effect
# Not run
\dontrun{
data("cptac_oca") # oca.set object
plot_pca(oca.set, phenotype = 'Batch')
oca.set.2 <- correct_batch_effect(oca.set, batch_name = "Batch")
plot_pca(oca.set.2, phenotype = 'Batch')
}
}
\seealso{
\code{\link[sva]{ComBat}} \code{\link[WGCNA]{empiricalBayesLM}}
}
