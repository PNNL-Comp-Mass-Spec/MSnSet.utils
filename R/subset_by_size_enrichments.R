#' Filter of ontology enrichment results by size
#'
#' This utility works to filter results of enrich* test
#' from \code{\link[clusterProfiler]{clusterProfiler-package}},
#' \code{\link[enrichit]{enrichit-package}} or
#' \code{\link[ReactomePA]{ReactomePA-package}}
#'
#' @param x is an \code{enrichResults} instance
#' @param minObsSize minimum number of genes must be observed per set
#' @param maxObsSize maximum number of genes must be observed per set
#' @param minSetSize minimum number of genes must be in the set
#' @param maxSetSize maximum number of genes must be in the set
#'
#' @return filtered \code{enrichResults} instance
#'
#' @importFrom qvalue qvalue
#' @importFrom enrichit geneInCategory
#' @importFrom stats p.adjust
#'
#' @export subset_by_size
#'
#' @examples
#' library(clusterProfiler)
#' library(enrichplot)
#' library(org.Hs.eg.db)
#' data(gcSample)
#' x <- enrichGO(gcSample[[6]],
#'               universe = unique(unlist((gcSample))),
#'               OrgDb = "org.Hs.eg.db",
#'               ont = "CC",
#'               pvalueCutoff = 1,
#'               qvalueCutoff = 1,
#'               minGSSize = 0)
#' res <- x@result
#' res$geneID <- NULL
#' rownames(res) <- NULL
#' head(res)
#' cnetplot(x)
#'
#' # E.g. the "cell", "cell part" & "cellular_component" seem too generic
#' x <- subset_by_size(x, maxObsSize=50)
#' res <- x@result
#' res$geneID <- NULL
#' rownames(res) <- NULL
#' # Now they are gone
#' head(res)
#' cnetplot(x)

subset_by_size <- function(x, minObsSize=0, maxObsSize=Inf,
                           minSetSize=0, maxSetSize=Inf){
    #
    szo <- sapply(geneInCategory(x), length)
    idxo <- szo >= minObsSize & szo <= maxObsSize
    szs <- sapply(x@geneSets, length)
    idxs <- szs >= minSetSize & szs <= maxSetSize
    selnames <- intersect(names(geneInCategory(x)[idxo]),
                          names(x@geneSets[idxs]))

    # the reason subsetting is done in this cumbersome way
    # is that I don't now if the order matters
    # x@geneInCategory <- x@geneInCategory[names(x@geneInCategory) %in% selnames]
    x@geneSets <- x@geneSets[names(x@geneSets) %in% selnames]
    x@result <- x@result[rownames(x@result) %in% selnames,]

    # stats re-adjustment
    x@result$p.adjust <-  p.adjust(x@result$pvalue, x@pAdjustMethod)
    qobj <- tryCatch(qvalue(p = Over$pvalue,
                            lambda = 0.05,
                            pi0.method = "bootstrap"),
                     error = function(e) NULL)
    if (is(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    }
    else {
        qvalues <- NA
    }
    x@result$qvalue <- qvalues

    #
    return(x)
}

utils::globalVariables("Over")


