#' Affinity Enrichment Statistics
#'
#' @description Models log2(Enriched/Non-Enriched) ratio as normal distribution.
#' The computes the p-values that the a particular point is part of the
#' distribution.
#'
#' @param m1 MSnSet with enriched data. Replicates will be averaged.
#' @param m2 MSnSet with non-enriched data. Replicates will be averaged.
#' @param prefix appends this string to the column names of the output
#' data.frame
#'
#' @return A data.frame containing:
#' \itemize{
#'   \item \code{diff}: log2(Enriched/Non-Enriched)
#'   \item \code{pval}: Probability of observing the ratio under the null normal distribution
#'   \item \code{padj}: Benjamini-Hochberg adjusted p-value
#' }
#'
#' @importFrom Biobase featureNames
#'
#' @export affinity_enrichment_stats

affinity_enrichment_stats <- function(m1, m2, prefix = NULL){

  stopifnot(all(featureNames(m1) == featureNames(m2)))

  mean_1 <- apply(exprs(m1), 1, mean, na.rm = TRUE)
  mean_2 <- apply(exprs(m2), 1, mean, na.rm = TRUE)
  diff <- mean_1 - mean_2
  pval <- pnorm(abs(diff), median(diff), mad(diff), lower.tail = FALSE)
  padj <- p.adjust(pval, method = "BH")
  df_out <- data.frame(diff, pval, padj)
  if(!is.null(prefix) && prefix != ""){
    colnames(df_out) <- paste(prefix, colnames(df_out), sep = ".")
  }
  rownames(df_out) <- featureNames(m1)
  return(df_out)

}

