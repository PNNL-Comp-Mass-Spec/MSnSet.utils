#' Calculating basic enrichment statistics.
#'
#' Calculates fold enrichment and Matthews Correlation Coefficient
#'
#'
#' @param enriched list of enriched proteins
#' @param target list of proteins targeted for enrichment
#' @param background list of background proteins
#'
#' @description
#' Type of identification does not matter as long as it is the same -- be
#' that gene symbols or UniProt accessions. Computes fold enrichment and
#' Matthews Correlation Coefficient (MCC), which unlike fold enrichment
#' also accounts for true negatives and is more robust to class imbalance.
#' Absolute values depend on the choice of `background`; comparisons
#' across enrichment lists should use the same `background`. Ideally, the
#' `background` should be the detectable proteome. However, if the detectable
#' proteome (e.g. from global samples) is unavailable, then it is OK to use
#' the entire genome. Use of entire genome will make the numbers biased, since
#' unexpressed proteins will inflate true negatives, but it is OK to use for
#' relative comparison of the enrichment methods.
#'
#' @return numeric
#' @export

compute_enrichment_stats <- function(enriched, target, background){

  TP <- length(intersect(enriched, target))
  FP <- length(enriched) - TP
  FN <- length(target) - TP
  TN <- length(background) - length(target) - FP

  TP <- as.numeric(TP)
  FP <- as.numeric(FP)
  TN <- as.numeric(TN)
  FN <- as.numeric(FN)

  fold_enrichment <- (TP/(TP+FP))/((TP+FN)/(FP+TN+TP+FN))
  MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  output <- c(fold_enrichment, MCC)
  names(output) <- c("fold enrichment", "MCC")
  return(output)
}

