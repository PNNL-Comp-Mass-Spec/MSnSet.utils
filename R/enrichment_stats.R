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
#' @return an object of class `enrichment_stats`: a named numeric vector
#'   with elements `"fold enrichment"` and `"MCC"`, with the underlying
#'   confusion matrix stored as a `"confusion_matrix"` attribute.
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

  confusion_matrix <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
                             dimnames = list(predicted = c("enriched", "not enriched"),
                                             actual = c("target", "not target")))

  output <- c("fold enrichment" = fold_enrichment, "MCC" = MCC)
  attr(output, "confusion_matrix") <- confusion_matrix
  class(output) <- "enrichment_stats"
  return(output)
}


#' Print method for enrichment_stats objects
#'
#' Prints the confusion matrix and summary statistics (fold enrichment,
#' MCC) for an `enrichment_stats` object in a readable format.
#'
#' @param x an object of class `enrichment_stats`
#' @param ... additional arguments (unused)
#'
#' @return `x`, invisibly.
#' @export
print.enrichment_stats <- function(x, ...) {
  cat("Enrichment Statistics\n")
  cat("=====================\n\n")
  cat("Confusion matrix:\n")
  print(attr(x, "confusion_matrix"))
  cat("\n")
  cat("Metrics:\n")
  cat(sprintf("  Fold enrichment: %.3f\n", x["fold enrichment"]))
  cat(sprintf("  MCC:             %.3f\n", x["MCC"]))
  invisible(x)
}


