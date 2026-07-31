#' Calculating basic enrichment statistics.
#'
#' Calculates fold enrichment and Matthews Correlation Coefficient
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

  confusion_matrix <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE,
                             dimnames = list(enriched = c("enriched", "not enriched"),
                                             annotated = c("target", "not target")))

  output <- c("fold enrichment" = fold_enrichment, "MCC" = MCC)
  attr(output, "confusion_matrix") <- confusion_matrix
  class(output) <- "enrichment_stats"
  return(output)
}


#' Print method for enrichment_stats objects
#'
#' Prints the confusion matrix with marginal totals and summary statistics
#' (fold enrichment, MCC) for an `enrichment_stats` object in a readable format.
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

  # Wrap confusion matrix in addmargins() to show row and column totals
  cm_with_margins <- addmargins(attr(x, "confusion_matrix"))
  print(cm_with_margins)

  cat("\n")
  cat("Metrics:\n")
  cat(sprintf("  Fold enrichment: %.3f\n", x["fold enrichment"]))
  cat(sprintf("  MCC:              %.3f\n", x["MCC"]))
  invisible(x)
}





#' Print method for enrichment_stats objects
#'
#' Prints the confusion matrix with marginal totals separated by table lines,
#' along with summary statistics (fold enrichment, MCC).
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

  cm <- attr(x, "confusion_matrix")
  cm_m <- addmargins(cm)

  # Convert matrix counts to formatted string matrix for uniform width alignment
  str_cm <- apply(cm_m, c(1, 2), function(val) format(val, big.mark = ","))

  # Determine width required for column alignment
  col_widths <- apply(rbind(colnames(str_cm), str_cm), 2, function(col) max(nchar(col)))
  col_widths <- pmax(col_widths, 8)  # ensure at least width 8

  row_names <- rownames(str_cm)
  max_row_name_len <- max(nchar(row_names))

  # Helper to format a single table row cleanly
  fmt_row <- function(label, vals) {
    lbl_pad <- sprintf(paste0("%-", max_row_name_len, "s"), label)
    val_pad <- mapply(function(v, w) sprintf(paste0("%", w, "s"), v), vals, col_widths)
    paste0(lbl_pad, " | ", paste(val_pad[1:2], collapse = " "), " | ", val_pad[3])
  }

  # 1. Header row
  lbl_header <- sprintf(paste0("%-", max_row_name_len, "s"), "enriched \\ target")
  col_headers <- mapply(function(v, w) sprintf(paste0("%", w, "s"), v), colnames(str_cm), col_widths)
  header_str <- paste0(lbl_header, " | ", paste(col_headers[1:2], collapse = " "), " | ", col_headers[3])

  # 2. Separator line matching total row character width
  sep_line <- paste0(
    paste(rep("-", max_row_name_len), collapse = ""),
    "-+-",
    paste(rep("-", col_widths[1] + col_widths[2] + 1), collapse = ""),
    "-+-",
    paste(rep("-", col_widths[3]), collapse = "")
  )

  # Print formatted ASCII table
  cat(header_str, "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(fmt_row(row_names[1], str_cm[1, ]), "\n", sep = "")
  cat(fmt_row(row_names[2], str_cm[2, ]), "\n", sep = "")
  cat(sep_line, "\n", sep = "")
  cat(fmt_row(row_names[3], str_cm[3, ]), "\n", sep = "")

  cat("\n")
  cat("Metrics:\n")
  cat(sprintf("  Fold enrichment: %.3f\n", x["fold enrichment"]))
  cat(sprintf("  MCC:              %.3f\n", x["MCC"]))

  invisible(x)
}





