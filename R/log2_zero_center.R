
#' Log2 transform and center around zero
#'
#' Log2 transform applied to expression crosstab of ExpressionSet or MSnSet
#' object followed by centering the median around zero. Converts one
#' MSnSet to another MSnSet.
#'
#' @param m ExpresionSet or MSnSet object
#'
#' @importFrom Biobase exprs<- exprs
#'
#' @return (MSnSet) MSnSet object
#'
#' @export log2_zero_center
#'
#' @examples
#' data(srm_msnset)
#' msnset2 <- log2_zero_center(msnset)


log2_zero_center <- function(m, record_med = FALSE){

  exprs(m) <- log2(exprs(m))

  if (any(is.infinite(exprs(m))) == TRUE) {
    stop("After transformation, infinite values are present.")
  }

  med <- apply(exprs(m), 1, median, na.rm = TRUE)

  if (record_med) {
     fData(m) <- fData(m) %>%
        rownames_to_column() %>%
        mutate(feat_median = med) %>%
        column_to_rownames()
  }

  exprs(m) <- sweep(exprs(m),
                    MARGIN = 1,
                    STATS = med,
                    FUN = "-")

  return(m)
}
