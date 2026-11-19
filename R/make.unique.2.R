#' Custom Make Unique
#'
#' Similar to \code{\link[base]{make.unique}}, but with some
#' customization. E.g. starts counting with ".1" instead of blank.
#'
#' @param names character vector with non-unique entries
#' @param sep character separator between the name and index
#' @return character vector where all non-unique entries have attached suffix
#'          resolving their non-uniqueness.
#'
#' @export make.unique.2
#'
#' @examples
#' make.unique(c('a','b','b'))
#' make.unique.2(c('a','b','b'))

# make.unique.2 <- function(x){
#     if(!any(duplicated(x)))
#         return(x)
#     xu <- unique(x)
#     x2 <- c(xu, x)
#     x2 <- make.unique(x2)
#     return(x2[-seq_along(xu)])
# }


make.unique.2 <- function(names, sep = ".") {
  if (!any(duplicated(names))) {
    return(names)
  }
  dups <- unique(names[duplicated(names)])
  res <- make.unique(c(dups, names), sep = sep)
  res[-seq_along(dups)]
}