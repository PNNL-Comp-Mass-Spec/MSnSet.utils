
#' Adjusting PTM profiles by global
#'
#'
#' Adjusting PTM relative abundance profiles by the global profiles of the
#' parent protein. The motivation behind this adjustment is that there are two
#' contributions to the observed PTM abundance changes - change in the
#' PTM stoichiometry and the change in the parent protein abundance.
#'
#'
#' @param m_ptm MSnSet with PTM data
#' @param m_global MSnSet with global data
#' @param method One of "subtraction", "compensation" or "regression"
#'
#' @details
#' \itemize{
#'    \item \bold{Subtraction:} This simple method works best if there is
#'    no dynamic range compression. Otherwise, if global data has higher
#'    compression, the PTM data will be undercorrected. It calculates
#'    the adjusted profile as: \eqn{PTM_{adj} = PTM_{obs} - Protein_{global}}.
#'
#'    \item \bold{Compensation:} Experimental! The algorithm tries to
#'    infer the overall dynamic range compression coefficient from the entire
#'    data by fitting \emph{total} least squares regression. Then this coefficient
#'    is use for the adjustment similar to the \bold{subtraction} algorithm:
#'    \eqn{PTM_{adj} = PTM_{obs} - coeff*Protein_{global}}
#'
#'    \item \bold{Regression:} This method fits a linear model for each feature:
#'    \eqn{PTM_{obs} \sim \beta_0 + \beta_1 \cdot Protein_{global}}. Then
#'    \eqn{PTM_{adj} = PTM_{obs} - \beta_0 + \beta_1 \cdot Protein_{global}}
#'    Works well only if there are lots of samples with varying change
#'    (e.g. cohort of human subjects). It probably won't work for A vs B
#'    comparison. Since this is basically two point and thus there is no
#'    room for regression.
#' }
#'
#' Features that do not have matching parent protein data or have insufficient
#' data points for regression are returned as \code{NA} and subsequently
#' filtered.
#'
#' @return MSnSet object
#' @importFrom Biobase fData featureNames exprs
#' @export adjust_for_parent_protein_profiles
#'
#'


adjust_for_parent_protein_profiles <- function(m_ptm, m_global, method = "subtraction"){

   # sample alignment
   common_samples <- intersect(sampleNames(m_ptm), sampleNames(m_global))
   stopifnot(length(common_samples) > 2)
   m_ptm <- m_ptm[,common_samples]
   m_global <- m_global[,common_samples]

   # gene matching
   stopifnot("Gene" %in% fvarLabels(m_ptm))

   gene_overlap <- intersect(fData(m_ptm)$Gene, featureNames(m_global))
   m_ptm <- m_ptm[fData(m_ptm)$Gene %in% gene_overlap,]
   x_global <- exprs(m_global)[fData(m_ptm)$Gene,]

   x_ptm <- exprs(m_ptm)

   if(method == "subtraction"){
      x_ptm_adj <- x_ptm - x_global
   }else if(method == "compensation"){
      #
      df <- data.frame(global = as.vector(x_global), ptm = as.vector(x_ptm))
      df_clean <- na.omit(df[, c("global", "ptm")])
      pca_mod  <- prcomp(df_clean, scale. = FALSE)
      loadings <- pca_mod$rotation
      slope <- loadings["ptm", 1] / loadings["global", 1]
      # bulletproof
      slope <- abs(slope) * sign(cor(df$global, df$ptm, use = "complete.obs"))
      x_ptm_adj <- x_ptm - slope*x_global
   }else if(method == "regression"){
      pf <- function(i){
         x <- as.numeric(x_global[i,])
         y <- as.numeric(x_ptm[i,])
         idx <- !is.na(x) & !is.na(y)
         if(sum(idx) > 1){
            mdl <- lm(y ~ x)
            ynot <- predict(mdl, newdata = data.frame(x))
            return(y - as.numeric(ynot))
         }else{
            return(rep(NA,ncol(x_ptm)))
         }
      }

      res <- lapply(1:nrow(m_ptm), pf)

      x_ptm_adj <- matrix(unlist(res), byrow=TRUE, nrow=nrow(x_ptm))
   } else{
      stop("unknown method")
   }

   exprs(m_ptm) <- x_ptm_adj

   idx_NAs <- apply(is.na(exprs(m_ptm)), 1, all)

   return(m_ptm[!idx_NAs,])

}


