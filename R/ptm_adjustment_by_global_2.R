
#' Adjusting PTM profiles by global. Version 2
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
#'    \item \bold{subtraction:} This simple method works best if there is
#'    no dynamic range compression. Otherwise, if global data has higher
#'    compression, the PTM data will be undercorrected. It calculates
#'    the adjusted profile as: \eqn{PTM_{adj} = PTM_{obs} - Protein_{global}}.
#'
#'    \item \bold{regression:} This method fits a linear model for each feature:
#'    \eqn{PTM_{obs} \sim \beta_0 + \beta_1 \cdot Protein_{global}}. Then
#'    \eqn{PTM_{adj} = PTM_{obs} - \beta_0 + \beta_1 \cdot Protein_{global}}
#'    Works well only if there are lots of samples with varying change
#'    (e.g. cohort of human subjects). It probably won't work for A vs B
#'    comparison. Since this is basically two points, there is no
#'    room for regression.
#'
#'    \item \bold{moderation:} DO NOT USE! The method is similar to `regression`.
#'    However, the \eqn{\beta_1} coefficients are modeled as a mixture of two
#'    normal distributions using the Expectation-Maximization approach.
#'    One, with smaller standard deviation is assumed
#'    composed of PTMs with non-changing stoichiometry and the \eqn{\beta_1}
#'    coefficient, probably less than 1, reflects the dynamic range compression.
#'    The other distribution, with larger standard deviation reflects PTMs with
#'    truly changing stoichiometry. Further, the probabilities of belonging to
#'    non-changing and changing distributions are derived. If \eqn{p} is the
#'    probability that the site is changing, then the coefficient
#'    is \eqn{\hat{\beta}_1 = p \cdot \mu_{nonchanging} + (1-p) \cdot \beta_1}.
#'    \eqn{PTM_{adj} = PTM_{obs} - \beta_0 + \hat{\beta}_1 \cdot Protein_{global}}
#'
#'    \item \bold{compensation:} This algorithm can be viewed as a continuation:
#'    `regression` -> `moderation` -> `compensation`. Essentially this is an
#'    extreme moderation. That is, all the \eqn{\beta_1} are shrunk on just one value -
#'    the mean of the narrower (smaller standard deviation) distribution (modeled via the global median).
#'    \eqn{PTM_{adj} = PTM_{obs} - \mu_{nonchanging} \cdot Protein_{global}}
#' }
#'
#' Features that do not have matching parent protein data or have insufficient
#' data points for regression are returned as \code{NA} and subsequently
#' filtered.
#'
#' @return MSnSet object
#' @importFrom Biobase fData featureNames exprs sampleNames fvarLabels
#' @importFrom ggplot2 ggplot geom_histogram geom_line scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual labs theme_minimal theme
#' @importFrom ggplot2 after_stat aes
#' @export

adjust_for_parent_protein_profiles_2 <- function(m_ptm, m_global, method = "subtraction"){

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

   # fit the LM models
   betas_unit <- rep(1, nrow(x_ptm))
   betas_reg <- calc_rowwise_betas(x_ptm, x_global)
   betas_reg_hat <- moderate_betas(betas_reg)
   median_beta <- median(betas_reg, na.rm = TRUE)
   betas_mean <- rep(median_beta, nrow(x_ptm))

   # select the according betas
   betas <- switch(method,
                   subtraction = betas_unit,
                   regression = betas_reg,
                   moderation = betas_reg_hat,
                   compensation = betas_mean)

   # adjust
   x_stoich <- calc_stoich_residuals(x_ptm, x_global, betas)
   exprs(m_ptm) <- x_stoich

   if(method == "compensation"){
     # derive p-values
     pvals <- calc_rowwise_pvals_for_beta_deviation(x_ptm, x_global, median_beta)
     beta_deviation <- betas_reg - median_beta
     #
     fData(m_ptm)$PTM_beta_pval <- pvals
     fData(m_ptm)$PTM_beta_dev <- beta_deviation
   }

   idx_NAs <- apply(is.na(exprs(m_ptm)), 1, all)
   m_ptm <- m_ptm[!idx_NAs,]

   return(m_ptm)

}



calc_rowwise_betas <- function(x_ptm, x_glob) {
  betas <- sapply(1:nrow(x_ptm), function(i) {
    y <- x_ptm[i, ]
    x <- x_glob[i, ]

    # 1. Identify which samples have valid, real numbers in BOTH variables
    keep <- is.finite(y) & is.finite(x)

    # 2. Guardrail: Make sure we have enough data points to fit a line
    if (sum(keep) < 3) return(NA)

    # 3. Filter out the NA samples dynamically for this specific iteration
    y_clean <- y[keep]
    X_clean <- cbind(Intercept = 1, x_glob = x[keep])

    # 4. Fit the cleaned matrix
    fit <- lm.fit(X_clean, y_clean)

    return(fit$coefficients["x_glob"])
  })

  names(betas) <- rownames(x_ptm)
  return(betas)
}



# calc_rowwise_betas <- function(x_ptm, x_glob) {
#   # Loop through each row index
#   betas <- sapply(1:nrow(x_ptm), function(i) {
#     y <- x_ptm[i, ]
#
#     # Create design matrix: column of 1s (intercept) and column of x (slope)
#     X <- cbind(Intercept = 1, x_glob = x_glob[i, ])
#
#     # lm.fit is the high-performance, matrix-backed engine inside lm()
#     fit <- lm.fit(X, y)
#
#     # Extract the slope coefficient (second element)
#     return(fit$coefficients["x_glob"])
#   })
#
#   # Name the output vector using the original row names
#   names(betas) <- rownames(x_ptm)
#   return(betas)
# }

# -------------------------------------------------------------------------
# Example Usage:
# -------------------------------------------------------------------------
# set.seed(42)
# mat_ptm  <- matrix(rnorm(50, 10, 2), nrow = 5)
# mat_glob <- matrix(rnorm(50, 12, 3), nrow = 5)
#
# row_betas <- calc_rowwise_betas(mat_ptm, mat_glob)
# print(row_betas)





# calc_rowwise_pvals_for_beta_deviation <- function(x_ptm, x_glob, beta_null = 0) {
#   pvals <- sapply(1:nrow(x_ptm), function(i) {
#     y <- x_ptm[i, ]
#     x <- x_glob[i, ]
#
#     keep <- is.finite(y) & is.finite(x)
#     if (sum(keep) < 3) return(NA)
#
#     y_clean <- y[keep]
#     X_clean <- cbind(Intercept = 1, x_glob = x[keep])
#
#     fit <- lm.fit(X_clean, y_clean)
#
#     rdf <- length(y_clean) - 2
#     resvar <- sum(fit$residuals^2) / rdf
#     se <- sqrt(diag(chol2inv(fit$qr$qr[1:2, 1:2])) * resvar)
#
#     # ADJUSTMENT MADE HERE: Subtract beta_null from the slope coefficient
#     slope_coef <- fit$coefficients[2]
#     t_stat <- (slope_coef - beta_null) / se[2]
#
#     return(2 * pt(abs(t_stat), df = rdf, lower.tail = FALSE))
#   })
#
#   names(pvals) <- rownames(x_ptm)
#   return(pvals)
# }



calc_rowwise_pvals_for_beta_deviation <- function(x_ptm, x_glob, beta_null = 0) {
  pvals <- sapply(1:nrow(x_ptm), function(i) {
    y <- x_ptm[i, ]
    x <- x_glob[i, ]

    # 1. Clean missing/infinite data
    keep <- is.finite(y) & is.finite(x)
    if (sum(keep) < 3) return(NA)

    y_clean <- y[keep]
    x_clean <- x[keep]

    # 2. Guardrail: If x has zero variance (all values identical), regression fails
    if (var(x_clean) == 0) return(NA)

    X_clean <- cbind(Intercept = 1, x_glob = x_clean)
    fit <- lm.fit(X_clean, y_clean)

    # 3. Handle degrees of freedom dynamically
    rdf <- fit$df.residual
    if (is.null(rdf) || rdf < 1) return(NA)

    # 4. Extract standard errors safely using dynamic rank
    # (Fixes the hardcoded [1:2, 1:2] index crash if a column drops)
    p <- fit$rank
    if (p < 2) return(NA) # Means the slope couldn't be estimated

    resvar <- sum(fit$residuals^2) / rdf
    X_inv  <- chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE])
    se     <- sqrt(diag(X_inv) * resvar)

    # 5. Calculate deviation from beta_null
    slope_coef <- fit$coefficients[2]

    # Guardrail: If SE is zero or NA, t-stat is undefined
    if (is.na(se[2]) || se[2] == 0) return(NA)

    t_stat <- (slope_coef - beta_null) / se[2]

    return(2 * pt(abs(t_stat), df = rdf, lower.tail = FALSE))
  })

  names(pvals) <- rownames(x_ptm)
  return(pvals)
}




moderate_betas <- function(x, visualize = F){

  # ==========================================
  # Expectation-Maximization (EM) Algorithm
  # ==========================================

  # Initial Guesses (Slightly separated to help convergence)
  est_mu1   <- quantile(x, 0.25)
  est_sd1   <- sd(x) / 2
  est_mu2   <- quantile(x, 0.75)
  est_sd2   <- sd(x) / 2
  est_pi    <- 0.5 # Initial mixing weight for Component 1

  max_iter  <- 200
  tol       <- 1e-6
  log_liks  <- c()

  for (iter in 1:max_iter) {

    # ---------- E-Step ----------
    # Calculate the likelihood of each point for both Normal components
    f_norm1 <- dnorm(x, mean = est_mu1, sd = est_sd1)
    f_norm2 <- dnorm(x, mean = est_mu2, sd = est_sd2)

    # Marginal density (denominator)
    marginal <- est_pi * f_norm1 + (1 - est_pi) * f_norm2
    marginal[marginal == 0] <- .Machine$double.xmin # error-proofing

    # Posterior probabilities (Responsibilities)
    gamma_1 <- (est_pi * f_norm1) / marginal
    gamma_2 <- 1 - gamma_1 # Because gamma_1 + gamma_2 = 1

    # Track Log-Likelihood to check for convergence
    log_lik <- sum(log(marginal))
    log_liks <- c(log_liks, log_lik)

    if (iter > 1 && abs(log_lik - log_liks[iter - 1]) < tol) {
      # cat("Converged in", iter, "iterations.\n")
      break
    }

    # ---------- M-Step ----------
    # Update mixing weight
    est_pi <- mean(gamma_1)

    # Update Component 1 parameters
    est_mu1 <- sum(gamma_1 * x) / sum(gamma_1)
    est_sd1 <- sqrt(sum(gamma_1 * (x - est_mu1)^2) / sum(gamma_1))

    # Update Component 2 parameters
    est_mu2 <- sum(gamma_2 * x) / sum(gamma_2)
    est_sd2 <- sqrt(sum(gamma_2 * (x - est_mu2)^2) / sum(gamma_2))
  }


  if(visualize) {

    # ==========================================
    # Visualization
    # ==========================================

    # Create a data frame for plotting the actual data
    df_data <- data.frame(x = x)

    # Create a grid for the fitted curve lines
    x_grid <- seq(min(x), max(x), length.out = 500)
    fit_norm1 <- est_pi * dnorm(x_grid, mean = est_mu1, sd = est_sd1)
    fit_norm2 <- (1 - est_pi) * dnorm(x_grid, mean = est_mu2, sd = est_sd2)
    fit_total <- fit_norm1 + fit_norm2

    df_curves <- data.frame(
      x = rep(x_grid, 3),
      Density = c(fit_norm1, fit_norm2, fit_total),
      Component = factor(rep(c("Fitted Normal 1",
                               "Fitted Normal 2",
                               "Total Mixture"),
                             each = length(x_grid)))
    )

    # Plot using ggplot2
    p <- ggplot() +
      # Histogram of raw data
      geom_histogram(data = df_data, aes(x = x, y = after_stat(density)),
                     bins = 80, fill = "gray85", color = "gray60") +
      # Fitted mixture component curves
      geom_line(data = df_curves, aes(x = x,
                                      y = Density,
                                      color = Component,
                                      linetype = Component),
                size = 1) +
      scale_color_manual(values = c("Fitted Normal 1" = "royalblue",
                                    "Fitted Normal 2" = "darkorange",
                                    "Total Mixture" = "red")) +
      scale_linetype_manual(values = c("Fitted Normal 1" = "dashed",
                                       "Fitted Normal 2" = "dashed",
                                       "Total Mixture" = "solid")) +
      labs(title = "EM Algorithm Fit: Two Normal Distributions Mixture",
           x = "Value",
           y = "Density",
           color = "Model Components",
           linetype = "Model Components") +
      theme_minimal() +
      theme(legend.position = "top")
    plot(p)
  }


  # OK. I'll assume the distribution with larger SD is the changing one.
  # distro with smaller SD is the non-regulated passenger one.
  nonchange_idx <- which.min(c(est_sd1, est_sd2))
  change_idx <- which.max(c(est_sd1, est_sd2))
  prob_changing <- list(gamma_1, gamma_2)[[change_idx]]
  mu_nonchanging <- list(est_mu1, est_mu2)[[nonchange_idx]]

  # if it is non-changing, use the mean coeff for correction
  # if it is changing,

  # if it non-changing, then regress out the parent protein
  # if it is changing, then subtract

  beta_hat <- prob_changing * mu_nonchanging + (1-prob_changing) * x

  return(beta_hat)

}





calc_stoich_residuals <- function(x_ptm, x_glob, betas) {
  # -------------------------------------------------------------------------
  # Row-wise Residual Calculation
  # -------------------------------------------------------------------------
  x_stoich <- t(sapply(1:nrow(x_ptm), function(i) {
    y <- x_ptm[i, ]
    x <- x_glob[i, ]
    b <- betas[i]

    # Remove the fixed slope effect
    y_adjusted <- y - (b * x)

    # The optimal intercept is simply the mean of the adjusted values
    intercept <- mean(y_adjusted, na.rm = TRUE)

    # Residuals = actual - predicted
    res <- y_adjusted - intercept
    return(res)
  }))

  # Keep the original row and column names if they exist
  rownames(x_stoich) <- rownames(x_ptm)
  colnames(x_stoich) <- colnames(x_ptm)

  return(x_stoich)
}

# -------------------------------------------------------------------------
# Example Usage:
# -------------------------------------------------------------------------
# set.seed(42)
# mat_ptm  <- matrix(rnorm(50, 10, 2), nrow = 5)
# mat_glob <- matrix(rnorm(50, 12, 3), nrow = 5)
# vec_b    <- c(0.5, 1.2, -0.3, 0.8, 1.5)
#
# result_matrix <- calc_stoich_residuals(mat_ptm, mat_glob, vec_b)
