#' @title Plot AUC
#'
#' @description Plot AUC after LOOCV model evaluation.
#'
#' @param modelingResult output of \code{lr_modeling} or \code{rf_modeling}.
#' @param CI (logical) whether to plot confidence intervals. Default is
#'   \code{FALSE}.
#' @param ... further arguments passed to the \code{plot} method of the
#'   \code{\link[pROC]{ci.se}} object of the \code{pROC} package.
#'
#' @seealso \code{\link{plotAUC_gg}}
#'
#' @importFrom ROCR performance
#' @importFrom graphics abline
#' @importFrom pROC roc ci.se
#'
#' @export plotAUC

plotAUC <- function(modelingResult,
                    CI = FALSE,
                    ...)
{
  if (!CI) {
    perf <- performance(modelingResult$pred, "tpr", "fpr")
    # x=1-spec, y=sens
    plot(perf,
         main = sprintf("AUC: %s", round(modelingResult$auc, 2)),
         col = 2, lwd = 2)
    abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
  } else {
    old_par <- par()
    par(pty = "s")
    pROC_obj <- roc(modelingResult$pred@labels[[1]],
                    modelingResult$pred@predictions[[1]],
                    direction = "<",
                    smoothed = T,
                    # arguments for ci
                    ci = TRUE,
                    ci.alpha = 0.95,
                    stratified = FALSE,
                    # arguments for plot
                    plot = TRUE,
                    auc.polygon = FALSE,
                    max.auc.polygon = FALSE,
                    grid = FALSE,
                    print.auc = TRUE,
                    print.auc.y = 0.1,
                    print.auc.x = 0.8,
                    show.thres = TRUE)

    sens.ci <- ci.se(pROC_obj)
    plot(sens.ci, type = "shape", col = "#FF8888", conf = 95, ...)
    abline(a = 1, b = -1, lty = 2, lwd = 2, col = "grey50")
    par(old_par)
  }
}


#' @title Plot AUC (with ggplot2)
#'
#' @description Plot AUC (with ggplot2) after LOOCV model evaluation.
#'
#' @inheritParams plotAUC
#' @param rectilinear (logical) whether to prevent diagonal lines being formed from jumps
#'   in TPR from CI boundaries. Default is \code{FALSE}.
#'   See technical, mathematical details of this operation under Details.
#' @param no_numeric_policy (character) either \code{"warning"},
#'   \code{"plot_blank"}, or \code{"error"}. Defaults to \code{"warning"}. If
#'   \code{modelingResult} does not contain any numeric values,
#'   \code{"plot_blank"} will plot a blank ROC curve, \code{"warning"} will emit
#'   a warning in addition to plotting a blank ROC curve, and \code{"error"}
#'   will throw an error.
#' @param seed (numeric) the random seed to use when bootstrapping for ROC confidence
#'   intervals. Passed to \code{\link[base]{set.seed}}.
#'
#' @details If \code{rectilinear = TRUE}, transforms segment \eqn{y = mx + b}
#'   between \eqn{fpr1} and \eqn{fpr2} to the line \eqn{x = avg(fpr1, fpr2)}
#'   between \eqn{y = m \times fpr1 + b} and \eqn{y = m \times fpr2 + b}.
#'   Surrounding horizontal segments are extended to this new vertical segment.
#'
#' @noMd
#'
#' @seealso \code{\link{plotAUC}}
#'
#' @importFrom ROCR performance
#' @importFrom pROC roc ci.se
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot geom_abline scale_linetype_manual labs theme
#'   scale_x_continuous scale_y_continuous geom_ribbon scale_fill_manual
#'   geom_line element_text
#' @importFrom tidyr fill
#' @importFrom dplyr %>% mutate select arrange relocate
#'
#' @export plotAUC_gg

plotAUC_gg <- function(modelingResult,
                       CI = FALSE,
                       rectilinear = FALSE,
                       no_numeric_policy = c("warning", "plot_blank", "error"),
                       seed = 0)
{
  no_numeric_policy <- match.arg(no_numeric_policy,
                                 choices = c("warning", "plot_blank", "error"))

  fpr <- lo <- hi <- NULL

  if (is.na(modelingResult$auc)) {
    switch(no_numeric_policy,
           warning = {
             warning("modelingResult$auc is NA. Plotting blank ROC Curve.")
             # plot_blank
           },
           plot_blank = {
             # plot_blank
           },
           error = {
             stop("modelingResult$auc is NA and `modelingResult` = \"error\".")
           })

    # plot_blank
    p <- ggplot() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "#888888", show.legend = FALSE) +
      labs(subtitle = "AUC: NA")
  } else {
    p <- plot_main(modelingResult = modelingResult,
                   seed = seed,
                   conf_int = CI,
                   rectilinear = rectilinear)
  }

  # Modifications common to all plots
  p <- p +
    scale_x_continuous(name = "False Positive Rate",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(name = "True Positive Rate",
                       limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggtitle("ROC Curve") +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5))

  return(p)
}


## plotAUC_gg helper functions ----
plot_main <- function(modelingResult,
                      seed = 0,
                      conf_int = FALSE,
                      rectilinear = FALSE) {
  perf <- performance(modelingResult$pred, "tpr", "fpr")

  set.seed(seed)

  pROC_obj <- suppressMessages({
    roc(modelingResult$pred@labels[[1]],
        modelingResult$pred@predictions[[1]],
        direction = "<",
        smoothed = TRUE,
        # arguments for ci
        ci = conf_int,
        ci.alpha = 0.95,
        stratified = FALSE,
        # arguments for plot
        plot = FALSE,
        auc.polygon = FALSE,
        max.auc.polygon = FALSE,
        grid = FALSE)
  })

  df <- data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]])

  if (conf_int) {
    the_ci <- ci.se(pROC_obj,
                    specificities = seq(0, 1, 0.025),
                    progress = "none") %>%
      as.data.frame() %>%
      `colnames<-`(c("lo", "mid", "hi")) %>%
      rownames_to_column("fpr") %>%
      mutate(fpr = 1 - as.double(sub("X", "", fpr))) %>%
      select(-c("mid")) %>%
      arrange(fpr)

    if (rectilinear) {
      the_ci <- make_rectilinear_roc_ci(the_ci)
    }
  }

  auc.ci.lo <- pROC_obj$ci[1]
  auc.ci.hi <- pROC_obj$ci[3]

  ribbon_col <- "#c96f6f"
  ci_col <- "red"
  diag_df <- data.frame(x = c(0, 1), y = c(0, 1))
  p <- ggplot(data = df) +
    geom_line(data = diag_df,
              mapping = aes(x = x, y = y, linetype = "Random model"),
              color = "#888888")

  if (conf_int) {
    p <- p +
      geom_ribbon(data = the_ci,
                  mapping = aes(x = fpr, ymin = lo, ymax = hi, fill = ""),
                  color = ci_col, alpha = 0.5) +
      scale_fill_manual(name = "95% bootstrap C.I.", values = ribbon_col) +
      labs(subtitle = sprintf("AUC: %.3f\nAUC C.I.: [%.3f, %.3f]",
                              modelingResult$auc, auc.ci.lo, auc.ci.hi))
  } else {
    p <- p +
      labs(subtitle = sprintf("AUC: %.3f", modelingResult$auc))
  }

  p <- p +
    geom_line(mapping = aes(x = x, y = y, linetype = "Model curve (LOOCV)"), lwd = 1) +
    scale_linetype_manual(
      name = NULL,
      values = c("Model curve (LOOCV)" = "solid", "Random model" = "dashed")
    )

  return(p)
}


make_rectilinear_roc_ci <- function(ci_df) {
  intermediate_points_lo <- intermediate_points_hi <- list()
  prev_y_lo <- prev_y_hi <- 0

  for (i in 2:nrow(ci_df)) {
    mid_val <- (ci_df$fpr[i - 1] + ci_df$fpr[i]) / 2

    if (ci_df$lo[i] != prev_y_lo) {
      intermediate_points_lo <- append(
        intermediate_points_lo,
        list(
          c(mid_val, ci_df$lo[i - 1]),
          c(mid_val, ci_df$lo[i])
        )
      )
    }

    if (ci_df$hi[i] != prev_y_hi) {
      intermediate_points_hi <- append(
        intermediate_points_hi,
        list(
          c(mid_val, ci_df$hi[i - 1]),
          c(mid_val, ci_df$hi[i])
        )
      )
    }

    prev_y_lo <- ci_df$lo[i]
    prev_y_hi <- ci_df$hi[i]
  }

  lo_df <- do.call(cbind, intermediate_points_lo) %>%
    t() %>%
    as.data.frame() %>%
    mutate(hi = NA) %>%
    `colnames<-`(c("fpr", "lo", "hi"))

  hi_df <- do.call(cbind, intermediate_points_hi) %>%
    t() %>%
    as.data.frame() %>%
    mutate(lo = NA) %>%
    `colnames<-`(c("fpr", "hi", "lo")) %>%
    relocate(lo, .before = hi)

  ci_df <- rbind(ci_df, lo_df, hi_df) %>%
    arrange(fpr) %>%
    fill(lo, hi, .direction = "down")

  return(ci_df)
}


utils::globalVariables(
  c("fpr", "hi", "lo", "modelingResult", "rectilinear")
)
