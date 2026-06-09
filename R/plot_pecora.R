#' @title Plot PeCorA Results
#'
#' @description Create boxplots (categorical \code{treatment}) or a scatterplot
#'   (numeric \code{treatment}) of PeCorA results for a given peptide.
#'
#' @param m \code{MSnSet} object; output of \code{\link[MSnSet.utils]{PeCorA}}.
#' @param peptide character; the peptide of interest to plot. One of
#'   \code{fData(m)[[peptide_column]]}.
#' @param protein character; the protein that contains \code{peptide} as one of
#'   its peptides. One of \code{fData(m)[[protein_column]]}.
#' @param treatment character; pattern used to search for columns of the form
#'   "pecora_adj_pval_{treatment}" in \code{fData(m)}. If \code{NULL} (default),
#'   this will be automatically determined if there is only one column that
#'   starts with "pecora_adj_pval_".
#' @inheritParams PeCorA
#'
#' @details This function was inspired by the boxplots presented in the PeCorA
#'   publication (\url{https://doi.org/10.1021/acs.jproteome.0c00602}).
#'
#' @seealso \code{\link[MSnSet.utils]{PeCorA}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972-1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom rlang !! sym
#' @importFrom dplyr %>% filter pull select mutate rename group_by distinct
#'   summarise
#' @importFrom tidyr pivot_longer
#' @importFrom stats lm update coef
#'
#' @export plot_PeCorA

plot_PeCorA <- function(m,
                        peptide,
                        protein,
                        treatment = NULL,
                        protein_column = "Protein",
                        peptide_column = "Peptide",
                        median_mod = TRUE)
{
  f_data <- fData(m)

  ## Must have "Protein" and "Peptide" columns
  for (col_i in c(protein_column, peptide_column)) {
    if (!col_i %in% colnames(f_data)) {
      stop(sprintf("fData missing '%s' column.", col_i))
    }
  }

  ## If empty treatment string, look for results in f_data.
  if (is.null(treatment)) {
    index <- grep("pecora_adj_pval", colnames(f_data))

    if (length(index) > 1) {
      stop(paste0("Found multiple treatment variables in the PeCorA results.",
                  "Must specify treatment variable."))
    } else if (length(index) == 0) {
      stop("Did not find any PeCorA results")
    }

    treatment <- sub("^pecora_adj_pval_", "", colnames(f_data)[[index]])
  }

  result_col <- paste0("pecora_adj_pval_", treatment)

  f_data <- f_data %>%
    filter(!!sym(protein_column) == protein) %>%
    mutate(feature_name = rownames(.))

  padj <- try(
    f_data %>%
      filter(!!sym(peptide_column) == peptide) %>%
      pull(result_col) %>%
      signif(3)
  )

  if (inherits(padj, "try-error")) {
    stop(paste("Did not find PeCorA results for", treatment))
  }

  if (is.na(padj)) {
    stop(sprintf(
      paste("Missing PeCorA p-value when treatment = %s,",
            "peptide = %s, and protein = %s"),
      treatment, peptide, protein))
  }

  if (length(padj) == 0) {
    stop(paste(feature, "not found within the results.")) ####
  }

  metadata <- pData(m) %>%
    select(!!sym(treatment)) %>%
    mutate(Sample = rownames(.))

  if (is.character(metadata[[treatment]])) {
    metadata[[treatment]] <- as.factor(metadata[[treatment]])
  }

  plot_df <- exprs(m)[rownames(f_data), ] %>%
    as.data.frame() %>%
    mutate(feature_name = rownames(.)) %>%
    pivot_longer(cols = -feature_name,
                 names_to = "Sample",
                 values_to = "value") %>%
    merge(f_data, by = "feature_name") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = !!sym(treatment)) %>%
    mutate(peptide_group = ifelse(!!sym(peptide_column) == peptide,
                                  peptide, "All other peptides")) %>%
    mutate(peptide_group = factor(peptide_group,
                                  levels = c("All other peptides",
                                             peptide))) %>%
    select(Sample, Condition, value, peptide_group)

  if (median_mod) {
    plot_df <- plot_df %>%
      group_by(Condition, Sample, peptide_group) %>%
      summarise(value = median(value, na.rm = TRUE))
  }

  ## Boxplot if condition is a factor. Otherwise, scatterplots when using
  ## numerical data.
  if (is.factor(plot_df$Condition)) {
    p <- ggplot(plot_df, aes(x = Condition, y = value,
                             fill = peptide_group)) +
      geom_boxplot(outlier.shape = NA, na.rm = TRUE)
  } else {
    plot_df <- plot_df %>%
      mutate(alpha = ifelse(peptide_group == "All other peptides",
                            0.005, 0.25))

    lm_df <- filter(plot_df, peptide_group == "All other peptides")
    allothers_lm <- lm(value ~ Condition, data = lm_df)

    lm_df <- filter(plot_df, peptide_group == peptide)
    chosen_lm <- update(allothers_lm, data = lm_df)

    p <- ggplot(plot_df, aes(x = Condition, y = value,
                             color = peptide_group, alpha = alpha)) +
      geom_point(na.rm = TRUE) +
      geom_abline(intercept = coef(allothers_lm)[1],
                  slope = coef(allothers_lm)[2],
                  color = "red", size = 1) +
      geom_abline(intercept = coef(chosen_lm)[1],
                  slope = coef(chosen_lm)[2],
                  color = scales::hue_pal()(2)[2], size = 1) +
      guides(alpha = "none")
  }

  p <- p +
    labs(title = sprintf("Protein: %s", protein),
         subtitle = paste("BH-adjusted p-value =", padj),
         x = treatment,
         y = "Log Intensity")

  return(p)
}


utils::globalVariables(
  c("pvalue", "peptide_group", "Sample", "Condition")
)

