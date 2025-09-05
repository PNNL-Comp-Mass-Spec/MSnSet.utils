#' @title Peptide Correlation Analysis (PeCorA)
#'
#' @description Perform Peptide Correlation Analysis (PeCorA) as described in
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/}.
#'
#' @param m \code{MSnSet} object with pData and fData tables. fData must have
#'   columns \code{protein_column} and \code{peptide_column}. The exprs data
#'   should be aggregated such that features are unique combinations of peptide
#'   and protein.
#' @param proteins character; one or more proteins in
#'   \code{fData(m)[[protein_column]]} with peptides to test. Only those
#'   proteins with at least 2 peptides will be considered.
#' @param treatment character; the name of column in \code{pData(m)} containing
#'   treatment information. The data can be numeric or categorical.
#'   \code{PeCorA} will test for an interaction between \code{treatment} and all
#'   peptides that map to the proteins in \code{proteins}.
#' @param control_group character; optional control group. Only applicable if
#'   \code{treatment} is categorical. See \code{center_peptides} parameter below
#'   for details.
#' @param protein_column character; the name of a column in \code{fData(m)}
#'   containing protein information. Default is "Protein".
#' @param peptide_column character; the name of a column in \code{fData(m)}
#'   containing peptide information. Default is "Peptide".
#' @param standardize_samples logical; whether to standardize samples (columns)
#'   in exprs data.
#' @param center_peptides logical; whether to mean-center each peptide. If the
#'   data in the \code{treatment} column is categorical and \code{control_group}
#'   is not \code{NULL}, peptides will be zero-centered on the control group
#'   mean. Otherwise, the mean of all samples are used for centering.
#' @param median_mod logical; whether to use a modified version of PeCorA. See
#'   details.
#' @param verbose logical; whether to produce progress messages.
#'
#' @details By default, the MSnSet provided by \code{m} will be normalized as
#'   described in the PeCorA publication: standardize columns and mean-center
#'   peptides using the mean of the \code{control_group} samples (if provided)
#'   or the mean of all samples.
#'
#'   If \code{median_mod = TRUE}, the intensities of "all other peptides" will
#'   be aggregated to their per-sample medians. This way, all peptides mapping
#'   to a given protein (aside from the peptide of interest as the function
#'   loops over all peptides) are not treated as 'independent' observations. The
#'   effect is that p-values are much higher (less significant) than in the
#'   original version of PeCorA.
#'
#'   This function (as well as any helper functions) was based on those from the
#'   PeCorA R package (\url{https://github.com/jessegmeyerlab/PeCorA}).
#'
#' @return An \code{MSnSet} where the fData table includes new or updated
#'   columns "pecora_pval_{treatment}" and "pecora_adj_pval_{treatment}". The 
#'   p-values are adjusted on a per protein basis, using the BH method.
#'
#' @seealso \code{\link[MSnSet.utils]{plot_PeCorA}}
#'
#' @references Dermit, M., Peters-Clarke, T. M., Shishkova, E., & Meyer, J. G.
#'   (2021). Peptide Correlation Analysis (PeCorA) Reveals Differential
#'   Proteoform Regulation. \emph{Journal of proteome research, 20(4)},
#'   1972–1980. \url{https://doi.org/10.1021/acs.jproteome.0c00602}
#'
#' @importFrom dplyr select mutate filter rename %>% distinct
#' @importFrom rlang sym
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect any_of all_of
#'
#' @export PeCorA

PeCorA <- function(m,
                   proteins,
                   treatment,
                   control_group = NULL,
                   protein_column = "Protein",
                   peptide_column = "Peptide",
                   standardize_samples = TRUE,
                   center_peptides = TRUE,
                   median_mod = TRUE,
                   verbose = FALSE)
{
  ## Preprocessing step
  if (verbose) {
    message("Preprocessing MSnSet for PeCorA")
  }

  m <- .PeCorA_preprocess(m = m,
                          treatment = treatment,
                          control_group = control_group,
                          standardize_samples = standardize_samples,
                          center_peptides = center_peptides)
  f_data <- fData(m)

  ## Must have "Protein" and "Peptide" columns
  for (col_i in c(protein_column, peptide_column)) {
    if (!col_i %in% colnames(f_data)) {
      stop(sprintf("fData missing '%s' column.", col_i))
    }
  }

  pval_label <- paste0("pecora_pval_", treatment)
  padj_label <- paste0("pecora_adj_pval_", treatment)

  ## Remove previous results for the same treatment variable.
  f_data <- select(f_data, -any_of(c(pval_label, padj_label)))

  ## Filter to just those proteins for which we want to run PeCorA
  peptide_mapping <- f_data %>%
    filter(!!sym(protein_column) %in% proteins) %>%
    select(all_of(c(protein_column, peptide_column))) %>%
    distinct()

  mat <- exprs(m)[rownames(peptide_mapping), ]

  metadata <- pData(m) %>%
    select(!!sym(treatment)) %>%
    mutate(Sample = rownames(.))

  if (is.character(metadata[[treatment]])) {
    metadata[[treatment]] <- as.factor(metadata[[treatment]])
  }

  # Prepare input and run PeCorA
  PeCorA_input <- mat %>%
    cbind(peptide_mapping) %>%
    as.data.frame() %>%
    pivot_longer(cols = all_of(colnames(mat)),
                 names_to = "Sample",
                 values_to = "LogRatio") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = !!sym(treatment))

  if (verbose) {
    message("Performing PeCorA")
  }

  PeCorA_result <- .PeCorA_all_proteins(x = PeCorA_input,
                                        protein_column = protein_column,
                                        peptide_column = peptide_column,
                                        median_mod = median_mod)

  # Rename p-value and adjusted p-value columns
  pval_index <- which(colnames(PeCorA_result) == "pvalue")
  adj_pval_index <- which(colnames(PeCorA_result) == "adj_pval")
  colnames(PeCorA_result)[[pval_index]] <- pval_label
  colnames(PeCorA_result)[[adj_pval_index]] <- padj_label

  # Add PeCorA results to fData
  if (verbose) {
    message("Updating fData(m).")
  }

  f_data <- f_data %>%
    mutate(feature_name = rownames(f_data)) %>%
    left_join(PeCorA_result, by = c(protein_column, peptide_column))
  rownames(f_data) <- f_data$feature_name

  f_data <- select(f_data, -feature_name)

  fData(m) <- f_data[rownames(exprs(m)), ]

  return(m)
}

