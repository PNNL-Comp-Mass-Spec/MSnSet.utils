#' @title Peptide significance map
#'
#' @description Mapping peptides onto the protein sequence with different ways
#'   to show statistical significance.
#'
#' @param x `data.frame` with required columns Last_AA, First_AA, and ProtLen,
#'   as well as any columns specified by `p_val_from`, `accession_col`,
#'   `title_from`, and `fill_by` (only if `type` is `"continuous"` or
#'   `"manhattan"`). Recommended to provide the `fData` table of an MSnSet
#'   object.
#' @param accession character; name of the protein (or accession, in general
#'   terms) for which to plot peptides.
#' @param p_val_from character. Name of the column in `x` with p-values.
#' @param type character. Different ways of representing p-values. Types are:
#'   "continuous" - continuous color scheme, "discrete" - discrete color scheme,
#'   "manhattan" - p-values shown on y-axis (akin to a Manhattan plot).
#' @param accession_col character. Column in `x` that `accession` is taken from.
#'   Typically, it is "accession".
#' @param title_from character. One or more columns in `x` used to construct the
#'   plot title. Defaults to `accession_col`.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide rectangle.
#' @param step_y numeric; step size to resolve peptide overlaps.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#' @param ... other arguments passed to functions determined by `type`. See
#'   details for more information.
#'
#' @details Parameters passed to helper functions determined by the value of
#'   `type`:
#'
#'   \describe{
#'     \item{fill_by}{character; (`type = c("continuous", "manhattan")`) name of
#'     a column in `x` used to fill the peptide rectangles.}
#'   }
#'
#' @return \code{ggplot} object.
#'
#' @importFrom dplyr filter mutate arrange distinct across %>% vars
#' @importFrom rlang !! sym
#'
#' @export peptide_significance_map
#'
#' @md
#'
#' @examples
#' data("peptide_count_data")
#'
#' library(msmsTests)
#'
#' alt.f <- "y ~ group + 1"
#' null.f <- "y ~ 1"
#' div <- colSums(exprs(peptide_count_data)) # normalization factor
#' res <- msms.glm.qlll(peptide_count_data, alt.f, null.f, div = div)
#' res$p.val.adj <- p.adjust(res$p.value, "BH")
#'
#' # adding significance metrics
#' f_data <- fData(peptide_count_data)
#' f_data$P.Value <- res[featureNames(peptide_count_data), "p.value"]
#' f_data$adj.P.Val <- res[featureNames(peptide_count_data), "p.val.adj"]
#' fData(peptide_count_data) <- f_data
#'
#' # Path to FASTA file
#' fasta <- system.file("extdata/FASTAs",
#'                      "H_sapiens_UniProt_SPROT_2021-06-20.fasta.gz",
#'                      package = "MSnSet.utils")
#'
#' peptide_count_data <- map_peptide_position(object = peptide_count_data,
#'                                            fasta = fasta,
#'                                            accession_col = "protein",
#'                                            peptide_col = "Sequence")
#'
#' peptide_significance_map(x = fData(peptide_count_data),
#'                          accession = "sp|O15240|VGF_HUMAN",
#'                          accession_col = "protein",
#'                          p_val_from = "adj.P.Val",
#'                          type = "manhattan")

peptide_significance_map <- function(x,
                                     accession,
                                     p_val_from,
                                     type = c("discrete",
                                              "continuous",
                                              "manhattan"),
                                     accession_col = "accession",
                                     title_from = accession_col,
                                     min_y = 0,
                                     prot_width_y = 0.025,
                                     pep_width_y = 0.025,
                                     step_y = 0.033,
                                     border_color = "white",
                                     border_size = NULL,
                                     aa_step = 20,
                                     ...)
{
  x <- x %>%
    filter(!!rlang::sym(accession_col) == !!accession) %>%
    mutate(Length = Last_AA - First_AA + 1) %>%
    arrange(First_AA, -Length)

  # === plot type selector ===
  type <- match.arg(type, several.ok = FALSE)
  fun_switch <- list(peptide_significance_map_discrete,
                     peptide_significance_map_continuous,
                     peptide_significance_map_manhattan)
  names(fun_switch) <- c("discrete", "continuous", "manhattan")
  FUN <- fun_switch[[type]]
  p <- FUN(x, p_val_from, ...)

  # === plot title ===
  prot_name_cols <- distinct(x, across(title_from))
  if (nrow(prot_name_cols) == 1) {
    prot_name <- paste0(prot_name_cols, collapse = ", ")
  } else {
    stop("title_from must be unique")
  }

  p <- p +
    labs(title = prot_name)

  return(p)
}



# ===== Helper functions =======================================================


#' @title Manhattan-style arrangement of peptides
#'
#' @param prot `data.frame` containing peptides from a single protein.
#' @param p_val_from character; name of a column in `prot` providing statistical
#'   significance of each peptide.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide rectangle.
#'
#' @details Plots peptides according to their statistical significance.
#'   Specifically, this function will vertically-center peptides according to
#'   their -log10(`p_val_from`) scaled to the largest value or a value of 5.

layout_manhattan <- function(prot,
                             p_val_from,
                             min_y = 0,
                             width_y = 0.025)
{
  prot$significance <- -log10(prot[[p_val_from]])
  sig_scaler <- max(c(ceiling(max(prot$significance)), 5))

  prot$ymin <- min_y + (prot$significance / sig_scaler) - (width_y / 2)
  prot$ymax <- prot$ymin + width_y

  return(prot)
}


#' @title Compact stacking of peptides
#'
#' @param prot `data.frame` containing peptides from a single protein with
#'   columns "First_AA" and "Last_AA".
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide. Must be less than
#'   `step_y`.
#' @param step_y numeric; step size to resolve peptide overlaps.
#'
#' @importFrom dplyr %>% slice filter pull

layout_compact_stack <- function(prot,
                                 min_y = 0,
                                 width_y = 0.025,
                                 step_y = 0.033)
{
  stopifnot(step_y > width_y)

  ## If there is more than one peptide:
  # Check if the second peptide overlaps the first. If so, increase
  # y position by step_y until the overlap is resolved. Reset the y position
  # to min_y and repeat this process for each of the remaining peptides.
  prot$ymin <- min_y

  if (nrow(prot) > 1) {
    for (i in 2:nrow(prot)) {
      current_y <- min_y

      while (TRUE) {
        max_last_residue <- prot %>%
          dplyr::slice(1:(i-1)) %>%
          dplyr::filter(ymin == current_y) %>%
          pull(Last_AA) %>%
          max()

        if (max_last_residue + 0 >= prot[i, "First_AA"]) {
          current_y <- current_y + step_y
        } else {
          break()
        }
      }

      prot[i, "ymin"] <- current_y
    }
  }

  prot$ymax <- prot$ymin + width_y

  return(prot)
}


#' @title Plot peptides on protein sequence
#'
#' @param prot `data.frame` for exactly one protein with columns "ProtLen",
#'   "ymax", "ymin", "First_AA", "Last_AA", and a column specified by `fill_by`.
#'   `prot` should first be modified with `layout_manhattan` or
#'   `layout_compact_stack`.
#' @param fill_by character; name of a column in `prot` used to fill the peptide
#'   rectangles.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#'
#' @importFrom dplyr %>% distinct
#' @importFrom ggplot2 ggplot geom_rect scale_x_continuous labs theme_classic
#'   theme element_text element_blank element_rect
#' @importFrom rlang !! sym

plot_peptide_on_protein_map <- function(prot,
                                        fill_by,
                                        border_color = "white",
                                        border_size = NULL,
                                        aa_step = 20,
                                        env = parent.frame())
{

  prot_width_y <- env$prot_width_y

  prot_len <- prot %>%
    distinct(ProtLen) %>%
    as.numeric()

  if (length(prot_len) != 1) {
    stop("prot must contain exactly one protein")
  }

  if (is.null(border_size)) {
    border_size <- 1.5 / log10(prot_len)
  }

  # width_y <- prot$ymax[[1]] - prot$ymin[[1]]

  p <-
    ggplot(data = prot) +
    # protein rectangle
    geom_rect(aes(xmin = 1 - 0.5,
                  xmax = prot_len + 0.5,
                  ymin = -prot_width_y * 2,
                  ymax = -prot_width_y)) +
    # peptide rectangles
    geom_rect(aes(xmin = First_AA - 0.5,
                  xmax = Last_AA + 0.5,
                  ymin = ymin,
                  ymax = ymax,
                  fill = !!rlang::sym(fill_by)),
              color = border_color,
              size = border_size) +
    scale_x_continuous(breaks = c(1, seq(aa_step, prot_len, aa_step)),
                       expand = rep(0.02, 2)) +
    labs(x = "residue",
         y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line.x = element_blank(),
          legend.key = element_rect(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 16))

  return(p)
}


#' @title Plot colored peptide map
#'
#' @param prot `data.frame` containing peptides from a single protein with
#'   columns "First_AA" and "Last_AA".
#' @param fill_by character; name of a column in `prot` used to fill the peptide
#'   rectangles.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide. Must be less than
#'   `step_y`.
#' @param step_y numeric; step size to resolve peptide overlaps.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#'
#' @importFrom ggplot2 scale_y_continuous theme

plot_colored_peptide_map <- function(prot,
                                     fill_by,
                                     min_y = 0,
                                     width_y = 0.025,
                                     step_y = 0.033,
                                     border_color = "white",
                                     border_size = NULL,
                                     aa_step = 20)
{
  # arranging peptide positions - adds ymin and ymax columns
  prot <- layout_compact_stack(prot,
                               min_y = min_y,
                               width_y = width_y,
                               step_y = step_y)

  y_lims <- if (max(prot$ymax) < 0.45) c(-width_y * 2, 0.45) else NULL

  p <- plot_peptide_on_protein_map(prot = prot,
                                   fill_by = fill_by,
                                   border_color = border_color,
                                   border_size = border_size,
                                   aa_step = aa_step)

  p <- p +
    scale_y_continuous(expand = rep(0, 2),
                       limits = y_lims) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())

  return(p)
}


#' @title Discrete Peptide Significance
#'
#' @description Plot peptides on their protein sequence and fill according to
#'   their discreetized p-value.
#'
#' @param x `data.frame` with columns
#' @param p_val_from character; name of a column in x containing p-values.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide rectangle.
#' @param step_y numeric; step size to resolve peptide overlaps.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#'
#' @details The function cuts p-values into intervals \[0, 0.0001\), \[0.0001,
#'   0.001\), \[0.001, 0.01\), \[0.01, 0.05\), and \[0.05, 1\]. Each interval is
#'   then assigned a color to fill the peptide rectangles.
#'
#' @importFrom ggplot2 scale_fill_manual

peptide_significance_map_discrete <- function(x,
                                              p_val_from,
                                              min_y = 0,
                                              width_y = 0.025,
                                              step_y = 0.033,
                                              border_color = "white",
                                              border_size = NULL,
                                              aa_step = 20)
{
  cuts <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  x[[p_val_from]] <- cut(x[[p_val_from]], cuts, ordered_result = TRUE,
                         right = FALSE, include.lowest = TRUE)

  # Maybe define these fill arguments somewhere else in the future
  fillings <-
    scale_fill_manual(
      breaks = levels(x[[p_val_from]]),
      values = c("#E53935", "#FDD835", "#43A047", "#1E88E5", "#BDBDBD"),
      limits = levels(x[[p_val_from]]),
      labels = c("< 1e-04", "[1e-04, 1e-03)", "[1e-03, 0.01)",
                 "[0.01, 0.05)", expression(paste("">=0.05)))
    )

  p <- plot_colored_peptide_map(prot = x,
                                fill_by = p_val_from,
                                min_y = min_y,
                                width_y = width_y,
                                step_y = step_y,
                                border_color = border_color,
                                border_size = border_size,
                                aa_step = aa_step)
  p <- p +
    fillings

  return(p)
}


#' @title Plot peptides on protein sequence with gradient fill
#'
#' @param x `data.frame`.
#' @param p_val_from character; name of a column in `x` containing p-values.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide rectangle.
#' @param step_y numeric; step size to resolve peptide overlaps.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#'
#' @details A gradient fill dependent on the values of `p_val_from` is used for
#'   the peptide rectangles.
#'
#' @importFrom scales pretty_breaks trans_new log_breaks
#' @importFrom ggplot2 scale_fill_gradientn

peptide_significance_map_continuous <- function(x,
                                                p_val_from,
                                                min_y = 0,
                                                width_y = 0.025,
                                                step_y = 0.033,
                                                border_color = "white",
                                                border_size = NULL,
                                                aa_step = 20)
{
  p <- plot_colored_peptide_map(prot = x,
                                fill_by = p_val_from,
                                min_y = min_y,
                                width_y = width_y,
                                step_y = step_y,
                                border_color = border_color,
                                border_size = border_size,
                                aa_step = aa_step)

  # Round lowest log10(p-value) to the nearest 0.5 (more negative)
  min_pval <- min(p$data[[p_val_from]], na.rm = TRUE)
  min_pval <- floor(log10(min_pval) / 0.5) * 0.5

  # Create pretty breaks on the log10 scale and undo transformation
  breaks <- 10 ^ (pretty_breaks()(0:min_pval))
  breaks <- signif(breaks, digits = 1)

  labels <- ifelse(breaks > 1e-3,
                   format(breaks, scientific = FALSE, drop0trailing = TRUE),
                   format(breaks))

  log10_rev_trans <- trans_new(name = "log10_rev",
                               transform = function(x) -log10(x),
                               inverse = function(x) 10 ^ (-x),
                               breaks = function(x) rev(log_breaks(10)(x)))

  p <- p +
    scale_fill_gradientn(trans = log10_rev_trans,
                         colors = MSnSet.utils::jet2.colors(21),
                         breaks = breaks,
                         labels = breaks,
                         limits = c(1, min(breaks)))

  return(p)
}


#' @title Manhattan arrangement of peptides with p-value fill
#'
#' @param prot `data.frame` with a column for p-values and all columns required
#'   by `layout_manhattan`.
#' @param p_val_from character; name of a column in `prot` containing p-values.
#' @param fill_by character; name of a column in `prot` used to fill the peptide
#'   rectangles. Default is `NULL`, which will use `p_val_from`.
#' @param min_y numeric; shifts peptides vertically.
#' @param width_y numeric; vertical width of each peptide rectangle.
#' @param border_color character; color used for the borders of the peptide
#'   rectangles.
#' @param border_size numeric or `NULL`; linewidth of the borders around the
#'   peptide rectangles. If `NULL` (default), the linewidth will be calculated
#'   from the protein length.
#' @param aa_step numeric; distance between x axis breaks, starting at `aa_step`
#'   and continuing along the length of the protein.
#'
#' @importFrom scales pretty_breaks scientific_format
#' @importFrom ggplot2 scale_y_continuous geom_hline labs scale_fill_gradient

peptide_significance_map_manhattan <- function(prot,
                                               p_val_from,
                                               fill_by = NULL,
                                               min_y = 0,
                                               # width_y = 0.025,
                                               border_color = "white",
                                               border_size = NULL,
                                               aa_step = 20,
                                               env = parent.frame())
{

  prot_width_y <- env$prot_width_y
  pep_width_y <- env$pep_width_y
  # arranging peptide positions
  prot <- layout_manhattan(prot = prot,
                           p_val_from = p_val_from,
                           min_y = min_y,
                           width_y = pep_width_y)

  if (is.null(fill_by)) {
    fill_by <- p_val_from ## Dummy variable
  }

  p <- plot_peptide_on_protein_map(prot = prot,
                                   fill_by = fill_by,
                                   border_color = border_color,
                                   border_size = border_size,
                                   aa_step = aa_step)

  prot$significance <- -log10(prot[[p_val_from]])
  sig_scaler <- num_labels <- max(c(ceiling(max(prot$significance)), 5))

  # some integer number between 5 and 10, inclusive
  if (num_labels > 10) num_labels <- 10
  if (num_labels < 5) num_labels <- 5

  lblz <- pretty_breaks(num_labels)(0:sig_scaler)
  # width_y <- prot$ymax[[1]] - prot$ymin[[1]]

  p <- p +
    scale_y_continuous(breaks = lblz / sig_scaler,
                       labels = scientific_format(digits = 1)(10 ^ (-lblz)),
                       limits = c(-2 * prot_width_y, 1 + prot_width_y),
                       expand = rep(0, 2)) +
    geom_hline(yintercept = -log10(0.05) / sig_scaler,
               color = "grey",
               linetype = "dashed") +
    labs(y = p_val_from)

  # handling peptide fill
  if (fill_by == p_val_from) {
    p <- p +
      scale_fill_gradient(low = "#212121", high = "#212121", guide = NULL)
  }

  return(p)
}



utils::globalVariables(
  c("accession", "cleanSeq", "First_AA_First", "Last_AA_First", "ProtLen",
    "First_AA", "Last_AA", "ymin", "Length", "ymax")
)

