#' Prepare Data for Density Plots
#'
#' Reshapes an \code{MSnSet}'s exprs data into a long-format data frame
#' ("long-format" being "each row is a single expression value."), suitable for density
#' and ridge plots using \code{ggplot2}.
#'
#' @param msnset MSnSet. The object containing expression data and pData.
#' @param facet_by Character. pData column name used for faceting.
#' @param fill_by Character. pData column name used for fill colors.
#' @param do_paste Logical. Whether to append the column name \code{facet_by}
#'    to its values. For example, if \code{facet_by = "tmt_plex"} and a value is
#'   \code{"Plex1"}, it becomes \code{"tmt_plex Plex1"}. Default is TRUE.
#'
#' @importFrom Biobase exprs pData sampleNames
#' @importFrom dplyr %>% across arrange filter inner_join mutate select where
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#'
#' @return A data.frame in long format.
.get_density_df <- function(msnset, facet_by, fill_by, do_paste = TRUE) {
   feature_id <- NULL
   df <- exprs(msnset) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("feature_id") %>%
      tidyr::pivot_longer(
         cols = -feature_id,
         names_to = "sample",
         values_to = "value"
      ) %>%
      dplyr::filter(
         !is.na(value)
      ) %>%
      inner_join(
         y = pData(msnset) %>%
            tibble::rownames_to_column("sample") %>%
            select(
               sample,
               !!sym(facet_by),
               !!sym(fill_by)
            ) %>%
            dplyr::filter(
               !is.na(!!sym(fill_by)),
            ),
         by = "sample"
      ) %>%
      arrange(
         !!sym(facet_by),
         # value
      ) %>%
      mutate(
         sample = factor(
            x = sample,
            levels = intersect(
               sampleNames(msnset),
               sample
            )
         ),
         !!sym(facet_by) := as.factor(!!sym(facet_by)),
         across(
            .cols = where(is.character),
            .fns = ~ factor(.x, levels = unique(.x))
         )
      )

   if (do_paste) {
      df <- df %>% mutate(!!sym(facet_by) := paste(facet_by, !!sym(facet_by)))
   }
   df
}


#' Create Faceted Density Ridge Plot
#'
#' Generates density ridge plots for each sample, chunked into faceted groups.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param facet_by Character. The pData column to group plots by.
#' @param fill_by Character. The pData column used to color the ridges.
#'
#'
#' @details Tip: to get a consistent/desired font, add to the output of this function the following:
#' \code{ggplot2::theme(
#'          text = element_text(family = "Font Name"),
#'          axis.text.x = element_text(family = "Font Name", size = 4
#'       )}
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes facet_wrap scale_fill_manual coord_cartesian rel element_blank element_rect element_text theme labs after_stat
#' @importFrom viridis scale_fill_viridis
#' @importFrom tools toTitleCase
#' @importFrom forcats fct_rev
#'
#' @return A ggplot object.
#' @export
create_faceted_density_plot <- function(msnset, facet_by, fill_by) {
   df <- .get_density_df(msnset, facet_by, fill_by)


   p <- ggplot(
      data = df,
      mapping = aes(
         x = value,
         y = forcats::fct_rev(sample),
         fill = .data[[fill_by]],
         height = after_stat(density)
      )
   ) +
      geom_density_ridges(
         color = "black",
         stat = "density",
         linewidth = 0.25
      ) +
      facet_wrap(
         ~ .data[[facet_by]],
         ncol = 4,
         scales = "free_y"
      ) +
      scale_fill_viridis(
         discrete = TRUE,
         begin = 0.3,
         end = 0.8,
         name = tools::toTitleCase(gsub("_", " ", fill_by))
      ) +
      coord_cartesian(
         clip = "off"
      ) +
      labs(
         x = NULL,
         y = NULL
      ) +
      theme_ridges(
         font_size = 10
      ) +
      theme(
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         strip.background = element_rect(
            fill = "white",
            color = NA
         ),
         strip.text = element_text(
            hjust = 0,
            color = "black",
            size = rel(0.8)
         ),
         plot.background = element_rect(
            fill = "white",
            color = NA
         ),
      )
   p
}


#' Create Layered Density Plot
#'
#' Overlays density distribution curves for all samples on a single set of axes.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param color_by Character. The pData column used to color the density lines.
#' @importFrom ggplot2 ggplot aes geom_density geom_vline labs theme_minimal theme element_line element_rect element_text
#' @importFrom viridis scale_color_viridis
#' @importFrom glue glue
#'
#' @return A ggplot object.
#' @export
create_layered_density_plot <- function(msnset, color_by) {
   df <- .get_density_df(msnset, facet_by = color_by, fill_by = color_by, do_paste = FALSE)

   n_features_label <- glue::glue("Zero-centered log2 values per sample.\n({nrow(msnset)} total features)")

   p <- ggplot(
      data = df,
      mapping = aes(
         x = value,
         group = sample,
         color = .data[[color_by]]
      )
   ) +
      geom_density() +
      geom_vline(
         xintercept = 0,
         color = "grey",
         lty = "dashed"
      ) +
      scale_color_viridis(
         discrete = TRUE,
         begin = 0.3,
         end = 0.8,
         name = tools::toTitleCase(gsub("_", " ", color_by))
      ) +
      labs(
         x = "Median-centered log2 Expression Level",
         y = "Density",
         subtitle = n_features_label
      ) +
      theme(
         # axis.ticks = element_line(
         #    color = "black"
         # ),
         # legend.position = "inside",
         # legend.justification = c(0.85, 0.75),
         # axis.line = element_line(
         #    color = "black"
         # ),
         # panel.grid = element_blank(),
         # plot.background = element_rect(
         #    fill = "white",
         #    color = NA
         # )
         text = element_text(family = "Fira Code")
      )

   return(p)
}

#' Create Overall Sample Correlation Heatmap
#'
#' Calculates and visualizes Pearson correlation between all sample pairs
#' in an MSnSet.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param annotation_by Character. The pData column name used to annotate samples.
#' @param width unit. The width of the heatmap.
#' @param height unit. The height of the heatmap.
#'
#' @details This function requires width and height to be specified because
#' the underlying plotting technology (ComplexHeatmap) is best used with explicit sizing.
#'
#' @importFrom Biobase exprs pData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom grid unit gpar
#' @importFrom circlize colorRamp2
#' @importFrom viridis viridis
#' @importFrom tools toTitleCase
#' @importFrom glue glue
#' @importFrom stats cor
#'
#' @return A HeatmapList object, returned invisibly after drawing.
#' @export
create_overall_correlation_heatmap <- function(msnset, annotation_by, width, height) {
   cor_mat <- cor(
      x = exprs(msnset),
      use = "pairwise.complete.obs",
      method = "pearson"
   )

   cor_mat_idxes <- seq_len(nrow(cor_mat))
   cor_mat <- cor_mat[rev(cor_mat_idxes), rev(cor_mat_idxes)]

   if (!annotation_by %in% colnames(pData(msnset))) {
      stop(glue::glue("annotation_by column '{annotation_by}' not found in pData(msnset). I only see: {paste(colnames(pData(msnset)), collapse=', ')}"))
   }

   # Build row/column annotation from pData(msnset)[[annotation_by]]
   ann_vec <- pData(msnset)[[annotation_by]]
   names(ann_vec) <- rownames(pData(msnset))

   # Ensure ordering matches the correlation matrix
   ann_vec <- ann_vec[colnames(cor_mat)]
   ann_factor <- factor(ann_vec)
   ann_levels <- levels(ann_factor)

   # Colors for annotation levels
   ann_cols <- viridis::viridis(n = length(ann_levels), begin = 0.3, end = 0.8)
   names(ann_cols) <- ann_levels

   # Common argument list for HeatmapAnnotation / rowAnnotation
   ann_args <- list()
   ann_args[[annotation_by]] <- ann_factor
   ann_args$col <- setNames(list(ann_cols), annotation_by)
   ann_args$annotation_name_gp <- gpar(fontfamily = "Fira Code")
   ann_args$annotation_legend_param <- list(
      title = tools::toTitleCase(gsub("_", " ", annotation_by)),
      title_gp = gpar(fontfamily = "Fira Code"),
      labels_gp = gpar(fontfamily = "Fira Code")
   )

   ann_args_left <- ann_args
   ann_args_left$show_legend <- FALSE

   # Create top (column) and left (row) annotations using the same arg list
   top_annotation <- do.call(ComplexHeatmap::HeatmapAnnotation, ann_args)
   left_annotation <- do.call(ComplexHeatmap::rowAnnotation, ann_args_left)

   ht <- Heatmap(
      matrix = cor_mat,
      col = circlize::colorRamp2(
         breaks = c(-1, 0, 1),
         colors = c("blue", "white", "red")
      ),
      heatmap_legend_param = list(
         title = "Pearson\nCorrelation",
         title_gp = gpar(fontfamily = "Fira Code"),
         labels_gp = gpar(fontfamily = "Fira Code"),
         legend_direction = "vertical"
      ),
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 2, fontfamily = "Fira Code"),
      height = height,
      width = width,
      row_title = "Samples",
      column_title = "Samples",
      row_title_gp = gpar(fontfamily = "Fira Code"),
      column_title_gp = gpar(fontfamily = "Fira Code"),
      top_annotation = top_annotation,
      left_annotation = left_annotation
   )

   d <- draw(ht,
      align_heatmap_legend = "heatmap_center", merge_legend = TRUE
   )

   return(d)
}


#' Create Detailed Sample Correlation Heatmap
#'
#' Generates a pairwise correlation heatmap for a specific set of samples with
#' numeric correlation values displayed in cells.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param the_samples Character. Vector of sample names to include.
#' @param the_samples_relabling_fn Function. A function applied to sample IDs
#'    for display on the heatmap axes.
#' @param width unit. The width of individual heatmap rectangles.
#' @param height unit. The height of individual heatmap rectangles.
#'
#' @details This function requires width and height to be specified because
#' the underlying plotting technology (ComplexHeatmap) is best used with explicit sizing.
#'
#' @importFrom Biobase exprs
#' @importFrom ComplexHeatmap Heatmap draw max_text_width
#' @importFrom grid unit gpar grid.rect grid.text
#' @importFrom circlize colorRamp2
#' @importFrom stats cor
#'
#' @return A HeatmapList object, returned invisibly after drawing.
#' @export
create_samples_correlation_heatmap <- function(msnset, the_samples, the_samples_relabling_fn, width, height) {
   cor_mat <- cor(
      exprs(msnset)[, the_samples],
      use = "pairwise.complete.obs"
   )

   new_the_samples <- the_samples_relabling_fn(the_samples)
   rownames(cor_mat) <- colnames(cor_mat) <- new_the_samples

   ht <- Heatmap(
      matrix = cor_mat,
      col = circlize::colorRamp2(
         breaks = c(-1, 0, 1),
         colors = c("#3366ff", "white", "darkred")
      ),
      rect_gp = gpar(
         type = "none"
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      row_names_max_width = max_text_width(
         text = rownames(cor_mat),
         gp = gpar(
            fontsize = 10
         )
      ),
      column_names_max_height = max_text_width(
         text = colnames(cor_mat),
         gp = gpar(
            # fontsize = 10
            fontfamily = "Fira Code"
         )
      ),
      column_names_gp = gpar(
         # fontsize = 10,
         fontfamily = "Fira Code"
      ),
      row_names_gp = gpar(
         # fontsize = 10,
         fontfamily = "Fira Code"
      ),
      heatmap_legend_param = list(
         title = "Pearson Correlation",
         at = -1:1,
         # direction = "horizontal",
         # title_position = "lefttop",
         title_gp = gpar(fontfamily = "Fira Code"),
         labels_gp = gpar(fontfamily = "Fira Code")
         # legend_width = unit(80, "pt")
      ),
      height = height,
      width = width,
      cell_fun = function(j, i, x, y, w, h, fill) {
         if (i >= j) {
            # Only fill lower triangular part
            grid.rect(
               x = x,
               y = y,
               width = w,
               height = h,
               gp = gpar(
                  fill = fill,
                  col = "grey70",
                  fontfamily = "Fira Code"
               )
            )

            # Label correlations
            grid.text(
               label = round(cor_mat[i, j], digits = 2L),
               x = x,
               y = y,
               gp = gpar(
                  col = "white",
                  # Programatically get the default font size, times 0.7
                  fontsize = 9,
                  fontfamily = "Fira Code"
               )
            )
         }
      }
   )

   d <- draw(ht)
   d
}

# Expression heatmap
# Boxplot
#' Create Sample Boxplots
#'
#' Generates boxplots of expression distributions across samples,
#' faceted and colored by pData columns.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param fill_by Character. The pData column name for color-coding boxes.
#' @param facet_by Character. The pData column name for faceting the plot.
#'
#' @importFrom Biobase exprs pData sampleNames
#' @importFrom glue glue
#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text xlab vars
#' @importFrom ggh4x facet_nested
#' @importFrom dplyr left_join mutate
#' @importFrom viridis scale_fill_viridis
#'
#' @return A ggplot object.
#' @export
create_boxplot <- function(msnset, fill_by, facet_by) {
   Value <- NULL
   boxplot_df <- as.data.frame(as.table(exprs(msnset)))
   names(boxplot_df) <- c("feature", "sample", "Value")

   if (!fill_by %in% colnames(pData(msnset))) {
      stop(glue("fill_by column '{fill_by}' not found in pData(msnset). I only see: {paste(colnames(pData(msnset)), collapse=', ')}"))
   }
   if (!facet_by %in% colnames(pData(msnset))) {
      stop(glue("facet_by column '{facet_by}' not found in pData(msnset). I only see: {paste(colnames(pData(msnset)), collapse=', ')}"))
   }

   boxplot_df <- boxplot_df %>%
      left_join(pData(msnset) %>% mutate(sample = rownames(.)), by = "sample") %>%
      mutate(sample = factor(sample, levels = sampleNames(msnset)))

   p <- ggplot(
      boxplot_df,
      aes(
         x = sample,
         y = Value
      )
   ) +
      geom_boxplot(aes(fill = .data[[fill_by]]), outliers = F, median.linewidth = .75, linewidth = 0.25) +
      facet_nested(
         cols = vars({{ facet_by }}, .data[[facet_by]]),
         scales = "free_x",
      ) +
      theme(
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
      ) +
      xlab("Samples") +
      scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, direction = -1)

   p
}

#' Create Feature Count Plot
#'
#' Visualizes the number of features (e.g., proteins or peptides)
#' divided into complete and non-complete observations across sample groups.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param x_by Character. The pData column name used for the x-axis groups.
#' @param feature_kind Character. Label for the type of features
#'    (e.g., "Feature", "Protein", "Site"). Only used for legend labeling.
#'
#' @importFrom Biobase exprs pData
#' @importFrom ggplot2 ggplot aes geom_col theme element_text xlab position_stack
#' @importFrom viridis scale_fill_viridis
#' @importFrom glue glue
#'
#' @return A ggplot object.
#' @export
create_feature_count <- function(msnset, x_by, feature_kind = c("Feature", "Protein", "Site")) {
   Count <- NULL
   Type <- NULL
   x_to_num_complete <- c()
   x_to_num_any <- c()
   x_by_uniques <- unique(pData(msnset)[[x_by]])
   n_complete_features_overall <- sum(rowSums(!is.na(exprs(msnset))) == ncol(msnset))
   for (x in x_by_uniques) {
      x_msnset <- msnset[, pData(msnset)[[x_by]] == x]
      num_incomplete <- sum(rowSums(!is.na(exprs(x_msnset))) != ncol(x_msnset))
      x_to_num_complete[[as.character(x)]] <- n_complete_features_overall
      x_to_num_any[[as.character(x)]] <- num_incomplete
   }

   barplot_df <- data.frame(
      x = rep(x_by_uniques, 2),
      Count = c(unlist(x_to_num_complete), unlist(x_to_num_any)),
      Type = factor(
         rep(c("Complete", "Non-complete"), each = length(x_by_uniques)),
         levels = c("Complete", "Non-complete")
      )
   )

   p <- ggplot(
      barplot_df,
      aes(
         x = x,
         y = Count,
         fill = Type,
         group = factor(Type, levels = c("Non-complete", "Complete")) # Control stacking order
      )
   ) +
      geom_col(linewidth = 0.25, position = "stack") +
      theme(
         text = element_text(family = "Fira Code")
      ) +
      xlab(x_by) +
      scale_fill_viridis(discrete = TRUE, begin = 0.3, end = 0.8, name = glue::glue("{feature_kind} Type"))
   p
}


#' Create Missingness Plot
#'
#' Generates a visualization of missing values across samples using
#' \code{MSnbase::plotNA} as an engine.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param font_family Character. The font family to use for text elements.
#' @param feature_kind Character. Label for the type of features (e.g., "Feature").
#'   Only used for legend labeling.
#'
#' @importFrom MSnbase plotNA
#' @importFrom ggplot2 theme element_text ylim
#' @importFrom viridis scale_color_viridis
#' @importFrom glue glue
#'
#' @return A ggplot object.
#' @export
create_missingness <- function(msnset, font_family = "Courier New", feature_kind = "Feature") {
   p <- plotNA(msnset) +
      scale_color_viridis(
         discrete = TRUE,
         begin = 0.3,
         end = 0.8,
         name = "Feature Count",
         labels = c("proteins" = glue::glue("Complete {feature_kind}s"), "data" = "Complete Data")
      ) +
      ylim(0, 1)

   text_layer_idx <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomText")))

   if (length(text_layer_idx) > 0) {
      for (i in text_layer_idx) {
         p$layers[[i]]$aes_params$family <- font_family
      }
   }

   p
}


#' Create Multi-Contrast PCA Plots
#'
#' Generates a series of Principal Component Analysis (PCA) plots, each
#' colored by a different pData contrast.
#'
#' @param msnset MSnSet. The object containing expression data.
#' @param contrasts_to_test Character. Vector of pData columns for which to make a PCA plot.
#' @param colors_fn Function. A function that returns a color palette (e.g., \code{pals::glasbey}).
#' @param ret_i Integer. Optional. If provided, returns only the i-th plot.
#' @param legend_position Character. Positioning of the legend ("right", "bottom", etc.).
#' @param as_list Logical. If TRUE, returns a list of individual ggplot objects.
#' @param ... Additional arguments passed to \code{cowplot::plot_grid}.
#'
#'
#' @importFrom Biobase pData
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggtitle theme element_text scale_color_manual
#'    scale_fill_manual unit
#' @importFrom viridis scale_color_viridis scale_fill_viridis
#' @importFrom glue glue
#' @importFrom pals glasbey
#'
#' @return A ggplot object (grid) or a list of ggplot objects.
#' @export
create_multi_pca <- function(
  msnset,
  contrasts_to_test,
  colors_fn = glasbey,
  ret_i = NULL,
  legend_position = "right",
  as_list = FALSE,
  ...
) {
   pcas <- lapply(contrasts_to_test, function(ct) {
      is_discrete <- is.factor(pData(msnset)[[ct]]) || is.character(pData(msnset)[[ct]])
      if (is_discrete) {
         colors <- unname(colors_fn(n = length(unique(pData(msnset)[[ct]]))))
         cscale <- scale_color_manual(values = colors)
         fscale <- scale_fill_manual(values = colors)
      } else {
         cscale <- scale_color_viridis(discrete = is_discrete, begin = 0.1, end = 0.9)
         fscale <- scale_fill_viridis(discrete = is_discrete, begin = 0.1, end = 0.9)
      }
      pca <- plot_pca(msnset, ct) +
         cscale +
         fscale +
         theme(
            text = element_text(family = "Fira Code"),
            legend.position = "right",
            legend.justification = "right",
            legend.box.just = "right",
            legend.box.margin = unit(c(0, 0, 0, 0), "cm")
         )
      pca
   })


   pcas <- lapply(seq_along(pcas), function(i) {
      pcas[[i]] + theme(legend.position = legend_position)
   })
   if (!is.null(ret_i)) {
      return(pcas[[ret_i]])
   }

   if (!as_list) {
      p <- plot_grid(
         plotlist = pcas,
         align = "v",
         scale = 0.9,
         # rel_heights = c(0.5, rep(1, length(contrasts_to_test))),
         ...
      )
      p
   } else {
      pcas
   }
}


#' Generate QC Report
#'
#' Streamlines the creation of multiple QC visualizations and compiles them
#' into a Quarto/Markdown report.
#'
#' @param msnset MSnSet. The processed version (e.g., after normalization if desired).
#' @param msnset_orig MSnSet. The original object (for missingness plots).
#' @param p_data_contrasts Character. Vector of pData columns to use as
#'    primary and secondary grouping variables
#' @param the_sch_samples Character. Vector of sample names for the detailed
#'    correlation heatmap. If NULL and that plot is requested, an error will be raised.
#' @param report_output_path Character. File path where the report (\code{.qmd})
#'    will be saved.
#' @param plot_cache_dir Character. Optional. Directory to cache plot objects
#'    as \code{.rds} files.
#' @param report_title Character. Title of the report.
#' @param self_contained Logical. Whether the HTML report should be self-contained.
#' @param feature_kind Character. Label for features (e.g., "Protein"). Used for legend labeling.
#' @param img_dims List. Named list of numeric vectors defining \code{c(width, height)}
#'    for each plot type. This is the primary interface for controlling which plots are generated.
#'    Supported names are: \code{boxplot}, \code{feature_count}, \code{missingness}, \code{multi_pca},
#'    \code{overall_correlation_heatmap}, \code{samples_correlation_heatmap}, \code{faceted_density_plot}, and \code{layered_density_plot}.
#'
#' @details ADDINS.PNNL is required.
#'
#' @importFrom Biobase pData
#' @importFrom glue glue
#' @importFrom tools toTitleCase
#' @importFrom grid unit
#' @importFrom utils readRDS saveRDS
#' @importFrom graphics png dev.off
#' @importFrom knitr include_graphics
#'
#' @return None. Saves a file to \code{report_output_path} and renders it.
#' @export
qc_report <- function(msnset, msnset_orig, p_data_contrasts, the_sch_samples, report_output_path, plot_cache_dir = NULL, report_title = "QC Report", self_contained = TRUE, feature_kind = "Feature", img_dims = list(
                         boxplot = c(8, 6),
                         feature_count = c(8, 6),
                         missingness = c(8, 6),
                         multi_pca = c(8, 12)
                      )) {
   prim <- p_data_contrasts[1]
   seco <- p_data_contrasts[2]

   creating <- function(name) {
      message(glue("Creating {name}..."))
   }
   cache <- function(name, expr) {
      if (!is.null(plot_cache_dir)) {
         cache_path <- file.path(as.character(plot_cache_dir), paste0(name, ".rds"))
         if (file.exists(cache_path)) {
            message(glue("Loading cached {name} from {cache_path}"))
            readRDS(cache_path)
         } else {
            message(glue("Creating and caching {name} to {cache_path}"))
            result <- expr
            saveRDS(result, cache_path)
            result
         }
      } else {
         expr
      }
   }

   if ("boxplot" %in% names(img_dims)) {
      creating("boxplot")
      boxplot <- create_boxplot(
         msnset,
         fill_by = seco,
         facet_by = prim
      )
      cache("boxplot", boxplot)
   }

   if ("feature_count" %in% names(img_dims)) {
      creating("feature_count")
      feature_count <- create_feature_count(
         msnset,
         x_by = prim,
         feature_kind = feature_kind
      )
      cache("feature_count", feature_count)
   }

   if ("missingness" %in% names(img_dims)) {
      creating("missingness")
      tf <- tempfile(fileext = ".png")
      png(tf)
      on.exit(
         {
            dev.off()
            file.remove(tf)
         },
         add = TRUE
      )
      missingness <- create_missingness(
         msnset_orig,
         font_family = "Fira Code",
         feature_kind = feature_kind
      )
      cache("missingness", missingness)
   }

   if ("multi_pca" %in% names(img_dims)) {
      creating("multi_pca")
      multi_pca <- create_multi_pca(
         msnset,
         contrasts_to_test = p_data_contrasts,
         legend_position = "right",
         nrow = 1
      )
      cache("multi_pca", multi_pca)
   }

   if ("overall_correlation_heatmap" %in% names(img_dims)) {
      creating("overall_correlation_heatmap")
      tf <- tempfile(fileext = ".png")
      png(tf)
      on.exit(
         {
            dev.off()
            file.remove(tf)
         },
         add = TRUE
      )
      overall_correlation_heatmap <- create_overall_correlation_heatmap(
         msnset,
         annotation_by = seco,
         width = unit(img_dims$overall_correlation_heatmap[1] * 0.6, "in"),
         height = unit(img_dims$overall_correlation_heatmap[2] * 0.6, "in")
      )
      cache("overall_correlation_heatmap", overall_correlation_heatmap)
   }

   if ("samples_correlation_heatmap" %in% names(img_dims)) {
      if (is.null(the_sch_samples)) {
         stop("the_sch_samples must be provided to create samples_correlation_heatmap")
      }
      creating("samples_correlation_heatmap")
      tf <- tempfile(fileext = ".png")
      png(tf)
      on.exit(
         {
            dev.off()
            file.remove(tf)
         },
         add = TRUE
      )
      samples_correlation_heatmap <- create_samples_correlation_heatmap(
         msnset,
         the_samples = the_sch_samples,
         the_samples_relabling_fn = function(x) {
            gsub("_", " ", gsub("^X", "", x))
         },
         width = unit(img_dims$samples_correlation_heatmap[1] * 0.6, "in"),
         height = unit(img_dims$samples_correlation_heatmap[2] * 0.6, "in")
      )
      cache("samples_correlation_heatmap", samples_correlation_heatmap)
   }

   if ("faceted_density_plot" %in% names(img_dims)) {
      creating("faceted_density_plot")
      faceted_density_plot <- create_faceted_density_plot(
         msnset,
         facet_by = prim,
         fill_by = seco
      )
      cache("faceted_density_plot", faceted_density_plot)
   }

   if ("layered_density_plot" %in% names(img_dims)) {
      creating("layered_density_plot")
      layered_density_plot1 <- create_layered_density_plot(
         msnset,
         color_by = seco
      )
      layered_density_plot2 <- create_layered_density_plot(
         msnset,
         color_by = prim
      )
      layered_density_plot <- cowplot::plot_grid(
         layered_density_plot1,
         layered_density_plot2,
         ncol = 1,
         align = "vh",
         scale = 0.9
      )
      cache("layered_density_plot", layered_density_plot)
   }

   temp_img_dir <- tempdir()
   for (img_name in names(img_dims)) {
      img_path <- file.path(temp_img_dir, paste0(img_name, ".png"))
      img_dim <- img_dims[[img_name]]
      if (grepl("heat", img_name)) {
         png(
            filename = img_path,
            width = img_dim[1],
            height = img_dim[2],
            units = "in",
            res = 300
         )
         print(get(img_name))
         dev.off()
      } else {
         ggplot2::ggsave(
            filename = img_path,
            plot = get(img_name),
            width = img_dim[1],
            height = img_dim[2],
            units = "in",
            dpi = 300,
            device = "png"
         )
      }
   }


   qmd_template <- glue::glue(
      '---
title: "{report_title}"
monofont: "Fira Code"
format:
  html:
    # css: "render_support/custom_styles.css"
    toc: true
    toc-depth: 10
    toc-location: left
    toc-expand: true
    number-sections: false
    pandoc-args: [
        "--lua-filter", "render_support/foldableCodeBlock.lua"
    ]
    self-contained: {tolower(as.character(self_contained))}
    self-contained-math: true
    code-fold: true
    embed-resources: {tolower(as.character(self_contained))}
---'
   )

   rmd_content <- qmd_template

   for (img_name in names(img_dims)) {
      img_path <- file.path(temp_img_dir, paste0(img_name, ".png"))
      header_text <- tools::toTitleCase(gsub("_", " ", img_name))

      rmd_content <- paste0(
         rmd_content,
         "\n\n## ", header_text, "\n\n",
         "```{r ", img_name, ", echo=FALSE}\n",
         "knitr::include_graphics('", img_path, "')\n",
         "```\n"
      )
   }

   writeLines(rmd_content, report_output_path)

   # Check if crochet is installed
   if (!requireNamespace("ADDINS.PNNL", quietly = TRUE)) {
      stop("The 'ADDINS.PNNL' package is required to render the QC report.")
   }
   if (!"crochet" %in% ls("package:ADDINS.PNNL")) {
      stop(glue(
         "The 'ADDINS.PNNL' package is installed but does not contain the ",
         "'crochet' function. As of last, this should be either on the ",
         "`feature/crochet-notebook-renderer` branch or perhaps `main` branch (in the future)."
      ))
   }

   ADDINS.PNNL::crochet(
      abs_input_files = report_output_path,
      abs_output_dir = dirname(report_output_path)
   )
}
