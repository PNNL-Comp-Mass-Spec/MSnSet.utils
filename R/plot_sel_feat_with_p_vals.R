adj.P.Val <- NULL # (For lintr)


#' @title MSnSet Feature Bar Chart with (Adjusted) P-Value Coloring
#'
#' @export plot_sel_feat_with_p_vals
#'
#' @param top_selected A named vector with names as top Boruta-selected features and
#'     values as the frequency at which said features were selected
#' @param msnset A \code{MSnSet} object that represented the input to an `rf_modeling`
#'     call (with Boruta turned on). This is required so P-values can be computed
#' @param response_colname The column name (in \code{pData(msnset)}) of the machine
#'     learning response variable
#' @param feat_name_colname If provided, a column in \code{fData(msnset)} by which
#'     to name the features by. If \code{NULL}, then
#'     \code{featureNames(msnset) will be used.}
#' @param highlight_feats If provided, put a star (U+2605) next to the label of these
#'     features. Note: this uses \code{feat_name_colnames} if passed, and not
#'     \code{msnset}'s native features.
#' @param highlight_reason Short string describing why certain features would be
#'     highlighted. Only applicable if \code{highlight_feats} is not \code{NULL}.
#' @param alpha The alpha-level that defines statistical significance (for purposes
#'     of the color scale)
#' @param title The title of the returned plot
plot_sel_feat_with_p_vals <- function(
    top_selected,
    msnset,
    response_colname,
    feat_name_colname = NULL,
    highlight_feats = NULL,
    highlight_reason = NULL,
    alpha = 0.05,
    title = "Frequency of\nBoruta-Selected Features") {
    lab <- as.vector(names(top_selected))
    cnt <- as.vector(unlist(unname(top_selected)))
    if (!is.null(feat_name_colname)) {
        display_lab <- Biobase::fData(msnset)[lab, ] %>% dplyr::pull(feat_name_colname)
    } else {
        display_lab <- lab
    }
    if (!is.null(highlight_feats)) {
        display_lab <- unlist(unname(
            lapply(
                display_lab,
                function(dl) ifelse(dl %in% highlight_feats, paste("\U2605 ", dl, sep = ""), dl)
            )
        ))
    }

    df_from_counts <- data.frame(
        lab = lab,
        display_lab = display_lab,
        cnt = cnt
    )

    limma <- .get_limma(top_selected, msnset, alpha, response_colname)
    df <- merge(df_from_counts, limma, by.x = "lab", by.y = 0)
    labels_breaks <- sort(append(
        c(1, alpha),
        unlist(lapply(seq(-2, -10, -2), function(x) 10^x))
    ), decreasing = TRUE)

    labels_breaks <- labels_breaks[labels_breaks <= alpha | labels_breaks == 1]
    hard_cut <- scales::rescale(x = log10(labels_breaks), to = c(0, 1))[2]
    if (!is.null(highlight_reason) && !is.null(highlight_feats)) {
        title_if_else <- labs(
            title = title,
            subtitle = glue::glue("\U2605 = {highlight_reason}")
        )
    } else {
        title_if_else <- labs(
            title = title
        )
    }

    colors <- viridis::viridis(10)
    p <- ggplot(data = df) +
        ggplot2::geom_col(
            mapping = ggplot2::aes(
                x = forcats::fct_reorder(display_lab, cnt, .desc = TRUE),
                y = cnt,
                fill = adj.P.Val
            ),
            width = 2 / 3
        ) +
        ggplot2::xlab("Feature Name") +
        ggplot2::ylab("Counts") +
        ggplot2::geom_hline(mapping = ggplot2::aes(
            yintercept = ncol(msnset),
            linetype = "Total Number of\nSubjects"
        )) +
        ggplot2::scale_linetype_manual(values = c("dashed"), name = "Legend") +
        title_if_else +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1),
            plot.background = NULL,
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            legend.key.height = ggplot2::unit(2.5, "char"),
            legend.key.width = ggplot2::unit(1, "char")
        ) +
        ggplot2::scale_fill_gradientn(
            name = "Adjusted\nP-value",
            colors = c(colors, "white", "white"),
            trans = "log",
            limits = c(labels_breaks[length(labels_breaks)], labels_breaks[1]),
            labels = append(
                labels_breaks[1:(length(labels_breaks) - 1)],
                c(sprintf("\U2264 %s", labels_breaks[length(labels_breaks)]))
            ),
            breaks = labels_breaks,
            oob = scales::squish,
            values = sort(append(
                seq(0, hard_cut, hard_cut / 10),
                c(hard_cut - 1e-6, 1)
            )),
            guide = ggplot2::guide_colorbar(
                frame.colour = "black", ticks.colour = "black"
            )
        )
    return(p)
}

.get_limma <- function(top_selected, msnset, alpha, response) {
    l <- MSnSet.utils::limma_a_b(
        msnset,
        model.str = glue::glue("~ {response}"),
        coef.str = response
    )
    ll <- l %>%
        merge(
            enframe(top_selected) %>%
                tibble::column_to_rownames("name"),
            by = 0, all = TRUE
        ) %>%
        tibble::column_to_rownames("Row.names") %>%
        mutate(pctg_subj_selected = value /
            (msnset %>%
                Biobase::sampleNames() %>%
                length()
            )) %>%
        select(-c("value")) %>%
        mutate(is.signif = adj.P.Val < alpha) %>%
        arrange(adj.P.Val)
    return(ll)
}
