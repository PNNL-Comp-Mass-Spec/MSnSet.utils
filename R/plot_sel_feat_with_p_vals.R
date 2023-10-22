adj.P.Val <- NULL # (For lintr)


#' @export plot_sel_feat_with_p_vals
#'
#' @param top_selected A named vector with names as top Boruta-selected features and \
#'     values as the frequency at which said features were selected
#' @param msnset A \code{MSnSet} object that represented the input to an `rf_modeling` \
#'     call (with Boruta turned on). This is required so P-values can be computed
#' @param response_colname The column name (in \code{pData(msnset)}) of the machine \
#'     learning response variable
#' @param feat_name_colname If provided, a column in \code{fData(msnset)} by which \
#'     to name the features by. If \code{NULL}, then \
#'     \code{featureNames(msnset) will be used.}
#' @param highlight_feats If provided, put an asterisk ❉ next to the label of these \
#'     features. Note: this uses \code{feat_name_colnames} if passed, and not \
#'     \code{msnset}'s native features.
#' @param highlight_reason Short string describing why certain features would be \
#'     highlighted. Only applicable if \code{highlight_feats} is not \code{NULL}.
#' @param alpha The alpha-level that defines statistical significance (for purposes \
#'     of the color scale)
plot_sel_feat_with_p_vals <- function(
    top_selected,
    msnset,
    response_colname,
    feat_name_colname = NULL,
    highlight_feats = NULL,
    highlight_reason = NULL,
    alpha = 0.05) {
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
                function(dl) ifelse(dl %in% highlight_feats, paste("** ", dl, sep = ""), dl)
            )
        ))
    }

    df_from_counts <- data.frame(
        lab = lab,
        display_lab = display_lab,
        cnt = cnt
    )

    limma <- MSnSet.utils:::.get_limma(top_selected, msnset, alpha, response_colname)
    df <- merge(df_from_counts, limma, by.x = "lab", by.y = 0)
    labels_breaks <- append(
        c(1, alpha),
        unlist(lapply(seq(-2, -10, -2), function(x) 10^x))
    )
    hard_cut <- scales::rescale(x = log10(labels_breaks), to = c(0, 1))[2]
    if (!is.null(highlight_reason) && !is.null(highlight_feats)) {
        title_if_else <- labs(
            title = "Frequency of\nBoruta-Selected Features",
            subtitle = glue::glue("❉ = {highlight_reason}")
        )
    } else {
        title_if_else <- labs(
            title = "Frequency of\nBoruta-Selected Features"
        )
    }

    p <- ggplot(data = df) +
        geom_col(
            mapping = aes(
                x = forcats::fct_reorder(display_lab, cnt, .desc = TRUE),
                y = cnt,
                fill = adj.P.Val
            ),
            width = 0.66
        ) +
        xlab("Feature Name") +
        ylab("Counts") +
        title_if_else +
        ggplot2::theme(
            axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1),
            # aspect.ratio = 1 / 2,
            plot.background = NULL,
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.key.height = unit(2.5, "char"),
            legend.key.width = unit(1, "char")
        ) +
        ggplot2::scale_fill_gradientn(
            name = "Adjusted\nP-value",
            colors = append(viridis::viridis(5), c("white", "white")),
            trans = "log",
            limits = c(labels_breaks[length(labels_breaks)], labels_breaks[1]),
            labels = append(
                labels_breaks[1:(length(labels_breaks) - 1)],
                c(sprintf("≤ %s", labels_breaks[length(labels_breaks)]))
            ),
            breaks = labels_breaks,
            oob = scales::squish,
            values = append(
                seq(0, hard_cut, hard_cut / 5),
                c(hard_cut - 1e-10, hard_cut, 1)
            ),
            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
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