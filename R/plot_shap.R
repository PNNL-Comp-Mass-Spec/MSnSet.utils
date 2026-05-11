#' Plot SHAP values for Random Forest model
#'
#' This function calculates and visualizes SHAP (SHapley Additive exPlanations) values
#' for a *Random Forest* model trained on an MSnSet object. It shows the impact of
#' each feature on the model's prediction of a specific class, along with the
#' underlying protein abundance.
#'
#' @param msnset MSnSet object containing the data for analysis.
#' @param features_or_modeling_out Either the output of a call to functions like
#'    \code{rf_modeling} (which contains an element \code{top}) OR
#'    a character vector of feature names to be used for SHAP value calculation.
#'    The former will be assumed if \code{features_or_modeling_out} has names
#'     "prob", "features", "top", "auc", "pred".
#' @param response Character. The name of the response variable in \code{pData(msnset)}.
#' @param pred.cls Character. The level of the response variable that corresponds to the positive class (the one being predicted).
#' @param interpretation_and_direction Integer. Either 1 or -1. Direction of interpretation.
#'    Default is 1. If 1, the interpretation is that feature abundances are indeed abundances.
#'    If -1, the interpretation is that feature abundances represent changes in abundance
#'    (e.g., log2 fold changes between conditions). This is for historical purposes. If this is not desired,
#'    you may change axis labels by modifying the returned ggplot object.
#'
#' @return A ggplot object showing SHAP values and protein abundances.
#'
#' @references
#'   SHAP paper:
#'   doi = 10.5555/3295222.3295230,
#'   author = Lundberg, Scott M. and Lee, Su-In,
#'   title = A unified approach to interpreting model predictions,
#'   year = 2017,
#'   isbn = 9781510860964,
#'   publisher = Curran Associates Inc.,
#'   address = Red Hook, NY, USA,
#'   abstract = Understanding why a model makes a certain prediction can be as crucial as the prediction's accuracy in many applications. However, the highest accuracy for large modern datasets is often achieved by complex models that even experts struggle to interpret, such as ensemble or deep learning models, creating a tension between accuracy and interpretability. In response, various methods have recently been proposed to help users interpret the predictions of complex models, but it is often unclear how these methods are related and when one method is preferable over another. To address this problem, we present a unified framework for interpreting predictions, SHAP (SHapley Additive exPlanations). SHAP assigns each feature an importance value for a particular prediction. Its novel components include: (1) the identification of a new class of additive feature importance measures, and (2) theoretical results showing there is a unique solution in this class with a set of desirable properties. The new class unifies six existing methods, notable because several recent methods in the class lack the proposed desirable properties. Based on insights from this unification, we present new methods that show improved computational performance and/or better consistency with human intuition than previous approaches.,
#'   booktitle = Proceedings of the 31st International Conference on Neural Information Processing Systems,
#'   pages = 4768–4777,
#'   numpages = 10,
#'   location = Long Beach, California, USA,
#'   series = NIPS'17


#'
#' @importFrom randomForest randomForest
#' @importFrom treeshap unify treeshap
#' @importFrom shapviz shapviz get_shap_values
#' @importFrom MSnbase exprs pData
#' @importFrom ggplot2 ggplot aes geom_jitter scale_fill_gradientn scale_shape_manual geom_vline scale_y_continuous xlab ggtitle expansion element_text theme
#' @importFrom dplyr %>% group_by summarise ungroup arrange left_join pull mutate filter
#' @importFrom glue glue
#' @importFrom stats setNames quantile sd median
#' @importFrom scales rescale label_number squish trans_new
#'
#' @export
plot_treeshap <- function(msnset, features_or_modeling_out, response, pred.cls, interpretation_and_direction = 1) {
    . <- Feature <- SHAP_Value <- mean_shap <- Response <- NULL
    Protein_Abundance <- y_feat <- mean_prot_abundance <- NULL
    if (is.list(features_or_modeling_out) &&
        all(c("prob", "features", "top", "auc", "pred") %in% names(features_or_modeling_out))
    ) {
        feature_subset <- names(features_or_modeling_out$top)
    } else if (is.character(features_or_modeling_out)) {
        feature_subset <- features_or_modeling_out
    } else {
        stop("features_or_modeling_out must be either a list output from rf_modeling or a character vector of feature names.")
    }

    if ("top" %in% feature_subset) {
        warning(glue(
            "It seems that 'top' is included in the internal variable feature_subset. ",
            "It's possible this is correct, but more likely the output ",
            "from a function like rf_modeling, that you passed to features_or_modeling_out, was malformed. ",
            "See how features_or_modeling_out is expected to be structured in documentation. ",
            sep = ""
        ))
    }

    if (length(feature_subset) == 0) {
        stop("No features provided for SHAP value calculation.")
    } else if (length(feature_subset) == 1) {
        stop("`treeshap::` requires at least two features for proper SHAP value calculation.")
    }

    if (!interpretation_and_direction %in% c(1, -1)) {
        stop("interpretation_and_direction must be either 1 or -1")
    }

    if (!response %in% colnames(pData(msnset))) {
        stop(glue("Response variable '{response}' not found in pData of the MSnSet."))
    }

    resp_vec <- pData(msnset)[[response]]
    if (!is.factor(resp_vec)) {
        resp_vec <- as.factor(resp_vec)
    }

    resp_levels <- levels(resp_vec)
    if (length(resp_levels) != 2) {
        stop(glue("Response variable '{response}' must have exactly two levels. Found: {paste(resp_levels, collapse = ', ')}"))
    }

    if (!pred.cls %in% resp_levels) {
        stop(glue("pred.cls '{pred.cls}' not found in levels of response variable '{response}'."))
    }

    ctrl.cls <- resp_levels[resp_levels != pred.cls]

    X <- msnset %>%
        exprs() %>%
        t() %>%
        as.data.frame() %>%
        .[, feature_subset, drop = FALSE]

    if (any(is.na(X))) {
        stop("The feature matrix contains missing values (NAs); `treeshap::` requires a complete dataset.")
    }

    feat_to_dir <- sapply(
        feature_subset,
        FUN = function(f) {
            print(f)
            sign(mean(X[pData(msnset)[[response]] == pred.cls, f]) -
                mean(X[pData(msnset)[[response]] == ctrl.cls, f]))
        }
    )

    y <- msnset %>%
        pData() %>%
        .[[response]] %>%
        factor(levels = c(ctrl.cls, pred.cls), labels = c("0", "1")) %>%
        as.character() %>%
        as.factor()
    X <- X * interpretation_and_direction
    model <- randomForest(x = X, y = y)
    message("Modeled with ", nrow(X), " samples and ", ncol(X), " feature(s).")
    model_unified <- unify(model, X)
    shap_values_obj <- treeshap(model_unified, X, verbose = TRUE)
    sv <- shapviz(shap_values_obj)

    act_sv <- get_shap_values(sv)

    combi_plot_df <- data.frame(
        Feature = rep(colnames(act_sv), each = nrow(act_sv)),
        SHAP_Value = as.vector(act_sv),
        Sample = rep(rownames(act_sv), times = ncol(act_sv)),
        Response = rep(y, times = ncol(act_sv)),
        Protein_Abundance = as.vector(as.matrix(X))
    )

    most_influential <- combi_plot_df %>%
        group_by(Feature, Response) %>%
        summarise(mean_prot_abundance = mean(Protein_Abundance)) %>%
        ungroup() %>%
        arrange(Feature, as.numeric(Response)) %>%
        group_by(Feature) %>%
        summarise(case_ctrl_diff = diff(mean_prot_abundance)) %>%
        left_join(
            combi_plot_df %>%
                group_by(Feature) %>%
                summarise(mean_shap = mean(abs(SHAP_Value))),
            by = "Feature"
        ) %>%
        arrange(desc(abs(mean_shap))) %>%
        pull("Feature")

    all_ordered_features <- rev(most_influential)

    combi_slim_df <- combi_plot_df %>%
        dplyr::filter(Feature %in% all_ordered_features) %>%
        mutate(Feature = factor(Feature, levels = all_ordered_features))

    y_to_feature <- levels(combi_slim_df$Feature)
    y_to_feature_show <- sapply(
        y_to_feature,
        FUN = function(f) {
            if (!f %in% names(feat_to_dir)) {
                stop(glue("Feature {f} not found in feat_to_dir"))
            }
            dir_sign <- feat_to_dir[f]
            if (dir_sign * interpretation_and_direction == 1) {
                paste0(f, " ( + )")
            } else if (dir_sign * interpretation_and_direction == -1) {
                paste0(f, " ( - )")
            } else {
                stop()
            }
        }
    )
    feature_to_y <- setNames(seq_along(y_to_feature), y_to_feature) %>% unlist()
    y_to_feature <- unname(y_to_feature)
    combi_slim_df$y_feat <- feature_to_y[as.character(combi_slim_df$Feature)]

    pctl_low <- quantile(combi_slim_df$Protein_Abundance, probs = 0.025)
    pctl_high <- quantile(combi_slim_df$Protein_Abundance, probs = 0.975)

    plt <- ggplot(
        data = combi_slim_df,
        mapping = aes(fill = Protein_Abundance, y = y_feat, x = SHAP_Value, shape = Response)
    ) +
        geom_jitter(alpha = 0.8, size = 1.5, stroke = .25, color = "black", height = 0.1, width = 0) +
        scale_fill_gradientn(
            colours = c("#0000aa", "white", "#aa0000"),
            values = rescale(c(pctl_low, 0, pctl_high)),
            name = ifelse(interpretation_and_direction == 1,
                "log2(Protein\nAbundance)",
                "log2(Change in\nProtein\nAbundance)"
            ),
            limits = c(pctl_low, pctl_high),
            oob = squish
        ) +
        scale_shape_manual(
            values = c(21, 24),
            labels = c(ctrl.cls, pred.cls),
            name = response
        ) +
        geom_vline(
            mapping = aes(xintercept = 0),
            linetype = "dotted",
            color = "grey40"
        ) +
        scale_y_continuous(
            breaks = seq_along(y_to_feature),
            labels = unname(y_to_feature_show),
            expand = expansion(add = c(0.25, 0.25)),
            name = ifelse(interpretation_and_direction == 1,
                "Protein (+/- abundance)",
                "Protein (+/- change in abundance)"
            )
        ) +
        xlab(glue("SHAP value: impact on model predicting \"{pred.cls}\" status")) +
        ggtitle(glue("SHAP and Abundance Values for Top Proteins"))

    return(plt)
}
