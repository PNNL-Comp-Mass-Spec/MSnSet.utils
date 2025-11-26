#' Export Expression Data from MSnSet
#'
#' Exports the expression data (\code{exprs(m)}) from an \code{MSnSet} object to a CSV file.
#' Sample names are cleaned (removing leading "X" and replacing "_" with " ") and
#' columns can be reordered based on \code{pData} variables.
#'
#' @param m An \code{MSnSet} object.
#' @param filename Character string. The path to the output CSV file.
#' @param arrange_by Character vector. Variable names in \code{pData(m)} to sort the columns by.
#'
#' @return No return value, called for side effect of writing a file.
#' @export
export_exprs <- function(m, filename, arrange_by = c()) {
    data_wd()
    new_colnames <- gsub("^X", "", sampleNames(m))
    new_colnames <- gsub("_", " ", new_colnames)
    Biobase::sampleNames(m) <- new_colnames

    new_colname_order <- pData(m) %>%
        arrange(!!!syms(arrange_by), rownames(.)) %>%
        rownames()

    new_row_order <- sort(rownames(m))

    m <- m[new_row_order, new_colname_order]
    data <- exprs(m)
    write.csv(
        data,
        file = filename,
        row.names = TRUE,
        quote = FALSE
    )
}

#' Process and Export Data Levels
#'
#' Iterates through a set of raw data files, performs filtering, normalization, and batch correction,
#' and exports the results.
#'
#' @param m_raw_level_files Character vector. Paths to RDS files containing raw \code{MSnSet} objects.
#'   \strong{Caveat:} Filenames must match the regex pattern \code{.*_(.*?)_level\\.RDS$}.
#' @param out_dir Character string. Output directory path. Defaults to current directory.
#' @param missingness_thresh_gt Numeric. Threshold for missingness filtering (proportion of non-NA values). Default is 0.3.
#' @param min_plex_thresh_gt Numeric. Minimum number of plexes a feature must be observed in. Default is 3.
#'
#' @return No return value, called for side effects (creating directories and files).
#' @export
export_data_levels <- function(
  m_raw_level_files,
  out_dir = ".",
  missingness_thresh_gt = 0.3,
  min_plex_thresh_gt = 3
) {
    for (file in m_raw_level_files) {
        if (!grepl(".*_(.*?)_level\\.RDS$", file)) {
            stop(
                glue(
                    "File {file} does not match expected regex pattern: ",
                    ".*_(.*?)_level\\.RDS$"
                )
            )
        }
        if (!dir.exists(out_dir)) {
            dir.create(out_dir, recursive = TRUE)
            message(glue::glue("Created output directory {out_dir}"))
        }
        setwd(out_dir)
        level <- gsub(".*_(.*?)_level\\.RDS$", "\\1", file)
        message(glue::glue("Processing {level} level"))

        m_raw_level <- readRDS(file)
        # m_raw_level$NCI7 <- grepl("NCI", sampleNames(m_raw_level))
        saveRDS(m_raw_level, file = glue::glue("levels/m_raw_{level}.RDS"))
        export_exprs(
            m_raw_level,
            filename = glue::glue("raw_global_tmt_{level}_level.csv")
        )


        m_nafilt_level <- m_raw_level[
            rowMeans(!is.na(exprs(m_raw_level))) >= missingness_thresh_gt,
        ]
        design_mat <- model.matrix(~ 0 + m_nafilt_level$PlexID)
        keep <- rowSums(
            !is.na(exprs(m_nafilt_level)) %*% design_mat > 0
        ) >= min_plex_thresh_gt
        m_nafilt_level <- m_nafilt_level[keep, ]
        message(glue::glue("Normalizing {level}"))
        m_nafilt_norm_level <- MSnbase::MSnSet(
            medpolish(exprs(m_nafilt_level),
                eps = .Machine$double.eps,
                maxiter = 1000,
                na.rm = T,
                trace.iter = F
            )$residuals,
            pData = pData(m_nafilt_level),
            fData = fData(m_nafilt_level)
        )
        message(glue::glue("Correcting batch effect {level}"))
        m_nafilt_norm_batchcorr_level <- suppressMessages(correct_batch_effect_NA(
            m_nafilt_norm_level,
            batch_name = "PlexID"
        ))

        export_exprs(
            m_nafilt_norm_batchcorr_level,
            filename = glue::glue(
                "global_tmt_filter_{THRESH*100}_minPlexFilter_",
                "{COUNT_THRESH}_normalized_{level}_level.csv"
            )
        )
        saveRDS(
            m_nafilt_norm_batchcorr_level,
            file = glue::glue("levels/m_final_{level}.RDS")
        )
    }
}
