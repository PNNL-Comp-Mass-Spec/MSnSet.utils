#' @export
.do_filter <- function(m, missingness_thresh_gt, min_plex_thresh_gt) {
    m <- m[rowMeans(!is.na(exprs(m))) >= missingness_thresh_gt, ]
    design_mat <- model.matrix(~ 0 + m$PlexID)
    keep <- rowSums(!is.na(exprs(m)) %*% design_mat > 0L) >= min_plex_thresh_gt
    m <- m[keep, ]
    m
}

#' @export
.do_normalize <- function(m) {
    m <- MSnbase::MSnSet(
        medpolish(exprs(m),
            eps = .Machine$double.eps,
            maxiter = 1000,
            na.rm = T,
            trace.iter = F
        )$residuals,
        pData = pData(m),
        fData = fData(m)
    )
    m
}

#' @export
.do_batch_correction <- function(m) {
    if (is.null(m$PlexID)) {
        stop("Batch correction requires the exact column name 'PlexID' in pData(m) to correct batch effect.")
    }
    m <- suppressMessages(correct_batch_effect_NA(
        m,
        batch_name = "PlexID"
    ))
    m
}

#' @export
.raw <- function(m) {
    m
}

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

export_msnset_to_rds <- function(m, filename) {
    saveRDS(
        object = m,
        file = filename
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
  funcs_to_apply_and_args = list(
      .do_filter = list(
          args = list(missingness_thresh_gt = 0.3, min_plex_thresh_gt = 3),
          exporting = c("csv", "msnset")
      ),
      .do_normalize = list(args = list(), exporting = c()),
      .do_batch_correction = list(args = list(), exporting = c("export_exprs", "export_msnset_to_rds"))
  ),
  out_dir = ".",
  dry_run = FALSE
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
        other <- gsub("(.*)_.*?_level\\.RDS$", "\\1", file)
        message(glue::glue("Processing {level} level ---"))
        running_m <- readRDS(file)
        for (func_name in names(funcs_to_apply_and_args)) {
            func <- get(func_name)
            message(glue::glue("Applying {func_name} --"))
            args <- c(list(m = running_m), funcs_to_apply_and_args[[func_name]]$args)
            if (!dry_run) {
                running_m <- do.call(func, args)
            }
            for (export_func_str in funcs_to_apply_and_args[[func_name]]$exporting) {
                export_func <- get(export_func_str)
                filename_export_func_str <- gsub("\\.", "", func_name)
                message(glue::glue("Exporting with {export_func_str} -"))
                export_args <- list(m = running_m)
                if (export_func_str == "export_exprs") {
                    export_args$filename <- glue::glue(
                        "{other}_{level}_level_after_{filename_export_func_str}.csv"
                    )
                } else if (export_func_str == "export_msnset_to_rds") {
                    export_args$filename <- glue::glue(
                        "{other}_{level}_level_after_{filename_export_func_str}.RDS"
                    )
                }
                if (!dry_run) {
                    do.call(export_func, export_args)
                }
            }
        }
    }
}
