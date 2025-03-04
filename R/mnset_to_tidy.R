#' Convert a proteomics MSnSet object to tidy data format
#'
#' This function converts a proteomics MSnSet object to a tidy data format,
#' optionally including feature data (fData) and/or phenotype data (pData).
#'
#' @param m An MSnSet object to convert to tidy format
#' @param include_fData Logical; whether to include feature data in the tidy object (default: FALSE)
#' @param include_pData Logical; whether to include phenotype data in the tidy object (default: TRUE)
#' @param intensity_colname Character; name for the intensity column in the tidy output (default: "intensity")
#' @param fData_id Character; column name to use for fData rownames (default: "protein_group")
#' @param pData_id Character; column name to use for pData rownames (default: "study_id")
#' @param quiet Logical; whether to suppress warnings and messages (default: FALSE)
#' @param remove_na Logical; whether to remove NA values from expression data (default: TRUE)
#'
#' @return A tibble in tidy format containing the expression data and optionally fData and/or pData
#'
#' @importFrom dplyr filter left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom rlang .data
#' @importFrom MSnbase exprs pData fData
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' tidy_data <- msnset_to_tidy(msnset_object)
#'
#' # Include feature data
#' tidy_data <- msnset_to_tidy(msnset_object, include_fData = TRUE)
#'
#' # Include feature data but exclude phenotype data
#' tidy_data <- msnset_to_tidy(msnset_object, include_fData = TRUE, include_pData = FALSE)
#'
#' # Custom column naming
#' tidy_data <- msnset_to_tidy(msnset_object,
#'                             intensity_colname = "abundance",
#'                             fData_id = "peptide_id",
#'                             pData_id = "patient_id")
#'
#' # Suppress all warnings and messages
#' tidy_data <- msnset_to_tidy(msnset_object, quiet = TRUE)
#'
#' # Keep NA values in the output
#' tidy_data <- msnset_to_tidy(msnset_object, remove_na = FALSE)
#' }
#'
#' @export
msnset_to_tidy <- function(m,
                           include_fData = FALSE,
                           include_pData = TRUE,
                           intensity_colname = "intensity",
                           fData_id = "protein_group",
                           pData_id = "study_id",
                           quiet = FALSE,
                           remove_na = TRUE) {

   # Check if input is an MSnSet object
   if (!inherits(m, "MSnSet")) {
      stop("Input must be an MSnSet object")
   }

   # Message about dimensions
   if (!quiet) {
      message("Converting MSnSet with ", nrow(exprs(m)), " features and ",
              ncol(exprs(m)), " samples to tidy format")
   }

   # Check if there are any NAs in expression data and warn if many
   na_count <- sum(is.na(exprs(m)))
   if (na_count > 0 && !quiet) {
      if (remove_na) {
         warning("Expression data contains ", na_count, " NA values (",
                 round(na_count / (nrow(exprs(m)) * ncol(exprs(m))) * 100, 2),
                 "%) which will be filtered out")
      } else {
         message("Expression data contains ", na_count, " NA values (",
                 round(na_count / (nrow(exprs(m)) * ncol(exprs(m))) * 100, 2),
                 "%) which will be retained")
      }
   }

   # Expression data
   eData <- as.data.frame(exprs(m)) |>
      rownames_to_column(var = fData_id) |>
      pivot_longer(cols = -all_of(fData_id),
                   names_to = "sample_id",
                   values_to = intensity_colname)

   # Filter NAs if requested
   if (remove_na) {
      eData <- eData |> filter(!is.na(.data[[intensity_colname]]))
   }

   # Initialize data as eData
   data <- eData

   # Add phenotype data if requested
   if (include_pData) {
      if (ncol(pData(m)) == 0) {
         if (!quiet) warning("pData is empty, no phenotype data to include")
      } else {
         pData <- pData(m) |>
            rownames_to_column(var = "sample_id")

         # Check if sample IDs match
         missing_samples <- setdiff(unique(eData$sample_id), pData$sample_id)
         if (length(missing_samples) > 0 && !quiet) {
            warning("Some samples in expression data are missing from pData (",
                    length(missing_samples), " samples)")
         }

         # Combine with phenotype data
         data <- left_join(data, pData, by = "sample_id")

         # Rename sample_id to pData_id if they differ
         if (pData_id != "sample_id" && !quiet) {
            message("Sample IDs kept as 'sample_id' for joining. Use rename() if you need to change this.")
         }
      }
   }

   # Add feature data if requested
   if (include_fData) {
      if (ncol(fData(m)) == 0) {
         if (!quiet) warning("fData is empty, no feature data to include")
      } else {
         fData <- fData(m) |>
            rownames_to_column(var = fData_id)

         # Check if feature IDs match
         missing_features <- setdiff(unique(eData[[fData_id]]), fData[[fData_id]])
         if (length(missing_features) > 0 && !quiet) {
            warning("Some features in expression data are missing from fData (",
                    length(missing_features), " features)")
         }

         # Combine with feature data
         data <- left_join(data, fData, by = fData_id)
      }
   }

   # Return the tidy data
   return(data)
}