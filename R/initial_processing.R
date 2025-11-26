#' Copy Raw Files from DMS Data Package
#'
#' Downloads raw mass spectrometry files from a DMS data package to a local
#' directory structure.
#'
#' @param data_package_num Integer. The DMS data package number.
#' @param base_path Character. Base directory path for storing files.
#' @param overwrite Logical. Whether to overwrite existing files. Default is FALSE.
#'
#' @return Invisible list containing paths to copied files and copy status.
#' @export
copy_raw_files_from_dms <- function(data_package_num, base_path, overwrite = FALSE) {
   if (!is.numeric(data_package_num) || data_package_num <= 0) {
      stop("'data_package_num' must be a positive integer.", call. = FALSE)
   }

   datasets <- tryCatch(
      PNNL.DMS.utils::get_datasets_by_data_package(data_package_num = data_package_num),
      error = function(e) {
         stop("Failed to retrieve datasets from DMS for data package ", data_package_num,
            ": ", conditionMessage(e),
            call. = FALSE
         )
      }
   )

   if (is.null(datasets) || nrow(datasets) == 0) {
      stop("No datasets found for data package ", data_package_num, call. = FALSE)
   }

   raw_files <- list.files(
      path = datasets$folder,
      pattern = "\\.raw$",
      recursive = FALSE,
      full.names = TRUE
   )

   if (length(raw_files) == 0) {
      warning("No .raw files found in DMS folders for data package ", data_package_num, call. = FALSE)
      return(invisible(list(files = character(0), status = logical(0))))
   }

   raw_dir <- file.path(base_path, "raw")

   if (!dir.exists(raw_dir)) {
      dir.create(raw_dir, recursive = TRUE)
   }

   raw_files_new <- file.path(raw_dir, basename(raw_files))

   if (!overwrite && any(file.exists(raw_files_new))) {
      existing <- raw_files_new[file.exists(raw_files_new)]
      warning("The following files already exist (use overwrite = TRUE to replace):\n  ",
         paste(basename(existing), collapse = "\n  "),
         call. = FALSE
      )
   }
   do_files_logical <- !file.exists(raw_files_new)
   raw_files_new <- raw_files_new[do_files_logical]
   raw_files <- raw_files[do_files_logical]
   stopifnot(length(raw_files) == length(raw_files_new))
   cli::cli_progress_bar(name = "Copying raw files from DMS", total = length(raw_files_new))
   all_copy_status <- logical(length(raw_files_new))
   for (i in seq_along(raw_files_new)) {
      cur_file_src <- raw_files[i]
      cur_file_dest <- raw_files_new[i]
      copy_status <- file.copy(
         from = cur_file_src,
         to = cur_file_dest,
         overwrite = overwrite
      )
      cli::cli_progress_update()
   }
   cli::cli_progress_done()
   all_copy_status[do_files_logical] <- copy_status
   copy_status <- all_copy_status
   if (!all(copy_status)) {
      failed <- raw_files[!copy_status]
      warning("Failed to copy the following files:\n  ",
         paste(basename(failed), collapse = "\n  "),
         call. = FALSE
      )
   }

   invisible(list(files = raw_files_new, status = copy_status))
}




#' Run msConvert CLI to Convert Raw Files to mzML
#'
#' Converts Thermo .raw files to .mzML format using the msConvert command-line tool.
#'
#' @param raw_files_dir Character. Directory containing .raw files to convert.
#' @param output_dir Character. Directory to save converted .mzML files.
#' @param msconvert_path Character. Path to msconvert executable. If NULL, will use msconvert from system PATH.
#' @param addl_cli_args Character vector. Additional command-line arguments for msconvert. Default includes zlib compression and peak picking.
#' @param overwrite Logical. Whether to overwrite existing files. Default is FALSE.
#'
#' @return Invisible list of msconvert command outputs for each file.
#' @export
#'
#' @examples
#' \dontrun{
#' run_msConvert_cli(
#'    raw_files_dir = "E:/Data_E/Tyler/5306_global/raw",
#'    output_dir = "E:/Data_E/Tyler/5306_global/mzML"
#' )
#' }
run_msConvert_cli <- function(
    raw_files_dir,
    output_dir,
    msconvert_path = NULL,
    addl_cli_args = c(
       "--zlib",
       "--filter",
       "\"peakPicking vendor msLevel=1-\"",
       "--filter",
       "\"titleMaker",
       "<RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>",
       "File:\"\"\"^<SourcePath^>\"\"\", NativeID:\"\"\"^<Id^>\"\"\"\""
    ),
    overwrite = FALSE) {
   if (is.null(msconvert_path)) {
      msconvert_path <- Sys.which("msconvert")
      if (msconvert_path == "") {
         stop("msconvert executable not found in system PATH. ",
            "Please provide the full path to msconvert.exe via the 'msconvert_path' argument.",
            call. = FALSE
         )
      }
   }

   if (!dir.exists(raw_files_dir)) {
      dir.create(raw_files_dir, recursive = TRUE)
   }

   input_files <- list.files(
      path = raw_files_dir,
      pattern = "\\.raw$",
      full.names = TRUE,
      recursive = FALSE
   )

   # Assemble command using glue
   addl_cli_args_str <- paste(addl_cli_args, collapse = " ")

   if (length(input_files) == 0) {
      stop("No .raw files found in the specified directory: ", raw_files_dir, call. = FALSE)
   }

   if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
   }

   cli::cli_progress_bar(
      name = "Raw -> mzML",
      total = length(input_files)
   )
   options(cli.progress_bar_style = "squares")
   all_outputs <- list()
   for (i in seq_along(input_files)) {
      input_file <- input_files[i]
      input_file_escaped <- gsub(" ", "\\ ", input_file)
      output_dir_escaped <- gsub(" ", "\\ ", output_dir)

      cmd <- glue::glue(
         "\"{msconvert_path}\" {addl_cli_args_str} \"{input_file_escaped}\" --outdir \"{output_dir_escaped}\""
      )
      cli::cli_progress_update()
      output <- system(cmd, intern = TRUE, show.output.on.console = TRUE)
      all_outputs[[basename(input_file)]] <- output
   }
   cli::cli_progress_done()
   saveRDS(output, file = file.path(output_dir, "msconvert_logs.RDS"))
}



#' Organize mzML Files by Plex
#'
#' Moves mzML files into plex-specific subdirectories based on filename regex.
#'
#' @param src_mzML_dir Character. Base directory containing the mzML folder.
#' @param dest_mzML_subdirs Character. Name of mzML subdirectory. Default is "plexes".
#' @param overwrite Logical. Whether to overwrite existing files. Default is FALSE.
#' @param sub_with_cap_grp_1 Character. Used to obtain plex ID through a
#'    call to `sub(sub_with_cap_grp_1, "\\1", <filename>)`.
#'    Default captures first two digits.
#'
#' @return Invisible list of new file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' organize_mzml_by_plex("E:/Data_E/Tyler/5306_global/")
#' }
organize_mzml_by_plex <- function(
    src_mzML_dir,
    dest_mzML_subdirs = "mzML/plexes/",
    overwrite = FALSE,
    sub_with_cap_grp_1 = "^(\\d{2}).*") {
   if (!dir.exists(src_mzML_dir)) {
      stop("mzML directory does not exist: ", src_mzML_dir, call. = FALSE)
   }

   mzml_files <- list.files(
      path = src_mzML_dir,
      pattern = "\\.mzML$",
      full.names = TRUE,
      recursive = FALSE
   )

   if (length(mzml_files) == 0) {
      warning("No mzML files found in ", src_mzML_dir, call. = FALSE)
      return(invisible(character(0)))
   }

   plex_prefix <- sub(sub_with_cap_grp_1, "\\1", basename(mzml_files))

   if (any(is.na(plex_prefix)) || any(!grepl(sub_with_cap_grp_1, plex_prefix))) {
      stop("Some mzML files do not have valid two-digit plex prefixes. ",
         "Expected format: '01_sample.mzML'",
         call. = FALSE
      )
   }

   plex <- paste0("Plex", as.numeric(plex_prefix))

   for (plex_i in unique(plex)) {
      plex_dir <- file.path(dest_mzML_subdirs, plex_i)
      if (!dir.exists(plex_dir)) {
         dir.create(plex_dir, recursive = TRUE)
      }
   }

   mzml_files_new <- file.path(dest_mzML_subdirs, plex, basename(mzml_files))

   if (!overwrite && any(file.exists(mzml_files_new))) {
      existing <- mzml_files_new[file.exists(mzml_files_new)]
      stop("The following destination files already exist:\n  ",
         paste(basename(existing), collapse = "\n  "),
         call. = FALSE
      )
   }

   move_status <- file.rename(from = mzml_files, to = mzml_files_new)

   if (!all(move_status)) {
      failed <- mzml_files[!move_status]
      stop("Failed to move the following files:\n  ",
         paste(basename(failed), collapse = "\n  "),
         call. = FALSE
      )
   }

   invisible(mzml_files_new)
}

#' Write Plex Annotation Files from Mapping
#'
#' For each plex subdirectory in mzML_plexes_dir, filter the mapping data.frame to matching plex,
#' remove the plex column, and write a tab-separated annotation file with no quotes, no column names, and no row names.
#'
#' @param mapping data.frame. Must contain columns 'plex', 'tmt', 'name'.
#' @param mzML_plexes_dir Character. Path to the directory containing plex subdirectories.
#' @importFrom dplyr `%>%`
#'
#' @return Invisible vector of written file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' write_plex_annotation_files(mapping, "E:/Data_E/Tyler/5306_global/mzML/plexes")
#' }
write_plex_annotation_files <- function(mapping, mzML_plexes_dir) {
   plex <- NULL # To appease lintr
   if (!is.data.frame(mapping)) {
      stop("'mapping' must be a data.frame.", call. = FALSE)
   }
   required_cols <- c("plex", "tmt", "name")
   if (!all(required_cols %in% colnames(mapping))) {
      stop("'mapping' must contain columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
   }
   if (length(required_cols) != ncol(mapping)) {
      warning("'mapping' contains extra columns beyond the required ones. Only 'plex', 'tmt', and 'name' will be used.", call. = FALSE)
   }
   if (!dir.exists(mzML_plexes_dir)) {
      stop("mzML_plexes_dir does not exist: ", mzML_plexes_dir, call. = FALSE)
   }
   plex_dirs <- list.dirs(mzML_plexes_dir, full.names = FALSE, recursive = FALSE)
   if (length(plex_dirs) == 0) {
      stop("No plex subdirectories found in ", mzML_plexes_dir, call. = FALSE)
   }
   written_files <- character(0)
   for (plex_dir in plex_dirs) {
      mapping_sub_dir <- mapping %>%
         dplyr::filter(plex == plex_dir) %>%
         dplyr::select(-plex)
      if (nrow(mapping_sub_dir) == 0) {
         stop("No mapping entries found for plex: \"", plex_dir, "\". Check for leading zeros in plex names.", call. = FALSE)
      }
      mapping_sub_dir$plex <- NULL
      out_file <- file.path(mzML_plexes_dir, plex_dir, paste0(plex_dir, "_annotation.txt"))
      write.table(
         mapping_sub_dir,
         file = out_file,
         sep = "\t",
         quote = FALSE,
         col.names = FALSE,
         row.names = FALSE
      )
      written_files <- c(written_files, out_file)
   }
   invisible(written_files)
}



#' Create MSnSet Object from Ratios and FASTA
#'
#' Generates an MSnSet object using ratio data and a FASTA file for feature data.
#'
#' @param path_to_ratios_or_ratios Character or data.frame. Path to a file containing ratios or a data.frame of ratios.
#' @param path_to_fasta_for_f_data Character. Path to the FASTA file used for feature data.
#' @param p_data_df data.frame. Metadata with row names matching sample names in ratios.
#'
#' @import dplyr
#' @importFrom Biostrings readAAStringSet
#' @importFrom tidyr separate
#' @importFrom fuzzyjoin regex_left_join
#'
#' @return An MSnSet object containing the processed data.
#' @export
create_msnset <- function(
    path_to_ratios_or_ratios,
    path_to_fasta_for_f_data,
    p_data_df,
    feature_name_col = "ensembl_protein",
    semicolon_warning = TRUE) {
   Index <- NULL
   MaxPepProb <- NULL
   feature_id <- NULL
   . <- NULL
   description <- NULL
   name <- NULL

   get_ratios <- function() {
      if (is.character(path_to_ratios_or_ratios)) {
         if (!file.exists(path_to_ratios_or_ratios)) {
            stop("Ratios file does not exist: ", path_to_ratios_or_ratios, call. = FALSE)
         }
         ratios <- read.delim(path_to_ratios_or_ratios, stringsAsFactors = FALSE, check.names = FALSE)
      } else if (is.data.frame(path_to_ratios_or_ratios)) {
         ratios <- path_to_ratios_or_ratios
      } else {
         stop("'path_to_ratios_or_ratios' must be either a file path or a data.frame.", call. = FALSE)
      }
      return(ratios)
   }

   ratios <- get_ratios()
   where_samples_start <- (which(colnames(ratios) == "ReferenceIntensity") + 1)

   get_f_data <- function() {
      fst <- readAAStringSet(path_to_fasta_for_f_data)
      check_flag <- FALSE
      withCallingHandlers(
         {
            feature_df <- data.frame(name = names(fst)) %>%
               filter(!grepl("^Cont", name)) %>%
               separate(
                  col = name,
                  into = c(
                     "ensembl_protein",
                     "ensembl_transcript",
                     "ensembl_gene",
                     "description"
                  ),
                  sep = "\\|",
                  remove = TRUE
               ) %>%
               dplyr::filter(!is.na(!!sym(feature_name_col))) %>%
               mutate(
                  gene = sub("(^[^ ]+) .*", "\\1", description),
                  description = sub("^[^ ]+ (.*$)", "\\1", description)
               ) %>%
               distinct(!!sym(feature_name_col), description,
                  .keep_all = TRUE
               )
         },
         warning = function(w) {
            if (grepl("Expected \\d+ pieces", conditionMessage(w))) {
               check_flag <<- TRUE
               invokeRestart("muffleWarning")
            } else {
               return()
            }
            if (check_flag) {
               missing_mappings <- setdiff(
                  unique(ratios$Index),
                  feature_df[[feature_name_col]]
               )
               if (length(missing_mappings) > 0) {
                  stop("The following IDs from the ratios data.frame could not ",
                     "be found in the provided to the FASTA file ",
                     normalizePath(path_to_fasta_for_f_data),
                     "':\n\n  * ",
                     paste(missing_mappings, collapse = "\n  * "),
                     call. = FALSE
                  )
               }
            }
         }
      )


      f_data <- ratios %>%
         mutate(feature_id = Index) %>%
         arrange(MaxPepProb) %>%
         filter(!duplicated(feature_id)) %>%
         as.data.frame() %>%
         `rownames<-`(.$feature_id) %>%
         select(
            feature_id,
            1:(where_samples_start - 1)
         ) %>%
         left_join(feature_df,
            by = c("Index" = feature_name_col)
         ) %>%
         `rownames<-`(.$feature_id)

      semicolon_rows <- which(grepl(";", f_data$ensembl_gene))
      if (length(semicolon_rows) > 0 && semicolon_warning) {
         warning("Removed ",
            length(semicolon_rows),
            " rows with IDs (semicolon ';' in 'Index' column).",
            call. = FALSE
         )
      }

      return(f_data)
   }

   get_p_data <- function() {
      if (!identical(sort(rownames(p_data_df)), sort(colnames(ratios)[where_samples_start:ncol(ratios)]))) {
         stop("Row names of 'p_data_df' must match the sample names in the ratios data.frame.", call. = FALSE)
      }
      p_data_df <- p_data_df[colnames(ratios)[where_samples_start:ncol(ratios)], ]
      return(p_data_df)
   }

   get_e_data <- function(f_data) {
      e_data <- ratios %>%
         `rownames<-`(.$Index) %>%
         .[
            rownames(f_data),
            where_samples_start:ncol(ratios)
         ] %>%
         as.matrix()
      return(e_data)
   }

   f_data <- get_f_data()
   p_data <- get_p_data()
   e_data <- get_e_data(f_data)

   m <- MSnSet(e_data, f_data, p_data)
   return(m)
}
