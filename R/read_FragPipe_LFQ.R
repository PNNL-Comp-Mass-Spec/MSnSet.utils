#' Reading MSFragger-generated LFQ-based MSstats from a file path as MSnSet object
#'
#' @description Function has only been tested with label-free intensity-based
#'   quantification data. MSstats.csv is
#'   an optional output file which needs to be specified in FP settings.
#'
#' @param path character; File path to the FragPipe-generated MSstats.csv file
#' @param type Either "CP" or "MSstats". In the first case reads the results
#' from MaxQuant-like output file "combined_protein.tsv". In the second the
#' reults are taken from "MSstats.csv".
#'
#' @return (MSnSet) MSnSet object of MSFragger LFQ results
#'
#' @importFrom MSnbase MSnSet
#' @importFrom data.table fread
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>% select filter distinct relocate everything mutate
#'
#' @examples
#' \dontrun{
#' file_path <- "C:/Users/fakeusr222/Desktop/MSF_LFQ_job/MSstats.csv"
#' msnset <- read_FragPipe_LFQ(file_path)
#' show(msnset)
#' }
#'
#' @export read_FragPipe_LFQ



read_FragPipe_LFQ <- function(path = NULL, type = c("CP", "MSstats")){

}


read_FragPipe_LFQ <- function(path = NULL, type = c("CP", "MSstats")) {

  type <- match.arg(type)

  FUN <- switch(type,
              "CP"      = read_FragPipe_LFQ.combined_protein,
              "MSstats" = read_FragPipe_LFQ.msstats)

  FUN(path = path)
}


read_FragPipe_LFQ.msstats <- function(path = NULL)
{
  path_to_file <- path

  if (!file.exists(path_to_file)) {
    stop(sprintf("MSstats.csv file not found in folder: %s", dirname(path_to_file)))
  }

  df <- fread(file = path_to_file, showProgress = FALSE, data.table = FALSE) %>%
    filter(!is.na(Intensity)) %>%
    # May add charge col later
    select(ProteinName, PeptideSequence, Run, Intensity) %>%
    mutate(featureName = paste0(ProteinName, "@", PeptideSequence)) %>%
    relocate(featureName, .before = everything())

  # Will sum intensity of unique features.
  x_data <- df %>%
    pivot_wider(id_cols = "featureName",
                names_from = "Run",
                values_from = "Intensity",
                values_fn = sum) %>%
    as.data.frame() %>%
    column_to_rownames(var = "featureName") %>%
    as.matrix()

  f_data <- df %>%
    distinct(featureName, ProteinName, PeptideSequence) %>%
    `rownames<-`(.[["featureName"]])

  p_data <- df %>%
    distinct(Run) %>%
    `rownames<-`(.[["Run"]])

  x_data <- x_data[rownames(f_data), rownames(p_data)]

  m <- MSnSet(exprs = x_data, fData = f_data, pData = p_data)

  return(m)
}


read_FragPipe_LFQ.combined_protein <- function(path){

  fpath <- list.files(path = path, pattern = "experiment_annotation.tsv", full.names = TRUE)
  expann <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)

  fpath <- list.files(path = path, pattern = "combined_protein.tsv", full.names = TRUE)
  x <- read.delim(fpath, check.names = FALSE, stringsAsFactors = FALSE)

  id.cols <- c("Protein", "Protein ID", "Entry Name", "Gene", "Protein Length",
               "Organism", "Protein Existence", "Description",
               "Protein Probability", "Top Peptide Probability",
               "Combined Total Peptides","Combined Spectral Count",
               "Combined Unique Spectral Count","Combined Total Spectral Count",
               "Indistinguishable Proteins")
  quantType <- "MaxLFQ Intensity"
  quant.cols <- paste(expann$sample_name, quantType, sep = " ")
  x.exprs <- as.matrix(x[, quant.cols])
  x.exprs[x.exprs == 0] <- NA
  rownames(x.exprs) <- x$Protein
  colnames(x.exprs) <- trimws(sub(quantType, "", colnames(x.exprs)))

  x.pdata <- data.frame(sample.name = colnames(x.exprs),
                        row.names = colnames(x.exprs),
                        stringsAsFactors = FALSE)

  x.fdata <- x[, c(id.cols)]
  rownames(x.fdata) <- x.fdata$Protein

  m <- MSnSet(exprs = x.exprs, fData = x.fdata, pData = x.pdata)
  validObject(m)

  # remove contaminants
  m <- m[!grepl("contam_", featureNames(m)),]

  # move feature names to genes
  # 1. remove missing value genes
  m <- m[fData(m)$Gene != "",]
  # 2. resolve duplicates
  proteins_of_unique_genes <- fData(m) %>%
    select(Protein, Gene, `Combined Unique Spectral Count`) %>%
    group_by(Gene) %>%
    slice_max(order_by = `Combined Unique Spectral Count`, with_ties = FALSE) %>%
    pull(Protein)
  m <- m[proteins_of_unique_genes,]
  featureNames(m) <- fData(m)$Gene
  validObject(m)
  return(m)
}


utils::globalVariables(
  c("ProteinName", "PeptideSequence", "Run", "Intensity", ".", "featureName")
)

