#' Reading MSFragger-generated tmt-report files from a file path as MSnSet object
#'
#' @description Function has only been tested with TMT intensity-based
#'   quantification data. Desired tmt-report output files (e.g., "ratio_multi-site_MD.tsv")
#'   must be properly selected in the FP settings.
#'
#' @param path character; File path to the desired FragPipe-generated tmt-report file.
#'   Any tmt-report file may be used.
#' @param org_to_retain character; Filtering out contaminants. The argument is the
#'   official organism name such as "Homo sapiens" or "Bos taurus" or
#'   "Sus scrofa". Default is NULL, meaning pass all.
#' @param use_gene_as_prot_id logical; Used only in case of `single-site` files. Switches
#'   notation from UniProt_Site to a more human-readable and conventional Gene-Site.
#'   Default it `TRUE`. In case there are duplicates in the new site IDs, returns
#'   the error message with the prompt to switch to `FALSE`.
#'
#' @return (MSnSet) MSnSet object of MSFragger TMT results
#'
#' @importFrom MSnbase MSnSet
#' @importFrom data.table fread
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>% select mutate
#'
#' @examples
#' \dontrun{
#' file_path <- "C:/Users/fakeusr222/Desktop/MSF_TMT_job/ratio_multi-site_MD.tsv"
#'   msnset <- read_FragPipe_TMT(file_path)
#'   show(msnset)
#' }
#'
#' @export read_FragPipe_TMT

read_FragPipe_TMT <- function(path = NULL, org_to_retain = NULL, use_gene_as_prot_id = TRUE)
{

   path_to_file <- path

   if (!file.exists(path_to_file)) {
      stop(sprintf("file not found in folder: %s", dirname(path_to_file)))
   }

   df <- fread(file = path_to_file, showProgress = FALSE, data.table = FALSE)

   if(!is.null(org_to_retain)){
      combined_protein_path <- file.path(dirname(dirname(path_to_file)), "combined_protein.tsv")
      retained_proteins <- fread(file = combined_protein_path,
                                 showProgress = FALSE, data.table = FALSE) %>%
         filter(Organism == org_to_retain) %>%
         distinct(`Protein ID`) %>%
         rename(ProteinID = `Protein ID`)

      df <- semi_join(df, retained_proteins, by = "ProteinID")

   }

   # make featureNames
   if (grepl("peptide", basename(path_to_file))) {
      df <- df %>%
         mutate(rowname = paste(Gene, ProteinID, Peptide, sep = "|"))
   }
   else if (grepl("multi-site", basename(path_to_file))) {
      df <- df %>% mutate(rowname = paste(Gene,
                                          ProteinID,
                                          str_split_i(Index, "_", i = 4),
                                          Peptide,
                                          sep = "|"))
      df <- df %>%
         dplyr::filter(
            str_split_i(Index, "_", i = 4) == str_split_i(Index, "_", i = 5))
   }
   else if (grepl("single-site", basename(path_to_file))) {
      if(use_gene_as_prot_id){
         df <- df %>%
            filter(Gene != "") %>%
            mutate(rowname = paste0(Gene,
                                    "-",
                                    sub("[^_]*_([A-Z]\\d+)","\\1",Index)))
         if(anyDuplicated(df$rowname)){
            # let's try to resolve by ReferenceIntensity
            if(!("ReferenceIntensity" %in% colnames(df))){
               stop("Duplicates in the gene-based site names. Can't resolve ambiguity.
                   Switch to use_gene_as_prot_id = FALSE.")
            } else {
               df <- df %>%
                  group_by(rowname) %>%
                  slice_max(ReferenceIntensity)
            }
         }
      }else{
         df <- df %>% mutate(rowname = Index)
      }
   }
   else if (grepl("gene", basename(path_to_file))) {
      df <- df %>%
         mutate(rowname = paste(Index, ProteinID, sep = "|"))
   }
   else if (grepl("protein", basename(path_to_file))) {
      df <- df %>%
         mutate(rowname = paste(Gene, Index, sep = "|"))
   }
   else{
      stop("unknown file")
   }

   df <- df %>%
      mutate(featureName = rowname, .before = colnames(.)[[1]]) %>%
      column_to_rownames(var = "rowname")

   x_data <- df %>%
      select(-c(colnames(.)[[1]]:ReferenceIntensity)) %>%
      as.matrix()

   p_data <- data.frame(sample_name = colnames(x_data),
                        row.names = colnames(x_data))

   f_data <- df %>%
      select(c(colnames(.)[[1]]:ReferenceIntensity))

   m <- MSnSet(exprs = x_data, fData = f_data, pData = p_data)

   return(m)
}



utils::globalVariables(
   c(".", "featureName")
)

