
#' Rolling-Up Features Abundance Measurements
#'
#' A primary purpose of this function is to roll-up peptide
#' measurements to protein level. Ultimately it can be used
#' for rolling-up or combining any features according to a grouping
#' factor. E.g. from peptides with charges to just the peptide sequences.
#'
#' @param msnset msnset (or most likely eset subclass) object
#' @param rollBy character. A column name in fData(msnset)
#' @param rollFun function for roll-up. "-" for log-transformed data,
#'                      "/" for normal scale
#' @param algorithm Either scaling to the most observed peptided ("reference")
#'                  or simply summing ("sum"). Summing makes sense only in
#'                  case of original label-free intensities.
#' @param verbose logical. Controls real-time output.
#'
#' @return MSnSet object
#' @importFrom Biobase exprs fData
#' @importFrom outliers grubbs.test outlier
#' @importFrom stats rnorm aggregate
#' @importFrom rlang as_function
#' @export rrollup
#'
#' @examples
#' data(srm_msnset)
#' dim(msnset)
#' msnset2 <- rrollup(msnset, rollBy="Gene", rollFun="-", verbose=TRUE)
#' dim(msnset2)


rrollup <- function(msnset, rollBy, rollFun,
                    algorithm = c("reference", "sum", "median", "mean", "min", "max"),
                    verbose=TRUE){

    algorithm <- match.arg(algorithm)

    if(algorithm == "reference"){
        summarisedFeatures <- list()
        unique_rollBy <- unique(fData(msnset)[[rollBy]])
        for (i in 1:length(unique_rollBy)) {
            # Subset msnset to each rollBy group
            msnset_sub <- msnset[fData(msnset)[[rollBy]] == unique_rollBy[i], ,
                                 drop = FALSE]
            summarisedFeatures[[i]] <- rrollup_a_feature_set(msnset_sub, rollBy,
                                                             rollFun, verbose)
        }
        exprs.new <- do.call(rbind, summarisedFeatures)
        rownames(exprs.new) <- unique_rollBy
    } else {
        temp <- data.frame(rollBy = fData(msnset)[[rollBy]],
                           exprs(msnset),
                           check.names = FALSE)
        temp[is.na(temp)] <- 0
        temp <- aggregate(. ~ rollBy, temp, as_function(algorithm))
        temp[temp == 0] <- NA
        exprs.new <- as.matrix(temp[,-1])
        rownames(exprs.new) <- temp[,1]
    }

    if(nrow(unique(fData(msnset))) == nrow(exprs.new)){
        fData.new <- unique(fData(msnset))
        rownames(fData.new) <- fData.new[, rollBy]
        fData.new <- fData.new[rownames(exprs.new),]
    }else{
        fData.new <- data.frame(rollBy = rownames(exprs.new),
                                stringsAsFactors = FALSE)
        colnames(fData.new) <- rollBy
        rownames(fData.new) <- fData.new[, rollBy]
    }
    msnset.new <- MSnSet(exprs = as.matrix(exprs.new),
                         fData = fData.new,
                         pData = pData(msnset))
    return(msnset.new)
}



rrollup_a_feature_set <- function(mat, rollBy, rollFun, verbose){
    # Get the group name (protein name)
    rollBy_name <- unique(fData(mat)[[rollBy]])

    # We don't want a named vector if nrow(exprs(mat)) == 1
    # We need a data frame
    mat <- exprs(mat) %>% as.data.frame()

    if (nrow(mat) > 1) {
        # Calculate reference index
        refIdx <- which.max(rowSums(!is.na(mat)))
        # Count number of times each row overlaps with reference
        overlap_count <-
            apply(mat, 1,
                  function(x)
                      sum(!is.na(as.numeric(x) == as.numeric(mat[refIdx, ]))))
        overlap_count <- as.numeric(overlap_count)

        # Subset to rows that overlap more than once with reference
        mat <- mat[overlap_count != 1, , drop = FALSE]

        # If rows were removed and verbose = TRUE, output a message
        # that includes the group name.
        if (any(overlap_count == 1) & verbose) {
            message(paste0("Some rows associated with group ",
                           rollBy_name, " have been removed."))
        }

        # Recalculate index using subset mat
        refIdx <- which.max(rowSums(!is.na(mat)))

        mat.ratios <- sweep(mat, 2, as.numeric(mat[refIdx,]), FUN = rollFun)
        scale.factors <- apply(mat.ratios, 1, median, na.rm = TRUE)
        mat <- sweep(mat, 1, scale.factors, FUN = rollFun)
        #
        # for each sample, check if any of the feature (peptide)
        # abundance measurements
        # representing the set (protein) are not in line with other.
        maxVals <- apply(mat, 1, max, na.rm = TRUE)
        # If all values in a row are NA, the max is -Inf
        # Change infinite values to NA
        maxVals[is.infinite(maxVals)] <- NA

        # use unique and non-NA values for outlier detection
        uniqueMaxVals <- unique(maxVals[!is.na(maxVals)])

        if(nrow(mat) > 2){
            while(grubbs.test(uniqueMaxVals)$p.value < 0.05 &
                  sum(!is.na(uniqueMaxVals)) > 2){
                i <- which(uniqueMaxVals == outliers::outlier(uniqueMaxVals))
                uniqueMaxVals[i] <- NA
            }
        }

        # drop NA value
        uniqueMaxVals <- uniqueMaxVals[!is.na(uniqueMaxVals)]

        retained_idx <- maxVals %in% uniqueMaxVals
        mat <- mat[retained_idx, , drop = FALSE]
        # protein abundance estimates
        protProfile <- apply(mat, 2, median, na.rm=TRUE)
    } else {
        # The output should be a named vector,
        # not a data frame with 1 row.
        protProfile <- mat[1, ]
    }

    return(protProfile)
}







