

#' @importFrom varSelRF varSelRF
select_features_varSelRF <- function(x, y, ...){
  rfvs <- varSelRF(xdata = x, Class = y,
                   verbose=FALSE,
                   recompute.var.imp=FALSE, # tweak
                   whole.range = TRUE,
                   ntree=10000, # tweak
                   ntreeIterat=10000, # tweak
                   vars.drop.num=NULL,
                   vars.drop.frac=0.2,
                   ...)
  hi <- rfvs$selec.history
  idx <- which.min(hi$OOB)
  idx <- min(idx + 5, nrow(hi)) # start from 5 steps back
  vrz <- hi[idx, "Vars.in.Forest"]
  vrz <- strsplit(as.character(vrz)," \\+ ")[[1]]
  #
  rfvs <- varSelRF(xdata = x[,vrz], Class = y,
                   verbose=FALSE,
                   recompute.var.imp=FALSE, # tweak
                   whole.range = TRUE,
                   ntree=10000, # tweak
                   ntreeIterat=10000, # tweak
                   vars.drop.num=1,
                   vars.drop.frac=NULL,
                   ...)
  return(rfvs$selected.vars)

}


#' @importFrom Boruta Boruta getSelectedAttributes
select_features_Boruta <- function(x, y, ...){

  #     # .. add on ..
  #     # .. preselection of features using plain random forest
  #     # .. perhaps I can skimp on Boruta settings
  #     #
  #     rf <- randomForest(x, y, data=dSet, ntree = 10000)
  #     f2 <- names(head(importance(rf)[order(importance(rf)[,1],decreasing=TRUE),],20))
  #     x <- x[,f2]
  #     #
  rfbo <- Boruta(x, y,
                 #maxRuns = 100, # 1000
                 doTrace=0,
                 #ntree=5000, # 10000
                 ...)
  # 5h with 10000
  # 10m with 1000
  # ~ 5 min with 500
  # features_i <- getSelectedAttributes(TentativeRoughFix(rfbo))
  features_i <- getSelectedAttributes(rfbo)
  if(identical(character(),features_i))
    features_i <- getSelectedAttributes(rfbo, withTentative = TRUE)
  if(identical(character(),features_i)){
    # return just the best one
    meds <- apply(rfbo$ImpHistory,2,median)
    meds <- meds[names(rfbo$finalDecision)]
    features_i <- names(rev(sort(meds)))[1]
  }
  return(features_i)
}


#' @importFrom utils head
#' @importFrom randomForest importance
select_features_top <- function(x, y, ...){
  r <- randomForest(x, y, ntree=1000, ...)
  ordr <- order(importance(r), decreasing = TRUE)
  imp <- importance(r)[ordr,,drop=FALSE]
  return(rownames(head(imp, 5)))
}



train_model_rf <- function(x, y, ...){
  r <- randomForest(x, y, ntree=1000, ...)
  return(r)
}



#' @title Random Forest Predictive Models
#'
#' @description The objective is - given the number of features, select the most
#'   informative ones and evaluate the predictive random forest model. The
#'   feature and model selection performed independently for each round of
#'   LOOCV.
#'
#' @param msnset MSnSet object
#' @param features character vector features to select from for building
#'   prediction model. The features can be either in featureNames(msnset) or in
#'   pData(msnset).
#' @param response factor to classify along. Must be only 2 levels.
#' @param pred.cls class to predict
#' @param K specifies the cross-validation type. Default NULL means LOOCV.
#'   Another typical value is 10.
#' @param sel.feat logical defining if to select features or use the entire set?
#' @param sel.alg character.
#'      \itemize{
#'          \item varSelRF
#'          \item Boruta
#'          \item top - just selecting top 3 features
#'      }
#' @param cores integer; the number of processes to use for Boruta.
#'   If \code{NULL}, one less than the maximum number of processes will be used.
#'   If this number is greater than 125, a maximum of 125 cores will be used.
#' @param seed numeric; the random seed to use during various stages of 
#'   \code{rf_modeling}. Passed to \code{\link[base]{set.seed}}.
#' @param ... Extra arguments. Currently passed only to Boruta algorithm.

#' @return list
#'      \describe{
#'          \item{\code{prob}}{is the probabilities (response) from LOOCV
#'          that the sample is "case"}
#'          \item{\code{features}}{list of selected features
#'          for each iteration of LOOCV}
#'          \item{\code{top}}{top features over all iterations}
#'          \item{\code{auc}}{AUC}
#'          \item{\code{pred}}{prediction perfomance obtained by
#'                  \code{ROCR::prediction}}
#'      }
#'
#' @importFrom randomForest randomForest
#' @importFrom ROCR prediction performance
#' @importFrom parallel mclapply detectCores parLapply clusterEvalQ
#' @importFrom stats predict
#' @importFrom Biobase featureNames
#'
#' @export rf_modeling
#'
#' @examples
#' # Not run
#' \dontrun{
#' data(srm_msnset)
#' head(varLabels(msnset))
#' head(msnset$subject.type)
#' # reduce to two classes
#' msnset <- msnset[,msnset$subject.type != "control.1"]
#' msnset$subject.type <- as.factor(msnset$subject.type)
#  out <- rf_modeling(msnset,
#                    features=featureNames(msnset)[1:5],
#                    response="subject.type",
#                    pred.cls="case")
#' plotAUC(out)
#' # top features consistently recurring in the models during LOOCV
#' print(out$top)
#' # the AUC
#' print(out$auc)
#' # probabilities of classifying the sample right, if the feature selection
#' # and model training was performed on other samples
#' plot(sort(out$prob))
#' abline(h=0.5, col='red')
#' }


rf_modeling <- function(msnset,
                        features,
                        response,
                        pred.cls,
                        K=NULL,
                        sel.feat=TRUE,
                        sel.alg = c("varSelRF","Boruta","top"),
                        cores = NULL,
                        seed = 0,
                        ...) {

  # to avoid running processes in case the function crashes
  on.exit(stopCluster(multiproc_cl))

  # prepare data
  dSet <- cbind(pData(msnset), t(exprs(msnset)))
  stopifnot(length(unique(dSet[,response])) == 2)
  stopifnot(pred.cls %in% dSet[,response] )
  if(unique(dSet[,response])[1] == pred.cls){
    lvlz <- rev(unique(dSet[,response]))
  }else{
    lvlz <- unique(dSet[,response])
  }
  dSet[,response] <- factor(dSet[,response], levels=lvlz)
  #---
  sel.alg <- match.arg(sel.alg)
  FUN <- switch(sel.alg,
                varSelRF = select_features_varSelRF,
                Boruta = function(x,y, ...) select_features_Boruta(x,y,...),
                top = select_features_top)

  # do K-fold split here
  if(is.null(K))
    K <- nrow(dSet)
  num_rep <- ceiling(nrow(dSet)/K)
  set.seed(seed)
  cv_idx <- sample(rep(seq_len(K), num_rep)[seq_len(nrow(dSet))])
  if(is.null(cores)){
    cores <- max(1, min(125, detectCores() - 1))
  }
  stopifnot(1 <= cores)
  if(cores > detectCores()){
    msg <- sub("\n", "", "The number of specified processes is greater than
                the number of cores available on this computer.
                This may lead to high computational overhead.")
    warning(msg)
  }
  multiproc_cl <- makeCluster(cores)

  clusterEvalQ(multiproc_cl, invisible(suppressWarnings({
    Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
    suppressWarnings(suppressPackageStartupMessages({
      library("MSnID")
      library("Biobase")
    }))
  })))
  set.seed(seed)
  seed_seq <- sample(1:.Machine$integer.max, size = K)
  invisible(clusterExport(multiproc_cl,
                          c("dSet","cv_idx","features",
                            "response", "seed_seq"),
                          envir = environment()))
  fn <- function(i){
    RNGkind("L'Ecuyer-CMRG")
    seed <- seed_seq[i]
    set.seed(seed)
    i <- cv_idx == i

    if(sel.feat){
      features.sel <- FUN(x=dSet[!i,features],
                          y=dSet[!i,response], seed = seed, ...)
    }else{
      features.sel <- features
    }

    # train model
    x=dSet[!i,features.sel,drop=FALSE]
    colnames(x) <- make.names(colnames(x))
    mdl <- train_model_rf(x=x, y=dSet[!i,response])
    # predict
    newdata <- dSet[i,features.sel,drop=FALSE]
    colnames(newdata) <- make.names(colnames(newdata))

    predProb <- predict(mdl, newdata=newdata, type='prob')[,pred.cls]
    names(predProb) <- rownames(newdata)
    # print(i)
    list(predProb, features.sel)
  }
  X <- mapply(function(i, seed) list(i, seed), 1:K, seed_seq)
  res <- parLapply(cl = multiproc_cl, X = 1:K, fun = fn)

  predProb <- unlist(sapply(res, '[[', 1, simplify = FALSE)) # unlist TODO
  predProb <- predProb[rownames(dSet)]
  names(predProb) <- NULL # for compatibility
  selected <- sapply(res, '[[', 2)
  #---
  pred <- prediction(predProb, dSet[,response] == pred.cls)
  return(list(prob = predProb,
              features = selected,
              top = rev(sort(table(unlist(selected)))),
              auc = performance(pred,"auc")@y.values[[1]],
              pred = pred))
}





rf_modeling_old <- function( msnset, features, response, pred.cls, sel.feat=TRUE,
                             sel.alg=c("varSelRF","Boruta","top"), ...){
  #
  if(.Platform$OS.type == "unix")
    mc.cores <- detectCores()
  else
    mc.cores <- 1
  # prepare data
  dSet <- cbind(pData(msnset), t(exprs(msnset)))
  #
  stopifnot(length(unique(dSet[,response])) == 2)
  stopifnot(pred.cls %in% dSet[,response] )
  if(unique(dSet[,response])[1] == pred.cls){
    lvlz <- rev(unique(dSet[,response]))
  }else{
    lvlz <- unique(dSet[,response])
  }
  dSet[,response] <- factor(dSet[,response], levels=lvlz)
  #---
  sel.alg <- match.arg(sel.alg)
  FUN <- switch(sel.alg,
                varSelRF = select_features_varSelRF,
                Boruta = function(x,y) select_features_Boruta(x,y,...),
                top = select_features_top)
  res <- mclapply(1:nrow(dSet), function(i){
    if(sel.feat){
      features.sel <- FUN(x=dSet[-i,features],
                          y=dSet[-i,response])
    }else{
      features.sel <- features
    }
    # train model
    x=dSet[-i,features.sel,drop=FALSE]
    colnames(x) <- make.names(colnames(x))
    mdl <- train_model_rf(x=x, y=dSet[-i,response])
    # predict
    newdata <- dSet[i,features.sel,drop=FALSE]
    colnames(newdata) <- make.names(colnames(newdata))
    predProb <- predict(mdl, newdata=newdata, type='prob')[1,pred.cls]
    print(i)
    list(predProb, features.sel)
  }, mc.cores=mc.cores, mc.set.seed = FALSE)
  #
  predProb <- sapply(res, '[[', 1)
  names(predProb) <- NULL # for compatibility
  selected <- sapply(res, '[[', 2)
  #---
  pred <- prediction(predProb, dSet[,response] == pred.cls)
  return(list(prob = predProb,
              features = selected,
              top = rev(sort(table(unlist(selected)))),
              auc = performance(pred,"auc")@y.values[[1]],
              pred = pred))
}


