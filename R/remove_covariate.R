
#' Removing Covariate Effect form Expression Data
#' 
#' The main purpose of this function is to remove batch effect from the data.
#' Batch can be associated with different days of sample processing (as factor)
#' or with run order (continuous). Can also be used to remove any unwanted
#' effects from the data.
#' 
#' @param x MSnSet or ExpressionSet object
#' @param covName covariate name. Must be in pData(x). At this point it can be
#' only one name.
#' 
#' @note The algorithm essentially uses an LM. The reason for re-inventing the 
#' wheel is presense of missing values in proteomics datasets more then usual.
#' @seealso \code{\link[sva]{ComBat}}
#' 
#' @importFrom Biobase exprs pData
#' @importClassesFrom Matrix dgCMatrix
#' @export remove_covariate
#' 
#' @examples
#' # example 1
#' set.seed(1)
#' means <- rep(c(1,2,3), each=3)
#' nrows <- 5
#' e <- matrix(rep(means, nrows), ncol=length(means), byrow=T)
#' e <- e + 
#'     matrix(rnorm(Reduce(`*`, dim(e)), sd=0.3), ncol=length(means), byrow=T)
#' 
#' # add missing values in increasing frequency
#' extreme <- 10 # controls how quickly increases propotion of NAs
#' # 1 - means it will reach 100% by the end
#' # 2 - means only 50% will be missing by the last row
#' # N - is 1/N-th
#' freqs <- (1:nrow(e)-1)/(extreme*(nrow(e)-1))
#' mis <- t(sapply(freqs, rbinom, n=ncol(e), size = 1))
#' 
#' mis[mis == 1] <- NA
#' e[5,8:9] <- NA
#' e[4,4:6] <- NA
#' e <- e + mis
#' image(e) 
#' library("ggplot2"); library("reshape2")
#' ggplot(melt(e), aes(x=Var1, y=Var2, fill=value)) + geom_raster()
#' 
#' # generating factors
#' facs <- gl(length(unique(means)),length(means)/length(unique(means)))
#' # alternative is correction for continuous variable
#' cova <- seq_along(means)
#' 
#' library("Biobase")
#' m <- ExpressionSet(e)
#' pData(m)$pesky <- facs
#' pData(m)$runorder <- cova
#' m2 <- remove_covariate(m, "pesky")
#' m3 <- remove_covariate(m, "runorder")
#' 
#' image(exprs(m))
#' image(exprs(m2))
#' image(exprs(m3))
#' 
#' # Example 2 (real-world)
#' data(cptac_oca)
#' # let's test for iTRAQ_Batch effect
#' res <- eset_lm(oca.set, "y ~ iTRAQ_Batch", "y ~ 1")
#' # not too strong, but there
#' hist(res$p.value, 50)
#' image_msnset(oca.set, facetBy="iTRAQ_Batch")
#' oca.fixed <- remove_covariate(oca.set, "iTRAQ_Batch")
#' res <- eset_lm(oca.fixed, "y ~ iTRAQ_Batch", "y ~ 1")
#' hist(res$p.value, 50)
#' image_msnset(oca.fixed, facetBy="iTRAQ_Batch")

remove_covariate <- function(x, covName){
    e <- exprs(x)
    rmns <- apply(e, 1, mean, na.rm=TRUE) # to add later
    cova <- pData(x)[[covName]]
    # the reason for splitting into factor vs continuous
    # is that I don't want to rely on (Intercept) reference group in factor
    if(is.factor(cova) || is.character(cova)){
        desmat <- model.matrix( ~ cova + 0)
        suppressWarnings(
            cfs <-
                Reduce(rbind,
                       lapply(1:nrow(e), function(i){
                           cf <- coefficients(lm(e[i,] ~ cova + 0))
                       }))
        )
    }else if(is.numeric(cova)){
        desmat <- model.matrix( ~ cova)
        cfs <-
            Reduce(rbind,
                   lapply(1:nrow(e), function(i){
                       cf <- coefficients(lm(e[i,] ~ cova))
                   }))
    }else{
        stop("unknown type of covariate")
    }
    btch <- t(as.matrix(as(desmat,'dgCMatrix') %*% as(t(cfs),'dgCMatrix')))
    e.nobatch <- e - btch
    # really necessary in case one factor is fully NA
    rmns <- rmns - apply(e.nobatch, 1, mean, na.rm=TRUE)
    #..
    e.backmeans <- sweep(e.nobatch, 1, rmns, '+')
    exprs(x) <- e.backmeans
    return(x)
}