#' @title Create some simulated data based on a real single cell dataset
#' 
#' @param J the number of features which is at least 2.
#' 
#' @param R the number of conditions which is at least 2.
#' 
#' @param effects a list consisting of the following components: 
#' \item{umat}{R x K matrix with each of the K (>=1) columns specifying a pattern of non-null effects.}
#' \item{prob}{Non-negative vector of length K specifying the relative frequency of each non-null effect among all J features.}
#' 
#' @param seed a random seed for reproducibility.
#' 
#' @examples 
#' pois_mash_sim_data(J=1000, R=5, effects=list(umat=cbind(c(1,1,0,0,0), c(0,0,0,1,0.5)), prob=c(0.05, 0.05))) 
#' 
#' @return a list containing the following components:
#' \item{Y}{J x N sparse matrix of counts with features (e.g., genes) as rows and observations (e.g., cells) as columns.}
#' \item{condition}{Factor of length N with R levels giving the condition for each observation.}


pois_mash_sim_data <- function(J=1000, R=5, effects=NULL, seed=1){
  library(Matrix)
  library(seqgendiff)
  set.seed(seed)
  
  scdata <- readRDS("data/sim_raw_count.rds")
  if(nrow(scdata) >= J){
    Y <- scdata[sample(1:nrow(scdata), J, replace = FALSE),]
  }
  else{
    Y <- scdata[sample(1:nrow(scdata), J, replace = TRUE),]
  }

  condition <- as.factor(sample(1:R, ncol(Y), replace = TRUE))
  
  # create null data if effects is not specified
  if(is.null(effects)){
    return(list(Y=Y, condition=condition))
  }
  
  else{
    umat <- effects$umat
    prob <- effects$prob
    
    if(nrow(umat)!=R){
      stop("The number of rows of umat should be equal to the number of conditions")
    }
    if(ncol(umat)!=length(prob)){
      stop("The number of columns of umat should be equal to the length of prob")
    }
    if(any(prob < 0) | sum(prob) > 1){
      stop("Each element of prob should be non-negative with its sum no greater than 1")
    }
    
    ### simulate non-null effects
    designmat <- model.matrix(~0+condition)
    beta <- matrix(0, nrow=J, ncol=R)
    non.null.idx <- sample(1:J, round(sum(prob)*J), replace = FALSE)
    prob.normalize <- prob/sum(prob)
    
    for(j in 1:J){
      if(j %in% non.null.idx){
        u.j <- umat[, which(as.numeric(rmultinom(1, 1, prob.normalize))==1)]
        beta[j,] <- rnorm(1, 1, 0.2)*ifelse(runif(1) > 0.5, 1, -1)*u.j
      }
    }
    
    ### create signals through thinning
    scdata.mod <- thin_diff(mat = as.matrix(Y), design_fixed = designmat, coef_fixed = beta)
    Y.mod <- as(scdata.mod$mat, "sparseMatrix") 
    rownames(Y.mod) <- rownames(Y)
    colnames(Y.mod) <- colnames(Y) 
    return(list(Y=Y.mod, condition=condition))
  }
}


