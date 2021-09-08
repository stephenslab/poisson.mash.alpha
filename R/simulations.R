#' @title Simulate Data Based on a Real Single-cell Dataset
#' 
#' @param J Number of features to simulate (should be 2 or greater).
#' 
#' @param R Number of conditions to simulate (should be 2 or greater).
#' 
#' @param effects List consisting of the following components:
#'   \dQuote{umat}, an R x K matrix with each of the K columns
#'   specifying a pattern of non-null effects; \dQuote{prob}, a
#'   non-negative vector of length K specifying the relative frequency
#'   of each non-null effect among all J features.
#' 
#' @param seed Random seed (for reproducing the output).
#' 
#' @return A list containing the following components:
#' 
#' \item{Y}{J x N sparse matrix of counts with features (e.g., genes)
#'   as rows and observations (e.g., cells) as columns.}
#' 
#' \item{condition}{Factor of length N with R levels giving the
#'   condition for each observation.}
#'
#' @examples
#' umat <- cbind(c(1,1,0,0,0),c(0,0,0,1,0.5))
#' pois_mash_sim_data(1000,5,effects = list(umat = umat,prob = rep(0.05,2)))
#'
#' @importFrom methods as
#' @importFrom stats model.matrix
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats rmultinom
#' @importFrom seqgendiff thin_diff
#' 
#' @export
#' 
pois_mash_sim_data <- function (J = 1000, R = 5, effects = NULL, seed = 1) {
  set.seed(seed)
  data(sim_raw)
  n <- nrow(sim_raw)
  if (n >= J) 
    Y <- sim_raw[sample(n,J,replace = FALSE),]
  else
    Y <- sim_raw[sample(n,J,replace = TRUE),]
  condition <- as.factor(sample(R,ncol(Y),replace = TRUE))
  
  if (is.null(effects)) {
      
    # Return null data if "effects" is not specified.
    return(list(Y = Y,condition = condition))
  } else {
    umat <- effects$umat
    prob <- effects$prob
    
    if (nrow(umat) != R)
      stop("The number of rows of umat should be equal to the number of ",
           "conditions")
    if (ncol(umat) != length(prob))
      stop("The number of columns of umat should be equal to the length ",
           "of prob")
    if (any(prob < 0) | sum(prob) > 1)
      stop("Each element of prob should be non-negative with its sum no ",
           "greater than 1")
    
    # Simulate non-null effects.
    designmat      <- model.matrix(~0 + condition)
    beta           <- matrix(0,J,R)
    non.null.idx   <- sample(J,round(sum(prob)*J),replace = FALSE)
    prob.normalize <- prob/sum(prob)
    
    for (j in 1:J) {
      if (j %in% non.null.idx) {
        u.j <- umat[,which(as.numeric(rmultinom(1,1,prob.normalize)) == 1)]
        beta[j,] <- rnorm(1,1,0.2) * ifelse(runif(1) > 0.5,1,-1)*u.j
      }
    }
    
    # Create signals through thinning.
    scdata.mod <- thin_diff(mat = as.matrix(Y),design_fixed = designmat,
                            coef_fixed = beta)
    Y.mod <- as(scdata.mod$mat,"dgCMatrix")
    rownames(Y.mod) <- rownames(Y)
    colnames(Y.mod) <- colnames(Y)
    return(list(Y = Y.mod,condition = condition))
  }
}
