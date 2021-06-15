#' @title Create Data Oject for Poisson Mash Analysis
#'
#' @description Add description here.
#' 
#' @param Y A J by N matrix of counts with features (e.g., genes) as
#'   rows and observations (e.g., cells) as columns
#' 
#' @param condition An N by 1 vector of factors with R levels denoting
#'   the condition status of each observation
#' 
#' @param si An N by 1 vector of size factors for each observation 
#' 
#' @param subgroup An R by 1 named vector of factors with M levels
#'   denoting the subgroup status of each of R conditions.  Default to
#'   no subgroup (M=1). If not set to \code{NULL}, the names of subgroup
#'   must match the levels of condition.
#' 
#' @return A pois.mash data object for poisson mash analysis,
#'   including the following components:
#' 
#' \item{X}{J x R matrix of count data collapsed over conditions, with
#'   features as rows and conditions as columns.}
#' 
#' \item{s}{R x 1 numeric vector adjusting for sequencing depth of
#'   each of R conditions.}
#' 
#' \item{subgroup}{R x 1 factor vector with M levels denoting the
#'   subgroup status of each of R conditions.}
#' 
#' @export
#' 
pois_mash_set_data <- function(Y, condition, si, subgroup=NULL) {
  
  if (ncol(Y) != length(condition))
    stop("The number of columns of Y and the length of condition do not match")
  if (ncol(Y) != length(si))
    stop("The number of columns of Y and the length of si do not match")
  
  J <- nrow(Y)
  trts <- sort(unique(condition))
  R <- length(trts)
  
  if(!is.null(subgroup)) {
    subgroup <- subgroup[order(names(subgroup))]
    if(sum(trts != names(subgroup)) > 0)
      stop("The levels of condition and the names of subgroup do not match")
  } else {
    subgroup <- rep(1, R)
    names(subgroup) <- trts
  }

  # Aggregate the cell level data into condition level data.
  X <- matrix(as.numeric(NA), nrow=J, ncol=R)
  rownames(X) <- rownames(Y)
  colnames(X) <- trts
  s <- rep(as.numeric(NA), R)
  names(s) <- trts
  
  for(r in 1:R) {
    Y.tmp <- Y[, condition==trts[r]]
    X[,r] <- rowSums(Y.tmp)
    s[r] <- sum(si[condition==trts[r]])
  }
  
  data <- list(X=X, s=s/min(s), subgroup=subgroup)
  class(data) <- "pois.mash"
  return(data)
}
