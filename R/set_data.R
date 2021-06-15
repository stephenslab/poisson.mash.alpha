#' @title Create Data Oject for Poisson Mash Analysis
#'
#' @description Add slightly more detailed description here.
#' 
#' @param Y J x N matrix of counts with features (e.g., genes) as
#'   rows and observations (e.g., cells) as columns.
#' 
#' @param condition Factor of length N with R levels denoting
#'   the condition status of each observation.
#' 
#' @param si Numeric vector of length N containing the size factor for
#'   each of the N observations.
#' 
#' @param subgroup Named vector of factors with M levels denoting the
#'   subgroup status of each of R conditions.  Default to no subgroup
#'   (M=1). If not set to \code{NULL}, the names of subgroup must match
#'   the levels of condition.
#' 
#' @return A pois.mash data object for poisson mash analysis. It is a list
#'   with the following components:
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
#' @seealso
#'
#' @examples
#' # Add examples here.
#' 
#' @import Matrix
#' 
#' @export
#' 
pois_mash_set_data <- function (Y, condition, si = colSums(Y),
                                subgroup = gl(1,nlevels(condition))) {

  # Check and process the inputs.
  if (ncol(Y) != length(condition))
    stop("The number of columns of Y and the length of condition do not match")
  if (ncol(Y) != length(si))
    stop("The number of columns of Y and the length of si do not match")
  
  J    <- nrow(Y)
  trts <- levels(condition)
  R    <- nlevels(condition)
  
  # subgroup <- subgroup[order(names(subgroup))]
  # if(sum(trts != names(subgroup)) > 0)
  #   stop("The levels of condition and the names of subgroup do not match")

  # Aggregate the individual-level data into condition-level data.
  X <- matrix(as.numeric(NA), nrow=J, ncol=R)
  rownames(X) <- rownames(Y)
  colnames(X) <- trts
  s <- tapply(si,condition,sum)
  names(s) <- trts
  for(r in 1:R) 
    X[,r] <- rowSums(Y[, condition == trts[r], drop = FALSE])

  dat <- list(X=X, s=s/min(s), subgroup=subgroup)
  class(dat) <- c("pois.mash","list")
  return(dat)
}
