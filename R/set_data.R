#' @title Create Data Object for Poisson MASH Analysis
#'
#' @param Y J x N matrix of counts with features (e.g., genes) as
#'   rows and observations (e.g., cells) as columns.
#' 
#' @param condition Factor with N entries and R levels in which
#'   condition[i] specifies the condition status of the ith observation.
#' 
#' @param si Vector of size factors for each observation.
#' 
#' @param subgroup Factor with N entries and M levels giving the
#'   subgroup status of each of R conditions. Default is no subgroups
#'   (that is, M = 1 level). If not set to \code{NULL}, the names of
#'   subgroup must match the levels of \code{condition}.
#' 
#' @return A \dQuote{pois.mash} data object for poisson mash analysis,
#'   including the following components:
#' 
#' \item{X}{J x R matrix of count data collapsed over conditions, with
#'   features as rows and conditions as columns.}
#' 
#' \item{s}{Vector of length R adjusting for sequencing depth of
#'   each of R conditions.}
#' 
#' \item{subgroup}{Vector with M levels specifying the subgroup status
#'   of each of R conditions.}
#'
#' @export
#' 
pois_mash_set_data <- function (Y, condition, si, subgroup) {
  J <- nrow(Y)
  N <- ncol(Y)
  if (N != length(condition))
    stop("The number of columns of Y and the length of condition do not ",
         "match!")
  if (N != length(si))
    stop("The number of columns of Y and the length of si do not match!")
  trts <- sort(unique(condition))
  R    <- length(trts)
  
  if (!is.missing(subgroup)) {
    subgroup <- subgroup[order(names(subgroup))]
    if (sum(trts != names(subgroup)) > 0)
      stop("The levels of condition and the names of subgroup do not match!")
  } else {
    subgroup        <- rep(1,R)
    names(subgroup) <- trts
  }

  # Aggregate the cell-level data into condition-level data.
  X <- matrix(as.numeric(NA),J,R)
  s <- rep(as.numeric(NA),R)
  rownames(X) <- rownames(Y)
  colnames(X) <- trts
  names(s)    <- trts
  
  for (r in 1:R) {
    i     <- which(condition == trts[r])
    X[,r] <- rowSums(Y[,i])
    s[r]  <- sum(si[i])
  }
  
  data <- list(X = X,s = s/min(s),subgroup = subgroup)
  class(data) <- c("pois.mash","list")
  return(data)
}
