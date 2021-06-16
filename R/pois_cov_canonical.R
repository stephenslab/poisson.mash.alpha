#' @title Generate Condition-specific Canonical Covariance Matrices
#' 
#' @param data A Poisson mash data object, typically created using
#'   \code{\link{pois_mash_set_data}}.
#' 
#' @return A list ii which each list element is a vector \code{u}
#' defining a rank-1 covariance matrix \code{tcrossprod(u)}
#'
#' @export
#' 
pois_cov_canonical <- function (dat) {
  R <- ncol(dat$X)
  ulist <- as.list(as.data.frame(diag(R)))
  names(ulist) <- paste("e", 1:R, sep="_")
  return(ulist)
}
