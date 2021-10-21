#' @title Get List of Canonical Covariance Matrices
#'
#' @description Compute a list of canonical covariance matrices
#'   modeling condition-specific effects.
#'
#' @param data \dQuote{pois.mash} data object, typically created by
#'   calling \code{\link{pois_mash_set_data}}.
#' 
#' @return A list of numerical vectors each of which forms a rank-1
#'   canonical covariance matrix.
#'
#' @export
#' 
pois_cov_canonical <- function (data) {
  R <- length(data$s)
  IR <- diag(R)
  ulist <- vector("list",R)
  for (r in 1:R)
    ulist[[r]] <- drop(IR[,r])
  names(ulist) <- paste0("e_",c(1:R))
  return(ulist)
}
