#' @title Compute a list of canonical covariance matrices modeling condition-specific effects
#' 
#' @param data A pois.mash data object containing the following
#'   components: X, an J x R matrix of count data collapsed over
#'   conditions, with features as rows and conditions as columns; s, an
#'   R x 1 numeric vector adjusting for sequencing depth of each of R
#'   conditions; and subgroup, an R x 1 factor vector with M levels
#'   denoting the subgroup status of each of R conditions.
#' 
#' @return a list of numerical vectors each of which forming a rank-1
#' canonical covariance matrix
#'
pois_cov_canonical <- function(data){
  R <- length(data$s)
  IR <- diag(R)
  ulist <- list(NULL)
  for(r in 1:R){
    ulist[[r]] <- as.numeric(IR[,r])
  }
  names(ulist) <- paste0("e_", c(1:R))
  return(ulist)
}
