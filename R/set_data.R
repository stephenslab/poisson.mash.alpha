#' @title Create Data Oject for Poisson Mash Analysis
#'
#' @description Prepare count data for an analysis with Poisson mash.
#' 
#' @param Y J x N matrix of counts (non-negative numbers) with
#'   features (e.g., genes) as rows and observations (e.g., cells) as
#'   columns. Y should have at least 2 rows and 2 columns.
#' 
#' @param condition Factor of length N wiith R levels giving the
#'   condition for each observation.
#' 
#' @param si Numeric vector of length N containing the size factors for
#'   the N observations.
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

  # Check and process input "Y".
  if (!(((is.matrix(Y) & is.numeric(Y)) | inherits(Y,"dgCMatrix")) &
        all(Y >= 0) & all(is.finite(Y)) & all(!is.na(Y))))
    stop("Input argument Y should be a non-negative numeric matrix ",
         "(a \"matrix\" or a \"dgCMatrix\") in which all elements are ",
         "finite, and not NA")
  if (nrow(Y) < 2 | ncol(Y) < 2)
    stop("Input matrix Y should have at least 2 rows and at least 2 columns")
  if (is.null(rownames(Y)))
    rownames(Y) <- paste0("j",1:row(Y))
  if (is.null(colnames(Y)))
    colnames(Y) <- paste0("r",1:ncol(Y))

  # Check and process input "condition".
  J    <- nrow(Y)
  R    <- nlevels(condition)
  trts <- levels(condition)
  if (!(is.factor(condition) &
        length(condition) == ncol(Y) &
        nlevels(condition) > 1 &
        all(table(condition) > 0)))
    stop("Input argument \"condition\" should be a factor with one entry ",
         "per column of Y, at least two levels, and all levels should occur ",
         "at least once")
  if (is.null(names(condition)))
    names(condition) <- colnames(Y)
  else if (colnames(Y) != names(condition))
    stop("names(condition) do not match colnames(Y)")
      
  # Check and process input "si".
  if (!(is.numeric(si) & all(si > 0) & ncol(Y) == length(si)))
    stop("Input \"si\" should be a numeric vector with one entry for each ",
         "column of Y")

  # Check and process input "subgroup".
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
