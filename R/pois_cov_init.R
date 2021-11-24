#' @title Initialize Data-driven Prior Covariance Matrices using PCA
#'
#' @description Initialize data-driven prior covariance matrices based
#'   on principal component analysis.
#' 
#' @param data \dQuote{pois.mash} data object, typically created by
#'   calling \code{\link{pois_mash_set_data}}.
#' 
#' @param ruv Logical scalar indicating whether to account for
#'   unwanted variation. Default is \code{FALSE}. If \code{ruv = TRUE},
#'   \code{Fuv} and \code{rho} must be provided.
#' 
#' @param Fuv J x D matrix of latent factors causing unwanted
#'   variation, with features as rows and latent factors as columns.
#' 
#' @param rho D x R matrix of effects corresponding to unwanted
#'   variation, such that \code{bias = Fuv \%*\% rho}.
#' 
#' @param prop The proportion by which to take a random subset of
#'   genes for prior covariance estimation (useful in case of many
#'   genes).
#' 
#' @param seed Useful for reproducibility when \code{prop} is less
#'   than 1.
#' 
#' @param npc The number of principal components to use.
#' 
#' @param cutoff The threshold for the maximum of absolute values of
#'    Z-scores taken across conditions to include as "strong" features
#'   used for prior covariance estimation.
#' 
#' @return A list with initial estimates of prior covariances, and
#'   indices of the features (j = 1,...,J) to include in the subsequent
#'   ED step.
#'
#' @importFrom stats sd
#' @importFrom stats cov
#' @importFrom stats pbinom
#' @importFrom stats dbinom
#' @importFrom stats qnorm
#' @importFrom stats chisq.test
#' 
#' @export
#' 
pois_cov_init <- function (data, ruv = FALSE, Fuv = NULL, rho = NULL, prop = 1,
                           seed = 1, npc = 5, cutoff = 3) {
  s        <- data$s
  subgroup <- data$subgroup
  data     <- as.matrix(data$X)
  
  # Take a random subset of genes.
  set.seed(seed)
  idx.select <- sort(sample(nrow(data),round(nrow(data) * prop)))
  data.test  <- data[idx.select,]
  J          <- nrow(data.test)
  R          <- ncol(data.test)
  M          <- length(unique(subgroup))
  subgroup   <- as.numeric(as.factor(subgroup))
  
  if (ruv) {
    Fuv  <- as.matrix(Fuv[match(rownames(data.test),rownames(data)),])
    bias <- Fuv %*% as.matrix(rho)
  }
  else
    bias <- matrix(0,J,R)
  
  # Get Z-scores for observing each X_jr assuming no differential
  # expression within each subgroup while accounting for unwanted
  # variation. Z-scores are obtained by inverting exact cdf.
  Z.score           <- matrix(as.numeric(NA),J,R)
  rownames(Z.score) <- rownames(data.test)
  colnames(Z.score) <- colnames(data.test)
  for (j in 1:J)
    for (i in 1:M) {
      idx.i <- which(subgroup == i)
      x     <- data.test[j,idx.i]
      s.j   <- s[idx.i] * exp(bias[j,idx.i])
      Z.score[j,idx.i] <- qnorm(pbinom(q=x,size=sum(x),prob=s.j/sum(s.j)) -
                                0.5*dbinom(x=x,size=sum(x),prob=s.j/sum(s.j)))
    }
  
  # Residuals obtained by normal approximation.
  chisq.test.res <- matrix(as.numeric(NA),J,R)
  for (j in 1:J)
    for (i in 1:M) {
      idx.i <- which(subgroup == i)
      x     <- data.test[j,idx.i]
      s.j   <- s[idx.i] * exp(bias[j,idx.i])
      chisq.test.res[j,idx.i] <-
        suppressWarnings(chisq.test(x,p = s.j/sum(s.j)))$residuals
    }
  
  # Replace inf Z-scores with those calculated based on normal
  # approximation.
  idx.inf <- which(apply(Z.score,1,function (x) sum(is.infinite(x))) > 0)
  Z.score[idx.inf,] <- chisq.test.res[idx.inf,]
  
  # Get the genes that are differentially expressed across conditions.
  idx.sig <- which(apply(abs(Z.score),1,max) >= cutoff)
  if (length(idx.sig) < 2)
    stop("Number of genes satisfying cutoff is less than 2; consider ",
         "increasing cutoff")
  
  # Get the Z-scores for these genes.
  Z.sig <- Z.score[idx.sig,]
  Z.cen <- scale(Z.sig,center = TRUE,scale = FALSE)
  
  # Get the empirical covariance matrix of Z.cen.
  Z.cov <- cov(Z.cen)
  
  # Get the initial data-driven covariances based on SVD.
  res.svd <- svd(Z.cen,nv = npc,nu = npc)
  f.svd   <- res.svd$v
  d.svd   <- res.svd$d[1:npc]
  
  # Create list of rank-1 data-driven covariance matrices U_g =
  # u_g*u'_g based on SVD
  ulist <- vector("list",length(d.svd))
  for(i in 1:length(d.svd))
    ulist[[i]] <- d.svd[i] * f.svd[,i] / sqrt(nrow(Z.cen))
  names(ulist) <- paste0("PC_",1:length(d.svd))
  
  # List of data-driven covariance matrices U_h that might have a rank
  # larger than 1.
  Ulist <- list(tPCA    = f.svd %*% diag(d.svd^2) %*% t(f.svd)/nrow(Z.cen),
                Emp_cov = Z.cov)
  
  # Avoid too large values in Ulist and ulist.
  upr_bd <- estimate_psi2_range(data,s)$upr_bd
  for (h in 1:length(Ulist)) {
    Uh <- Ulist[[h]]
    if (max(diag(Uh)) > upr_bd)
      Uh <- upr_bd * Uh / max(diag(Uh))
    Ulist[[h]] <- Uh
  }
  for (g in 1:length(ulist)) {
    ug <- ulist[[g]]
    if(max(abs(ug)) > sqrt(upr_bd))
      ug <- sqrt(upr_bd) * ug / max(abs(ug))
    ulist[[g]] <- ug
  }
  
  # Add in the zero vector.
  G                 <- length(ulist)
  ulist[[G+1]]      <- rep(0,R)
  names(ulist)[G+1] <- "u_0"
  
  # Get the index of the genes used to compute data-driven covariances
  # in the original dataset.
  idx <- which(rownames(data) %in% rownames(Z.sig))
  return(list(ulist  = ulist,
              Ulist  = Ulist,
              subset = idx))
}
