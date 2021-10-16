# Update the D x R matrix of rho for all conditions given current rho,
# where X is a J x R matrix of counts, s is a vector of sequencing
# depths of length R, mu is J x R matrix of means, Fuv is a J x D
# matrix of latent factors causing unwanted variation, rho is a D x R
# matrix of current rho, L is a J x R matrix, and tol is a small
# positive number to assess convergence.
update_rho_all <- function (X, s, mu, Fuv, rho, L, maxiter = 100,
                            tol = 1e-6, maxrho = 100/max(abs(Fuv)),
                            version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R")
    rho.new <- update_rho_all_r(X,s,mu,Fuv,rho,L,maxiter,tol,maxrho)
  else if (version == "Rcpp")
    rho.new <- update_rho_all_rcpp(X,Fuv,s,mu,L,rho,maxiter,tol,maxrho)
  return(rho.new)
}

# This implements update_rho_all with version = "R".
update_rho_all_r <- function (X, s, mu, Fuv, rho, L, maxiter, tol, maxrho) {
  R <- ncol(X)
  rho.new <- matrix(as.numeric(NA),nrow(rho),ncol(rho))
  for (r in 1:R)
    rho.new[,r] <- update_rho(X[,r],Fuv,s[r],mu[,r],L[,r],init = rho[,r], 
                              maxiter,tol,maxrho)$rho
  return(rho.new)
}

# Update vector rho (of length D) for a given condition r.
#
#' @importFrom utils modifyList
update_rho <- function (Xr, Fuv, sr, mu, Lr, init, maxiter, tol, maxrho) {
  rho <- init
  Fr  <- rep(as.numeric(NA),maxiter)
  D   <- length(rho)
  d1F <- rep(as.numeric(NA),D)
  d2F <- matrix(as.numeric(NA),D,D)
  
  for (iter in 1:maxiter) {
    bias <- Fuv %*% rho
    Fr[iter] <- sum(Xr * bias) - sr * sum(exp(mu + bias) * Lr)
    for (d in 1:D)
      d1F[d] <- sum(Xr * Fuv[,d]) - sr * sum(Fuv[,d] * Lr * exp(mu + bias))
    for (d in 1:D)
      for (t in 1:D)
        d2F[d,t] <- -sr * sum(Fuv[,d] * Fuv[,t] *  Lr * exp(mu + bias))
    rho_new <- rho - solve(d2F,d1F)
    rho_new <- pmin(pmax(rho_new,-maxrho),maxrho)
    if (max(abs(rho_new - rho)) < tol)
      break
    rho <- rho_new
  }
  
  return(list(rho = rho,Fr = Fr))
}
