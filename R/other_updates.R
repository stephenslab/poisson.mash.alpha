# Update the J x R matrix of means mu, where X is J x R matrix of
# counts, subgroup is factor with R entries and M levels, zeta is J x K
# matrix of posterior weights, and tmp.mu is J x K x M array.
update_mu <- function (X, subgroup, zeta, tmp.mu) {
  M  <- length(unique(subgroup))
  mu <- matrix(as.numeric(NA),nrow(X),ncol(X))
  for (i in 1:M) {
    k        <- which(subgroup == i)
    mu.i.new <- log(rowSums(X[,k])) - log(rowSums(zeta * tmp.mu[,,i]))
    mu[,k]   <- mu.i.new
  }
  return(mu)
}

# Update the vector of dispersion parameter psi2 of length K, where
# zeta is J x K matrix of posterior weights, tmp.psi2 is J x K matrix,
# R is the number of conditions, minpsi2 and maxpsi2 are respectively
# positive scalars giving the lower and upper bound of psi2.
update_psi2 <- function (zeta, tmp.psi2, R, minpsi2, maxpsi2) {
  psi2.new <- rowSums(zeta * tmp.psi2)/R
  return(pmin(pmax(psi2.new,minpsi2),maxpsi2))
}

# Update the vector of prior weights pi of length K, where zeta is J x
# K matrix of posterior weights.
update_pi <- function (zeta) {
  pi <- colMeans(zeta)
  return(pmax(pi,1e-8))
}

# Update the J x K matrix of posterior weights zeta, where ELBOs is J
# x K matrix of local ELBO, pi is the vector of prior weights of length K.
update_zeta <- function (ELBOs, pi) {
  ELBOs.cen <- ELBOs - apply(ELBOs,1,max)
  zeta <- scale.cols(exp(ELBOs.cen),pi)
  zeta <- normalize.rows(zeta)
  return(pmax(zeta,1e-15))
}

# Update the J x R matrix Ruv needed to update rho:
# Ruv[j,r] = sum_k zeta[j,k] * exp(A[j,k,r]).
update_ruv <- function (zeta, A) {
  J <- nrow(zeta)
  R <- dim(A)[3]
  Ruv <- matrix(as.numeric(NA),J,R)
  for (r in 1:R)
    Ruv[,r] <- rowSums(zeta * exp(A[,,r]))
  return(Ruv)
}

# Calculate solve(solve(U) + diag(a)) in a numerically more stable manner, where a is R x 1 vector
update_V <- function (U, a)
  return(U %*% solve(U*a + diag(nrow=nrow(U)), tol = 1e-50))

update_m <- function (U, V, a, x, m)
  drop(m - V %*% (a - x + solve(U,m)))
