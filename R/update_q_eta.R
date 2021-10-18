# Update q(eta) for a given unit and prior covariance w*U, where U has
# full rank.
update_q_eta_general <- function (theta_m, theta_V, c2, psi2, w=1, U) {
    
  # Make sure eigenvalues of U are all strictly positive..
  U       <- w*U
  eig.U   <- eigen(U)
  eig.val <- pmax(eig.U$values,1e-8)
  U       <- tcrossprod(eig.U$vectors %*% diag(sqrt(eig.val)))
  R       <- length(theta_m)
  S_inv   <- 1/(psi2*c2)
  tmp1    <- solve(solve(U) + diag(S_inv))
  tmp2    <- solve(diag(R) + t(t(U)*S_inv))
  eta2_m  <- diag(tmp1 + tmp2 %*% (theta_m %*% t(theta_m) + theta_V) %*%
                  t(tmp2))
  return(eta2_m)
}

# Update q(eta) for a given unit and prior covariance w*U, where U = u*u'.
update_q_eta_rank1 <- function (theta_m, theta_V, a2_m, a_theta_m, u)
  diag(theta_V) + theta_m^2 + a2_m*u^2 - 2*u*a_theta_m

# Uupdate q(eta) for a given unit and prior covariance w*U, where
# U = u*u' + epsilon2
update_q_eta_rank1_robust <- function (theta_m, theta_V, c2, psi2, w = 1,
                                       u, epsilon2) {
  tmp1   <- mat_inv_rank1(1/(w*epsilon2) + 1/(psi2*c2),-u/w/epsilon2,
                          u/epsilon2/(1 + sum(u^2)/epsilon2))
  tmp2   <- mat_inv_rank1(1 + w*epsilon2/(psi2*c2),w*u,u/(psi2*c2))
  eta2_m <- diag(tmp1 + tmp2 %*% (theta_m %*% t(theta_m) + theta_V) %*%
              t(tmp2))
  return(eta2_m)
}
