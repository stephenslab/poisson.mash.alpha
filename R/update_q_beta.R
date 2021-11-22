# Update q(beta) for a given unit and prior covariance w*U, where U
# has full rank.
update_q_beta_general <- function (theta_m, theta_V, c2, psi2, w = 1, U) {
    
  # Make sure eigenvalues of U are all strictly positive.
  U       <- w*U
  eig.U   <- eigen(U)
  eig.val <- pmax(eig.U$values,1e-8)
  U       <- tcrossprod(eig.U$vectors %*% diag(sqrt(eig.val)))
  S_inv   <- 1/(psi2*c2)
  tmp     <- update_V(U, S_inv)
  beta_m  <- tmp %*% (theta_m * S_inv)
  beta_V  <- tmp + tmp %*% t(t(theta_V*S_inv)*S_inv) %*% tmp
  beta2_m <- beta_m %*% t(beta_m) + beta_V
  return(list(beta_m = beta_m,beta_V = beta_V,beta2_m = beta2_m))
}

# Update q(a,theta) for a given unit and prior covariance w*U, where
# U = u*u'
update_q_beta_rank1 <- function (theta_m, theta_V, c2, psi2, w = 1, u) {
  S_inv     <- 1/(psi2*c2)
  tmp       <- (u*S_inv)/(sum(u^2*S_inv) + 1/w)
  a_m       <- sum(tmp*theta_m)
  a_sigma2  <- 1/(sum(u^2*S_inv)+1/w) + as.numeric(t(tmp) %*% theta_V %*% tmp)
  a2_m      <- a_m^2 + a_sigma2
  a_theta_m <- (theta_m%*%t(theta_m) + theta_V) %*% tmp
  return(list(a_m = a_m,a_sigma2 = a_sigma2,a2_m = a2_m,a_theta_m = a_theta_m))
}

# Update q(beta), for a given unit and prior covariance w*U, where
# U = u*u'+epsilon2
update_q_beta_rank1_robust <- function (theta_m, theta_V, c2, psi2, w = 1, u,
                                        epsilon2) {
  tmp1    <- mat_inv_rank1(1/(w*epsilon2) + 1/(psi2*c2),-u/w/epsilon2,
                           u/epsilon2/(1 + sum(u^2)/epsilon2))
  tmp2    <- mat_inv_rank1(1 + psi2*c2/(w*epsilon2),-psi2*c2*u/w/epsilon2,
                           u/epsilon2/(1 + sum(u^2)/epsilon2))
  beta_m  <- tmp2 %*% theta_m
  beta_V  <- tmp1 + tmp2 %*% theta_V %*% tmp2
  beta2_m <- beta_m %*% t(beta_m) + beta_V
  return(list(beta_m = beta_m,beta_V = beta_V,beta2_m = beta2_m))
}
