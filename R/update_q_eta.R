# Update q(theta) = q(eta) for a given unit without fixed effects,
# i.e., with prior covariance U = 0.
update_q_eta_only <- function (x, s, mu, bias, c2, psi2, init = list(),
                               maxiter = 25, tol = 0.01, lwr = -10, upr = 10) {
  R <- length(x)
  m <- init$m
  V <- init$V
  Utilde <- psi2 * c2
  if (is.null(m))
    m <- rep(0, R)
  if(is.null(V))
    V <- Utilde
  
  bias <- drop(bias)
  a    <- compute_poisson_rates(s,mu,bias,m,V)

  for (iter in 1:maxiter) {
    V_new <- 1/(1/Utilde + a)
    a     <- compute_poisson_rates(s,mu,bias,m,V_new)
    m_new <- m - V_new*(a - x + m/Utilde)
    
    # Make sure the updated posterior mean is not unreasonably large
    # or small.
    m_new[m_new < lwr] <- lwr
    m_new[m_new > upr] <- upr
    
    # Decide whether to stop based on change in mu and V.
    m_tmp <- m
    m_tmp[abs(m_tmp) < 1e-15] <- 1e-15
    V_tmp <- V
    V_tmp[abs(V_tmp) < 1e-15] <- 1e-15
    idx.mu <- (max(abs(m_new-m)) < tol/100 | max(abs(m_new/m_tmp-1)) < tol)
    idx.V  <- (max(abs(V_new-V)) < tol/100 | max(abs(V_new/V_tmp-1)) < tol)
    
    if (idx.mu & idx.V)
      break
    
    m <- m_new
    V <- V_new
    a <- compute_poisson_rates(s,mu,bias,m,V)
  }
  
  # Calculate "local" ELBO F_j.
  ELBO <- sum(x*(mu + bias + m)) - sum(a) -
          0.5*(sum(V/Utilde) + sum(m^2/Utilde) - R
               + sum(log(Utilde)) - sum(log(V)))
  return(list(m = m,V = V,a = a,ELBO = ELBO))
}

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
