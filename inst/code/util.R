#' @title Utility functions to update model and variational parameters and speed up matrix inversions
#' 
#' @description These are internal functions which users might not want to call directly. 

### function to update q(theta)=q(eta) for a given unit without fixed effects, i.e., prior covariance U=0
update_q_eta_only <- function(x, s, mu, bias, c2, psi2, init=list(NULL), control=list(maxiter=25, tol=1e-2, lwr=-10, upr=10)){
  R <- length(x)
  maxiter <- control$maxiter
  tol <- control$tol
  lwr <- control$lwr
  upr <- control$upr
  m <- init$m
  V <- init$V
  
  Utilde <- psi2*c2
  
  if(is.null(maxiter))
  {
    maxiter <- 25
  }
  
  if(is.null(tol)){
    tol <- 1e-2
  }
  
  if(is.null(lwr)){
    lwr <- -10
  }
  
  if(is.null(upr)){
    upr <- 10
  }
  
  if(is.null(m)){
    m <- rep(0, R)
  }
  
  if(is.null(V)){
    V <- Utilde
  }
  
  bias <- as.numeric(bias)
  a <- s*exp(mu + bias + m + V/2)
  
  for(iter in 1:maxiter){
    V_new <- 1/(1/Utilde + a)
    a <- s*exp(mu + bias + m + V_new/2)
    m_new <- m - V_new*(a - x + m/Utilde)
    
    # make sure the updated posterior mean is not unreasonably large or small
    m_new[m_new < lwr] <- lwr
    m_new[m_new > upr] <- upr
    
    # decide whether to stop based on change in mu and V
    m_tmp <- m
    m_tmp[abs(m_tmp) < 1e-15] <- 1e-15
    V_tmp <- V
    V_tmp[abs(V_tmp) < 1e-15] <- 1e-15
    idx.mu <- (max(abs(m_new-m)) < tol/1e2 | max(abs(m_new/m_tmp-1)) < tol)
    idx.V <- (max(abs(V_new-V)) < tol/1e2 | max(abs(V_new/V_tmp-1)) < tol)
    
    if(idx.mu & idx.V) break
    
    m <- m_new
    V <- V_new
    a <- s*exp(mu + bias + m + V/2)
  }
  
  # calculate local ELBO F_j
  ELBO <- sum(x*(mu+bias+m)) - sum(a) - 0.5*(sum(V/Utilde) + sum(m^2/Utilde) - R + sum(log(Utilde)) - sum(log(V)))
  
  return(list(m=m, V=V, a=a, ELBO=ELBO))
}



### function to update q(theta), for a given unit and prior covariance w*U, where U has full rank
update_q_theta_general <- function(x, s, mu, bias, c2, psi2, w=1, U, init=list(NULL), control=list(maxiter=25, tol=1e-2, lwr=-10, upr=10)){
  R <- length(x)
  maxiter <- control$maxiter
  tol <- control$tol
  lwr <- control$lwr
  upr <- control$upr
  m <- init$m
  V <- init$V
  
  Utilde <- w*U + psi2*diag(c2) 
  
  if(is.null(maxiter))
  {
    maxiter <- 25
  }
  
  if(is.null(tol)){
    tol <- 1e-2
  }
  
  if(is.null(lwr)){
    lwr <- -10
  }
  
  if(is.null(upr)){
    upr <- 10
  }
  
  if(is.null(m)){
    m <- rep(0, R)
  }
  
  if(is.null(V)){
    V <- Utilde
  }
  
  bias <- as.numeric(bias)
  a <- s*exp(mu + bias + m + diag(V)/2)
  
  for(iter in 1:maxiter){
    V_new <- solve(solve(Utilde) + diag(a), tol=1e-50)
    a <- s*exp(mu + bias + m + diag(V_new)/2)
    m_new <- as.numeric(m - V_new %*% (a - x + solve(Utilde, m)))
    
    # make sure the updated posterior mean is not unreasonably large or small
    m_new[m_new < lwr] <- lwr
    m_new[m_new > upr] <- upr

    # decide whether to stop based on change in mu and V
    m_tmp <- m
    m_tmp[abs(m_tmp) < 1e-15] <- 1e-15
    V_tmp <- V
    V_tmp[abs(V_tmp) < 1e-15] <- 1e-15
    idx.mu <- (max(abs(m_new-m)) < tol/1e2 | max(abs(m_new/m_tmp-1)) < tol)
    idx.V <- (max(abs(V_new-V)) < tol/1e2 | max(abs(V_new/V_tmp-1)) < tol)
  
    if(idx.mu & idx.V) break
     
    m <- m_new
    V <- V_new
    a <- s*exp(mu + bias + m + diag(V)/2)
  }
  
  # calculate local ELBO F_jhl
  ELBO <- sum(x*(mu+bias+m)) - sum(a) - 0.5*(tr(solve(Utilde, V)) + as.numeric(t(m)%*%solve(Utilde, m)) - R + 
                                               determinant(Utilde, log=TRUE)$modulus - determinant(V, log=TRUE)$modulus)   
  
  return(list(Utilde=Utilde, m=m, V=V, a=a, ELBO=ELBO))
}



### function to update q(theta), for a given unit and prior covariance w*U, where U = u*u'
update_q_theta_rank1 <- function(x, s, mu, bias, c2, psi2, w=1, u, init=list(NULL), control=list(maxiter=25, tol=1e-2, lwr=-10, upr=10)){
  R <- length(x)
  maxiter <- control$maxiter
  tol <- control$tol
  lwr <- control$lwr
  upr <- control$upr
  m <- init$m
  V <- init$V
  
  Utilde <- w*u%*%t(u) + psi2*diag(c2) 
  S_inv <- 1/(psi2*c2)
  
  if(is.null(maxiter))
  {
    maxiter <- 25
  }
  
  if(is.null(tol)){
    tol <- 1e-2
  }
  
  if(is.null(lwr)){
    lwr <- -10
  }
  
  if(is.null(upr)){
    upr <- 10
  }
  
  if(is.null(m)){
    m <- rep(0, R)
  }
  
  if(is.null(V)){
    V <- Utilde
  }
  
  bias <- as.numeric(bias)
  a <- s*exp(mu + bias + m + diag(V)/2)
  
  for(iter in 1:maxiter){
    V_new <- mat_inv_rank1(v1=a+S_inv, v2=-w*u*S_inv, v3=(u*S_inv)/(1+w*sum(u^2*S_inv)))
    
    a <- s*exp(mu + bias + m + diag(V_new)/2)
    Utilde_inv_m <- m*S_inv - w*u*S_inv*sum(u*m*S_inv)/(1+w*sum(u^2*S_inv))
    m_new <- as.numeric(m - V_new %*% (a - x + Utilde_inv_m))
    
    # make sure the updated posterior mean is not unreasonably large or small
    m_new[m_new < lwr] <- lwr
    m_new[m_new > upr] <- upr  
    
    # decide whether to stop based on change in mu and V
    m_tmp <- m
    m_tmp[abs(m_tmp) < 1e-15] <- 1e-15
    V_tmp <- V
    V_tmp[abs(V_tmp) < 1e-15] <- 1e-15
    idx.mu <- (max(abs(m_new-m)) < tol/1e2 | max(abs(m_new/m_tmp-1)) < tol)
    idx.V <- (max(abs(V_new-V)) < tol/1e2 | max(abs(V_new/V_tmp-1)) < tol)
    
    if(idx.mu & idx.V) break
    
    m <- m_new
    V <- V_new
    a <- s*exp(mu + bias + m + diag(V)/2)
  }
  
  # calculate local ELBO F_jhl
  ELBO <- sum(x*(mu+bias+m)) - sum(a) - 0.5*(tr(solve(Utilde, V)) + as.numeric(t(m)%*%solve(Utilde, m)) - R + 
                                               determinant(Utilde, log=TRUE)$modulus - determinant(V, log=TRUE)$modulus) 
  
  return(list(Utilde=Utilde, m=m, V=V, a=a, ELBO=ELBO))
}



### function to update q(beta), for a given unit and prior covariance w*U, where U has full rank
update_q_beta_general <- function(theta_m, theta_V, c2, psi2, w=1, U){
  # make sure eigenvalues of U are all strictly positive
  U <- w*U
  eig.U <- eigen(U)
  eig.val <- pmax(eig.U$values, 1e-8)
  U <- tcrossprod(eig.U$vectors %*% diag(sqrt(eig.val)))
  
  S_inv <- 1/(psi2*c2)
  tmp <- solve(solve(U) + diag(S_inv))
  beta_m <- tmp%*%(theta_m*S_inv)
  beta_V <- tmp + tmp%*%t(t(theta_V*S_inv)*S_inv)%*%tmp
  beta2_m <- beta_m%*%t(beta_m) + beta_V
  return(list(beta_m=beta_m, beta_V=beta_V, beta2_m=beta2_m))
}



### function to update q(a, theta), for a given unit and prior covariance w*U, where U = u*u'
update_q_beta_rank1 <- function(theta_m, theta_V, c2, psi2, w=1, u){
  S_inv <- 1/(psi2*c2)
  tmp <- (u*S_inv)/(sum(u^2*S_inv) + 1/w)
  a_m <- sum(tmp*theta_m)
  a_sigma2 <- 1/(sum(u^2*S_inv) + 1/w) + as.numeric(t(tmp)%*%theta_V%*%tmp)
  a2_m <- a_m^2 + a_sigma2
  a_theta_m <- (theta_m%*%t(theta_m) + theta_V)%*%tmp
  return(list(a_m=a_m, a_sigma2=a_sigma2, a2_m=a2_m, a_theta_m=a_theta_m))
}



### function to update q(beta), for a given unit and prior covariance w*U, where U=u*u'+epsilon2
update_q_beta_rank1_robust <- function(theta_m, theta_V, c2, psi2, w=1, u, epsilon2){
  tmp1 <- mat_inv_rank1(v1=1/(w*epsilon2)+1/(psi2*c2), v2=-u/w/epsilon2, v3=u/epsilon2/(1+sum(u^2)/epsilon2))
  tmp2 <- mat_inv_rank1(v1=1+psi2*c2/(w*epsilon2), v2=-psi2*c2*u/w/epsilon2, v3=u/epsilon2/(1+sum(u^2)/epsilon2))
  beta_m <- tmp2%*%theta_m
  beta_V <- tmp1 + tmp2%*%theta_V%*%tmp2
  beta2_m <- beta_m%*%t(beta_m) + beta_V
  return(list(beta_m=beta_m, beta_V=beta_V, beta2_m=beta2_m))
}



### function to update q(eta), for a given unit and prior covariance w*U, where U has full rank
update_q_eta_general <- function(theta_m, theta_V, c2, psi2, w=1, U){
  # make sure eigenvalues of U are all strictly positive
  U <- w*U
  eig.U <- eigen(U)
  eig.val <- pmax(eig.U$values, 1e-8)
  U <- tcrossprod(eig.U$vectors %*% diag(sqrt(eig.val)))
  
  R <- length(theta_m)
  S_inv <- 1/(psi2*c2)
  tmp1 <- solve(solve(U) + diag(S_inv))
  tmp2 <- solve(diag(R) + t(t(U)*S_inv))
  eta2_m <- diag(tmp1 + tmp2%*%(theta_m%*%t(theta_m)+theta_V)%*%t(tmp2))
  return(eta2_m)
}



### function to update q(eta), for a given unit and prior covariance w*U, where U = u*u'
update_q_eta_rank1 <- function(theta_m, theta_V, a2_m, a_theta_m, u){
  eta2_m <- diag(theta_V) + theta_m^2 + a2_m*u^2 - 2*u*a_theta_m
  return(eta2_m)
}
  


### function to update q(eta), for a given unit and prior covariance w*U, where U=u*u'+epsilon2
update_q_eta_rank1_robust <- function(theta_m, theta_V, c2, psi2, w=1, u, epsilon2){
  tmp1 <- mat_inv_rank1(v1=1/(w*epsilon2)+1/(psi2*c2), v2=-u/w/epsilon2, v3=u/epsilon2/(1+sum(u^2)/epsilon2))
  tmp2 <- mat_inv_rank1(v1=1+w*epsilon2/(psi2*c2), v2=w*u, v3=u/(psi2*c2))
  eta2_m <- diag(tmp1 + tmp2%*%(theta_m%*%t(theta_m)+theta_V)%*%t(tmp2))
  return(eta2_m)
}



### function to take inverse of the matrix diag(v1) + v2%*%t(v3), where v1, v2, v3 are vectors
mat_inv_rank1 <- function(v1, v2, v3){
  V_inv <- diag(1/v1) - (v2/v1)%*%t(v3/v1)/(1+sum(v2*v3/v1))
  return(V_inv)
}



### function to update the D x 1 vector rho for a given condition r
update_rho <- function(Xr, Fuv, sr, mu, Lr, init, control=list(maxiter=100, tol=1e-8, maxrho=1e2)){
  maxiter <- control$maxiter
  tol <- control$tol
  maxrho <- control$maxrho
  
  if(is.null(maxiter))
  {
    maxiter <- 100
  }
  
  if(is.null(tol)){
    tol <- 1e-8
  }
  
  if(is.null(maxrho)){
    maxrho <- 1e2
  }
  
  rho <- init
  Fr <- rep(NA, maxiter)
  
  D <- length(rho)
  d1F <- rep(NA, D)
  d2F <- matrix(NA, D, D)
  
  for(iter in 1:maxiter){
    bias <- Fuv %*% rho
    Fr[iter] <- sum(Xr*bias) - sr*sum(exp(mu+bias)*Lr)
    
    for(d in 1:D){
      d1F[d] <- sum(Xr*Fuv[,d]) - sr*sum(Fuv[,d]*Lr*exp(mu+bias))
    }
    
    for(d in 1:D){
      for(t in 1:D){
        d2F[d,t] <- -sr*sum(Fuv[,d]*Fuv[,t]*Lr*exp(mu+bias))
      }
    }
     
    rho_new <- rho - solve(d2F, d1F)
    rho_new <- pmin(pmax(rho_new, -maxrho), maxrho)
    
    if(max(abs(rho_new - rho)) < tol) break
    rho <- rho_new
  }
  
  return(list(rho=rho, Fr=Fr))
}



### function to scale the J x R matrix of bias to avoid too large values in bias
scale_bias <- function(bias, maxbias){
  range.bias <- apply(bias, 1, function(x){max(x)-min(x)})
  idx.bias <- which(range.bias > maxbias)
  
  if(length(idx.bias) > 0){
    bias[idx.bias,] <- maxbias*bias[idx.bias,]/range.bias[idx.bias]
  }
  
  return(bias)
}