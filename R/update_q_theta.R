# Update q(theta) for a given unit and for all prior mixture
# components w_l * U_k.
update_q_theta_all <- function (x, s, mu, bias, c2 = rep(1,length(x)),
                                psi2, wlist = 1, Ulist, ulist, init = list(),
                                maxiter.q = 25, tol.q = 0.01) {
  R <- length(x)
  H <- length(Ulist)
  G <- length(ulist)
  L <- length(wlist)
  K <- L*(H+G)
  
  gamma <- matrix(as.numeric(NA),K,R)
  Sigma <- array(as.numeric(NA),c(K,R,R))
  A     <- matrix(as.numeric(NA),K,R)
  ELBOs <- rep(as.numeric(NA),K)
  
  if (H > 0) {
    hl <- 0
    for (h in 1:H) {
      for(l in 1:L) {
        hl <- hl + 1
        if (!is.null(init$gamma) & !is.null(init$Sigma)) {
          init.m <- init$gamma[hl,]
          init.V <- init$Sigma[hl,,]
        }
        else {
          init.m <- NULL
          init.V <- NULL
        }
        out <-
          update_q_theta_general(x,s,mu,bias,c2,psi2,wlist[l],Ulist[[h]],
                                 list(m = init.m,V = init.V),
                                 maxiter = maxiter.q,tol = tol.q)
        gamma[hl,]  <- out$m
        Sigma[hl,,] <- out$V
        ELBOs[hl]   <- out$ELBO
        A[hl,]      <- with(out,m + diag(V)/2)
      }
    }
  }
  
  gl <- 0
  for (g in 1:G) {
    for (l in 1:L) {
      gl <- gl + 1
      if (!is.null(init$gamma) & !is.null(init$Sigma)) {
        init.m <- init$gamma[H*L+gl,]
        init.V <- init$Sigma[H*L+gl,,]
      }
      else {
        init.m <- NULL
        init.V <- NULL
      }
      out <- update_q_theta_rank1(x,s,mu,bias,c2,psi2,wlist[l],ulist[[g]], 
                                  list(m = init.m,V = init.V),
                                  maxiter = maxiter.q,tol = tol.q)
      gamma[H*L+gl,]  <- out$m
      Sigma[H*L+gl,,] <- out$V
      ELBOs[H*L+gl]   <- out$ELBO
      A[H*L+gl,]      <- with(out,m + diag(V)/2)
    }
  }
  
  return(list(ELBOs = ELBOs,A = A,gamma = gamma,Sigma = Sigma))
}

# Update q(theta) for a given unit and prior covariance w*U, where U
# has full rank.
update_q_theta_general <- function (x, s, mu, bias, c2, psi2, w = 1, U,
                                    init = list(), maxiter = 25, tol = 0.01,
                                    lwr = -10, upr = 10) {
  R       <- length(x)
  m       <- init$m
  V       <- init$V
  Utilde  <- w*U + psi2*diag(c2) 
  bias    <- drop(bias)
  if (is.null(m))
    m <- rep(0,R)
  if (is.null(V))
    V <- Utilde
  
  a <- compute_poisson_rates(s,mu,bias,m,diag(V))
  
  for (iter in 1:maxiter) {
    V_new <- update_V(Utilde,a)
    a     <- compute_poisson_rates(s,mu,bias,m,diag(V_new))
    m_new <- update_m(Utilde,V_new,a,x,m)
    
    # Make sure the updated posterior mean is not unreasonably large
    # or small.
    m_new[m_new < lwr] <- lwr
    m_new[m_new > upr] <- upr

    # Decide whether to stop based on change in mu and V.
    m_tmp <- m
    m_tmp[abs(m_tmp) < 1e-15] <- 1e-15
    V_tmp <- V
    V_tmp[abs(V_tmp) < 1e-15] <- 1e-15
    idx.mu <- (max(abs(m_new - m)) < tol/100 | max(abs(m_new/m_tmp - 1)) < tol)
    idx.V  <- (max(abs(V_new - V)) < tol/100 | max(abs(V_new/V_tmp - 1)) < tol)
    if (idx.mu & idx.V)
      break
    
    m <- m_new
    V <- V_new
    a <- compute_poisson_rates(s,mu,bias,m,diag(V))
  }
  
  # Calculate "local" ELBO F_jhl.
  ELBO <- sum(x*(mu + bias + m)) - sum(a) -
          0.5*(tr(solve(Utilde,V)) + drop(t(m) %*% solve(Utilde,m)) 
               - R + logdet(Utilde) - logdet(V))
  return(list(Utilde = Utilde,m = m,V = V,a = a,ELBO = ELBO))
}

# Update q(theta) for a given unit and prior covariance w*U, where U = u*u'.
update_q_theta_rank1 <- function (x, s, mu, bias, c2, psi2, w = 1, u,
                                  init = list(), maxiter = 25, tol = 0.01,
                                  lwr = -10, upr = 10) {
  R <- length(x)
  m <- init$m
  V <- init$V
  
  Utilde <- w*u %*% t(u) + psi2 * diag(c2) 
  S_inv  <- 1/(psi2 * c2)
  
  if (is.null(m))
    m <- rep(0,R)
  if (is.null(V))
    V <- Utilde
  
  bias <- drop(bias)
  a    <- compute_poisson_rates(s,mu,bias,m,diag(V))
  
  for (iter in 1:maxiter) {
    V_new <- mat_inv_rank1(a+S_inv,-w*u*S_inv,(u*S_inv)/(1+w*sum(u^2*S_inv)))
    a     <- compute_poisson_rates(s,mu,bias,m,diag(V_new))
    Utilde_inv_m <- m*S_inv - w*u*S_inv*sum(u*m*S_inv)/(1 + w*sum(u^2*S_inv))
    m_new        <- drop(m - V_new %*% (a - x + Utilde_inv_m))
    
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
    a <- compute_poisson_rates(s,mu,bias,m,diag(V))
  }
  
  # Calculate "local" ELBO F_jhl
  ELBO <- sum(x*(mu + bias + m)) - sum(a) -
          0.5*(tr(solve(Utilde, V)) + drop(t(m) %*% solve(Utilde,m))
               - R + logdet(Utilde) - logdet(V))
  
  return(list(Utilde = Utilde,m = m,V = V,a = a,ELBO = ELBO))
}

