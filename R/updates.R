# Update q(theta) = q(eta) for a given unit without fixed effects,
# i.e., with prior covariance U = 0.
update_q_eta_only <- function (x, s, mu, bias, c2, psi2, init = list(NULL),
                               control = list(maxiter = 25, tol = 0.01,
                                              lwr = -10, upr = 10)) {
  R       <- length(x)
  maxiter <- control$maxiter
  tol     <- control$tol
  lwr     <- control$lwr
  upr     <- control$upr
  m       <- init$m
  V       <- init$V
  
  Utilde <- psi2 * c2
  
  if (is.null(maxiter))
    maxiter <- 25
  if (is.null(tol))
    tol <- 0.01
  if (is.null(lwr))
    lwr <- -10
  if (is.null(upr))
    upr <- 10
  if (is.null(m))
    m <- rep(0, R)
  if(is.null(V))
    V <- Utilde
  
  bias <- as.numeric(bias)
  a    <- s*exp(mu + bias + m + V/2)
  
  for (iter in 1:maxiter) {
    V_new <- 1/(1/Utilde + a)
    a     <- s * exp(mu + bias + m + V_new/2)
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
    a <- s*exp(mu + bias + m + V/2)
  }
  
  # Calculate "local" ELBO F_j.
  ELBO <- sum(x*(mu + bias + m)) - sum(a) -
          0.5*(sum(V/Utilde) + sum(m^2/Utilde) - R +
               sum(log(Utilde)) - sum(log(V)))
  return(list(m = m,V = V,a = a,ELBO = ELBO))
}

# Update q(theta) for a given unit and prior covariance w*U, where U
# has full rank.
update_q_theta_general <- function (x, s, mu, bias, c2, psi2, w = 1,
                                    U, init = list(NULL),
                                    control = list(maxiter = 25, tol = 0.01,
                                                   lwr = -10, upr = 10)) {
  R       <- length(x)
  maxiter <- control$maxiter
  tol     <- control$tol
  lwr     <- control$lwr
  upr     <- control$upr
  m       <- init$m
  V       <- init$V
  Utilde  <- w*U + psi2*diag(c2) 
  
  if (is.null(maxiter))
    maxiter <- 25
  if (is.null(tol))
    tol <- 0.01
  if (is.null(lwr))
    lwr <- -10
  if (is.null(upr))
    upr <- 10
  if (is.null(m))
    m <- rep(0,R)
  if (is.null(V))
    V <- Utilde
  
  bias <- as.numeric(bias)
  a    <- s*exp(mu + bias + m + diag(V)/2)
  
  for (iter in 1:maxiter) {
    V_new <- solve(solve(Utilde) + diag(a),tol = 1e-50)
    a     <- s*exp(mu + bias + m + diag(V_new)/2)
    m_new <- as.numeric(m - V_new %*% (a - x + solve(Utilde,m)))
    
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
    a <- s*exp(mu + bias + m + diag(V)/2)
  }
  
  # Calculate "local" ELBO F_jhl.
  ELBO <- sum(x*(mu+bias+m)) - sum(a) -
          0.5*(tr(solve(Utilde, V)) + as.numeric(t(m) %*% solve(Utilde,m)) -
          R + determinant(Utilde,logarithm = TRUE)$modulus -
          determinant(V,logarithm = TRUE)$modulus)   
  return(list(Utilde = Utilde,m = m,V = V,a = a,ELBO = ELBO))
}

# Update q(theta) for a given unit and prior covariance w*U, where U = u*u'.
update_q_theta_rank1 <- function (x, s, mu, bias, c2, psi2, w = 1, u,
                                  init = list(NULL),
                                  control=list(maxiter = 25, tol = 0.01,
                                               lwr = -10, upr = 10)) {
  R       <- length(x)
  maxiter <- control$maxiter
  tol     <- control$tol
  lwr     <- control$lwr
  upr     <- control$upr
  m       <- init$m
  V       <- init$V
  
  Utilde <- w*u %*% t(u) + psi2 * diag(c2) 
  S_inv  <- 1/(psi2 * c2)
  
  if (is.null(maxiter))
    maxiter <- 25
  if (is.null(tol))
    tol <- 0.01
  if (is.null(lwr))
    lwr <- -10
  if (is.null(upr))
    upr <- 10
  if (is.null(m))
    m <- rep(0,R)
  if (is.null(V))
    V <- Utilde
  
  bias <- as.numeric(bias)
  a    <- s*exp(mu + bias + m + diag(V)/2)
  
  for (iter in 1:maxiter) {
    V_new <- mat_inv_rank1(a + S_inv,-w*u*S_inv,
                           (u*S_inv)/(1 + w*sum(u^2*S_inv)))
    a            <- s*exp(mu + bias + m + diag(V_new)/2)
    Utilde_inv_m <- m*S_inv - w*u*S_inv*sum(u*m*S_inv)/(1 + w*sum(u^2*S_inv))
    m_new        <- as.numeric(m - V_new %*% (a - x + Utilde_inv_m))
    
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
    a <- s*exp(mu + bias + m + diag(V)/2)
  }
  
  # Calculate "local" ELBO F_jhl
  ELBO <- sum(x*(mu +  bias + m)) - sum(a) -
          0.5*(tr(solve(Utilde, V)) + as.numeric(t(m) %*% solve(Utilde, m)) -
          R + determinant(Utilde,logarithm = TRUE)$modulus -
          determinant(V,logarithm = TRUE)$modulus) 
  return(list(Utilde = Utilde,m = m,V = V,a = a,ELBO = ELBO))
}

# Update q(beta) for a given unit and prior covariance w*U, where U
# has full rank.
update_q_beta_general <- function (theta_m, theta_V, c2, psi2, w = 1, U) {
    
  # Make sure eigenvalues of U are all strictly positive.
  U       <- w*U
  eig.U   <- eigen(U)
  eig.val <- pmax(eig.U$values,1e-8)
  U       <- tcrossprod(eig.U$vectors %*% diag(sqrt(eig.val)))
  S_inv   <- 1/(psi2*c2)
  tmp     <- solve(solve(U) + diag(S_inv))
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

# Update the D x 1 vector rho for a given condition r.
update_rho <- function (Xr, Fuv, sr, mu, Lr, init,
                        control = list(maxiter = 100, tol = 1e-8,
                                       maxrho = 100)) {
  maxiter <- control$maxiter
  tol     <- control$tol
  maxrho  <- control$maxrho
  
  if (is.null(maxiter))
    maxiter <- 100
  if (is.null(tol))
    tol <- 1e-8
  if (is.null(maxrho))
    maxrho <- 100
  
  rho <- init
  Fr  <- rep(as.numeric(NA),maxiter)
  D   <- length(rho)
  d1F <- rep(as.numeric(NA),D)
  d2F <- matrix(as.numeric(NA),D,D)
  
  for (iter in 1:maxiter) {
    bias <- Fuv %*% rho
    Fr[iter] <- sum(Xr*bias) - sr*sum(exp(mu + bias)*Lr)
    for (d in 1:D)
      d1F[d] <- sum(Xr*Fuv[,d]) - sr*sum(Fuv[,d]*Lr*exp(mu + bias))
    for (d in 1:D)
      for (t in 1:D)
        d2F[d,t] <- -sr*sum(Fuv[,d]*Fuv[,t]*Lr*exp(mu + bias))
    rho_new <- rho - solve(d2F,d1F)
    rho_new <- pmin(pmax(rho_new,-maxrho),maxrho)
    if (max(abs(rho_new - rho)) < tol)
      break
    rho <- rho_new
  }
  
  return(list(rho = rho,Fr = Fr))
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


# Update the J x R matrix of means mu, where X is J x R matrix of counts, subgroup is R x 1 factor vector with M levels, 
# zeta is J x K matrix of posterior weights, tmp.mu is J x K x M array.
update_mu <- function(X, subgroup, zeta, tmp.mu){    
  M <- length(unique(subgroup))
  mu <- matrix(as.numeric(NA), nrow(X), ncol(X))
  for (i in 1:M) {
    mu.i.new <- log(rowSums(X[,subgroup == i])) - log(rowSums(zeta * tmp.mu[,,i]))
    mu[,subgroup == i] <- mu.i.new
  }
  return(mu)
}


# Update the J x 1 vector of dispersion parameter psi2, where zeta is J x K matrix of posterior weights, tmp.psi2 is J x K matrix,
# R is the number of conditions, minpsi2 and maxpsi2 are respectively positive scalars giving the lower and upper bound of psi2.
update_psi2 <- function (zeta, tmp.psi2, R, minpsi2, maxpsi2){
  psi2.new <- rowSums(zeta * tmp.psi2)/R
  return(pmin(pmax(psi2.new,minpsi2),maxpsi2))
}

# Update the K x 1 vector of prior weights pi, where zeta is J x K
# matrix of posterior weights.
update_pi <- function (zeta) {
  pi <- colMeans(zeta)
  return(pmax(pi,1e-8))
}

# Update the D x R matrix of rho for all conditions given current rho, where X is J x R matrix of counts, s is R x 1 vector of sequencing depths, 
# mu is J x R matrix of means, Fuv is J x D matrix of latent factors causing unwanted variation, rho is D x R matrix of current rho,
# tmp.ruv is J x R matrix, tol.rho is a small positive number to assess convergence. 
update_rho_all <- function(X, s, mu, Fuv, rho, tmp.ruv, tol.rho=1e-6){
  rho.new <- matrix(as.numeric(NA),nrow(rho),ncol(rho))
  for (r in 1:ncol(X)) 
    rho.new[,r] <- update_rho(Xr = X[,r], Fuv = Fuv, sr = s[r], mu = mu[,r], Lr = tmp.ruv[,r], init = rho[,r], 
                              control = list(maxiter = 100, tol = tol.rho, maxrho=100/max(abs(Fuv))))$rho
  return(rho.new)
}


# Update the J x K matrix of posterior weights zeta, where ELBOs is J x K matrix of local ELBO, pi is K x 1 vector of prior weights.
update_zeta <- function(ELBOs, pi){
  ELBOs.cen <- ELBOs - apply(ELBOs,1,max)
  zeta <- t(t(exp(ELBOs.cen)) * pi)
  zeta <- zeta*(1/rowSums(zeta)) 
  zeta <- pmax(zeta, 1e-15)
}


# Update q(theta) for a given unit and all prior mixtures w_l*U_k for beta.
update_q_theta_all <- function(x, s, mu, bias, c2=rep(1,length(x)), psi2, wlist=1, Ulist, ulist, init = list(NULL), maxiter.q=25, tol.q=1e-2){
  R <- length(x)
  H <- length(Ulist)
  G <- length(ulist)
  L <- length(wlist)
  K <- L*(H+G)
  
  gamma <- matrix(as.numeric(NA), nrow=K, ncol=R)
  Sigma <- array(as.numeric(NA), c(K,R,R))
  A <- matrix(as.numeric(NA), nrow=K, ncol=R)
  ELBOs <- rep(as.numeric(NA), K)
  
  if (H > 0) {
    hl <- 0
    for (h in 1:H) {
      for(l in 1:L) {
        hl <- hl + 1
        if(!is.null(init$gamma) & !is.null(init$Sigma)){
          init.m <- init$gamma[hl,]
          init.V <- init$Sigma[hl,,]
        }
        else{
          init.m <- NULL
          init.V <- NULL
        }
        theta.qjhl <- update_q_theta_general(x = x, s = s, mu = mu, bias = bias, c2 = c2, psi2 = psi2, w = wlist[l], U = Ulist[[h]],
                                             init = list(m = init.m,V = init.V), control = list(maxiter = maxiter.q, tol = tol.q))
        gamma[hl,] <- theta.qjhl$m
        Sigma[hl,,] <- theta.qjhl$V
        ELBOs[hl] <- theta.qjhl$ELBO
        A[hl,] <- theta.qjhl$m + diag(theta.qjhl$V)/2        
      }
    }
  }
  
  gl <- 0
  for (g in 1:G) {
    for(l in 1:L) {
      gl <- gl + 1
      if(!is.null(init$gamma) & !is.null(init$Sigma)){
        init.m <- init$gamma[H*L+gl,]
        init.V <- init$Sigma[H*L+gl,,]
      }
      else{
        init.m <- NULL
        init.V <- NULL
      }
      theta.qjgl <- update_q_theta_rank1(x = x, s = s, mu = mu, bias = bias, c2 = c2, psi2 = psi2, w = wlist[l], u = ulist[[g]], 
                                         init = list(m = init.m,V = init.V), control = list(maxiter = maxiter.q, tol = tol.q))
      gamma[H*L+gl,]  <- theta.qjgl$m
      Sigma[H*L+gl,,] <- theta.qjgl$V
      ELBOs[H*L+gl] <- theta.qjgl$ELBO
      A[H*L+gl,] <- theta.qjgl$m + diag(theta.qjgl$V)/2
    }
  }
  
  return(list(ELBOs=ELBOs, A=A, gamma=gamma, Sigma=Sigma))
}


# Update posterior distribution of theta, beta and eta and accompanying quantities (ELBOs, gamma, A, tmp.mu, tmp.psi2) for selected units.
update_q_by_j <- function(X, s, subgroup, idx.update, mu, bias, psi2, wlist=1, Ulist, ulist, ulist.epsilon2,
                          gamma, A, ELBOs, tmp.mu, tmp.psi2, maxiter.q=25, tol.q=1e-2){
  J <- nrow(X)
  R <- ncol(X)
  M <- length(unique(subgroup))
  H <- length(Ulist)
  G <- length(ulist)
  L <- length(wlist)
  
  for (j in 1:length(idx.update)) {
    j.idx <- idx.update[j]
    
    if (H > 0) {
      hl <- 0
      for (h in 1:H) {
        for (l in 1:L) {
          hl <- hl + 1
          if(any(is.na(gamma[j.idx,hl,])) | any(is.na(A[j.idx,hl,]))){
            init.m <- NULL
            init.V <- NULL
          }
          else{
            init.m <- gamma[j.idx,hl,]
            Utilde <- wlist[l] * Ulist[[h]] + psi2[j.idx] * diag(R)
            a.tmp <- s * exp(mu[j.idx,] + bias[j.idx,] + A[j.idx,hl,])
            init.V <- solve(solve(Utilde) + diag(a.tmp), tol = 1e-50)
          }
          theta.qjhl <- update_q_theta_general(x = X[j.idx,], s = s, mu = mu[j.idx,], bias = bias[j.idx,], c2 = rep(1,R),
                                               psi2 = psi2[j.idx], w = wlist[l], U = Ulist[[h]],
                                               init = list(m = init.m,V = init.V), control = list(maxiter=maxiter.q,tol=tol.q))
          ELBOs[j.idx,hl] <- theta.qjhl$ELBO
          gamma.tmp <- theta.qjhl$m
          Sigma.tmp <- theta.qjhl$V
          for (i in 1:M)
            tmp.mu[j.idx,hl,i] <- sum(s[subgroup == i] * exp(bias[j.idx,subgroup == i] + gamma.tmp[subgroup == i] + diag(Sigma.tmp)[subgroup == i]/2))
          eta.qjhl <- update_q_eta_general(theta_m = gamma.tmp, theta_V = Sigma.tmp, c2 = rep(1,R), psi2 = psi2[j.idx], w = wlist[l], U = Ulist[[h]])
          tmp.psi2[j.idx,hl] <- sum(eta.qjhl)
          gamma[j.idx,hl,] <- gamma.tmp
          A[j.idx,hl,] <- gamma.tmp + diag(Sigma.tmp)/2 
        }
      }        
    }
    
    gl <- 0
    for (g in 1:G) {
      ug <- ulist[[g]]
      epsilon2.g <- ulist.epsilon2[g]
      for (l in 1:L) {
        gl <- gl + 1
        if(any(is.na(gamma[j.idx,H*L+gl,])) | any(is.na(A[j.idx,H*L+gl,]))){
          init.m <- NULL
          init.V <- NULL
        }
        else{
          init.m <- gamma[j.idx,H*L+gl,]
          a.tmp  <- s * exp(mu[j.idx,] + bias[j.idx,] + A[j.idx,H*L+gl,])
          S_inv  <- 1/(psi2[j.idx] * rep(1,R) + wlist[l] * epsilon2.g)
          init.V <- mat_inv_rank1(a.tmp + S_inv,-wlist[l] * ug * S_inv, (ug*S_inv)/(1 + wlist[l]*sum(ug^2*S_inv)))      
        }
        theta.qjgl <- update_q_theta_rank1(x = X[j.idx,], s = s, mu = mu[j.idx,], bias = bias[j.idx,], c2 = rep(1,R),
                                           psi2 = psi2[j.idx] + wlist[l] * epsilon2.g, w = wlist[l], u = ug,
                                           init = list(m = init.m,V = init.V), control = list(maxiter = maxiter.q,tol = tol.q))
        ELBOs[j.idx,H*L+gl] <- theta.qjgl$ELBO
        gamma.tmp <- theta.qjgl$m
        Sigma.tmp <- theta.qjgl$V
        for (i in 1:M)
          tmp.mu[j.idx,H*L+gl,i] <- sum(s[subgroup == i] * exp(bias[j.idx,subgroup == i] + gamma.tmp[subgroup == i] + diag(Sigma.tmp)[subgroup == i]/2))
        gamma[j.idx,H*L+gl,] <- gamma.tmp
        A[j.idx,H*L+gl,] <- gamma.tmp + diag(Sigma.tmp)/2
        
        # If ug is zero vector.
        if (sum(ug != 0) == 0) {
          eta.qjgl <- gamma.tmp^2 + diag(Sigma.tmp) 
          tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
        }
        else if (epsilon2.g > 1e-4) {
          eta.qjgl <- update_q_eta_rank1_robust(theta_m = gamma.tmp, theta_V = Sigma.tmp, c2 = rep(1,R), psi2 = psi2[j.idx],
                                                w = wlist[l], u = ug, epsilon2 = epsilon2.g)
          tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
        }
        else {
          beta.qjgl <- update_q_beta_rank1(theta_m = gamma.tmp, theta_V = Sigma.tmp, c2 = rep(1,R), psi2 = psi2[j.idx], w = wlist[l], u = ug)
          eta.qjgl <- update_q_eta_rank1(theta_m = gamma.tmp, theta_V = Sigma.tmp, a2_m = beta.qjgl$a2_m, a_theta_m = beta.qjgl$a_theta_m, u = ug)
          tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
        }
      }
    }      
  }
  
  return(list(ELBOs=ELBOs, A=A, gamma=gamma, tmp.mu=tmp.mu, tmp.psi2=tmp.psi2))
}
