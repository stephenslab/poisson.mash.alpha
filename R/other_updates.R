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

# Update posterior distribution of theta, beta and eta and
# accompanying quantities (ELBOs, gamma, A, tmp.mu, tmp.psi2) for
# selected units.
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
