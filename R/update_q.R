# Update posterior distribution of theta, beta and eta and
# accompanying quantities (ELBOs, gamma, A, tmp.mu, tmp.psi2) for
# selected units.
update_q_by_j <- function (X, s, subgroup, idx.update, mu, bias, psi2,
                           wlist = 1, Ulist, ulist, ulist.epsilon2, gamma,
                           A, ELBOs, tmp.mu, tmp.psi2, maxiter.q = 25,
                           tol.q = 0.01, version = c("Rcpp","R")) {
  version <- match.arg(version)
  out <- update_q_by_j_r(X,s,subgroup,idx.update,mu,bias,psi2,wlist,Ulist,
                         ulist,ulist.epsilon2,gamma,A,ELBOs,tmp.mu,tmp.psi2,
                         maxiter.q,tol.q,version)
  return(out)
}
  
# This implements update_q_by_j with version = "R".
update_q_by_j_r <- function (X, s, subgroup, idx.update, mu, bias, psi2,
                             wlist, Ulist, ulist, ulist.epsilon2, gamma,
                             A, ELBOs, tmp.mu, tmp.psi2, maxiter.q, tol.q,
                             version) {
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
          if(any(is.na(gamma[j.idx,hl,])) | any(is.na(A[j.idx,hl,]))) {
            init.m <- NULL
            init.V <- NULL
          }
          else {
            init.m <- gamma[j.idx,hl,]
            Utilde <- wlist[l] * Ulist[[h]] + psi2[j.idx] * diag(R)
            a.tmp  <- compute_poisson_rates(s,mu[j.idx,],bias[j.idx,],
                                            A[j.idx,hl,])
            init.V <- update_V(Utilde,a.tmp)
          }
          out <- update_q_theta_general(X[j.idx,],s,mu[j.idx,],bias[j.idx,],
                                        rep(1,R),psi2[j.idx],wlist[l],
                                        Ulist[[h]],list(m = init.m,V = init.V),
                                        maxiter = maxiter.q,tol = tol.q)
          ELBOs[j.idx,hl] <- out$ELBO
          gamma.tmp       <- out$m
          Sigma.tmp       <- out$V
          for (i in 1:M) {
            k <- which(subgroup == i)
            tmp.mu[j.idx,hl,i] <-
              sum(compute_poisson_rates(s[k],bias = bias[j.idx,k],
                                        gamma = gamma.tmp[k],
                                        V = diag(Sigma.tmp)[k]))
          }
          tmp.psi2[j.idx,hl] <-
            sum(update_q_eta_general(theta_m = gamma.tmp,theta_V = Sigma.tmp,
                                     c2 = rep(1,R),psi2 = psi2[j.idx],
                                     w = wlist[l],U = Ulist[[h]]))
          gamma[j.idx,hl,]   <- gamma.tmp
          A[j.idx,hl,]       <- gamma.tmp + diag(Sigma.tmp)/2 
        }
      }        
    }
    
    gl <- 0
    for (g in 1:G) {
      ug         <- ulist[[g]]
      epsilon2.g <- ulist.epsilon2[g]
      for (l in 1:L) {
        gl <- gl + 1
        if (any(is.na(gamma[j.idx,H*L+gl,])) | any(is.na(A[j.idx,H*L+gl,]))) {
          init.m <- NULL
          init.V <- NULL
        }
        else {
          init.m <- gamma[j.idx,H*L+gl,]
          a.tmp  <- compute_poisson_rates(s,mu[j.idx,],bias[j.idx,],
                                          A[j.idx,H*L+gl,])
          S_inv  <- 1/(psi2[j.idx] * rep(1,R) + wlist[l] * epsilon2.g)
          init.V <- mat_inv_rank1(a.tmp + S_inv,-wlist[l] * ug * S_inv,
                                  (ug*S_inv)/(1 + wlist[l]*sum(ug^2*S_inv)))
        }
        out <- update_q_theta_rank1(X[j.idx,],s, mu[j.idx,],bias[j.idx,],
                                    rep(1,R),psi2[j.idx]+wlist[l]*epsilon2.g,
                                    wlist[l],ug,list(m = init.m,V = init.V),
                                    maxiter = maxiter.q,tol = tol.q)
        ELBOs[j.idx,H*L+gl] <- out$ELBO
        gamma.tmp           <- out$m
        Sigma.tmp           <- out$V
        for (i in 1:M) {
          k <- which(subgroup == i)
          tmp.mu[j.idx,H*L+gl,i] <-
            sum(compute_poisson_rates(s[k],bias = bias[j.idx,k],
                                      gamma = gamma.tmp[k],
                                      V = diag(Sigma.tmp)[k]))
        }
        gamma[j.idx,H*L+gl,] <- gamma.tmp
        A[j.idx,H*L+gl,]     <- gamma.tmp + diag(Sigma.tmp)/2
       
        # If ug is zero vector.
        if (sum(ug != 0) == 0)
          tmp.psi2[j.idx,H*L+gl] <- sum(gamma.tmp^2 + diag(Sigma.tmp))
        else if (epsilon2.g > 1e-4)
          tmp.psi2[j.idx,H*L+gl] <-
            sum(update_q_eta_rank1_robust(theta_m = gamma.tmp,
                                          theta_V = Sigma.tmp,c2 = rep(1,R),
                                          psi2 = psi2[j.idx],w = wlist[l],
                                          u = ug,epsilon2 = epsilon2.g))
        else {
          out <- update_q_beta_rank1(theta_m = gamma.tmp,theta_V = Sigma.tmp,
                                     c2 = rep(1,R),psi2 = psi2[j.idx],
                                     w = wlist[l],u = ug)
          tmp.psi2[j.idx,H*L+gl] <-
            sum(update_q_eta_rank1(theta_m = gamma.tmp,theta_V = Sigma.tmp,
                                   a2_m = out$a2_m,a_theta_m = out$a_theta_m,
                                   u = ug))
        }
      }
    }
  }
  
  return(list(ELBOs = ELBOs,A = A,gamma = gamma,
              tmp.mu = tmp.mu,tmp.psi2 = tmp.psi2))
}
