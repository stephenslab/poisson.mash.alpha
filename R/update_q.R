# Update posterior distribution of theta, beta and eta and
# accompanying quantities (ELBOs, gamma, A, tmp.mu, tmp.psi2) for
# selected units.
update_q_by_j <- function (X, s, subgroup, idx.update, mu, bias, psi2,
                           wlist = 1, Ulist, ulist, ulist.epsilon2, gamma,
                           A, ELBOs, tmp.mu, tmp.psi2, maxiter.q = 25,
                           tol.q = 0.01) {
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
                                        Ulist[[h]],list(m=init.m,V=init.V),
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

  return(list(ELBOs = ELBOs,A = A,gamma = gamma,tmp.mu = tmp.mu,
              tmp.psi2 = tmp.psi2))
}

# This is the same as update_q_by_j, except that "dat" is a list
# containing all the arguments that contain the gene-wise data; that
# is, data that is indexed by the indices specified in idx.update.
update_q_by_j_with_dat <- function (dat, s, subgroup, wlist, Ulist, ulist,
                                    ulist.epsilon2, maxiter.q = 25, tol.q)
  update_q_by_j(dat$X,s,subgroup,1:nrow(dat$X),dat$mu,dat$bias,dat$psi2,
                wlist,Ulist,ulist,ulist.epsilon2,dat$gamma,dat$A,dat$ELBOs,
                dat$tmp.mu,dat$tmp.psi2,maxiter.q,tol.q)

# This is the multithreading variant of update_q_by_j. It could
# produce the same result as update_q_by_j, but possibly faster when
# nc >= 2.
#
#' @importFrom parallel splitIndices
#' @importFrom parallel mclapply
update_q_by_j_multicore <- function (X, s, subgroup, idx.update, mu, bias,
                                     psi2, wlist = 1, Ulist, ulist,
                                     ulist.epsilon2, gamma, A, ELBOs, tmp.mu,
                                     tmp.psi2, maxiter.q = 25, tol.q = 0.01,
                                     nc = 1) {
  if (nc == 1)
    return(update_q_by_j(X,s,subgroup,idx.update,mu,bias,psi2,wlist,Ulist,
                         ulist,ulist.epsilon2,gamma,A,ELBOs,tmp.mu,
                         tmp.psi2,maxiter.q,tol.q))
  else {
      
    # Split the data.
    n       <- length(idx.update)
    indices <- splitIndices(n,nc)
    indices <- lapply(indices,function (x) idx.update[x])
    dat     <- vector("list",nc)
    for (i in 1:nc) {
      js       <- indices[[i]]
      dat[[i]] <- list(X        = X[js,,drop=FALSE],
                       mu       = mu[js,,drop=FALSE],
                       bias     = bias[js,,drop=FALSE],
                       psi2     = psi2[js],
                       gamma    = gamma[js,,,drop=FALSE],
                       A        = A[js,,,drop=FALSE],
                       ELBOs    = ELBOs[js,,drop=FALSE],
                       tmp.mu   = tmp.mu[js,,,drop=FALSE],
                       tmp.psi2 = tmp.psi2[js,,drop=FALSE])
    }

    # Distribute the calculations using mclapply.
    ans <- mclapply(mc.cores = nc,dat,update_q_by_j_with_dat,s,subgroup,
                    wlist,Ulist,ulist,ulist.epsilon2,maxiter.q,tol.q)

    # Combine the individual update_q_by_j outputs, then output the
    # combined result.
    for (i in 1:nc) {
      js            <- indices[[i]]
      ELBOs[js,]    <- ans[[i]]$ELBOs
      A[js,,]       <- ans[[i]]$A
      gamma[js,,]   <- ans[[i]]$gamma
      tmp.mu[js,,]  <- ans[[i]]$tmp.mu
      tmp.psi2[js,] <- ans[[i]]$tmp.psi2
    }
  }

  return(list(ELBOs = ELBOs,A = A,gamma = gamma,tmp.mu = tmp.mu,
              tmp.psi2 = tmp.psi2))
}
