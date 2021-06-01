#' @title Perform extreme deconvolution (ED) to estimate data-driven prior covariance matrices
#'
#' @param data A pois.mash data object containing the following components:
#' \item{X}{J x R matrix of count data collapsed over conditions, with features as rows and conditions as columns.} 
#' \item{s}{R x 1 numeric vector adjusting for sequencing depth of each of R conditions.}
#' \item{subgroup}{R x 1 factor vector with M levels denoting the subgroup status of each of R conditions.}
#' 
#' @param subset The indices of features to be used. Defaults to all of them.
#' 
#' @param Ulist A list of H full-rank covariance matrices (e.g., initialized by \code{pois_cov_init})
#' 
#' @param ulist A list of G numeric vectors each of which forming a rank-1 covariance matrix  
#' 
#' @param ulist.dd A G by 1 logical vector denoting whether each element in ulist is data-driven (\code{TRUE}) or canonical (\code{FALSE}). 
#' Defaults to data-driven for all elements. For canonical covariances, the spanned space is not updated.
#' 
#' @param ruv A logical scalar indicating whether to account for unwanted variation. Default is \code{FALSE}. If \code{TRUE}, Fuv must be provided. 
#' 
#' @param Fuv A J by D matrix of latent factors causing unwanted variation, with features as rows and latent factors as columns
#' 
#' @param verbose A logical scalar indicating whether to print ELBO at each iteration. Default is \code{FALSE}.
#' 
#' @param init A list of initial values for model parameters (e.g., returned by \code{pois_mash_ruv_prefit}). Could be empty.
#' 
#' @param control A list of control parameters with the following elements:
#' \item{maxiter}{Maximum number of ED iterations. Default is 500.}
#' \item{maxiter.q}{Maximum number of inner loop iterations to update variational parameters at each ED iteration. Default is 25.}
#' \item{tol.stop}{Tolerance for assessing convergence of ED, as measured by relative change in ELBO. Default is 1e-6.} 
#' \item{tol.q}{Relative tolerance for assessing convergence of variational parameters at each ED iteration. Default is 1e-2.}
#' \item{tol.rho}{Tolerance for assessing convergence of effects corresponding to unwanted variation. Default is 1e-6.}
#' 
#' @return A list with the following elements: 
#' \item{Ulist}{A list of H full-rank covariance matrices.}
#' \item{ulist}{A list of G numeric vectors each of which forming a rank-1 covariance matrix.}
#' \item{pi}{(H+G) by 1 numeric vector of mixture proportions for Ulist and ulist.}


pois_cov_ed <- function(data, subset=NULL, Ulist, ulist, ulist.dd=NULL, ruv=FALSE, Fuv=NULL, verbose=FALSE, init=list(NULL),  
                        control=list(maxiter=500, maxiter.q=25, tol.q=1e-2, tol.rho=1e-6, tol.stop=1e-6)){
  X <- data$X
  s <- data$s
  subgroup <- data$subgroup
  
  if(is.null(ulist)){
    stop("ulist cannot be empty!")
  }
  
  if(is.null(subset)){
    subset <- 1:nrow(X)
  }
  
  data.ed <- as.matrix(X[subset,])
  J <- nrow(data.ed)
  R <- ncol(data.ed)
  M <- length(unique(subgroup))
  subgroup <- as.numeric(as.factor(subgroup))
  maxiter <- control$maxiter
  maxiter.q <- control$maxiter.q
  tol.q <- control$tol.q
  tol.rho <- control$tol.rho
  tol.stop <- control$tol.stop
  
  if(is.null(maxiter)){
    maxiter <- 500
  }
  
  if(is.null(maxiter.q)){
    maxiter.q <- 25
  }
  
  if(is.null(tol.q)){
    tol.q <- 1e-2
  }
  
  if(is.null(tol.rho)){
    tol.rho <- 1e-6
  }
  
  if(is.null(tol.stop)){
    tol.stop <- 1e-6
  }
  
  H <- length(Ulist)
  G <- length(ulist)
  K <- H+G
  
  if(is.null(ulist.dd)){
    ulist.dd <- rep(TRUE, G)
    # set to FALSE if zero vector
    for(g in 1:G){
      if(sum(ulist[[g]]!=0)==0){
        ulist.dd[g] <- FALSE
      }
    }
  }
  
  # initialize pi_h and pi_g
  pi <- init$pi
  if(is.null(pi)){
    pi <- rep(1/K, K)
  }
  
  # initialize rho and bias
  rho <- init$rho
  diff.rho <- NULL
  bias <- matrix(0, nrow=J, ncol=R)
  
  if(ruv){
    if(is.null(Fuv)){
      stop("The matrix Fuv must be provided if ruv is set to TRUE")
    }
    F.ed <- as.matrix(Fuv[subset,])
    D <- ncol(F.ed)
    if(is.null(rho)){
      rho <- matrix(0, nrow=D, ncol=R)      
    }
    else{
      rho <- as.matrix(rho)
    }
    bias <- F.ed %*% rho
  }
  
  # initialize mu by ignoring condition-specific effects (i.e., theta) and unwanted variation
  mu <- init$mu
  if(is.null(mu)){
    mu <- matrix(NA, nrow=J, ncol=R)
    for(i in 1:M){
      mu[,subgroup==i] <- log(rowSums(data.ed[,subgroup==i])) - log(sum(s[subgroup==i]))
    }
  }
  else{
    mu <- mu[subset,]
  }
  
  # get a rough estimate of log lambda, which is useful for estimating the range of psi2, Ulist, ulist
  s.mat <- rep(1,J) %*% t(s)
  loglambda <- log((data.ed+0.1)/s.mat)
  upr_bd <- 4*max(apply(loglambda, 1, sd)^2) 
  minpsi2 <- pmax(min(apply(loglambda, 1, sd)^2)/1e2, 1e-8)
  maxpsi2 <- max(apply(loglambda, 1, sd)^2)
  
  # use grid search to initialize psi^2 by fitting a poisson-log-normal model while ignoring fixed effects (i.e., beta_j) and unwanted variation
  psi2 <- init$psi2
  if(is.null(psi2)){
    psi2 <- rep(NA, J)
    
    for(j in 1:J){
      psi2_max <- pmax(sd(loglambda[j,])^2, 1)
      log2_psi2_grid <- seq(log2(1e-4), log2(psi2_max), length.out = 25)
      psi2_grid <- 2^log2_psi2_grid
      
      logdens <- rep(0, length(psi2_grid))
      
      for(l in 1:length(psi2_grid)){
        for(r in 1:R){
          logdens[l] <- logdens[l] + log(dpoilog(data.ed[j,r], mu[j,r] + log(s[r]), sqrt(psi2_grid[l])))
        }
      }
      
      psi2[j] <- psi2_grid[which.max(logdens)]
    }
  }
  else{
    psi2 <- psi2[subset]
  }
  
  
  # matrices and arrays to store the posterior mean and covariance of theta, i.e., gamma_jk, Sigma_jk
  gamma_jk <- list(NULL)
  Sigma_jk <- list(NULL)
  for(k in 1:K){
    gamma_jk[[k]] <- matrix(NA, nrow=J, ncol=R)
    Sigma_jk[[k]] <- array(NA, c(J, R, R))
  }
  
  # J x K x R array A to store the quantities related to q_jk, s.t. A[j,k,r] = gamma_jkr + 0.5*Sigma_jk,rr
  A <- array(NA, c(J, K, R))
  
  # matrices to store local ELBO F_jk
  ELBOs <- matrix(0, nrow=J, ncol=K)     
  
  # update posterior mean and covariance of theta and local ELBO
  for(j in 1:J){
    if(H > 0){
      for(h in 1:H){
        theta.qjh <- update_q_theta_general(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], U=Ulist[[h]],
                                            control=list(maxiter=maxiter.q, tol=tol.q))
        gamma_jk[[h]][j,] <- theta.qjh$m
        Sigma_jk[[h]][j,,] <- theta.qjh$V
        ELBOs[j,h] <- theta.qjh$ELBO
        A[j,h,] <- theta.qjh$m + diag(theta.qjh$V)/2
      }      
    }
    
    for(g in 1:G){
      theta.qjg <- update_q_theta_rank1(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], u=ulist[[g]],
                                        control=list(maxiter=maxiter.q, tol=tol.q))
      gamma_jk[[H+g]][j,] <- theta.qjg$m
      Sigma_jk[[H+g]][j,,] <- theta.qjg$V
      ELBOs[j,H+g] <- theta.qjg$ELBO
      A[j,H+g,] <- theta.qjg$m + diag(theta.qjg$V)/2
    }
  }
  
  
  # update J x K matrix zeta of posterior weights
  ELBOs.cen <- ELBOs - apply(ELBOs, 1, max)
  zeta <- t(t(exp(ELBOs.cen)) * pi)
  zeta <- zeta*(1/rowSums(zeta)) 
  zeta <- pmax(zeta, 1e-15)
  
  # update J x R matrix tmp.ruv needed to update rho, s.t. tmp.ruv[j,r] = sum_k zeta[j,k]*exp(A[j,k,r])
  tmp.ruv <- matrix(NA, nrow=J, ncol=R)
  for(r in 1:R){
    tmp.ruv[,r] <- rowSums(zeta*exp(A[,,r]))
  }
  
  const <- sum(data.ed%*%log(s)) - sum(lgamma(data.ed+1))
  
  
  # overall ELBO after updating all parameters at each iteration
  ELBOs.overall <- c()
  
  if(verbose){
    cat("Start running extreme deconvolution to estimate prior covariance matrices.\n")
  }
  
  for(iter in 1:maxiter){
    # calculate overall ELBO at the current iteration
    ELBO.overall <- sum(zeta*(log(rep(1,J)%*%t(pi)) + ELBOs - log(zeta))) + const
    ELBOs.overall <- c(ELBOs.overall, ELBO.overall)
    
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO.overall))
    }
    
    if(iter >= 50){
      if(is.finite(ELBOs.overall[iter]) & is.finite(ELBOs.overall[iter-1])){
        if(abs(ELBOs.overall[iter]-ELBOs.overall[iter-1])/abs(ELBOs.overall[iter-1]) < tol.stop)
          break
      }
    }
    
    # calculate or update quantities related to model parameters mu, psi2, U_h, u_g, rho
    tmp.mu <- array(0, c(J, K, M))
    tmp.psi2 <- matrix(0, nrow=J, ncol=K)
    diff.U <- rep(0, K)
    
    if(H > 0){
      for(h in 1:H){
        tmp.U <- matrix(0, nrow=R, ncol=R)
        for(j in 1:J){
          gamma.tmp <- gamma_jk[[h]][j,]
          Sigma.tmp <- Sigma_jk[[h]][j,,]
          beta.qjh <- update_q_beta_general(theta_m=gamma.tmp, theta_V=Sigma.tmp, c2=rep(1,R), psi2=psi2[j], U=Ulist[[h]])$beta2_m
          tmp.U <- tmp.U + zeta[j,h]*beta.qjh
          for(i in 1:M){
            tmp.mu[j,h,i] <- sum(s[subgroup==i]*exp(bias[j,subgroup==i] + gamma.tmp[subgroup==i] + diag(Sigma.tmp)[subgroup==i]/2))
          }
          eta.qjh <- update_q_eta_general(theta_m=gamma.tmp, theta_V=Sigma.tmp, c2=rep(1,R), psi2=psi2[j], U=Ulist[[h]])
          tmp.psi2[j,h] <- sum(eta.qjh)
        }
        
        # update U_h
        Uh.new <- tmp.U/sum(zeta[,h])
        Uh.new <- (Uh.new + t(Uh.new))/2
        
        # avoid too large values in U_h
        if(max(diag(Uh.new)) > upr_bd){
          Uh.new <- upr_bd*Uh.new/max(diag(Uh.new))
        }
        
        diff.U[h] <- max(abs(Uh.new - Ulist[[h]]))
        Ulist[[h]] <- Uh.new
      }
    }
    
    for(g in 1:G){
      # if u is zero vector
      if(sum(ulist[[g]]!=0)==0){
        for(j in 1:J){
          gamma.tmp <- gamma_jk[[H+g]][j,]
          Sigma.tmp <- Sigma_jk[[H+g]][j,,]
          for(i in 1:M){
            tmp.mu[j,H+g,i] <- sum(s[subgroup==i]*exp(bias[j,subgroup==i] + gamma.tmp[subgroup==i] + diag(Sigma.tmp)[subgroup==i]/2))
          }
          eta.qjg <- gamma.tmp^2 + diag(Sigma.tmp) 
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
      }
      
      # if u is data-driven 
      else if (ulist.dd[g]){
        tmp1.u <- 0
        tmp2.u <- rep(0, R)
        for(j in 1:J){
          gamma.tmp <- gamma_jk[[H+g]][j,]
          Sigma.tmp <- Sigma_jk[[H+g]][j,,]
          beta.qjg <- update_q_beta_rank1(theta_m=gamma.tmp, theta_V=Sigma.tmp, c2=rep(1,R), psi2=psi2[j], u=ulist[[g]])
          tmp1.u <- tmp1.u + zeta[j,H+g]*beta.qjg$a2_m/psi2[j]
          tmp2.u <- tmp2.u + zeta[j,H+g]*beta.qjg$a_theta_m/psi2[j]
          for(i in 1:M){
            tmp.mu[j,H+g,i] <- sum(s[subgroup==i]*exp(bias[j,subgroup==i] + gamma.tmp[subgroup==i] + diag(Sigma.tmp)[subgroup==i]/2))
          }
          eta.qjg <- update_q_eta_rank1(theta_m=gamma.tmp, theta_V=Sigma.tmp, a2_m=beta.qjg$a2_m, a_theta_m=beta.qjg$a_theta_m, u=ulist[[g]])
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
        
        # update u_g
        ug.new <- tmp2.u/pmax(tmp1.u, 1e-8)
        
        # avoid too large values in u_g
        if(max(abs(ug.new)) > sqrt(upr_bd)){
          ug.new <- sqrt(upr_bd)*ug.new/max(abs(ug.new))
        }
        
        diff.U[H+g] <- max(abs(ug.new - ulist[[g]]))
        ulist[[g]] <- ug.new
      }
      
      # if u is canonical
      else{
        ug <- ulist[[g]]
        ug <- ug/ug[which.max(abs(ug))]
        tmp1.u <- 0
        tmp2.u <- 0
        for(j in 1:J){
          gamma.tmp <- gamma_jk[[H+g]][j,]
          Sigma.tmp <- Sigma_jk[[H+g]][j,,]
          beta.qjg <- update_q_beta_rank1(theta_m=gamma.tmp, theta_V=Sigma.tmp, c2=rep(1,R), psi2=psi2[j], u=ulist[[g]])
          tmp1.u <- tmp1.u + zeta[j,H+g]*sum(ug^2)*beta.qjg$a2_m/psi2[j]
          tmp2.u <- tmp2.u + zeta[j,H+g]*sum(ug*beta.qjg$a_theta_m)/psi2[j]
          for(i in 1:M){
            tmp.mu[j,H+g,i] <- sum(s[subgroup==i]*exp(bias[j,subgroup==i] + gamma.tmp[subgroup==i] + diag(Sigma.tmp)[subgroup==i]/2))
          }
          eta.qjg <- update_q_eta_rank1(theta_m=gamma.tmp, theta_V=Sigma.tmp, a2_m=beta.qjg$a2_m, a_theta_m=beta.qjg$a_theta_m, u=ulist[[g]])
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
        
        # update u_g
        ug.new <- pmin(pmax(tmp2.u/tmp1.u, 1e-2), sqrt(upr_bd))*ug
        diff.U[H+g] <- max(abs(ug.new - ulist[[g]]))
        ulist[[g]] <- ug.new
      }
    }
    
    # update mu
    for(i in 1:M){
      mu.i.new <- log(rowSums(data.ed[,subgroup==i])) - log(rowSums(zeta*tmp.mu[,,i]))
      mu[,subgroup==i] <- mu.i.new
    }
    
    # update psi2 
    psi2.new <- rowSums(zeta*tmp.psi2)/R
    psi2 <- pmin(pmax(psi2.new, minpsi2), maxpsi2)
    
    # update pi
    pi.new <- colMeans(zeta)
    pi.new <- pmax(pi.new, 1e-6)
    diff.pi <- pi.new - pi
    pi <- pi.new
    
    # update rho and bias if ruv=TRUE
    if(ruv){
      rho.new <- matrix(NA, nrow=nrow(rho), ncol=ncol(rho))
      for(r in 1:R){
        rho.new[,r] <- update_rho(Xr=data.ed[,r], Fuv=F.ed, sr=s[r], mu=mu[,r], Lr=tmp.ruv[,r], init=rho[,r], 
                                  control=list(maxiter=100, tol=tol.rho, maxrho=1e2/max(abs(F.ed))))$rho 
      }
      diff.rho <- rho.new - rho
      rho <- rho.new
      bias <- F.ed %*% rho
    }
    
    # update posterior mean and covariance of theta and local ELBO F_jk
    for(j in 1:J){
      if(H > 0){
        for(h in 1:H){
          theta.qjh <- update_q_theta_general(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], U=Ulist[[h]],
                                              init=list(m=gamma_jk[[h]][j,], V=Sigma_jk[[h]][j,,]), control=list(maxiter=maxiter.q, tol=tol.q))
          gamma_jk[[h]][j,] <- theta.qjh$m
          Sigma_jk[[h]][j,,] <- theta.qjh$V
          ELBOs[j,h] <- theta.qjh$ELBO
          A[j,h,] <- theta.qjh$m + diag(theta.qjh$V)/2
        }
      }
      
      for(g in 1:G){
        theta.qjg <- update_q_theta_rank1(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], u=ulist[[g]],
                                          init=list(m=gamma_jk[[H+g]][j,], V=Sigma_jk[[H+g]][j,,]), control=list(maxiter=maxiter.q, tol=tol.q))
        gamma_jk[[H+g]][j,] <- theta.qjg$m
        Sigma_jk[[H+g]][j,,] <- theta.qjg$V
        ELBOs[j,H+g] <- theta.qjg$ELBO
        A[j,H+g,] <- theta.qjg$m + diag(theta.qjg$V)/2
      }
    }
    
    # update J x K matrix zeta of posterior weights
    ELBOs.cen <- ELBOs - apply(ELBOs, 1, max)
    zeta <- t(t(exp(ELBOs.cen)) * pi)
    zeta <- zeta*(1/rowSums(zeta)) 
    zeta <- pmax(zeta, 1e-15)
    
    # update J x R matrix tmp.ruv needed to update rho, s.t. tmp.ruv[j,r] = sum_k zeta[j,k]*exp(A[j,k,r])
    for(r in 1:R){
      tmp.ruv[,r] <- rowSums(zeta*exp(A[,,r]))
    }
  }
  
  # name the model paramter estimates
  rownames(mu) <- rownames(data.ed)
  colnames(mu) <- names(s)
  names(psi2) <- rownames(data.ed)
  colnames(rho) <- colnames(data.ed)
  names(pi) <- c(names(Ulist), names(ulist))
  
  if(verbose){
    cat("Finish running extreme deconvolution to estimate prior covariance matrices.\n")
  }
  
  return(list(mu=mu, psi2=psi2, rho=rho, Ulist=Ulist, ulist=ulist, pi=pi, zeta=zeta, ELBO=ELBOs.overall, 
              diff.U=diff.U, diff.pi=diff.pi, diff.rho=diff.rho))
}
