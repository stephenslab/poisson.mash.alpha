#' @title Prefit the poisson mash model to get an initial estimate of model parameters if accounting for unwanted variation
#' 
#' @description There is no need to call this function if unwanted variation is not accounted for.
#' 
#' @param data A pois.mash data object containing the following components:
#' \item{X}{J x R matrix of count data collapsed over conditions, with features as rows and conditions as columns.}
#' \item{s}{R x 1 numeric vector adjusting for sequencing depth of each of R conditions.}
#' \item{subgroup}{R x 1 factor vector with M levels denoting the subgroup status of each of R conditions.}
#' 
#' @param Fuv J by D matrix of latent factors causing unwanted variation, with features as rows and latent factors as columns
#' 
#' @param verbose A logical scalar indicating whether to print ELBO at each iteration. Default is \code{FALSE}.
#' 
#' @param init A list of initial values for model parameters which could be empty
#' 
#' @param control A list of control parameters with the following elements:
#' \item{maxiter}{Maximum number of iterations. Default is 100.}
#' \item{maxiter.q}{Maximum number of inner loop iterations to update variational parameters at each iteration. Default is 25.}
#' \item{tol.stop}{Tolerance for assessing convergence of prefit, as measured by relative change in ELBO. Default is 1e-6.} 
#' \item{tol.q}{Relative tolerance for assessing convergence of variational parameters at each iteration. Default is 1e-2.}
#' \item{tol.rho}{Tolerance for assessing convergence of effects corresponding to unwanted variation. Default is 1e-4.}
#' 
#' @return A list with initial estimates of model parameters


pois_mash_ruv_prefit <- function(data, Fuv, verbose=FALSE, init=list(NULL), 
                                 control=list(maxiter=100, maxiter.q=25, tol.q=1e-2, tol.rho=1e-4, tol.stop=1e-6)){
  s <- data$s
  subgroup <- data$subgroup
  data <- as.matrix(data$X)
  J <- nrow(data)
  R <- ncol(data)
  M <- length(unique(subgroup))
  subgroup <- as.numeric(as.factor(subgroup))
  Fuv <- as.matrix(Fuv)
  maxiter <- control$maxiter
  maxiter.q <- control$maxiter.q
  tol.q <- control$tol.q
  tol.rho <- control$tol.rho
  tol.stop <- control$tol.stop
  
  if(is.null(maxiter)){
    maxiter <- 100
  }
  
  if(is.null(maxiter.q)){
    maxiter.q <- 25
  }
  
  if(is.null(tol.q)){
    tol.q <- 1e-2
  }
  
  if(is.null(tol.rho)){
    tol.rho <- 1e-4
  }
  
  if(is.null(tol.stop)){
    tol.stop <- 1e-6
  }
  
  # initialize mu by ignoring random effects and unwanted variation
  mu <- init$mu
  if(is.null(mu)){
    mu <- matrix(NA, nrow=J, ncol=R)
    for(i in 1:M){
      mu[,subgroup==i] <- log(rowSums(data[,subgroup==i])) - log(sum(s[subgroup==i]))
    }
  }
  
  # get a rough estimate of log lambda, which is useful for estimating the range of psi2
  s.mat <- rep(1,J) %*% t(s)
  loglambda <- log((data+0.1)/s.mat)
  minpsi2 <- pmax(min(apply(loglambda, 1, sd)^2)/1e2, 1e-4)
  maxpsi2 <- max(apply(loglambda, 1, sd)^2) 
  
  # use grid search to initialize psi2 by fitting a poisson-log-normal model while ignoring the unwanted variation
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
          logdens[l] <- logdens[l] + log(dpoilog(data[j,r], mu[j,r] + log(s[r]), sqrt(psi2_grid[l])))
        }
      }
      
      psi2[j] <- psi2_grid[which.max(logdens)]
    }
  }
  
  # initialize rho and bias
  D <- ncol(Fuv)
  rho <- init$rho
  if(is.null(rho)){
    rho <- matrix(0, nrow=D, ncol=R)      
  }
  else{
    rho <- as.matrix(rho)
  }
  bias <- Fuv %*% rho
    
  # matrices to store the posterior mean and covariance of theta_j (equal to eta_j in this case), i.e., gamma_j, diag(Sigma_j)
  gamma <- matrix(NA, nrow=J, ncol=R)
  Sigma <- matrix(NA, nrow=J, ncol=R)
  
  # matrix to store the quantities related to q_j, s.t. A[j,r] = gamma_jr + 0.5*Sigma_j,rr
  A <- matrix(NA, nrow=J, ncol=R)
  
  # update posterior mean and covariance of theta
  for(j in 1:J){
    eta.qj <- update_q_eta_only(x=data[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], control=list(maxiter=maxiter.q, tol=tol.q))
    gamma[j,] <- eta.qj$m
    Sigma[j,] <- eta.qj$V
    A[j,] <- eta.qj$m + eta.qj$V/2
  }
  
  # vector to store local ELBO for each j
  ELBOs <- rep(NA, J)  
  
  # overall ELBO after updating all parameters at each iteration
  ELBOs.overall <- c()
  
  const <- sum(data%*%log(s)) - sum(lgamma(data+1))
  
  if(verbose){
    cat("Start prefitting Poisson mash model to initialize model parameters.\n")
  }
  
  for(iter in 1:maxiter){
    # calculate or update quantities related to model parameters mu, psi2, rho
    tmp.mu <- matrix(NA, nrow=J, ncol=M)
    tmp.psi2 <- rep(NA, J)
    
    for(j in 1:J){
      gamma.tmp <- gamma[j,]
      Sigma.tmp <- Sigma[j,]
      for(i in 1:M){
        tmp.mu[j,i] <- sum(s[subgroup==i]*exp(bias[j,subgroup==i] + gamma.tmp[subgroup==i] + Sigma.tmp[subgroup==i]/2))
      }
      tmp.psi2[j] <- sum(gamma.tmp^2 + Sigma.tmp)
    }
    
    # update mu
    for(i in 1:M){
      mu.i.new <- log(rowSums(data[,subgroup==i])) - log(tmp.mu[,i])
      mu[,subgroup==i] <- mu.i.new
    }
    
    # update psi2 
    psi2.new <- tmp.psi2/R
    psi2 <- pmin(pmax(psi2.new, minpsi2), maxpsi2)
     
    # update rho and bias
    rho.new <- matrix(NA, nrow=nrow(rho), ncol=ncol(rho))
    for(r in 1:R){
      rho.new[,r] <- update_rho(Xr=data[,r], Fuv=Fuv, sr=s[r], mu=mu[,r], Lr=exp(A[,r]), init=rho[,r], 
                                control=list(maxiter=100, tol=tol.rho, maxrho=1e2/max(abs(Fuv))))$rho 
    }
    diff.rho <- rho.new - rho
    rho <- rho.new
    bias <- Fuv %*% rho
    
    # update posterior mean and covariance of theta and local ELBO F_j
    for(j in 1:J){
      eta.qj <- update_q_eta_only(x=data[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], init=list(m=gamma[j,], V=Sigma[j,]),
                                  control=list(maxiter=maxiter.q, tol=tol.q))
      gamma[j,] <- eta.qj$m
      Sigma[j,] <- eta.qj$V
      A[j,] <- eta.qj$m + eta.qj$V/2
      ELBOs[j] <- eta.qj$ELBO
    }
    
    # calculate overall ELBO at the current iteration
    ELBO.overall <- sum(ELBOs) + const
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
  }
  
  # name the model paramter estimates
  rownames(mu) <- rownames(data)
  colnames(mu) <- names(s)
  names(psi2) <- rownames(data)
  colnames(rho) <- colnames(data)
  
  if(verbose){
    cat("Finish prefitting Poisson mash model to initialize model parameters.\n")
  }
  
  return(list(mu=mu, psi2=psi2, rho=rho, ELBO=ELBOs.overall, diff.rho=diff.rho))
}
