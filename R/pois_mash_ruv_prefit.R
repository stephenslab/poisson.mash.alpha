#' @title Prefit Poisson MASH Model Accounting for Unwanted Variation
#' 
#' @description Prefit the poisson mash model to get an initial
#'   estimate of model parameters if accounting for unwanted
#'   variation. Don't call this function if unwanted variation is not
#'   accounted for.
#' 
#' @param data \dQuote{pois.mash} data object, typically created by
#'   calling \code{\link{pois_mash_set_data}}.
#'
#' @param Fuv J x D matrix of latent factors causing unwanted
#'   variation, with features as rows and latent factors as columns.
#' 
#' @param verbose Logical scalar indicating whether to print ELBO at
#'   each iteration.
#' 
#' @param init List of initial values for model parameters which
#'   could be empty.
#'
#' @param version R (slower) and C++ (faster) implementations of the
#'   model fitting algorithm are provided; these are selected with
#'   \code{version = "R"} and \code{version = "Rcpp"}.
#' 
#' @param control List of control parameters with the following
#'   elements: \dQuote{maxiter}, maximum number of iterations;
#'   \dQuote{maxiter.q}, maximum number of inner-loop iterations to
#'   update variational parameters at each iteration; \dQuote{tol.stop},
#'   tolerance for assessing convergence of prefit, as measured by
#'   relative change in ELBO; \dQuote{tol.q}, relative tolerance for
#'   assessing convergence of variational parameters at each iteration;
#'   \dQuote{tol.rho}, tolerance for assessing convergence of effects
#'   corresponding to unwanted variation.
#' 
#' @return A list containing initial estimates of model parameters.
#'
#' @importFrom stats sd
#' @importFrom poilog dpoilog
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib poisson.mash.alpha
#' 
#' @export
#' 
pois_mash_ruv_prefit <- function (data, Fuv, verbose = FALSE,
                                  init = list(NULL),
                                  version = c("Rcpp","R"),
                                  control = list(maxiter   = 100,
                                                 maxiter.q = 25,
                                                 tol.q     = 0.01,
                                                 tol.rho   = 1e-4,
                                                 tol.stop  = 1e-6)) {
  s         <- data$s
  subgroup  <- data$subgroup
  data      <- as.matrix(data$X)
  J         <- nrow(data)
  R         <- ncol(data)
  M         <- length(unique(subgroup))
  subgroup  <- as.numeric(as.factor(subgroup))
  Fuv       <- as.matrix(Fuv)
  version   <- match.arg(version)
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  tol.stop  <- control$tol.stop
  
  if (is.null(maxiter))
    maxiter <- 100
  if(is.null(maxiter.q))
    maxiter.q <- 25
  if(is.null(tol.q))
    tol.q <- 0.01
  if(is.null(tol.rho))
    tol.rho <- 1e-4
  if(is.null(tol.stop))
    tol.stop <- 1e-6

  t_eta   <- proc.time()
  t_rho   <- proc.time()
  t_eta[] <- 0
  t_rho[] <- 0
  
  # Initialize mu by ignoring random effects and unwanted variation.
  mu <- init$mu
  if(is.null(mu)) {

    # CAN THIS BE A FUNCTION? e.g., initialize_mu
    # (start of function)
    mu <- matrix(as.numeric(NA),J,R)
    for(i in 1:M)
      mu[,subgroup == i] <- log(rowSums(data[,subgroup == i])) -
                            log(sum(s[subgroup == i]))
    # (end of function)
  }
    
  # Get a rough estimate of log-lambda, which is useful for estimating
  # the range of psi2.
  #
  # CAN THIS BE A FUNCTION? e.g., estimate_psi2_range
  # (start of function)
  s.mat     <- rep(1,J) %*% t(s)
  loglambda <- log((data + 0.1)/s.mat)
  minpsi2   <- pmax(min(apply(loglambda,1,sd)^2)/100,1e-4)
  maxpsi2   <- max(apply(loglambda,1,sd)^2) 
  # (end of function)
  
  # Use grid search to initialize psi2 by fitting a poisson-log-normal
  # model while ignoring the unwanted variation.
  cat("Initializing psi2 via grid search.\n")
  t0 <- proc.time()
  psi2 <- init$psi2
  if (is.null(psi2)) {
      
    # CAN THIS BE A FUNCTION? e.g., initialize_psi2
    # (start of function)
    psi2 <- rep(as.numeric(NA),J)
    for (j in 1:J) {
      psi2_max       <- pmax(sd(loglambda[j,])^2,1)
      log2_psi2_grid <- seq(log2(1e-4),log2(psi2_max),length.out = 25)
      psi2_grid      <- 2^log2_psi2_grid
      logdens        <- rep(0,length(psi2_grid))
      for(l in 1:length(psi2_grid))
        for(r in 1:R)
          logdens[l] <- logdens[l] + log(dpoilog(data[j,r],mu[j,r] + log(s[r]),
                                                 sqrt(psi2_grid[l])))
      psi2[j] <- psi2_grid[which.max(logdens)]
    }
    # (end of function)
  }
  t1 <- proc.time()
  print(t1 - t0)
  
  # Initialize rho and bias.
  D   <- ncol(Fuv)
  rho <- init$rho
  if (is.null(rho))
    rho <- matrix(0,D,R)
  else
    rho <- as.matrix(rho)
  bias <- Fuv %*% rho
    
  # matrices to store the posterior mean and covariance of theta_j
  # (equal to eta_j in this case), i.e., gamma_j, diag(Sigma_j).
  gamma <- matrix(as.numeric(NA),J,R)
  Sigma <- matrix(as.numeric(NA),J,R)
  
  # Matrix to store the quantities related to q_j, s.t. A[j,r] =
  # gamma_jr + 0.5*Sigma_j,rr.
  A <- matrix(as.numeric(NA),J,R)
  
  # Update posterior mean and covariance of theta.
  t0 <- proc.time()
  for(j in 1:J) {
    eta.qj <- update_q_eta_only(x = data[j,],s = s,mu = mu[j,],bias = bias[j,],
                                c2 = rep(1,R),psi2 = psi2[j],
                                control = list(maxiter = maxiter.q,
                                               tol = tol.q))
    gamma[j,] <- eta.qj$m
    Sigma[j,] <- eta.qj$V
    A[j,]     <- eta.qj$m + eta.qj$V/2
  }
  t1    <- proc.time()
  t_eta <- t_eta + (t1 - t0)
  
  # Vector to store local ELBO for each j.
  ELBOs <- rep(as.numeric(NA),J)  
  
  # Overall ELBO after updating all parameters at each iteration.
  ELBOs.overall <- c()
  
  # CAN THIS BE A FUNCTION? e.g., compute_elbo_const
  # (start of function)
  const <- sum(data %*% log(s)) - sum(lgamma(data + 1))
  # (end of function)
  
  if (verbose)
    cat("Start prefitting Poisson mash model to initialize model",
        "parameters.\n")
  for (iter in 1:maxiter) {
      
    # Calculate or update quantities related to model parameters mu,
    # psi2, rho.
    tmp.mu   <- matrix(as.numeric(NA),J,M)
    tmp.psi2 <- rep(as.numeric(NA),J)
    
    for (j in 1:J) {
      gamma.tmp <- gamma[j,]
      Sigma.tmp <- Sigma[j,]
      for (i in 1:M) {
        tmp.mu[j,i] <- sum(s[subgroup == i] * exp(bias[j,subgroup == i] +
                           gamma.tmp[subgroup == i] +
                           Sigma.tmp[subgroup == i]/2))
      }
      tmp.psi2[j] <- sum(gamma.tmp^2 + Sigma.tmp)
    }
    
    # CAN THIS BE A FUNCTION? e.g., update_mu
    # (start of function)
    for(i in 1:M){
      mu.i.new <- log(rowSums(data[,subgroup == i])) - log(tmp.mu[,i])
      mu[,subgroup == i] <- mu.i.new
    }
    # (end of function)
    
    # CAN THIS BE A FUNCTION? e.g., update_psi2
    # (start of function)
    psi2.new <- tmp.psi2/R
    psi2     <- pmin(pmax(psi2.new,minpsi2),maxpsi2)
    # (end of function)
    
    # Update rho and bias.
    t0 <- proc.time()
    if (version == "R") {

      # CAN THIS BE A FUNCTION? e.g., update_rhos
      # (start of function)
      rho.new <- matrix(as.numeric(NA),nrow(rho),ncol(rho))
      for (r in 1:R)
        rho.new[,r] <- update_rho(Xr = data[,r],Fuv = Fuv,sr = s[r],
                                  mu = mu[,r],Lr = exp(A[,r]),init = rho[,r], 
                                  control = list(maxiter = 100,tol = tol.rho,
                                    maxrho = 100/max(abs(Fuv))))$rho
      # (end of function)
    } else 
      rho.new <- update_rho_rcpp(data,Fuv,s,mu,exp(A),rho,maxiter = 100,
                                 tol = tol.rho,maxrho = 100/max(abs(Fuv)))
    t1       <- proc.time()
    t_rho    <- t_rho + (t1 - t0)
    diff.rho <- rho.new - rho
    rho      <- rho.new
    bias     <- Fuv %*% rho
    
    # Update posterior mean and covariance of theta and local ELBO F_j.
    t0 <- proc.time()
    for (j in 1:J) {
      eta.qj <- update_q_eta_only(x = data[j,],s = s,mu = mu[j,],
                                  bias = bias[j,],c2 = rep(1,R),psi2 = psi2[j],
                                  init = list(m = gamma[j,],V = Sigma[j,]),
                                  control = list(maxiter=maxiter.q,tol=tol.q))
      gamma[j,] <- eta.qj$m
      Sigma[j,] <- eta.qj$V
      A[j,]     <- eta.qj$m + eta.qj$V/2
      ELBOs[j]  <- eta.qj$ELBO
    }
    t1    <- proc.time()
    t_eta <- t_eta + (t1 - t0)
    
    # Calculate overall ELBO at the current iteration.
    ELBO.overall  <- sum(ELBOs) + const
    ELBOs.overall <- c(ELBOs.overall, ELBO.overall)
    
    if (verbose) {
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO.overall))
    }
    if (iter >= 50)
      if (is.finite(ELBOs.overall[iter]) & is.finite(ELBOs.overall[iter - 1]))
        if (abs(ELBOs.overall[iter] -
                ELBOs.overall[iter-1])/abs(ELBOs.overall[iter-1]) < tol.stop)
          break
  }
  
  # Name the model paramter estimates.
  rownames(mu)  <- rownames(data)
  colnames(mu)  <- names(s)
  names(psi2)   <- rownames(data)
  colnames(rho) <- colnames(data)
  
  if(verbose)
    cat("Finish prefitting Poisson mash model to initialize model",
        "parameters.\n")
  print(t_rho)
  print(t_eta)
  return(list(mu = mu,psi2 = psi2,rho = rho,ELBO = ELBOs.overall,
              diff.rho = diff.rho))
}
