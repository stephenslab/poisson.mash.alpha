#' @rdname pois_mash_ruv_prefit
#' 
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
#' @param control List of control parameters with one or more of the
#'   following elements: \dQuote{maxiter}, maximum number of iterations;
#'   \dQuote{maxiter.q}, maximum number of inner-loop iterations to
#'   update variational parameters at each iteration; \dQuote{tol.stop},
#'   tolerance for assessing convergence of prefit, as measured by
#'   relative change in ELBO; \dQuote{tol.q}, relative tolerance for
#'   assessing convergence of variational parameters at each iteration;
#'   \dQuote{tol.rho}, tolerance for assessing convergence of effects
#'   corresponding to unwanted variation. Any named components will
#'   override the default optimization algorithm settings (as they are
#'   defined by \code{pois_mash_ruv_prefit_control_default}).
#' 
#' @return A list containing initial estimates of model parameters.
#'
#' @importFrom utils modifyList
#' @importFrom stats sd
#' @importFrom poilog dpoilog
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib poisson.mash.alpha
#' 
#' @export
#' 
pois_mash_ruv_prefit <- function (data, Fuv, verbose = FALSE,
                                  init = list(),
                                  control = list()) {

  s         <- data$s
  subgroup  <- data$subgroup
  data      <- as.matrix(data$X)
  J         <- nrow(data)
  R         <- ncol(data)
  M         <- length(unique(subgroup))
  subgroup  <- as.numeric(as.factor(subgroup))
  Fuv       <- as.matrix(Fuv)
  control   <- modifyList(pois_mash_ruv_prefit_control_default(),control,
                          keep.null = TRUE)

  # Get the optimization settings.
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  tol.stop  <- control$tol.stop

  # Initialize mu by ignoring random effects and unwanted variation.
  mu <- init$mu
  if (is.null(mu))
    mu <- initialize_mu(data,s,subgroup)
    
  # Get a rough estimate of log-lambda, which is useful for estimating
  # the range of psi2.
  out     <- estimate_psi2_range(data,s,epsilon = 1e-4)
  minpsi2 <- out$minpsi2
  maxpsi2 <- out$maxpsi2
  
  # Use grid search to initialize psi2 by fitting a poisson-log-normal
  # model while ignoring the unwanted variation.
  cat("Initializing psi2 via grid search.\n")
  psi2 <- init$psi2
  if (is.null(psi2))
    psi2 <- initialize_psi2(data,s,mu)
  
  # Initialize rho and bias.
  D   <- ncol(Fuv)
  rho <- init$rho
  if (is.null(rho))
    rho <- matrix(0,D,R)
  else
    rho <- as.matrix(rho)
  bias <- Fuv %*% rho
    
  # Matrices to store the posterior mean and covariance of theta_j
  # (equal to eta_j in this case), i.e., gamma_j, diag(Sigma_j).
  gamma <- matrix(as.numeric(NA),J,R)
  Sigma <- matrix(as.numeric(NA),J,R)
  
  # Matrix to store the quantities related to q_j,
  # s.t. A[j,r] = gamma_jr + 0.5 * Sigma_j,rr.
  A <- matrix(as.numeric(NA),J,R)
  
  # Update posterior mean and covariance of theta.
  for (j in 1:J) {
    out <- update_q_eta_only(data[j,],s,mu[j,],bias[j,],rep(1,R),psi2[j],
                             maxiter = maxiter.q,tol = tol.q)
    gamma[j,] <- out$m
    Sigma[j,] <- out$V
    A[j,]     <- out$m + out$V/2
  }
  
  # Vector to store local ELBO for each j.
  ELBOs         <- rep(as.numeric(NA),J)  
  const         <- compute_elbo_const(data,s)
  ELBOs.overall <- rep(0,maxiter)
  
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
        k <- which(subgroup == i)
        tmp.mu[j,i] <- sum(compute_poisson_rates(s[k],0,bias[j,k],gamma.tmp[k],
                                                 Sigma.tmp[k]))
      }
      tmp.psi2[j] <- sum(gamma.tmp^2 + Sigma.tmp)
    }
    
    # Update mu.
    for (i in 1:M) {
      k        <- which(subgroup == i)
      mu.i.new <- log(rowSums(data[,k])) - log(tmp.mu[,i])
      mu[,k]   <- mu.i.new
    }
    
    # Update the dispersion parameter psi2.
    psi2.new <- tmp.psi2/R
    psi2     <- pmin(pmax(psi2.new,minpsi2),maxpsi2)
    
    # Update rho and bias.
    rho.new  <- update_rho_all(data,s,mu,Fuv,rho,exp(A),tol = tol.rho)
    diff.rho <- rho.new - rho
    rho      <- rho.new
    bias     <- Fuv %*% rho
    
    # Update posterior mean and covariance of theta and local ELBO F_j.
    for (j in 1:J) {
      out <- update_q_eta_only(data[j,],s,mu[j,],bias[j,],rep(1,R),psi2[j],
                               init = list(m = gamma[j,],V = Sigma[j,]),
                               maxiter = maxiter.q,tol = tol.q)
      gamma[j,] <- out$m
      Sigma[j,] <- out$V
      A[j,]     <- out$m + out$V/2
      ELBOs[j]  <- out$ELBO
    }
    
    # Calculate overall ELBO at the current iteration.
    ELBOs.overall[iter] <- sum(ELBOs) + const
    
    if (verbose) {
      print("iter         ELBO")
      print(sprintf("%d:    %f",iter,ELBOs.overall[iter]))
    }

    # Check the convergence criteria.
    #
    # NOTE: Consider fixing this to (1) check for increases in ELBO;
    # and (2) compare absolute change (not relative change).
    # 
    if (iter >= 50)
      if (is.finite(ELBOs.overall[iter]) & is.finite(ELBOs.overall[iter-1]))
        if (abs(ELBOs.overall[iter] -
                ELBOs.overall[iter-1])/abs(ELBOs.overall[iter-1]) < tol.stop)
          break
  }
  
  # Add names to the parameter estimates.
  rownames(mu)  <- rownames(data)
  colnames(mu)  <- names(s)
  names(psi2)   <- rownames(data)
  colnames(rho) <- colnames(data)
  
  if (verbose)
    cat("Finish prefitting Poisson mash model to initialize model",
        "parameters.\n")
  return(list(mu       = mu,
              psi2     = psi2,
              rho      = rho,
              ELBO     = ELBOs.overall[1:iter],
              diff.rho = diff.rho))
}

#' @rdname pois_mash_ruv_prefit
#'
#' @export
#' 
pois_mash_ruv_prefit_control_default <- function()
  list(maxiter   = 100,
       maxiter.q = 25,
       tol.q     = 0.01,
       tol.rho   = 1e-4,
       tol.stop  = 1e-6)
