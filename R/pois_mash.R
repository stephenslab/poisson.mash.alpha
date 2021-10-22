#' @rdname pois_mash
#' 
#' @title Fit poisson mash to data
#' 
#' @param data \dQuote{pois.mash} data object, typically created by
#'   calling \code{\link{pois_mash_set_data}}.
#' 
#' @param Ulist List of H full-rank covariance matrices, such as the
#'   list of covariance matrices \code{Ulist} returned by
#'   \code{\link{pois_cov_ed}}).
#' 
#' @param ulist List of G numeric vectors each of which forms a
#'   rank-1 covariance matrix, such as the list of vectors \code{ulist}
#'   returned by \code{\link{pois_cov_ed}}.
#' 
#' @param ulist.epsilon2 Numeric vector of length G used to add a
#'   small positive value to the diagonals of each rank-1 prior
#'   covariance matrix to avoid tight error bars.
#' 
#' @param normalizeU Logical scalar indicating whether to normalize
#'   the prior covariances to have a maximum of 1 on diagonal.
#' 
#' @param gridmult Numeric scalar indicating factor by which
#'   adjacent grid values should differ; use a number close to 1 for
#'   fine grid.
#' 
#' @param wlist Numeric vector of length L giving the scaling
#'   factors for the prior covariance matrices
#' 
#' @param ruv Logical scalar indicating whether to account for
#'   unwanted variation. When \code{ruv = TRUE}, \code{Fuv} must be
#'   provided.
#' 
#' @param Fuv J x D matrix of latent factors causing unwanted
#'   variation, with features as rows and latent factors as columns.
#' 
#' @param rho D x R matrix of effects corresponding to unwanted
#'   variation, such that \code{bias = Fuv \%*\% rho}.
#' 
#' @param update.rho A logical scalar indicating whether to update
#'   effects corresponding to unwanted variation. Ignored if \code{ruv =
#'   FALSE}.
#' 
#' @param verbose A logical scalar indicating whether to print ELBO at
#'   each iteration.
#' 
#' @param C Q x R matrix of contrasts for effects. The default
#'   contrasts matrix is an matrix of condition-wise differences
#'   relative to the mean across all conditions.
#' 
#' @param res.colnames Character vector of length Q giving the names
#'   of the contrasts.
#'
#' @param init List of initial values for model parameters, such as an
#'   output from \code{\link{pois_mash_ruv_prefit}}).
#' 
#' @param version R (slower) and C++ (faster) implementations of the
#'   model fitting algorithm are provided; these are selected with
#'   \code{version = "R"} and \code{version = "Rcpp"}.
#' 
#' @param control List of control parameters with the following
#'   elements: \dQuote{maxiter}, maximum number of outer loop
#'   iterations; \dQuote{maxiter.q}, maximum number of inner loop
#'   iterations for updating the variational parameters at each outer
#'   loop iteration; \dQuote{maxpsi2}, maximum value for the
#'   gene-specific dispersion parameter psi2.; \dQuote{maxbias}, maximum
#'   value for the gene-specific range of bias caused by unwanted
#'   variation; \dQuote{tol.mu}, threshold for mu (gene-specific,
#'   subgroup-specific means on the log scale) to skip update;
#'   \dQuote{tol.psi2}, relative threshold for psi2 (gene-specific
#'   dispersion parameter) to skip update; \dQuote{tol.bias}, threshold
#'   for bias caused by unwanted variation to skip update;
#'   \dQuote{tol.q}, relative tolerance for assessing convergence of
#'   variational parameters at each iteration; \dQuote{tol.rho},
#'   tolerance for assessing convergence of effects corresponding to
#'   unwanted variation. Any named components will override the
#'   default optimization algorithm settings (as they are defined by
#'   \code{pois_mash_control_default}).
#' 
#' @return List with the following elements:
#'
#' \item{result}{List containing the posterior summaries of the J x Q
#'   matrix of effects.}
#'
#' \item{pois.mash.fit}{List containing the parameter estimates of the
#'   poisson mash model.}
#'
#' @importFrom poilog dpoilog
#' 
#' @export
#' 
pois_mash <- function (data, Ulist, ulist,
                       ulist.epsilon2 = rep(1e-8,length(ulist)),
                       normalizeU = TRUE, gridmult = 2, wlist, ruv = FALSE,
                       Fuv, rho, update.rho = TRUE, verbose = FALSE, C,
                       res.colnames, init = list(), version = c("Rcpp","R"),
                       control = list()) {
  
  s        <- data$s
  subgroup <- data$subgroup
  data     <- as.matrix(data$X)
  J        <- nrow(data)
  R        <- ncol(data)
  M        <- length(unique(subgroup))
  subgroup <- as.numeric(as.factor(subgroup))
  version   <- match.arg(version)
  control   <- modifyList(pois_mash_control_default(),control,keep.null = TRUE)
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  maxpsi2   <- control$maxpsi2
  maxbias   <- control$maxbias
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  tol.mu    <- control$tol.mu
  tol.psi2  <- control$tol.psi2
  tol.bias  <- control$tol.bias
  
  # Initialize rho and bias.
  if (ruv) {
    if (missing(Fuv))
      stop("The matrix Fuv must be provided if ruv is set to TRUE")
    Fuv <- as.matrix(Fuv)
    D   <- ncol(Fuv)
    if (missing(rho))
      rho <- matrix(0,D,R)
    else
      rho <- as.matrix(rho)
    diff.rho <- matrix(0,D,R)
    bias     <- Fuv %*% rho
    bias     <- scale_bias(bias,maxbias)
  }
  else {
    rho      <- NULL
    diff.rho <- NULL
    bias     <- matrix(0,J,R)
  }
  
  # Initialize mu by ignoring condition-specific effects (i.e.,
  # theta) and unwanted variation.
  mu <- init$mu
  if (is.null(mu))
    mu <- initialize_mu(data,s,subgroup,bias)
  
  # Get a rough estimate of log lambda, which is useful for estimating
  # the range of psi2.
  out       <- estimate_psi2_range(data,s,maxpsi2)
  loglambda <- out$loglambda
  minpsi2   <- out$minpsi2
  maxpsi2   <- out$maxpsi2
  
  # Use grid search to initialize psi2 by fitting a poisson-log-normal
  # model while ignoring fixed effects (i.e., beta) and unwanted
  # variation.
  psi2 <- init$psi2
  if (is.null(psi2))
    psi2 <- initialize_psi2(data,s,mu,bias)
  else
    psi2 <- pmin(psi2,maxpsi2)
  
  # Calculate wlist if not provided.
  if (missing(wlist)) {
    w_max      <- 4*max(apply(loglambda,1,sd)^2) 
    w_min      <- pmax(min(apply(loglambda,1,sd)^2)/100,1e-8)
    log2_wlist <- seq(log2(w_min),log2(w_max),by = log2(gridmult))
    wlist      <- 2^log2_wlist
  }
  
  H <- length(Ulist)
  G <- length(ulist)
  L <- length(wlist)
  K <- (H + G)*L
  
  # Normalize prior covariance matrices if normalizeU = TRUE.
  if (normalizeU) {
    if (H > 0) {
      for (h in 1:H) {
        Uh <- Ulist[[h]]
        Uh <- Uh/max(diag(Uh))
        Ulist[[h]] <- Uh
      }
    }
    
    for (g in 1:G) {
      ug <- ulist[[g]]
      if (sum(ug != 0) != 0) {
        ug <- ug/ug[which.max(abs(ug))]
        ulist[[g]] <- ug
      }
    }
  }
  
  # Initialize pi.
  pi <- init$pi
  if (is.null(pi))
    pi <- rep(1/K,K)
  
  const <- compute_elbo_const(data,s)
  
  if (verbose)
    cat("Start fitting Poisson mash model.\n")
  
  # J x K x R arrays to store the posterior mean gamma_jklr.
  gamma <- array(as.numeric(NA),c(J,K,R))
  
  # J x K x R arrays to store the quantities related to q_jkl,
  # s.t. A[j,kl,r] = gamma_jklr + 0.5 * Sigma_jkl,rr.
  A <- array(as.numeric(NA),c(J,K,R))
  
  # J x K matrix of local ELBO and temporary quantities needed to
  # update mu and psi2.
  ELBOs    <- matrix(0,J,K)
  tmp.mu   <- array(0,c(J,K,M))
  tmp.psi2 <- matrix(0,J,K)
  
  # Update posterior mean and covariance of theta and local ELBO.
  out <- update_q_by_j(data,s,subgroup,1:J,mu,bias,psi2,wlist,Ulist,
                       ulist,ulist.epsilon2,gamma,A,ELBOs,tmp.mu,tmp.psi2,
                       maxiter.q,tol.q)
  gamma    <- out$gamma
  A        <- out$A
  ELBOs    <- out$ELBOs
  tmp.mu   <- out$tmp.mu
  tmp.psi2 <- out$tmp.psi2
  
  # Update zeta.
  zeta <- update_zeta(ELBOs,pi)
  
  # Update J x R matrix tmp.ruv needed to update rho,
  # s.t. tmp.ruv[j,r] = sum_kl zeta[j,kl] * exp(A[j,kl,r]).
  tmp.ruv <- update_ruv(zeta,A)
  
  # Store the overall ELBO at each iteration.
  ELBOs.overall <- rep(0,maxiter)
  
  # Store the number of j to be updated at each iteration.
  j.update <- c()
  
  for (iter in 1:maxiter) {
      
    # Calculate overall ELBO at the current iteration.
    ELBOs.overall[iter] <- compute_overall_elbo(ELBOs,pi,zeta,const)
    
    # Update pi.
    pi.new  <- update_pi(zeta)
    diff.pi <- pi.new - pi
    pi      <- pi.new
    
    # Calculate the new mu.
    mu.new        <- update_mu(data,subgroup,zeta,tmp.mu)
    idx.update.mu <- apply(abs(mu.new - mu),1,max) > tol.mu
    diff.mu       <- mu.new - mu
    
    # Calculate the new psi2.
    psi2.new        <- update_psi2(zeta,tmp.psi2,R,minpsi2,maxpsi2)
    diff.psi2       <- psi2.new/psi2
    idx.update.psi2 <- abs(diff.psi2 - 1) > tol.psi2
    
    # Calculate the new rho and bias.
    if (ruv & update.rho) {
      rho.new <- update_rho_all(data,s,mu,Fuv,rho,tmp.ruv,tol.rho = tol.rho)
      diff.rho        <- rho.new - rho
      bias.new        <- Fuv %*% rho.new
      bias.new        <- scale_bias(bias.new,maxbias)
      idx.update.bias <- apply(abs(bias.new - bias),1,max) > tol.bias
      rho             <- rho.new
    }
    else {
      bias.new        <- bias
      idx.update.bias <- rep(FALSE,J)
    }
    
    # Update the set of indices j that need update, and update mu,
    # psi2, bias for these j.
    idx.update <- which(idx.update.mu | idx.update.psi2 | idx.update.bias)
    mu[idx.update,]   <- mu.new[idx.update,]
    psi2[idx.update]  <- psi2.new[idx.update]
    bias[idx.update,] <- bias.new[idx.update,]
    j.update          <- c(j.update,length(idx.update))
    
    if (verbose) {
      print("iter         ELBO")
      print(sprintf("%d:    %f",iter,ELBOs.overall[iter]))
      print("iter         number_of_j_to_update")
      print(sprintf("%d:    %d",iter,length(idx.update)))
    }
    
    if (length(idx.update) == 0)
      break
    
    # Update posterior mean and covariance of theta and local ELBO for
    # these j.
    out <- update_q_by_j(X=data, s=s, subgroup=subgroup, idx.update=idx.update, mu=mu, bias=bias, psi2=psi2,
                                     wlist=wlist, Ulist=Ulist, ulist=ulist, ulist.epsilon2=ulist.epsilon2,
                                     gamma=gamma, A=A, ELBOs=ELBOs, tmp.mu=tmp.mu, tmp.psi2=tmp.psi2, maxiter.q=maxiter.q, tol.q=tol.q)
    gamma    <- out$gamma
    A        <- out$A
    ELBOs    <- out$ELBOs
    tmp.mu   <- out$tmp.mu
    tmp.psi2 <- out$tmp.psi2    
    
    # Update zeta.
    zeta <- update_zeta(ELBOs,pi)
    
    # Update J x R matrix tmp.ruv needed to update rho,
    # s.t. tmp.ruv[j,r] = sum_kl zeta[j,kl] * exp(A[j,kl,r]).
    tmp.ruv <- update_ruv(zeta,A)
  }
  
  # Name the model paramter estimates.
  rownames(mu)   <- rownames(data)
  colnames(mu)   <- colnames(data)
  names(psi2)    <- rownames(data)
  rownames(zeta) <- rownames(data)
  if (ruv)
    colnames(rho) <- colnames(data)
  
  # Create the list with model parameters.
  pois.mash.fit <- list(mu = mu,psi2 = psi2,pi = pi,ulist = ulist,
                        ulist.epsilon2 = ulist.epsilon2,Ulist = Ulist,
                        wlist = wlist,zeta = zeta,Fuv = Fuv,rho = rho,
                        bias = bias,ELBO = ELBOs.overall,j.update = j.update)
  if (verbose) {
    cat("Finished fitting Poisson mash model.\n")
    cat("Start calculating posterior summary.\n")
  }
  
  # Calculate posterior summaries for the matrix of effects.
  result <- pois_mash_posterior(data = data,s = s,mu = mu,psi2 = psi2,
                                bias = bias,wlist = wlist,Ulist = Ulist,
                                ulist = ulist,ulist.epsilon2 = ulist.epsilon2, 
                                zeta = zeta,C = C,res.colnames = res.colnames)
  if (verbose)
    cat("Finish calculating posterior summary.\n")
  return(list(pois.mash.fit = pois.mash.fit,result = result))
}

#' @rdname pois_mash
#'
#' @export
#' 
pois_mash_control_default <- function()
  list(maxiter   = 500,
       maxiter.q = 25,
       maxbias   = 10,
       maxpsi2   = NULL,
       tol.mu    = 0.01,
       tol.psi2  = 0.02,
       tol.bias  = 0.01,
       tol.q     = 0.01,
       tol.rho   = 1e-6)
