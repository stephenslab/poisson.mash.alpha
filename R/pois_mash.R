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
#'   unwanted variation.
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
pois_mash <- function (data, Ulist, ulist, ulist.epsilon2 = NULL,
                       normalizeU = TRUE, gridmult = 2, wlist = NULL, 
                       ruv = FALSE, Fuv = NULL, rho = NULL, update.rho = TRUE,
                       verbose = FALSE, C = NULL, res.colnames = NULL,
                       init = list(NULL),
                       control = list(maxiter = 500, maxiter.q = 25,
                                      maxpsi2 = NULL, maxbias = 10,
                                      tol.mu = 0.01, tol.psi2 = 0.02,
                                      tol.bias = 0.01, tol.q = 0.01,
                                      tol.rho = 1e-6)) {
  
  s        <- data$s
  subgroup <- data$subgroup
  data     <- as.matrix(data$X)
  J        <- nrow(data)
  R        <- ncol(data)
  M        <- length(unique(subgroup))
  subgroup <- as.numeric(as.factor(subgroup))
  
  # Check if ulist is empty.
  if (is.null(ulist))
    stop("ulist cannot be empty!")
  
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  maxpsi2   <- control$maxpsi2
  maxbias   <- control$maxbias
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  tol.mu    <- control$tol.mu
  tol.psi2  <- control$tol.psi2
  tol.bias  <- control$tol.bias
  
  if (is.null(maxiter))
    maxiter <- 500
  if (is.null(maxiter.q))
    maxiter.q <- 25
  if (is.null(maxbias))
    maxbias <- 10
  if (is.null(tol.q))
    tol.q <- 0.01
  if (is.null(tol.rho))
    tol.rho <- 1e-6
   if (is.null(tol.mu))
    tol.mu <- 0.01
  if (is.null(tol.psi2))
    tol.psi2 <- 0.02
  if (is.null(tol.bias))
    tol.bias <- 0.01
  
  # Initialize rho and bias.
  if (ruv) {
    if (is.null(Fuv))
      stop("The matrix Fuv must be provided if ruv is set to TRUE")
    Fuv <- as.matrix(Fuv)
    D   <- ncol(Fuv)
    if (is.null(rho))
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
  if (is.null(wlist)) {
    if (is.null(gridmult))
      gridmult <- 2
    w_max      <- 4*max(apply(loglambda,1,sd)^2) 
    w_min      <- pmax(min(apply(loglambda,1,sd)^2)/100,1e-8)
    log2_wlist <- seq(log2(w_min),log2(w_max),by = log2(gridmult))
    wlist      <- 2^log2_wlist
  }
  
  H <- length(Ulist)
  G <- length(ulist)
  L <- length(wlist)
  K <- (H + G)*L
  
  # Specify ulist.epsilon2 if not provided.
  if (is.null(ulist.epsilon2))
    ulist.epsilon2 <- rep(1e-8,G)
  
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
  
  # J x K matrix of local ELBO and temporary quantities needed to update mu and psi2.
  ELBOs    <- matrix(0,J,K)
  tmp.mu   <- array(0,c(J,K,M))
  tmp.psi2 <- matrix(0,J,K)
  
  # Update posterior mean and covariance of theta and local ELBO.
  # CAN THIS BE A FUNCTION? e.g., update_q_by_j.
  res.update.by.j <- update_q_by_j(X=data, s=s, subgroup=subgroup, idx.update=1:J, mu=mu, bias=bias, psi2=psi2,
                                   wlist=wlist, Ulist=Ulist, ulist=ulist, ulist.epsilon2=ulist.epsilon2,
                                   gamma=gamma, A=A, ELBOs=ELBOs, tmp.mu=tmp.mu, tmp.psi2=tmp.psi2, maxiter.q=maxiter.q, tol.q=tol.q)
  gamma <- res.update.by.j$gamma
  A <- res.update.by.j$A
  ELBOs <- res.update.by.j$ELBOs
  tmp.mu <- res.update.by.j$tmp.mu
  tmp.psi2 <- res.update.by.j$tmp.psi2
  
  # Update zeta.
  # CAN THIS BE A FUNCTION? e.g., update_zeta.
  zeta <- update_zeta(ELBOs=ELBOs, pi=pi)
  
  # Update J x R matrix tmp.ruv needed to update rho,
  # s.t. tmp.ruv[j,r] = sum_kl zeta[j,kl] * exp(A[j,kl,r]).
  tmp.ruv <- matrix(as.numeric(NA),J,R)
  for (r in 1:R)
    tmp.ruv[,r] <- rowSums(zeta*exp(A[,,r]))
  
  # Store the overall ELBO at each iteration.
  ELBOs.overall <- c()
  
  # Store the number of j to be updated at each iteration.
  j.update <- c()
  
  for (iter in 1:maxiter) {
      
    # Calculate overall ELBO at the current iteration.
    ELBO.overall  <- compute_overall_elbo(ELBOs,pi,zeta,const)
    ELBOs.overall <- c(ELBOs.overall,ELBO.overall)
    
    # Update pi.
    pi.new  <- update_pi(zeta)
    diff.pi <- pi.new - pi
    pi      <- pi.new
    
    # Calculate the new mu.
    # CAN THIS BE A FUNCTION? e.g., update_mu.
    mu.new <- update_mu(X=data, subgroup=subgroup, zeta=zeta, tmp.mu=tmp.mu)
    idx.update.mu <- apply(abs(mu.new - mu),1,max) > tol.mu
    diff.mu       <- mu.new - mu
    
    # Calculate the new psi2.
    psi2.new        <- update_psi2(zeta,tmp.psi2,R,minpsi2,maxpsi2)
    diff.psi2       <- psi2.new/psi2
    idx.update.psi2 <- abs(diff.psi2 - 1) > tol.psi2
    
    # Calculate the new rho and bias.
    if (ruv & update.rho) {
      # CAN THIS BE A FUNCTION? e.g., update_rho_all.
      rho.new <- update_rho_all(X=data, s=s, mu=mu, Fuv=Fuv, rho=rho, tmp.ruv=tmp.ruv, tol.rho=tol.rho)
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
      print(sprintf("%d:    %f", iter, ELBO.overall))
      print("iter         number_of_j_to_update")
      print(sprintf("%d:    %d", iter, length(idx.update)))
    }
    
    if (length(idx.update) == 0)
      break
    
    # Update posterior mean and covariance of theta and local ELBO for
    # these j.
    # CAN THIS BE A FUNCTION? e.g., update_q_by_j.
    res.update.by.j <- update_q_by_j(X=data, s=s, subgroup=subgroup, idx.update=idx.update, mu=mu, bias=bias, psi2=psi2,
                                     wlist=wlist, Ulist=Ulist, ulist=ulist, ulist.epsilon2=ulist.epsilon2,
                                     gamma=gamma, A=A, ELBOs=ELBOs, tmp.mu=tmp.mu, tmp.psi2=tmp.psi2, maxiter.q=maxiter.q, tol.q=tol.q)
    gamma <- res.update.by.j$gamma
    A <- res.update.by.j$A
    ELBOs <- res.update.by.j$ELBOs
    tmp.mu <- res.update.by.j$tmp.mu
    tmp.psi2 <- res.update.by.j$tmp.psi2    
    
    # Update zeta.
    # CAN THIS BE A FUNCTION? e.g., update_zeta.
    zeta <- update_zeta(ELBOs=ELBOs, pi=pi)
    
    # Update J x R matrix tmp.ruv needed to update rho,
    # s.t. tmp.ruv[j,r] = sum_kl zeta[j,kl] * exp(A[j,kl,r]).
    for (r in 1:R)
      tmp.ruv[,r] <- rowSums(zeta*exp(A[,,r]))
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
