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
  if (is.null(mu)) {
    mu <- matrix(as.numeric(NA),J,R)
    for (i in 1:M)
      mu[,subgroup == i] <- log(rowSums(data[,subgroup == i])) -
          log(exp(bias[,subgroup == i]) %*% s[subgroup == i])
  }
  
  # Get a rough estimate of log lambda, which is useful for estimating
  # the range of psi2.
  s.mat     <- rep(1,J) %*% t(s)
  loglambda <- log((data + 0.1)/s.mat)
  minpsi2   <- pmax(min(apply(loglambda,1,sd)^2)/100,1e-8)
  if (is.null(maxpsi2))
    maxpsi2 <- max(apply(loglambda,1,sd)^2)
  
  # Use grid search to initialize psi2 by fitting a poisson-log-normal
  # model while ignoring fixed effects (i.e., beta) and unwanted
  # variation.
  psi2 <- init$psi2
  if (is.null(psi2)) {
    psi2 <- rep(as.numeric(NA),J)
    for (j in 1:J) {
      psi2_max       <- pmax(sd(loglambda[j,])^2,1)
      log2_psi2_grid <- seq(log2(1e-4),log2(psi2_max),length.out = 25)
      psi2_grid      <- 2^log2_psi2_grid
      logdens        <- rep(0,length(psi2_grid))
      for (l in 1:length(psi2_grid))
        for (r in 1:R)
          logdens[l] <- logdens[l] +
            log(dpoilog(data[j,r],mu[j,r] + bias[j,r] + log(s[r]),
                        sqrt(psi2_grid[l])))
      psi2[j] <- psi2_grid[which.max(logdens)]
    }
  }
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
  const <- sum(data %*% log(s)) - sum(lgamma(1 + data))
  if (verbose)
    cat("Start fitting Poisson mash model.\n")
  
  # J x K x R arrays to store the posterior mean gamma_jklr.
  gamma <- array(as.numeric(NA),c(J,K,R))
  
  # J x K x R arrays to store the quantities related to q_jkl,
  # s.t. A[j,kl,r] = gamma_jklr + 0.5 * Sigma_jkl,rr.
  A <- array(as.numeric(NA),c(J,K,R))
  
  # Update posterior mean and covariance of theta and local ELBO.
  ELBOs    <- matrix(0,J,K)
  tmp.mu   <- array(0,c(J,K,M))
  tmp.psi2 <- matrix(0,J,K)
  
  for (j in 1:J) {
    if (H > 0) {
      hl <- 0
      for (h in 1:H) {
        for (l in 1:L) {
          hl <- hl + 1
          theta.qjhl <- update_q_theta_general(x = data[j,],s = s,mu = mu[j,],
            bias = bias[j,],c2 = rep(1,R),psi2 = psi2[j],w = wlist[l],
            U = Ulist[[h]],control = list(maxiter = maxiter.q,tol = tol.q))
          ELBOs[j,hl] <- theta.qjhl$ELBO
          gamma.tmp   <- theta.qjhl$m
          Sigma.tmp   <- theta.qjhl$V
          for (i in 1:M)
            tmp.mu[j,hl,i] <- sum(s[subgroup == i]*exp(bias[j,subgroup == i] +
                              gamma.tmp[subgroup == i] +
                              diag(Sigma.tmp)[subgroup == i]/2))
          eta.qjhl <- update_q_eta_general(theta_m = gamma.tmp,
                                           theta_V = Sigma.tmp,c2 = rep(1,R),
                                           psi2 = psi2[j],w = wlist[l],
                                           U = Ulist[[h]])
          tmp.psi2[j,hl] <- sum(eta.qjhl)
          gamma[j,hl,]   <- gamma.tmp
          A[j,hl,]       <- gamma.tmp + diag(Sigma.tmp)/2 
        }
      }      
    }
    
    gl <- 0
    for (g in 1:G) {
      ug         <- ulist[[g]]
      epsilon2.g <- ulist.epsilon2[g]
      for (l in 1:L) {
        gl <- gl + 1
        theta.qjgl <- update_q_theta_rank1(x = data[j,],s = s,mu = mu[j,],
                        bias = bias[j,],c2 = rep(1,R),
                        psi2 = psi2[j] + wlist[l] * epsilon2.g, 
                        w = wlist[l],u = ug,
                        control = list(maxiter = maxiter.q,tol = tol.q))
        ELBOs[j,H*L+gl] <- theta.qjgl$ELBO
        gamma.tmp       <- theta.qjgl$m
        Sigma.tmp       <- theta.qjgl$V
        for (i in 1:M)
          tmp.mu[j,H*L+gl,i] <-
            sum(s[subgroup == i] * exp(bias[j,subgroup == i] +
            gamma.tmp[subgroup == i] + diag(Sigma.tmp)[subgroup == i]/2))
        gamma[j,H*L+gl,] <- gamma.tmp
        A[j,H*L+gl,]     <- gamma.tmp + diag(Sigma.tmp)/2
        
        # If ug is zero vector.
        if (sum(ug != 0) == 0) {
          eta.qjgl <- gamma.tmp^2 + diag(Sigma.tmp) 
          tmp.psi2[j,H*L+gl] <- sum(eta.qjgl)
        }
        else if (epsilon2.g > 1e-4) {
          eta.qjgl <- update_q_eta_rank1_robust(theta_m = gamma.tmp,
                                                theta_V = Sigma.tmp,
                                                c2 = rep(1,R),psi2 = psi2[j],
                                                w = wlist[l],u = ug,
                                                epsilon2 = epsilon2.g)
          tmp.psi2[j,H*L+gl] <- sum(eta.qjgl)
        }
        else {
          beta.qjgl <- update_q_beta_rank1(theta_m = gamma.tmp,
                                           theta_V = Sigma.tmp,c2 = rep(1,R),
                                           psi2 = psi2[j],w = wlist[l],u = ug)
          eta.qjgl <- update_q_eta_rank1(theta_m = gamma.tmp,
                                         theta_V = Sigma.tmp,
                                         a2_m = beta.qjgl$a2_m,
                                         a_theta_m = beta.qjgl$a_theta_m,
                                         u = ug)
          tmp.psi2[j,H*L+gl] <- sum(eta.qjgl)
        }
      }
    }
  }
  
  # Update zeta.
  ELBOs.cen <- ELBOs - apply(ELBOs,1,max)
  zeta      <- t(t(exp(ELBOs.cen)) * pi)
  zeta      <- zeta*(1/rowSums(zeta))  
  zeta      <- pmax(zeta,1e-15)
  
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
    ELBO.overall <- sum(zeta*(log(rep(1,J) %*% t(pi)) + ELBOs - log(zeta))) +
                    const
    ELBOs.overall <- c(ELBOs.overall,ELBO.overall)
    
    # Update pi.
    pi.new  <- colMeans(zeta)
    pi.new  <- pmax(pi.new,1e-8)
    diff.pi <- pi.new - pi
    pi      <- pi.new
    
    # Calculate the new mu.
    mu.new <- matrix(as.numeric(NA),J,R)
    for (i in 1:M) {
      mu.i.new <- log(rowSums(data[,subgroup == i])) -
                  log(rowSums(zeta * tmp.mu[,,i]))
      mu.new[,subgroup == i] <- mu.i.new
    }
    idx.update.mu <- apply(abs(mu.new - mu),1,max) > tol.mu
    diff.mu       <- mu.new - mu
    
    # Calculate the new psi2.
    psi2.new        <- rowSums(zeta * tmp.psi2)/R
    psi2.new        <- pmin(pmax(psi2.new,minpsi2),maxpsi2)
    diff.psi2       <- psi2.new/psi2
    idx.update.psi2 <- abs(diff.psi2 - 1) > tol.psi2
    
    # Calculate the new rho and bias.
    if (ruv & update.rho) {
      rho.new <- matrix(as.numeric(NA),nrow(rho),ncol(rho))
      for (r in 1:R)
        rho.new[,r] <- update_rho(Xr = data[,r],Fuv = Fuv,sr = s[r],
                                  mu = mu[,r],Lr = tmp.ruv[,r],init = rho[,r],
                                  control = list(maxiter = 100,tol = tol.rho,
                                    maxrho = 100/max(abs(Fuv))))$rho 
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
    for (j in 1:length(idx.update)) {
      j.idx <- idx.update[j]
      
      if (H > 0) {
        hl <- 0
        for (h in 1:H) {
          for (l in 1:L) {
            hl         <- hl + 1
            Utilde     <- wlist[l] * Ulist[[h]] + psi2[j.idx] * diag(R)
            a.tmp      <- s * exp(mu[j.idx,] + bias[j.idx,] + A[j.idx,hl,])
            theta.qjhl <- update_q_theta_general(x = data[j.idx,],s = s,
                            mu = mu[j.idx,],bias = bias[j.idx,],c2 = rep(1,R),
                            psi2 = psi2[j.idx],w = wlist[l],U = Ulist[[h]],
                            init = list(m = gamma[j.idx,hl,],
                                        V = solve(solve(Utilde) + diag(a.tmp),
                                                  tol = 1e-50)), 
                                control = list(maxiter=maxiter.q,tol=tol.q))
            ELBOs[j.idx,hl] <- theta.qjhl$ELBO
            gamma.tmp       <- theta.qjhl$m
            Sigma.tmp       <- theta.qjhl$V
            for (i in 1:M)
              tmp.mu[j.idx,hl,i] <- sum(s[subgroup == i] *
                                      exp(bias[j.idx,subgroup == i] +
                                          gamma.tmp[subgroup == i] +
                                          diag(Sigma.tmp)[subgroup == i]/2))
            eta.qjhl <- update_q_eta_general(theta_m = gamma.tmp,
                                             theta_V = Sigma.tmp,c2 = rep(1,R),
                                             psi2 = psi2[j.idx],w = wlist[l],
                                             U = Ulist[[h]])
            tmp.psi2[j.idx,hl] <- sum(eta.qjhl)
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
          gl     <- gl + 1
          a.tmp  <- s * exp(mu[j.idx,] + bias[j.idx,] + A[j.idx,H*L+gl,])
          S_inv  <- 1/(psi2[j.idx] * rep(1,R) + wlist[l] * epsilon2.g)
          init.V <- mat_inv_rank1(a.tmp + S_inv,-wlist[l] * ug * S_inv,
                                  (ug*S_inv)/(1 + wlist[l]*sum(ug^2*S_inv)))
          theta.qjgl <- update_q_theta_rank1(x = data[j.idx,],s = s,
                          mu = mu[j.idx,],bias = bias[j.idx,],c2 = rep(1,R),
                          psi2 = psi2[j.idx] + wlist[l] * epsilon2.g, 
                          w = wlist[l],u = ug,
                          init = list(m = gamma[j.idx,H*L+gl,],V = init.V),
                          control = list(maxiter = maxiter.q,tol = tol.q))
          ELBOs[j.idx,H*L+gl] <- theta.qjgl$ELBO
          gamma.tmp           <- theta.qjgl$m
          Sigma.tmp           <- theta.qjgl$V
          for (i in 1:M)
            tmp.mu[j.idx,H*L+gl,i] <-
              sum(s[subgroup == i] * exp(bias[j.idx,subgroup == i] +
              gamma.tmp[subgroup == i] + diag(Sigma.tmp)[subgroup == i]/2))
          gamma[j.idx,H*L+gl,] <- gamma.tmp
          A[j.idx,H*L+gl,]     <- gamma.tmp + diag(Sigma.tmp)/2
          
          # If ug is zero vector.
          if (sum(ug != 0) == 0) {
            eta.qjgl <- gamma.tmp^2 + diag(Sigma.tmp) 
            tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
          }
          else if (epsilon2.g > 1e-4) {
            eta.qjgl <- update_q_eta_rank1_robust(theta_m = gamma.tmp,
                          theta_V = Sigma.tmp,c2 = rep(1,R),psi2 = psi2[j.idx],
                          w = wlist[l],u = ug,epsilon2 = epsilon2.g)
            tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
          }
          else {
            beta.qjgl <- update_q_beta_rank1(theta_m = gamma.tmp,
                           theta_V = Sigma.tmp,c2 = rep(1,R),
                           psi2 = psi2[j.idx],w = wlist[l],u = ug)
            eta.qjgl <- update_q_eta_rank1(theta_m = gamma.tmp,
                           theta_V = Sigma.tmp,a2_m = beta.qjgl$a2_m,
                           a_theta_m = beta.qjgl$a_theta_m,u = ug)
            tmp.psi2[j.idx,H*L+gl] <- sum(eta.qjgl)
          }
        }
      }      
    }
    
    # Update zeta.
    ELBOs.cen <- ELBOs - apply(ELBOs, 1, max)
    zeta      <- t(t(exp(ELBOs.cen)) * pi)
    zeta      <- zeta*(1/rowSums(zeta))  
    zeta      <- pmax(zeta, 1e-15)
    
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
