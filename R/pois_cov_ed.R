#' @rdname pois_cov_ed
#' 
#' @title Estimate Prior Covariance Matrices using Extreme Deconvolution
#'
#' @description Perform extreme deconvolution (ED) to estimate
#'   data-driven prior covariance matrices.
#' 
#' @param data \dQuote{pois.mash} data object, typically created by
#'   calling \code{\link{pois_mash_set_data}}.
#' 
#' @param subset The indices of features to be used. Defaults to using
#'   all features.
#' 
#' @param Ulist A list of H full-rank covariance matrices (e.g.,
#'   initialized by \code{pois_cov_init}).
#' 
#' @param ulist A list of G numeric vectors each of which defines a
#'   rank-1 covariance matrix.
#' 
#' @param ulist.dd Logical vector of length G denoting whether each
#'   element in \code{ulist} is data-driven (\code{TRUE}) or canonical
#'   (\code{FALSE}). Defaults to data-driven for all elements. For
#'   canonical covariances, the spanned space is not updated.
#' 
#' @param ruv Logical scalar indicating whether to account for
#'   unwanted variation. If \code{ruv = TRUE}, \code{Fuv} must be
#'   provided.
#' 
#' @param Fuv J x D matrix of latent factors causing unwanted
#'   variation, with features as rows and latent factors as columns.
#' 
#' @param verbose Logical scalar indicating whether to print ELBO at
#'   each iteration.
#' 
#' @param init Optional list of initial values for model parameters
#'   (e.g., returned by \code{pois_mash_ruv_prefit}).
#' 
#' @param version R (slower) and C++ (faster) implementations of the
#'   model fitting algorithm are provided; these are selected with
#'   \code{version = "R"} and \code{version = "Rcpp"}.
#' 
#' @param control A list of control parameters with the following
#'   elements: \dQuote{maxiter}, maximum number of ED iterations;
#'   \dQuote{maxiter.q}, maximum number of inner loop iterations to
#'   update variational parameters at each ED iteration;
#'   \dQuote{maxpsi2}, maximum for the gene-specific dispersion
#'   parameter \code{psi2}; \dQuote{maxbias}, maximum for the
#'   gene-specific range of bias caused by unwanted variation;
#'   \dQuote{tol.stop}, tolerance for assessing convergence of ED, as
#'   measured by relative change in ELBO; \dQuote{tol.q}, relative
#'   tolerance for assessing convergence of variational parameters at
#'   each ED iteration; and \dQuote{tol.rho}{tolerance for assessing
#'   convergence of effects corresponding to unwanted variation.} Any
#'   named components will override the default optimization algorithm
#'   settings (as they are defined by
#'   \code{pois_cov_ed_control_default}).
#' 
#' @return A list including the following elements:
#' 
#' \item{Ulist}{List of H full-rank covariance matrices.}
#' 
#' \item{ulist}{List of G numeric vectors each of which forms a
#'   rank-1 covariance matrix.}
#' 
#' \item{pi}{Numeric vector of length H + G containing the mixture
#'   proportions for Ulist and ulist.}
#'
#' @importFrom utils modifyList
#' @importFrom stats sd
#' @importFrom poilog dpoilog
#' 
#' @export
#' 
pois_cov_ed <- function (data, subset, Ulist, ulist, ulist.dd,
                         ruv = FALSE, Fuv, verbose = FALSE,
                         init = list(), version = c("Rcpp","R"),
                         control = list()) {
  X        <- data$X
  s        <- data$s
  subgroup <- data$subgroup
  if (missing(subset))
    subset <- 1:nrow(X)
  data.ed   <- as.matrix(X[subset,])
  J         <- nrow(data.ed)
  R         <- ncol(data.ed)
  M         <- length(unique(subgroup))
  H         <- length(Ulist)
  G         <- length(ulist)
  K         <- H + G
  subgroup  <- as.numeric(as.factor(subgroup))
  version   <- match.arg(version)
  control   <- modifyList(pois_cov_ed_control_default(),control,
                          keep.null = TRUE)
  
  # Get the optimization settings.
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  maxpsi2   <- control$maxpsi2
  maxbias   <- control$maxbias
  tol.stop  <- control$tol.stop
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  
  if (missing(ulist.dd)) {
    ulist.dd <- rep(TRUE,G)
    for (g in 1:G) 
      if (sum(ulist[[g]] != 0) == 0)
        ulist.dd[g] <- FALSE
  }
  
  # Initialize pi_h and pi_g.
  pi <- init$pi
  if (is.null(pi))
    pi <- rep(1/K,K)
  
  # Initialize rho and bias.
  rho      <- init$rho
  diff.rho <- NULL
  bias     <- matrix(0,J,R)
  
  if (ruv) {
    if (missing(Fuv))
      stop("The matrix Fuv must be provided if ruv is set to TRUE")
    F.ed <- as.matrix(Fuv[subset,])
    D    <- ncol(F.ed)
    if (is.null(rho))
      rho <- matrix(0,D,R)
    else
      rho <- as.matrix(rho)
    bias <- F.ed %*% rho
    bias <- scale_bias(bias,maxbias)
  }
  
  # Initialize mu by ignoring condition-specific effects (i.e.,
  # theta) and unwanted variation.
  mu <- init$mu
  if (is.null(mu))
    mu <- initialize_mu(data.ed,s,subgroup)
  else
    mu <- mu[subset,]
  
  # Get a rough estimate of log-lambda, which is useful for estimating
  # the range of psi2, Ulist, ulist.
  out     <- estimate_psi2_range(data.ed,s,maxpsi2)
  minpsi2 <- out$minpsi2
  maxpsi2 <- out$maxpsi2
  upr_bd  <- out$upr_bd
  
  # Use grid search to initialize psi^2 by fitting a
  # poisson-log-normal model while ignoring fixed effects (i.e.,
  # beta_j) and unwanted variation.
  psi2 <- init$psi2
  if (is.null(psi2))
    psi2 <- initialize_psi2(data.ed,s,mu)
  else
    psi2 <- pmin(psi2[subset],maxpsi2)
  
  # Create matrices and arrays to store the posterior mean and
  # covariance of theta, i.e., gamma_jk, Sigma_jk.
  gamma_jk <- array(as.numeric(NA),c(J,K,R))
  Sigma_jk <- list(NULL)
  for (j in 1:J)
    Sigma_jk[[j]] <- array(as.numeric(NA),c(K,R,R))
  
  # Create a J x K x R array to store the quantities related to
  # q_jk, s.t. A[j,k,r] = gamma_jkr + 0.5*Sigma_jk,rr
  A <- array(as.numeric(NA),c(J,K,R))
  
  # Matrix to store "local" ELBOs F_jk.
  ELBOs <- matrix(0,J,K)
  
  # Update posterior mean and covariance of theta, and the "local" ELBO.
  for (j in 1:J) {
    out <- update_q_theta_all(data.ed[j,],s,mu[j,],bias[j,],rep(1,R),
                              psi2[j],1,Ulist,ulist,list(),maxiter.q,tol.q)
    gamma_jk[j,,] <- out$gamma
    Sigma_jk[[j]] <- out$Sigma
    A[j,,]        <- out$A
    ELBOs[j,]     <- out$ELBOs
  }

  # Update J x K matrix zeta of posterior weights.
  zeta <- update_zeta(ELBOs,pi)
  
  # Update J x R matrix tmp.ruv needed to update rho.
  tmp.ruv <- update_ruv(zeta,A)

  const <- compute_elbo_const(data.ed,s)
  
  # Overall ELBO after updating all parameters at each iteration.
  ELBOs.overall <- rep(0,maxiter)
  
  if (verbose)
    cat("Start running extreme deconvolution to estimate prior covariance",
        "matrices.\n")

  for (iter in 1:maxiter) {
      
    # Calculate overall ELBO at the current iteration.
    ELBOs.overall[iter] <- compute_overall_elbo(ELBOs,pi,zeta,const)
    if (verbose) {
      print("iter         ELBO")
      print(sprintf("%d:    %f",iter,ELBOs.overall[iter]))
    }

    if (iter >= 50) 
      if (is.finite(ELBOs.overall[iter]) & is.finite(ELBOs.overall[iter - 1]))
        if (abs(ELBOs.overall[iter] -
               ELBOs.overall[iter-1])/abs(ELBOs.overall[iter-1]) < tol.stop)
          break

    # Calculate or update quantities related to model parameters mu,
    # psi2, U_h, u_g and, rho.
    tmp.mu   <- array(0,c(J,K,M))
    tmp.psi2 <- matrix(0,J,K)
    diff.U   <- rep(0,K)
    if (H > 0) {
      for (h in 1:H) {
        tmp.U <- matrix(0,R,R)
        for (j in 1:J) {
          gamma.tmp <- gamma_jk[j,h,]
          Sigma.tmp <- Sigma_jk[[j]][h,,]
          beta.qjh  <- update_q_beta_general(theta_m = gamma.tmp,
                                             theta_V = Sigma.tmp,
                                             c2 = rep(1,R),psi2 = psi2[j],
                                             U = Ulist[[h]])$beta2_m
          tmp.U <- tmp.U + zeta[j,h] * beta.qjh
          for (i in 1:M) {
            k <- which(subgroup == i)
            tmp.mu[j,h,i] <-
              sum(s[k] * exp(bias[j,k] + gamma.tmp[k] + diag(Sigma.tmp)[k]/2))
          }
          eta.qjh <- update_q_eta_general(theta_m = gamma.tmp,
                                          theta_V = Sigma.tmp,
                                          c2 = rep(1,R),psi2 = psi2[j],
                                          U = Ulist[[h]])
          tmp.psi2[j,h] <- sum(eta.qjh)
        }
        
        # Update U_h.
        Uh.new <- tmp.U/sum(zeta[,h])
        Uh.new <- (Uh.new + t(Uh.new))/2
        
        # Avoid too large values in U_h.
        if (max(diag(Uh.new)) > upr_bd)
          Uh.new <- upr_bd * Uh.new / max(diag(Uh.new))
        diff.U[h]  <- max(abs(Uh.new - Ulist[[h]]))
        Ulist[[h]] <- Uh.new
      }
    }
    
    for (g in 1:G) {
        
      # If u is zero vector.
      if (sum(ulist[[g]] != 0) == 0) {
        for (j in 1:J) {
          gamma.tmp <- gamma_jk[j,H+g,]
          Sigma.tmp <- Sigma_jk[[j]][H+g,,]
          for (i in 1:M) {
            k <- which(subgroup == i)
            tmp.mu[j,H+g,i] <- sum(s[k] * exp(bias[j,k] + gamma.tmp[k] +
                                              diag(Sigma.tmp)[k]/2))
          }
          eta.qjg <- gamma.tmp^2 + diag(Sigma.tmp) 
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
      }
      
      # If u is data-driven. 
      else if (ulist.dd[g]) {
        tmp1.u <- 0
        tmp2.u <- rep(0, R)
        for (j in 1:J) {
          gamma.tmp <- gamma_jk[j,H+g,]
          Sigma.tmp <- Sigma_jk[[j]][H+g,,]
          beta.qjg  <- update_q_beta_rank1(theta_m = gamma.tmp,
                                           theta_V = Sigma.tmp,
                                           c2 = rep(1,R),psi2 = psi2[j],
                                           u = ulist[[g]])
          tmp1.u    <- tmp1.u + zeta[j,H+g] * beta.qjg$a2_m/psi2[j]
          tmp2.u    <- tmp2.u + zeta[j,H+g] * beta.qjg$a_theta_m/psi2[j]
          for (i in 1:M) {
            k <- which(subgroup == i)
            tmp.mu[j,H+g,i] <- sum(s[k] * exp(bias[j,k] + gamma.tmp[k] +
                                              diag(Sigma.tmp)[k]/2))
          }
          eta.qjg <- update_q_eta_rank1(theta_m = gamma.tmp,
                                        theta_V = Sigma.tmp,
                                        a2_m = beta.qjg$a2_m,
                                        a_theta_m = beta.qjg$a_theta_m,
                                        u = ulist[[g]])
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
        
        # Update u_g.
        ug.new <- tmp2.u/pmax(tmp1.u,1e-8)
        
        # Avoid too large values in u_g.
        if (max(abs(ug.new)) > sqrt(upr_bd)) 
          ug.new <- sqrt(upr_bd)*ug.new/max(abs(ug.new))
        diff.U[H + g] <- max(abs(ug.new - ulist[[g]]))
        ulist[[g]] <- ug.new
      }
      
      # If u is canonical.
      else {
        ug     <- ulist[[g]]
        ug     <- ug / ug[which.max(abs(ug))]
        tmp1.u <- 0
        tmp2.u <- 0
        for (j in 1:J) {
          gamma.tmp <- gamma_jk[j,H+g,]
          Sigma.tmp <- Sigma_jk[[j]][H+g,,]
          beta.qjg  <- update_q_beta_rank1(theta_m = gamma.tmp,
                                           theta_V = Sigma.tmp,
                                           c2 = rep(1,R),psi2 = psi2[j],
                                           u = ulist[[g]])
          tmp1.u <- tmp1.u + zeta[j,H+g] * sum(ug^2) * beta.qjg$a2_m/psi2[j]
          tmp2.u <- tmp2.u + zeta[j,H+g] * sum(ug * beta.qjg$a_theta_m)/psi2[j]
          for (i in 1:M) {
            k <- which(subgroup == i)
            tmp.mu[j,H+g,i] <- sum(s[k] * exp(bias[j,k] + gamma.tmp[k] +
                                              diag(Sigma.tmp)[k]/2))
          }
          eta.qjg <- update_q_eta_rank1(theta_m = gamma.tmp,
                                        theta_V = Sigma.tmp,
                                        a2_m = beta.qjg$a2_m,
                                        a_theta_m = beta.qjg$a_theta_m,
                                        u = ulist[[g]])
          tmp.psi2[j,H+g] <- sum(eta.qjg)
        }
        
        # Update u_g.
        ug.new      <- pmin(pmax(tmp2.u/tmp1.u,0.01),sqrt(upr_bd)) * ug
        diff.U[H+g] <- max(abs(ug.new - ulist[[g]]))
        ulist[[g]]  <- ug.new
      }
    }
    
    # Update the matrix of means mu.
    mu <- update_mu(data.ed,subgroup,zeta,tmp.mu)

    # Update the dispersion parameter psi2.
    psi2 <- update_psi2(zeta,tmp.psi2,R,minpsi2,maxpsi2)
    
    pi.new  <- update_pi(zeta)
    diff.pi <- pi.new - pi
    pi      <- pi.new
    
    # Update rho and bias if ruv = TRUE.
    if (ruv) {
      rho.new  <- update_rho_all(data.ed,s,mu,F.ed,rho,tmp.ruv,tol = tol.rho,
                                 version = version)
      diff.rho <- rho.new - rho
      rho      <- rho.new
      bias     <- F.ed %*% rho
      bias     <- scale_bias(bias, maxbias)
    }
    
    # Update posterior mean and covariance of theta and local ELBO F_jk.
    for (j in 1:J) {
      out <- update_q_theta_all(data.ed[j,],s,mu[j,],bias[j,],rep(1,R),
                                psi2[j],1,Ulist,ulist,
                                list(gamma = gamma_jk[j,,],
                                     Sigma = Sigma_jk[[j]]),
                                maxiter.q,tol.q)
      gamma_jk[j,,] <- out$gamma
      Sigma_jk[[j]] <- out$Sigma
      A[j,,]        <- out$A
      ELBOs[j,]     <- out$ELBOs
    }
    
    # Update J x K matrix zeta of posterior weights.
    zeta <- update_zeta(ELBOs,pi)
    
    # Update J x R matrix tmp.ruv needed to update rho,
    tmp.ruv <- update_ruv(zeta,A)
  }
  
  # Name the model paramter estimates.
  rownames(mu) <- rownames(data.ed)
  colnames(mu) <- colnames(data.ed)
  names(psi2)  <- rownames(data.ed)
  if (ruv)
    colnames(rho) <- colnames(data.ed)
  names(pi) <- c(names(Ulist),names(ulist))
  
  if (verbose)
    cat("Finish running extreme deconvolution to estimate prior covariance",
        "matrices.\n")
  
  return(list(subset = subset,mu = mu,psi2 = psi2,rho = rho,Ulist = Ulist,
              ulist = ulist,ulist.dd = ulist.dd,pi = pi,zeta = zeta,
              ELBO = ELBOs.overall[1:iter],iff.U = diff.U,
              diff.pi = diff.pi,diff.rho = diff.rho))
}

#' @rdname pois_cov_ed
#'
#' @export
#' 
pois_cov_ed_control_default <- function()
  list(maxiter   = 500,
       maxiter.q = 25,
       maxbias   = 10,
       maxpsi2   = NULL,
       tol.q     = 0.01,
       tol.rho   = 1e-6,
       tol.stop  = 1e-6)
