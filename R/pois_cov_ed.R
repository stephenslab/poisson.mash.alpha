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
#'   each ED iteration; and \dQuote{tol.rho}{tolerance for
#'   assessing convergence of effects corresponding to unwanted
#'   variation.}
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
#' @importFrom stats sd
#' @importFrom poilog dpoilog
#' 
#' @export
#' 
pois_cov_ed <- function (data, subset = NULL, Ulist, ulist, ulist.dd = NULL,
                         ruv = FALSE, Fuv = NULL, verbose = FALSE,
                         init = list(NULL),  
                         control = list(maxiter = 500, maxiter.q = 25,
                                        maxpsi2 = NULL, maxbias = 10,
                                        tol.stop = 1e-6, tol.q = 0.01,
                                        tol.rho=1e-6)) {
  X        <- data$X
  s        <- data$s
  subgroup <- data$subgroup
  if (is.null(ulist))
    stop("ulist cannot be empty!")
  if (is.null(subset))
    subset <- 1:nrow(X)
  data.ed   <- as.matrix(X[subset,])
  J         <- nrow(data.ed)
  R         <- ncol(data.ed)
  M         <- length(unique(subgroup))
  subgroup  <- as.numeric(as.factor(subgroup))
  maxiter   <- control$maxiter
  maxiter.q <- control$maxiter.q
  maxpsi2   <- control$maxpsi2
  maxbias   <- control$maxbias
  tol.stop  <- control$tol.stop
  tol.q     <- control$tol.q
  tol.rho   <- control$tol.rho
  
  if (is.null(maxiter))
    maxiter <- 500
  if (is.null(maxiter.q))
    maxiter.q <- 25
  if (is.null(maxbias))
    maxbias <- 10
  if (is.null(tol.stop))
    tol.stop <- 1e-6
  if (is.null(tol.q))
    tol.q <- 0.01
  if (is.null(tol.rho))
    tol.rho <- 1e-6
  
  H <- length(Ulist)
  G <- length(ulist)
  K <- H + G
  
  if (is.null(ulist.dd)) {
    ulist.dd <- rep(TRUE,G)
    
    # Set to FALSE if zero vector.
    for (g in 1:G) 
      if (sum(ulist[[g]]!=0)==0)
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
    if (is.null(Fuv))
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
  # CAN THIS BE A FUNCTION? e.g., estimate_psi2_range.
  psi2_range <- estimate_psi2_range(X=data.ed, s=s, maxpsi2=maxpsi2)
  minpsi2 <- psi2_range$minpsi2
  maxpsi2 <- psi2_range$maxpsi2
  upr_bd <- psi2_range$upr_bd
  
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
  gamma_jk <- array(as.numeric(NA), c(J,K,R))
  Sigma_jk <- list(NULL)
  for (j in 1:J) {
    Sigma_jk[[j]] <- array(as.numeric(NA), c(K,R,R))
  }
  
  # Create a J x K x R array to store the quantities related to
  # q_jk, s.t. A[j,k,r] = gamma_jkr + 0.5*Sigma_jk,rr
  A <- array(as.numeric(NA),c(J,K,R))
  
  # Matrix to store "local" ELBOs F_jk.
  ELBOs <- matrix(0,J,K)
  
  # Update posterior mean and covariance of theta and local ELBO.
  for (j in 1:J) {
    # CAN THIS BE A FUNCTION? e.g., update_q_theta_all.
    theta.q.all <- update_q_theta_all(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], 
                                      wlist=1, Ulist=Ulist, ulist=ulist, maxiter.q=maxiter.q, tol.q=tol.q)
    gamma_jk[j,,] <- theta.q.all$gamma
    Sigma_jk[[j]] <- theta.q.all$Sigma
    A[j,,] <- theta.q.all$A
    ELBOs[j,] <- theta.q.all$ELBOs
  }
  
  # Update J x K matrix zeta of posterior weights.
  # CAN THIS BE A FUNCTION? e.g., update_zeta.
  zeta <- update_zeta(ELBOs=ELBOs, pi=pi)
  
  # Update J x R matrix tmp.ruv needed to update rho,
  # s.t. tmp.ruv[j,r] = sum_k zeta[j,k] * exp(A[j,k,r]).
  tmp.ruv <- matrix(as.numeric(NA),J,R)
  for (r in 1:R)
    tmp.ruv[,r] <- rowSums(zeta*exp(A[,,r]))

  const <- compute_elbo_const(data.ed,s)
  
  # Overall ELBO after updating all parameters at each iteration.
  ELBOs.overall <- c()
  
  if (verbose)
    cat("Start running extreme deconvolution to estimate prior covariance",
        "matrices.\n")
  
  for (iter in 1:maxiter) {
      
    # Calculate overall ELBO at the current iteration.
    ELBO.overall  <- compute_overall_elbo(ELBOs,pi,zeta,const)
    ELBOs.overall <- c(ELBOs.overall,ELBO.overall)
    
    if (verbose) {
      print("iter         ELBO")
      print(sprintf("%d:    %f",iter,ELBO.overall))
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
          for (i in 1:M)
            tmp.mu[j,h,i] <-
              sum(s[subgroup == i] * exp(bias[j,subgroup == i] +
              gamma.tmp[subgroup == i] + diag(Sigma.tmp)[subgroup == i]/2))
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
            tmp.mu[j,H+g,i] <- sum(s[subgroup == i] *
                                   exp(bias[j,subgroup == i] +
                                       gamma.tmp[subgroup == i] +
                                       diag(Sigma.tmp)[subgroup == i]/2))
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
          for (i in 1:M)
            tmp.mu[j,H+g,i] <- sum(s[subgroup == i] *
                                   exp(bias[j,subgroup == i] +
                                       gamma.tmp[subgroup == i] +
                                       diag(Sigma.tmp)[subgroup == i]/2))
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
          for (i in 1:M)
            tmp.mu[j,H+g,i] <- sum(s[subgroup == i] *
                                   exp(bias[j,subgroup == i] +
                                       gamma.tmp[subgroup == i] +
                                       diag(Sigma.tmp)[subgroup == i]/2))
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
    
    # CAN THIS BE A FUNCTION? e.g., update_mu.
    mu <- update_mu(X=data.ed, subgroup=subgroup, zeta=zeta, tmp.mu=tmp.mu)

    # Update the dispersion parameter psi2.
    psi2 <- update_psi2(zeta,tmp.psi2,R,minpsi2,maxpsi2)
    
    pi.new  <- update_pi(zeta)
    diff.pi <- pi.new - pi
    pi      <- pi.new
    
    # Update rho and bias if ruv = TRUE.
    if (ruv) {
      # CAN THIS BE A FUNCTION? e.g., update_rho_all.
      rho.new <- update_rho_all(X=data.ed, s=s, mu=mu, Fuv=F.ed, rho=rho, tmp.ruv=tmp.ruv, tol.rho=tol.rho)
      diff.rho <- rho.new - rho
      rho      <- rho.new
      bias     <- F.ed %*% rho
      bias     <- scale_bias(bias, maxbias)
    }
    
    # Update posterior mean and covariance of theta and local ELBO F_jk.
    for (j in 1:J) {
      # CAN THIS BE A FUNCTION? e.g., update_q_theta_all.
      theta.q.all <- update_q_theta_all(x=data.ed[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], wlist=1, Ulist=Ulist, ulist=ulist,
                                        init=list(gamma=gamma_jk[j,,], Sigma=Sigma_jk[[j]]), maxiter.q=maxiter.q, tol.q=tol.q)
      gamma_jk[j,,] <- theta.q.all$gamma
      Sigma_jk[[j]] <- theta.q.all$Sigma
      A[j,,] <- theta.q.all$A
      ELBOs[j,] <- theta.q.all$ELBOs
    }
    
    # Update J x K matrix zeta of posterior weights.
    # CAN THIS BE A FUNCTION? e.g., update_zeta.
    zeta <- update_zeta(ELBOs=ELBOs, pi=pi)
    
    # Update J x R matrix tmp.ruv needed to update rho,
    # s.t. tmp.ruv[j,r] = sum_k zeta[j,k] * exp(A[j,k,r])
    for (r in 1:R)
      tmp.ruv[,r] <- rowSums(zeta * exp(A[,,r]))
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
              ELBO = ELBOs.overall,iff.U = diff.U,diff.pi = diff.pi,
              diff.rho = diff.rho))
}
