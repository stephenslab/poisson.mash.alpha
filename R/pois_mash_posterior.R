#' @title Compute Posterior Summaries From Poisson Mash Fit
#' 
#' @description Compute posterior summaries of the given contrasts of
#'   the matrix of effects based on the poisson mash fit. Mixture
#'   components with very small posterior weights are ignored in the
#'   posterior calculations.
#'
#'   Note this is an internal function which users might not want to
#'   call directly.
#' 
#' @param data J x R matrix of counts collapsed over conditions, with
#' features as rows and conditions as columns.
#' 
#' @param s Numeric vector of length R adjusting for sequencing depth of
#'   each of R conditions.
#' 
#' @param mu J x R matrix of gene-specific means (R conditions are
#'   assumed to belong to M subgroups, so each row should have at most M
#'   distinct values).
#' 
#' @param psi2 Vector of length J giving the gene-specific dispersion
#'   parameters.
#' 
#' @param bias J x R matrix of bias caused by unwanted
#'   variation. Default is a matrix of all zeros.
#' 
#' @param wlist Numeric vector of length L specifying the scaling
#'   factors for the prior covariance matrices.
#' 
#' @param Ulist List of H full-rank covariance matrices
#' 
#' @param ulist List of G numeric vectors each of which forms a
#'   rank-1 covariance matrix.
#' 
#' @param ulist.epsilon2 Numeric vector of length G added to the
#'   diagonals of each rank-1 prior covariance matrix.
#' 
#' @param zeta J x K matrix of posterior weights, where K = L*(H + G) is
#'   the number of prior mixture components.
#' 
#' @param thresh Posterior weights thresholds. Below this threshold, the
#'    mixture components are ignored in the posterior calculations.
#' 
#' @param C Q x R matrix of contrasts for effects.
#' 
#' @param res.colnames Character vector of length Q giving the names
#' of the contrasts.
#' 
#' @param posterior_samples The number of samples to be drawn from the
#' posterior distribution of each effect.
#' 
#' @param median_deviations Logical scalar indicating whether to calculate
#' posterior summary of deviation of condition-specific effects 
#' relative to the median over all conditions.
#' 
#' @param seed a random number seed to use when sampling from the posteriors.
#' It is used when \code{posterior_samples > 0} 
#' or \code{median_deviations = TRUE}.
#' 
#' @return The return value is a list with the following components:
#' 
#' \item{PosteriorMean}{J x Q matrix of posterior means.}
#' 
#' \item{PosteriorSD}{J x Q matrix of posterior standard deviations.}
#' 
#' \item{ZeroProb}{J x Q matrix in which each entry is the posterior
#'   probability of the true effect being zero.}
#' 
#' \item{NegativeProb}{J x Q matrix in which each entry is the
#'   posterior probability of the true effect being negative.}
#' 
#' \item{lfsr}{J x Q matrix of local false sign rate estimates.}
#' 
#' \item{PosteriorSamples}{J x R x posterior_samples array of samples of effects, 
#' if \code{posterior_samples > 0}.}
#' 
#' \item{beta_median_dev_post}{a list containing posterior summary of 
#' deviation of effects relative to the median, 
#' if \code{median_deviations = TRUE}.}
#'
#' @keywords internal
#'
#' @importFrom stats pnorm rmultinom
#' @importFrom ashr compute_lfsr
#' @importFrom mvtnorm rmvnorm
#' @importFrom abind abind
#' 
#' @export
#' 
pois_mash_posterior <- function (data, s, mu, psi2,
                                 bias = matrix(0,nrow(data),ncol(data)),
                                 wlist, Ulist, ulist,
                                 ulist.epsilon2 = rep(1e-8,length(ulist)),
                                 zeta, thresh = 1/(500*ncol(zeta)),
                                 C = diag(ncol(data)) - 1/ncol(data),
                                 res.colnames=paste0(colnames(data),"-mean"),
                                 posterior_samples = 0, median_deviations = FALSE, seed = 1) {
  data <- as.matrix(data)
  J    <- nrow(data)
  R    <- ncol(data)
  L    <- length(wlist)
  H    <- length(Ulist)
  G    <- length(ulist)
  K    <- ncol(zeta)

  # Matrices to store returned results
  Q <- nrow(C)
  res_post_mean  <- matrix(as.numeric(NA),J,Q)
  res_post_mean2 <- matrix(as.numeric(NA),J,Q)
  res_post_neg   <- matrix(as.numeric(NA),J,Q)
  res_post_zero  <- matrix(as.numeric(NA),J,Q)
  
  # Whether to draw posterior samples of effects
  PosteriorSamples <- NULL
  beta_median_dev_post <- NULL
  
  if(posterior_samples > 0 | median_deviations){
    N <- ifelse(posterior_samples > 0, posterior_samples, 5e3)
    set.seed(seed)
    
    # List to store posterior samples of effects
    if(posterior_samples > 0)
      res_post_samples <- vector("list", J)
    
    # posterior summary of deviations of beta relative to the median
    if(median_deviations){
      median_post_mean <- matrix(as.numeric(NA), J, R)
      median_post_sd <- matrix(as.numeric(NA), J, R)
      median_post_neg <- matrix(as.numeric(NA), J, R)
      median_post_zero <- matrix(as.numeric(NA), J, R)      
    }
  }
  
  # Calculate the posterior summary for each j.
  for (j in 1:J) {
      
    # Matrices to temporarily store results for a given j.
    tmp_post_mean  <- matrix(0,K,Q)
    tmp_post_mean2 <- matrix(0,K,Q)
    tmp_post_neg   <- matrix(0,K,Q)
    tmp_post_zero  <- matrix(1,K,Q)
    
    # Draw posterior samples of effects for a given j.
    if(posterior_samples > 0 | median_deviations){
      samples_j <- rep(list(matrix(0, 0, R)), K)
      z <- rowSums(rmultinom(N, 1, zeta[j,]))
    }
    
    # full-rank prior covariances
    if (H > 0) {
      hl <- 0
      for (h in 1:H) {
        for (l in 1:L) {
          hl <- hl + 1
            
          # If posterior weight exceeds the specified threshold.
          if (zeta[j,hl] > thresh) {
                
            # Calculate posterior mean and covariance of theta.
            out <- update_q_theta_general(data[j,],s,mu[j,],bias[j,],rep(1,R),
                                          psi2[j],wlist[l],Ulist[[h]])
              
            # Calculate posterior mean and covariance of beta.
            out <- update_q_beta_general(out$m,out$V,rep(1,R),psi2[j],
                                         wlist[l],Ulist[[h]])

            # Calculate posterior mean and variance of C %*% beta.
            m.qjhl              <- drop(C %*% out$beta_m)
            sigma2.qjhl         <- pmax(0,diag(C %*% out$beta_V %*% t(C)))
            tmp_post_mean[hl,]  <- m.qjhl
            tmp_post_mean2[hl,] <- m.qjhl^2 + sigma2.qjhl
            tmp_post_neg[hl,]   <- ifelse(sigma2.qjhl == 0,0,
                                          pnorm(0,m.qjhl,sqrt(sigma2.qjhl)))
            tmp_post_zero[hl,]  <- ifelse(sigma2.qjhl == 0,1,0)
            
            # Draw poterior samples of beta
            if(posterior_samples > 0 | median_deviations){
              if(z[hl] > 0)
                samples_j[[hl]] <- rmvnorm(z[hl], mean = as.numeric(out$beta_m), sigma = out$beta_V)
            }
          }
        }
      }
    }
    
    # Rank-1 prior covariances.
    gl <- 0
    for (g in 1:G) {
      ug         <- ulist[[g]]
      utildeg    <- drop(C %*% ug)
      epsilon2.g <- ulist.epsilon2[g]
      for (l in 1:L) {
        gl <- gl + 1
          
        # If posterior weight exceeds the specified threshold.
        if (zeta[j,H*L+gl] > thresh) {
          # If ug is zero vector
          if(sum(utildeg != 0) == 0){
            # Draw posterior samples of beta
            if(posterior_samples > 0 | median_deviations)
              samples_j[[H*L+gl]] <- matrix(0, z[H*L+gl], R)
          }
          
          # If ug is not zero vector and epsilon2.g is not zero
          else if (epsilon2.g > 1e-4) {
                
            # Calculate posterior mean and covariance of theta.
            out <- update_q_theta_rank1(data[j,],s,mu[j,],bias[j,],rep(1,R),
                                        psi2[j] + wlist[l] * epsilon2.g,
                                        wlist[l],ug)
                    
            # Calculate posterior mean and covariance of beta.
            out <- update_q_beta_rank1_robust(out$m,out$V,rep(1,R),psi2[j],
                                              wlist[l],ug,epsilon2.g)
                
            # Calculate posterior mean and variance of C %*% beta.
            m.qjgl      <- drop(C %*% out$beta_m)
            sigma2.qjgl <- pmax(0,diag(C %*% out$beta_V %*% t(C)))
            tmp_post_mean[H*L+gl,]  <- m.qjgl
            tmp_post_mean2[H*L+gl,] <- m.qjgl^2 + sigma2.qjgl
            tmp_post_neg[H*L+gl,]   <-
              ifelse(sigma2.qjgl == 0,0,
                     pnorm(0,mean = m.qjgl,sd = sqrt(sigma2.qjgl),
                           lower.tail = TRUE))
            tmp_post_zero[H*L+gl,]  <- ifelse(sigma2.qjgl == 0,1,0)
            
            # Draw poterior samples of beta
            if(posterior_samples > 0 | median_deviations){
              if(z[H*L+gl] > 0)
                samples_j[[H*L+gl]] <- rmvnorm(z[H*L+gl], mean = as.numeric(out$beta_m), sigma = out$beta_V)
            }
          }
          
          # If ug is not zero vector and epsilon2.g is zero
          else {
              
            # Calculate posterior mean and covariance of theta.
            out <- update_q_theta_rank1(data[j,],s,mu[j,],bias[j,],rep(1,R),
                                        psi2[j],wlist[l],ug)

            # Calculate posterior mean and covariance of beta.
            out <- update_q_beta_rank1(out$m,out$V,rep(1,R),psi2[j],
                                       wlist[l],ug)
            
            # Draw posterior samples of beta
            if(posterior_samples > 0 | median_deviations){
              if(z[H*L+gl] > 0)
                samples_j[[H*L+gl]] <- rnorm(z[H*L+gl], mean=out$a_m, sd=sqrt(out$a_sigma2)) %*% t(ug)
            }

            # Calculate posterior mean and variance of C %*% beta.
            out <- pois_mash_compute_posterior_rank1(out$a_m,out$a_sigma2,
                                                     utildeg)
            tmp_post_mean[H*L+gl,]  <- out$post_mean
            tmp_post_mean2[H*L+gl,] <- out$post_mean2
            tmp_post_neg[H*L+gl,]   <- out$post_neg
            tmp_post_zero[H*L+gl,]  <- out$post_zero
          }
        }
      }
    }
    
    res_post_mean[j,]  <- zeta[j,] %*% tmp_post_mean
    res_post_mean2[j,] <- zeta[j,] %*% tmp_post_mean2
    res_post_neg[j,]   <- zeta[j,] %*% tmp_post_neg
    res_post_zero[j,]  <- zeta[j,] %*% tmp_post_zero
    
    if(posterior_samples > 0 | median_deviations){
      samples_j_full <- do.call(rbind, samples_j)
      
      # Store posterior samples of effects for a given j.
      if(posterior_samples > 0){
        if(nrow(samples_j_full) >= posterior_samples)
          res_post_samples[[j]] <- samples_j_full[sample(nrow(samples_j_full), posterior_samples),]
        else
          res_post_samples[[j]] <- samples_j_full[sample(nrow(samples_j_full), posterior_samples, replace = TRUE),]
      }
      
      # Calculate posterior summary of deviations of beta relative to the median
      if(median_deviations){
        samples_j_full <- samples_j_full - apply(samples_j_full, 1, median) %*% t(rep(1,R))
        median_post_mean[j,] <- apply(samples_j_full, 2, mean)
        median_post_sd[j,] <- apply(samples_j_full, 2, sd)
        median_post_neg[j,] <- apply(samples_j_full, 2, function(x){mean(x<0)})
        median_post_zero[j,] <- apply(samples_j_full, 2, function(x){mean(x==0)})
      }
    }
  }
  
  # Calculate posterior standard deviation.
  res_post_sd <- sqrt(res_post_mean2 - res_post_mean^2)
  
  # Calculate local false sign rate.
  lfsr <- compute_lfsr(res_post_neg,res_post_zero)
  
  rownames(res_post_mean) <- rownames(data)
  colnames(res_post_mean) <- res.colnames
  rownames(res_post_sd)   <- rownames(data)
  colnames(res_post_sd)   <- res.colnames
  rownames(res_post_zero) <- rownames(data)
  colnames(res_post_zero) <- res.colnames
  rownames(res_post_neg)  <- rownames(data)
  colnames(res_post_neg)  <- res.colnames
  rownames(lfsr)          <- rownames(data)
  colnames(lfsr)          <- res.colnames
  
  if(posterior_samples > 0){
    res_post_samples = abind(res_post_samples, along = 0, force.array=TRUE) # dim J x posterior_samples x R
    res_post_samples = array(res_post_samples, dim=c(J,posterior_samples,R))
    res_post_samples = aperm(res_post_samples, c(1,3,2)) # dim J x R x posterior_samples
    dimnames(res_post_samples) <- list(rownames(data), colnames(data), paste0("sample_", (1:posterior_samples)))
    PosteriorSamples = res_post_samples
  }
  
  if(median_deviations){
    rownames(median_post_mean) <- rownames(data)
    colnames(median_post_mean) <- colnames(data)
    rownames(median_post_sd) <- rownames(data)
    colnames(median_post_sd) <- colnames(data)
    rownames(median_post_neg) <- rownames(data)
    colnames(median_post_neg) <- colnames(data)
    rownames(median_post_zero) <- rownames(data)
    colnames(median_post_zero) <- colnames(data)
    beta_median_dev_post <- list(PosteriorMean=median_post_mean, PosteriorSD=median_post_sd, ZeroProb=median_post_zero, NegativeProb=median_post_neg)
  }
  
  return(list(PosteriorMean = res_post_mean,PosteriorSD = res_post_sd,
              ZeroProb = res_post_zero,NegativeProb = res_post_neg,lfsr = lfsr,
              PosteriorSamples = PosteriorSamples, beta_median_dev_post = beta_median_dev_post))
}

#' @importFrom stats pnorm
pois_mash_compute_posterior_rank1 <- function (m, sigma2, u) {
  R <- length(u)
  
  post_mean   <- m*u
  post_mean2  <- (m^2 + sigma2) * (u^2)
  post_neg    <- rep(0,R)
  post_zero   <- rep(0,R)
  post_neg_v  <- ifelse(sigma2 == 0,0,pnorm(0,m,sqrt(sigma2)))
  post_zero_v <- ifelse(sigma2 == 0,1,0)
  post_pos_v  <- pmax(0,1 - post_neg_v - post_zero_v)
  
  for (r in 1:R) {
    ur <- u[r]
    if (ur > 0) {
      post_neg[r]  <- post_neg_v
      post_zero[r] <- post_zero_v
    }
    else if (ur < 0) {
      post_neg[r]  <- post_pos_v
      post_zero[r] <- post_zero_v
    }
    else  
      post_zero[r] <- 1
  }
  
  return(list(post_mean  = post_mean,
              post_mean2 = post_mean2,
              post_neg   = post_neg,
              post_zero  = post_zero))
}
