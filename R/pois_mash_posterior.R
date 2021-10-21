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
#' @keywords internal
#'
#' @importFrom stats pnorm
#' @importFrom ashr compute_lfsr
#' 
#' @export
#' 
pois_mash_posterior <- function (data, s, mu, psi2,
                                 bias = matrix(0,nrow(data),ncol(data)),
                                 wlist, Ulist, ulist,
                                 ulist.epsilon2 = rep(1e-8,length(ulist)),
                                 zeta, thresh = 1/(500*ncol(zeta)),
                                 C = diag(ncol(data)) - 1/ncol(data),
                                 res.colnames=paste0(colnames(data),"-mean")) {
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
  
  # Calculate the posterior summary for each j.
  for (j in 1:J) {
      
    # Matrices to temporarily store results for a given j.
    tmp_post_mean  <- matrix(0,K,Q)
    tmp_post_mean2 <- matrix(0,K,Q)
    tmp_post_neg   <- matrix(0,K,Q)
    tmp_post_zero  <- matrix(1,K,Q)
    
    # full-rank prior covariances
    if (H > 0) {
      hl <- 0
      for (h in 1:H) {
        for (l in 1:L) {
          hl <- hl + 1
            
          # If posterior weight exceeds the specified threshold.
          if (zeta[j,hl] > thresh) {
                
            # Calculate posterior mean and covariance of theta.
            theta.qjhl <- update_q_theta_general(x = data[j,],s = s,
                                                 mu = mu[j,],bias = bias[j,],
                                                 c2 = rep(1,R),psi2 = psi2[j],
                                                 w = wlist[l],U = Ulist[[h]])
              
            # Calculate posterior mean and covariance of beta.
            beta.qjhl <- update_q_beta_general(theta_m = theta.qjhl$m,
                                               theta_V = theta.qjhl$V,
                                               c2 = rep(1,R),psi2 = psi2[j],
                                               w = wlist[l],U = Ulist[[h]])
              
            # Calculate posterior mean and variance of C %*% beta.
            m.qjhl <- as.numeric(C %*% beta.qjhl$beta_m)
            sigma2.qjhl <- pmax(0,diag(C %*% beta.qjhl$beta_V %*% t(C)))
            tmp_post_mean[hl,]  <- m.qjhl
            tmp_post_mean2[hl,] <- m.qjhl^2 + sigma2.qjhl
            tmp_post_neg[hl,]   <-
              ifelse(sigma2.qjhl == 0,0,
                     pnorm(0,mean = m.qjhl,sd = sqrt(sigma2.qjhl),
                           lower.tail=TRUE))
            tmp_post_zero[hl,]  <- ifelse(sigma2.qjhl == 0,1,0)            
          }
        }
      }
    }
    
    # Rank-1 prior covariances.
    gl <- 0
    for (g in 1:G) {
      ug         <- ulist[[g]]
      utildeg    <- as.numeric(C %*% ug)
      epsilon2.g <- ulist.epsilon2[g]
      for (l in 1:L) {
        gl <- gl + 1
          
        # Check if posterior weight exceeds the specified threshold
        # and ug is not zero vector.
        if (zeta[j,H*L+gl] > thresh & sum(utildeg != 0) > 0) {
          if (epsilon2.g > 1e-4) {
                
            # Calculate posterior mean and covariance of theta.
            theta.qjgl <- update_q_theta_rank1(x = data[j,],s = s,mu = mu[j,],
                            bias = bias[j,],c2 = rep(1,R),
                            psi2 = psi2[j] + wlist[l] * epsilon2.g,
                            w = wlist[l],u = ug)
                
            # Calculate posterior mean and covariance of beta.
            beta.qjgl <- update_q_beta_rank1_robust(theta_m = theta.qjgl$m,
                           theta_V = theta.qjgl$V,c2 = rep(1,R),psi2 = psi2[j],
                           w = wlist[l],u = ug,epsilon2 = epsilon2.g)
                
            # Calculate posterior mean and variance of C %*% beta.
            m.qjgl      <- as.numeric(C %*% beta.qjgl$beta_m)
            sigma2.qjgl <- pmax(0,diag(C %*% beta.qjgl$beta_V %*% t(C)))
            tmp_post_mean[H*L+gl,]  <- m.qjgl
            tmp_post_mean2[H*L+gl,] <- m.qjgl^2 + sigma2.qjgl
            tmp_post_neg[H*L+gl,]   <-
              ifelse(sigma2.qjgl == 0,0,
                     pnorm(0,mean = m.qjgl,sd = sqrt(sigma2.qjgl),
                           lower.tail = TRUE))
            tmp_post_zero[H*L+gl,]  <- ifelse(sigma2.qjgl == 0,1,0)
          }
          else {
              
            # Calculate posterior mean and covariance of theta.
            theta.qjgl <- update_q_theta_rank1(x = data[j,],s = s,mu = mu[j,],
                                               bias = bias[j,],c2 = rep(1,R),
                                               psi2 = psi2[j],w = wlist[l],
                                               u = ug)
            
            # Calculate posterior mean and covariance of beta.
            beta.qjgl <- update_q_beta_rank1(theta_m = theta.qjgl$m,
                                             theta_V = theta.qjgl$V,
                                             c2 = rep(1,R),psi2 = psi2[j],
                                             w = wlist[l],u = ug)
            
            # Calculate posterior mean and variance of C %*% beta.
            post.qjgl <- pois_mash_compute_posterior_rank1(m = beta.qjgl$a_m,
                           sigma2 = beta.qjgl$a_sigma2,u = utildeg)
            tmp_post_mean[H*L+gl,]  <- post.qjgl$post_mean
            tmp_post_mean2[H*L+gl,] <- post.qjgl$post_mean2
            tmp_post_neg[H*L+gl,]   <- post.qjgl$post_neg
            tmp_post_zero[H*L+gl,]  <- post.qjgl$post_zero
          }
        }
      }
    }
    
    res_post_mean[j,] <- zeta[j,] %*% tmp_post_mean
    res_post_mean2[j,] <- zeta[j,] %*% tmp_post_mean2
    res_post_neg[j,] <- zeta[j,] %*% tmp_post_neg
    res_post_zero[j,] <- zeta[j,] %*% tmp_post_zero
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
  
  return(list(PosteriorMean = res_post_mean,PosteriorSD = res_post_sd,
              ZeroProb = res_post_zero,NegativeProb = res_post_neg,
              lfsr = lfsr))
}

#' @importFrom stats pnorm
pois_mash_compute_posterior_rank1 <- function (m, sigma2, u) {
  R <- length(u)
  
  post_mean   <- m*u
  post_mean2  <- (m^2 + sigma2) * (u^2)
  post_neg    <- rep(0,R)
  post_zero   <- rep(0,R)
  post_neg_v  <- ifelse(sigma2 == 0,0,
                        pnorm(0,mean = m,sd = sqrt(sigma2),lower.tail = TRUE))
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
  
  return(list(post_mean = post_mean,post_mean2 = post_mean2,
              post_neg = post_neg,post_zero = post_zero))
}
