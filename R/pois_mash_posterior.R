#' @title Compute posterior summaries of the given contrasts of the matrix of effects based on the poisson mash fit
#' 
#' @description This is an internal function which users might not want to call directly. 
#' Mixtures components with very small posterior weights are ignored in the calculation of effect posterior summaries. 
#' 
#' @param data J by R matrix of counts collapsed over conditions, with features as rows and conditions as columns
#' 
#' @param s R by 1 numeric vector s adjusting for sequencing depth of each of R conditions
#' 
#' @param mu J by R matrix of gene-specific means (R conditions are assumed to belong to M subgroups, so each row should have at most M distinct values)
#' 
#' @param psi2 J by 1 vector of gene-specific dispersion parameters
#' 
#' @param bias J by R matrix of bias caused by unwanted variation. Default to matrix of all zeros.
#' 
#' @param wlist L by 1 numeric vector of scaling factors for the prior covariance matrices
#' 
#' @param Ulist A list of H full-rank covariance matrices
#' 
#' @param ulist A list of G numeric vectors each of which forming a rank-1 covariance matrix
#' 
#' @param ulist.epsilon2 G by 1 numeric vector that adds a small positive epsilon2 to the diagonals of each rank-1 prior covariance matrix
#' 
#' @param zeta J by K matrix of posterior weights, where K=L*(H+G) is the number of mixture components in the prior
#' 
#' @param thresh the threshold for posterior weights below which the corresponding mixture components are ignored in posterior summary calculation
#' 
#' @param C Q by R matrix of contrasts for effects
#' 
#' @param res.colnames Q by 1 character vector giving the names of the contrasts
#' 
#' @return A list with the following components:
#' \item{PosteriorMean}{J x Q matrix of posterior means.}
#' \item{PosteriorSD}{J x Q matrix of posterior standard deviations.}
#' \item{ZeroProb}{J x Q matrix of posterior probability of being zero.}
#' \item{NegativeProb}{J x Q matrix of posterior probability of being negative.}
#' \item{lfsr}{J x Q matrix of local false sign rate estimates.}
#'
#' @keywords internal
#' 
#' @export
#' 
pois_mash_posterior <- function(data, s, mu, psi2, bias=NULL, wlist, Ulist, ulist, ulist.epsilon2=NULL, zeta, thresh=NULL,
                                C=NULL, res.colnames=NULL){
  data <- as.matrix(data)
  J <- nrow(data)
  R <- ncol(data)
  L <- length(wlist)
  H <- length(Ulist)
  G <- length(ulist)
  K <- ncol(zeta)

  # calculate bias caused by unwanted variation  
  if(is.null(bias)){
    bias <- matrix(0, nrow=J, ncol=R)
  }
  
  # specify ulist.epsilon2 if not provided
  if(is.null(ulist.epsilon2)){
    ulist.epsilon2 <- rep(1e-8, G)
  }
  
  # specify the threshold for posterior weights to ignore mixture components
  if(is.null(thresh)){
    thresh <- 1/K/500
  }
  
  # specify the default contrast matrix
  if(is.null(C)){
    C <- diag(R) - 1/R
  }

  # specify the colnames of the contrasts
  if(is.null(res.colnames)){
    res.colnames <- paste0(colnames(data), "-mean")
  }
  
  # matrices to store returned results
  Q <- nrow(C)
  res_post_mean <- matrix(NA, nrow=J, ncol=Q)
  res_post_mean2 <- matrix(NA, nrow=J, ncol=Q)
  res_post_neg <- matrix(NA, nrow=J, ncol=Q)
  res_post_zero <- matrix(NA, nrow=J, ncol=Q)
  
  
  # calculate the posterior summary for each j
  for(j in 1:J){
    # matrices to temporarily store results for a given j
    tmp_post_mean <- matrix(0, nrow=K, ncol=Q)
    tmp_post_mean2 <- matrix(0, nrow=K, ncol=Q)
    tmp_post_neg <- matrix(0, nrow=K, ncol=Q)
    tmp_post_zero <- matrix(1, nrow=K, ncol=Q)
    
    # full-rank prior covariances
    if(H > 0){
      hl <- 0
      for(h in 1:H){
        for(l in 1:L){
          hl <- hl + 1
          # if posterior weight exceeds the specified threshold
          if(zeta[j, hl] > thresh){
            # calculate posterior mean and covariance of theta
            theta.qjhl <- update_q_theta_general(x=data[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], w=wlist[l], U=Ulist[[h]])
            # calculate posterior mean and covariance of beta
            beta.qjhl <- update_q_beta_general(theta_m=theta.qjhl$m, theta_V=theta.qjhl$V, c2=rep(1,R), psi2=psi2[j], w=wlist[l], U=Ulist[[h]])
            # calculate posterior mean and variance of C%*%beta
            m.qjhl <- as.numeric(C %*% beta.qjhl$beta_m)
            sigma2.qjhl <- pmax(0, diag(C %*% beta.qjhl$beta_V %*% t(C)))
            tmp_post_mean[hl,] <- m.qjhl
            tmp_post_mean2[hl,] <- m.qjhl^2 + sigma2.qjhl
            tmp_post_neg[hl,] <- ifelse(sigma2.qjhl==0, 0, pnorm(0, mean=m.qjhl, sqrt(sigma2.qjhl), lower.tail=TRUE))
            tmp_post_zero[hl,] <- ifelse(sigma2.qjhl==0, 1, 0)            
          }
        }
      }
    }
    
    # rank-1 prior covariances  
    gl <- 0
    for(g in 1:G){
      ug <- ulist[[g]]
      utildeg <- as.numeric(C %*% ug)
      epsilon2.g <- ulist.epsilon2[g]
      for(l in 1:L){
        gl <- gl + 1
        # if posterior weight exceeds the specified threshold and ug is not zero vector
        if(zeta[j, H*L+gl] > thresh & sum(utildeg!=0) > 0){
          if(epsilon2.g > 1e-4){
            # calculate posterior mean and covariance of theta
            theta.qjgl <- update_q_theta_rank1(x=data[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j]+wlist[l]*epsilon2.g, w=wlist[l], u=ug)
            # calculate posterior mean and covariance of beta
            beta.qjgl <- update_q_beta_rank1_robust(theta_m=theta.qjgl$m, theta_V=theta.qjgl$V, c2=rep(1,R), psi2=psi2[j], w=wlist[l], u=ug, epsilon2=epsilon2.g)
            # calculate posterior mean and variance of C%*%beta
            m.qjgl <- as.numeric(C %*% beta.qjgl$beta_m)
            sigma2.qjgl <- pmax(0, diag(C %*% beta.qjgl$beta_V %*% t(C)))
            tmp_post_mean[H*L+gl,] <- m.qjgl
            tmp_post_mean2[H*L+gl,] <- m.qjgl^2 + sigma2.qjgl
            tmp_post_neg[H*L+gl,] <- ifelse(sigma2.qjgl==0, 0, pnorm(0, mean=m.qjgl, sqrt(sigma2.qjgl), lower.tail=TRUE))
            tmp_post_zero[H*L+gl,] <- ifelse(sigma2.qjgl==0, 1, 0)
          }
          else{
            # calculate posterior mean and covariance of theta
            theta.qjgl <- update_q_theta_rank1(x=data[j,], s=s, mu=mu[j,], bias=bias[j,], c2=rep(1,R), psi2=psi2[j], w=wlist[l], u=ug)     
            # calculate posterior mean and covariance of beta
            beta.qjgl <- update_q_beta_rank1(theta_m=theta.qjgl$m, theta_V=theta.qjgl$V, c2=rep(1,R), psi2=psi2[j], w=wlist[l], u=ug)
            # calculate posterior mean and variance of C%*%beta
            post.qjgl <- pois_mash_compute_posterior_rank1(m=beta.qjgl$a_m, sigma2=beta.qjgl$a_sigma2, u=utildeg)
            tmp_post_mean[H*L+gl,] <- post.qjgl$post_mean
            tmp_post_mean2[H*L+gl,] <- post.qjgl$post_mean2
            tmp_post_neg[H*L+gl,] <- post.qjgl$post_neg
            tmp_post_zero[H*L+gl,] <- post.qjgl$post_zero
          }
        }
      }
    }
    
    res_post_mean[j,] <- zeta[j,] %*% tmp_post_mean
    res_post_mean2[j,] <- zeta[j,] %*% tmp_post_mean2
    res_post_neg[j,] <- zeta[j,] %*% tmp_post_neg
    res_post_zero[j,] <- zeta[j,] %*% tmp_post_zero
  }
  
  # calculate posterior standard deviation
  res_post_sd <- sqrt(res_post_mean2 - res_post_mean^2)
  
  # calculate local false sign rate
  lfsr <- ashr::compute_lfsr(res_post_neg, res_post_zero)
  
  rownames(res_post_mean) <- rownames(data)
  colnames(res_post_mean) <- res.colnames
  rownames(res_post_sd) <- rownames(data)
  colnames(res_post_sd) <- res.colnames
  rownames(res_post_zero) <- rownames(data)
  colnames(res_post_zero) <- res.colnames
  rownames(res_post_neg) <- rownames(data)
  colnames(res_post_neg) <- res.colnames
  rownames(lfsr) <- rownames(data)
  colnames(lfsr) <- res.colnames
  
  
  return(list(PosteriorMean=res_post_mean, PosteriorSD=res_post_sd, ZeroProb=res_post_zero, NegativeProb=res_post_neg, lfsr=lfsr))
}

pois_mash_compute_posterior_rank1 <- function(m, sigma2, u){
  R <- length(u)
  
  post_mean <- m*u
  post_mean2 <- (m^2 + sigma2)*(u^2)
  post_neg <- rep(0, R)
  post_zero <- rep(0, R)
  
  post_neg_v <- ifelse(sigma2==0, 0, pnorm(0, mean=m, sqrt(sigma2), lower.tail = TRUE))
  post_zero_v <- ifelse(sigma2==0, 1, 0)
  post_pos_v <- pmax(0, 1-post_neg_v-post_zero_v)
  
  for(r in 1:R){
    ur <- u[r]
    if(ur > 0){
      post_neg[r] <- post_neg_v
      post_zero[r] <- post_zero_v
    }
    else if(ur < 0){
      post_neg[r] <- post_pos_v
      post_zero[r] <- post_zero_v
    }
    else{
      post_zero[r] <- 1
    }
  }
  
  return(list(post_mean=post_mean, post_mean2=post_mean2, post_neg=post_neg, post_zero=post_zero))
}
