# Compute the trace of matrix m.
tr <- function (m)
  sum(diag(m),na.rm = TRUE)


# Take inverse of the matrix diag(v1) + v2 %*% t(v3), where
# v1, v2, v3 are vectors.
mat_inv_rank1 <- function (v1, v2, v3)
  diag(1/v1) - (v2/v1) %*% t(v3/v1) / (1 + sum(v2*v3/v1))


# Scale the J x R matrix of bias to avoid too large values in bias.
scale_bias <- function (bias, maxbias) {
  range.bias <- apply(bias,1,function(x) max(x) - min(x))
  idx.bias <- which(range.bias > maxbias)
  if(length(idx.bias) > 0)
    bias[idx.bias,] <- maxbias*bias[idx.bias,]/range.bias[idx.bias]
  return(bias)
}


# Initialize the J x R matrix of means mu, where X is J x R matrix of counts, s is R x 1 vector of sequencing depths, 
# subgroup is R x 1 factor vector with M levels, bias is J x R matrix of bias caused by unwanted variation.
initialize_mu <- function(X, s, subgroup, bias=matrix(0, nrow(X), ncol(X))){
  M <- length(unique(subgroup))
  mu <- matrix(as.numeric(NA), nrow(X), ncol(X))
  for (i in 1:M)
    mu[,subgroup == i] <- log(rowSums(X[,subgroup == i])) - log(exp(bias[,subgroup == i]) %*% s[subgroup == i])
  return(mu)
}


# Estimate the range of dispersion parameter psi2, where X is J x R matrix of counts, s is R x 1 vector of sequencing depths, 
# maxpsi2 is a positive scalar specifying the upper bound for psi2, epsilon is a small positive number to avoid psi2 being exactly 0.
estimate_psi2_range <- function (X, s, maxpsi2=NULL, epsilon=1e-8){
  s.mat <- rep(1, nrow(X)) %*% t(s)
  loglambda <- log((X + 0.1)/s.mat)
  upr_bd  <- 4*max(apply(loglambda,1,sd)^2) 
  minpsi2 <- pmax(min(apply(loglambda,1,sd)^2)/100, epsilon)
  if (is.null(maxpsi2))
    maxpsi2 <- max(apply(loglambda,1,sd)^2)
  return(list(loglambda=loglambda, upr_bd=upr_bd, minpsi2=minpsi2, maxpsi2=maxpsi2))
}


# Initialize the J x 1 vector of dispersion parameter psi2, where X is
# J x R matrix of counts, s is R x 1 vector of sequencing depths, mu
# is J x R matrix of means, bias is J x R matrix of bias caused by
# unwanted variation.
initialize_psi2 <- function (X, s, mu, bias = matrix(0,nrow(X),ncol(X))) {
    
  # Get a rough estimate of log-lambda, which is useful for estimating
  # the range of psi2.
  J         <- nrow(X)
  R         <- ncol(X)
  s.mat     <- rep(1,J) %*% t(s)
  loglambda <- log((X + 0.1)/s.mat)
  
  # Estimate psi2 for each j
  psi2 <- rep(as.numeric(NA),J)
  for (j in 1:J) {
    psi2_max       <- pmax(sd(loglambda[j,])^2,1)
    log2_psi2_grid <- seq(log2(1e-4),log2(psi2_max),length.out = 25)
    psi2_grid      <- 2^log2_psi2_grid
    logdens        <- rep(0,length(psi2_grid))
    for(l in 1:length(psi2_grid))
      for(r in 1:R)
        logdens[l] <- logdens[l] + log(dpoilog(X[j,r],mu[j,r] + bias[j,r] +
                                               log(s[r]),sqrt(psi2_grid[l])))
    psi2[j] <- psi2_grid[which.max(logdens)]
  }
  
  return(psi2)
}


# Calculate the constant term in overall ELBO, where X is J x R matrix of counts, s is R x 1 vector of sequencing depths.
compute_elbo_const <- function(X, s){
  sum(X %*% log(s)) - sum(lgamma(1 + X))
}


# Calculate the overall ELBO, where ELBOs is J x K matrix of local ELBO, 
# pi is K x 1 vector of prior weights, zeta is J x K matrix of posterior weights, const is a scalar. 
compute_overall_elbo <- function(ELBOs, pi, zeta, const){
  sum(zeta*(log(rep(1,nrow(zeta)) %*% t(pi)) + ELBOs - log(zeta))) + const
}
