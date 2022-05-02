# Initialize the J x R matrix of means mu, where X is J x R matrix of
# counts, s is R x 1 vector of sequencing depths, subgroup is R x 1
# factor vector with M levels, bias is J x R matrix of bias caused by
# unwanted variation.
initialize_mu <- function (X, s, subgroup, bias = matrix(0,nrow(X),ncol(X))) {
  M  <- length(unique(subgroup))
  mu <- matrix(as.numeric(NA),nrow(X),ncol(X))
  for (i in 1:M) {
    cols <- which(subgroup == i)
    mu[,cols] <- log(rowSums(X[,cols])) - log(exp(bias[,cols]) %*% s[cols])
  }
  return(mu)
}

# Initialize the J x 1 vector of dispersion parameter psi2, where X is
# J x R matrix of counts, s is R x 1 vector of sequencing depths, subgroup
# is R x 1 factor vector with M levels, mu is J x R matrix of means, 
# bias is J x R matrix of bias caused by unwanted variation.
initialize_psi2 <- function (X, s, subgroup = rep(1, length(s)), mu, bias = matrix(0,nrow(X),ncol(X))) {
    
  # Get a rough estimate of log-lambda, which is useful for estimating
  # the range of psi2.
  J         <- nrow(X)
  R         <- ncol(X)
  s.mat     <- rep(1,J) %*% t(s)
  loglambda <- log((X + 0.1)/s.mat)
  M  <- length(unique(subgroup))
  
  # Estimate psi2 for each j
  psi2 <- rep(as.numeric(NA),J)
  for (j in 1:J) {
    psi2_tmp <- rep(0, M)
    for(i in 1:M){
      cols <- which(subgroup == i)
      psi2_tmp[i] <- sd(loglambda[j,cols])^2
    }
    psi2_max       <- pmax(max(psi2_tmp),1)
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
