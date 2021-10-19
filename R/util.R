# Compute the trace of matrix m.
tr <- function (m)
  sum(diag(m),na.rm = TRUE)

# Compute log(det(x)).
logdet <- function (x)
  determinant(x,logarithm = TRUE)$modulus

# scale.cols(A,b) scales each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Scale each row of A so that the entries of each row sum to 1.
normalize.rows <- function (A)
  A / rowSums(A)

# Take inverse of the matrix diag(v1) + v2 %*% t(v3), where
# v1, v2, v3 are vectors.
mat_inv_rank1 <- function (v1, v2, v3)
  diag(1/v1) - (v2/v1) %*% t(v3/v1) / (1 + sum(v2*v3/v1))

# Compute the expected Poisson rates (under the variational
# approximation) from the parameters of the Poisson mash model.
compute_poisson_rates <- function (s, mu, bias, gamma, V) {
  return(s*exp(mu + bias + gamma + V/2))
}

# Scale the J x R matrix of bias to avoid too large values in bias.
scale_bias <- function (bias, maxbias) {
  range.bias <- apply(bias,1,function (x) max(x) - min(x))
  i          <- which(range.bias > maxbias)
  if(length(i) > 0)
    bias[i,] <- maxbias * bias[i,]/range.bias[i]
  return(bias)
}

# Estimate the range of dispersion parameter psi2, where X is J x R
# matrix of counts, s is R x 1 vector of sequencing depths, maxpsi2 is
# a positive scalar specifying the upper bound for psi2, epsilon is a
# small positive number to avoid psi2 being exactly 0.
estimate_psi2_range <- function (X, s, maxpsi2 = NULL, epsilon = 1e-8) {
  J         <- nrow(X)
  s.mat     <- outer(rep(1,J),s)
  loglambda <- log((X + 0.1)/s.mat)
  v         <- apply(loglambda,1,sd)^2
  upr_bd    <- 4*max(v) 
  minpsi2   <- pmax(min(v)/100,epsilon)
  if (is.null(maxpsi2))
    maxpsi2 <- max(v)
  return(list(loglambda = loglambda,
              upr_bd    = upr_bd,
              minpsi2   = minpsi2,
              maxpsi2   = maxpsi2))
}
