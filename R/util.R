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

