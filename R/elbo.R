# Calculate the constant term in overall ELBO, where X is J x R matrix
# of counts, s is R x 1 vector of sequencing depths.
compute_elbo_const <- function (X, s)
  sum(X %*% log(s)) - sum(lgamma(1 + X))

# Calculate the overall ELBO, where ELBOs is J x K matrix of local
# ELBO, pi is K x 1 vector of prior weights, zeta is J x K matrix of
# posterior weights, const is a scalar.
compute_overall_elbo <- function (ELBOs, pi, zeta, const)
  sum(zeta * (log(rep(1,nrow(zeta)) %*% t(pi)) + ELBOs - log(zeta))) + const

# Calculate the ELBO for each gene, where ELBOs is J x K matrix of local
# ELBO, pi is K x 1 vector of prior weights, zeta is J x K matrix of
# posterior weights, X is J x R matrix of counts, 
# s is R x 1 vector of sequencing depths.
compute_genewise_elbo <- function (ELBOs, pi, zeta, X, s)
  rowSums(zeta * (log(rep(1,nrow(zeta)) %*% t(pi)) + ELBOs - log(zeta))) + X %*% log(s) - rowSums(lgamma(1 + X))