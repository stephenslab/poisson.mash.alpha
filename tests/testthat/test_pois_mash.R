context("pois_mash")

test_that("Describe test here",{
  library(scran)
  library(glmpca)
  
  # Simulate a toy single-cell data for 100 genes and 4 conditions.
  dat <- pois_mash_sim_data(J = 100,R = 4,seed = 1)
  Y   <- dat$Y
  condition <- dat$condition
  
  # Compute cell-specific size factors using scran.
  clusters <- quickCluster(Y)
  s <- calculateSumFactors(Y,clusters = clusters)

  # Create a data object for poisson mash analysis.
  dat <- pois_mash_set_data(Y,condition,s)

  # Run glmpca to estimate matrix of latent factors.
  design <- model.matrix(~condition)
  fit.glmpca <- glmpca(Y = as.matrix(Y),X = design[,-1],L = 10,
                       fam = "nb2",sz = s,
                       ctl = list(maxIter = 100,tol = 1e-5))
  Fuv <- as.matrix(fit.glmpca$loadings)

  # Prefit the model to initialize parameters. The ELBO should
  # increase (or at least not decrease) over two successive
  # iterations.
  prefit <- pois_mash_ruv_prefit(dat,Fuv)
  expect_gte(min(diff(prefit$ELBO)),0)
})
