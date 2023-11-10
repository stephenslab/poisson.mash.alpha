context("pois_mash")

test_that("Perform basic checks of the model fitting",{
  library(scran)
  library(glmpca)
  
  # Simulate a toy single-cell data for 100 genes and 4 conditions.
  dat <- pois_mash_sim_data(J = 200,R = 4,seed = 1)
  Y   <- dat$Y
  condition <- dat$condition
  
  # Compute cell-specific size factors using scran.
  clusters <- quickCluster(Y)
  s <- calculateSumFactors(Y,clusters = clusters)

  # Create a data object for poisson mash analysis.
  dat <- pois_mash_set_data(Y,condition,s)

  # Run glmpca to estimate matrix of latent factors.
  design <- model.matrix(~condition)
  fit.glmpca <- glmpca(Y = as.matrix(Y),X = design[,-1],L = 2,
                       fam = "nb2",sz = s,
                       ctl = list(maxIter = 100,tol = 1e-5))
  Fuv <- as.matrix(fit.glmpca$loadings)

  # Prefit the model to initialize parameters. The ELBO should
  # increase (or at least not decrease) over two successive
  # iterations.
  capture.output(prefit <- pois_mash_prefit(dat,ruv = TRUE,Fuv))
  expect_gte(min(diff(prefit$ELBO.overall)),0)

  # Initialize the data-driven covariance matrices.
  res.pca <- pois_cov_init(dat,ruv = TRUE,Fuv = Fuv,rho = prefit$rho,
                           npc = 2,cutoff = 2)

  # Run the ED step.
  capture.output(
    fit.ed <- pois_cov_ed(dat,subset = res.pca$subset,Ulist = res.pca$Ulist,
                          ulist = res.pca$ulist,ruv = TRUE,Fuv = Fuv,
                          verbose = TRUE,init = prefit,
                          control = list(maxiter = 100)))

  # The ELBO should increase (or at least not decrease) over two
  # successive iterations.
  expect_gte(min(diff(fit.ed$ELBO)),0)
 
  # Construct the canonical prior covariance matrices.
  ulist.c <- pois_cov_canonical(dat)

  # Combine all the rank-1 prior covariance matrices.
  ulist <- c(fit.ed$ulist,ulist.c)

  # Run Poisson mash RUV.
  capture.output(
    res <- pois_mash(data = dat,Ulist = fit.ed$Ulist,ulist = ulist,
                     normalizeU = TRUE,gridmult = 2.5,ruv = TRUE,Fuv = Fuv,
                     rho = prefit$rho,verbose = TRUE,
                     init = list(mu = prefit$mu,psi2 = prefit$psi2),
                     control = list(maxiter = 100,nc = 2)))

  # The ELBO should increase (or at least not decrease) over two
  # successive iterations.
  expect_gte(min(diff(res$pois.mash.fit$ELBO.overall)),0)
})
