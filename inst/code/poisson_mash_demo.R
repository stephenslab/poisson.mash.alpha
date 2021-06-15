suppressMessages(suppressWarnings({
  library(Matrix)
  library(glmpca)
  library(scran)
  library(psych)
  library(poilog)
}))

# Read in functions and data
# --------------------------
# load in functions to run poisson mash
# source("code/util.R")
# source("code/set_data.R")
# source("code/pois_mash_ruv_prefit.R")
# source("code/pois_cov_init.R")
# source("code/pois_cov_canonical.R")
# source("code/pois_cov_ed.R")
# source("code/pois_mash_ruv.R")
# source("code/pois_mash_posterior.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# read in a simulated single cell dataset
scdata <- readRDS("../datafiles/scdata.Rds")

# Take a random subset of the rows (genes) and columns (samples).
i <- sample(8358, 2000)
j <- sample(2096, 800)
scdata$Y <- scdata$Y[i,j]
scdata$condition <- scdata$condition[j]

# compute cell-specific size factors using scrna
clusters <- quickCluster(scdata$Y)
si <- calculateSumFactors(scdata$Y, clusters=clusters)

# create a data object for poisson mash analysis
dat <- pois_mash_set_data(scdata$Y, scdata$condition, si)

# Estimate matrix of latent factors causing unwanted variation
# ------------------------------------------------------------
# Run glmpca to estimate matrix of latent factors while adjusting for
# cell-specific size factors and gene-specific, condition-specific
# means.
start_time <- proc.time()
cat("start fitting GLMPCA to estimate matrix of latent factors\n")
design <- model.matrix(~scdata$condition)
fit.glmpca <- glmpca(Y = as.matrix(scdata$Y), X=design[,-1], L=10, fam="nb2",
                     sz=si, ctl=list(verbose=TRUE, maxIter=100, tol=1e-5))
Fuv <- as.matrix(fit.glmpca$loadings)
cat("finish fitting GLMPCA to estimate matrix of latent factors\n")
runtime <- proc.time() - start_time
fit.glmpca$runtime <- runtime
saveRDS(fit.glmpca, file="fit_glmpca.Rds")

# Prefit the model to initialize parameters
# -----------------------------------------
cat("start prefit for parameter initialization\n")
prefit <- pois_mash_ruv_prefit(dat, Fuv, verbose=TRUE)
cat("finish prefit for parameter initialization\n")

# Estimate data-driven prior covariance matrices
# ----------------------------------------------
# Initialize the data-driven prior covariance matrices.
res.pca <- pois_cov_init(dat, ruv=TRUE, Fuv=Fuv, rho=prefit$rho, npc=5)

# Run the ED step.
start_time <- proc.time()
cat("start fitting ED step\n")
fit.ed <- pois_cov_ed(dat, subset=res.pca$subset, Ulist=res.pca$Ulist,
                      ulist=res.pca$ulist, ruv=TRUE, Fuv=Fuv, verbose=TRUE,
                      init=prefit, control=list(maxiter=100))
cat("finish fitting ED step\n")
runtime <- proc.time() - start_time
fit.ed$runtime <- runtime
saveRDS(fit.ed, file = "pois_mash_ruv_ed.Rds")

# Run Poisson mash ruv
# --------------------
# Construct the canonical prior covariance matrices.
ulist.c <- pois_cov_canonical(dat)

# Aombine all the rank-1 prior covariance matrice.
ulist <- c(fit.ed$ulist, ulist.c)

start_time <- proc.time()
cat("start fitting poisson mash with ruv\n")
res <- pois_mash(data=dat, Ulist=fit.ed$Ulist, ulist=ulist, normalizeU=TRUE,
                 gridmult=2.5, ruv=TRUE, Fuv=Fuv, rho=prefit$rho, verbose=TRUE,
                 init=list(mu=prefit$mu, psi2=prefit$psi2)) 
cat("finish fitting poisson mash with ruv\n")
runtime <- proc.time() - start_time
res$runtime <- runtime
saveRDS(res, file = "pois_mash_ruv_fit.Rds")
print(sessionInfo())
