library(Matrix)
library(scran)
library(glmpca)
library(poisson.mash.alpha)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Take a random subset of the rows (genes) and columns (samples).
scdata <- readRDS("../datafiles/scdata.Rds")
i <- sample(8358,2000)
j <- sample(2096,800)
scdata$Y <- scdata$Y[i,j]
scdata$condition <- scdata$condition[j]

# Compute cell-specific size factors using scran.
clusters <- quickCluster(scdata$Y)
si <- calculateSumFactors(scdata$Y, clusters=clusters)

# Create a data object for poisson mash analysis.
dat <- pois_mash_set_data(scdata$Y,scdata$condition,si)

# Estimate matrix of latent factors causing unwanted variation
# ------------------------------------------------------------
# Run glmpca to estimate matrix of latent factors while adjusting for
# cell-specific size factors and gene-specific, condition-specific
# means.
start_time <- proc.time()
cat("start fitting GLMPCA to estimate matrix of latent factors\n")
design <- model.matrix(~scdata$condition)
fit.glmpca <- glmpca(Y = as.matrix(scdata$Y),X = design[,-1],L = 10,
                     fam = "nb2",sz = si,
                     ctl = list(verbose = TRUE,maxIter = 100,tol = 1e-5))
Fuv <- as.matrix(fit.glmpca$loadings)
cat("finish fitting GLMPCA to estimate matrix of latent factors\n")
runtime <- proc.time() - start_time
fit.glmpca$runtime <- runtime
saveRDS(fit.glmpca,"fit_glmpca.Rds")

# Prefit the model to initialize parameters
# -----------------------------------------
cat("start prefit for parameter initialization\n")
prefit <- pois_mash_ruv_prefit(dat,Fuv,verbose = TRUE)
cat("finish prefit for parameter initialization\n")

# Estimate data-driven prior covariance matrices
# ----------------------------------------------
# Initialize the data-driven prior covariance matrices.
res.pca <- pois_cov_init(dat,ruv = TRUE,Fuv = Fuv,rho = prefit$rho,npc = 5)

# Run the ED step.
start_time <- proc.time()
cat("start fitting ED step\n")
fit.ed <- pois_cov_ed(dat,subset = res.pca$subset,Ulist = res.pca$Ulist,
                      ulist = res.pca$ulist,ruv = TRUE,Fuv = Fuv,
                      verbose = TRUE,init = prefit,
                      control = list(maxiter = 100))
cat("finish fitting ED step\n")
runtime <- proc.time() - start_time
fit.ed$runtime <- runtime
saveRDS(fit.ed,"pois_mash_ruv_ed.Rds")

# Run Poisson mash ruv
# --------------------
# Construct the canonical prior covariance matrices.
ulist.c <- pois_cov_canonical(dat)

# Combine all the rank-1 prior covariance matrices.
ulist <- c(fit.ed$ulist,ulist.c)

start_time <- proc.time()
cat("start fitting poisson mash with ruv\n")
res <- pois_mash(data = dat,Ulist = fit.ed$Ulist,ulist = ulist,
                 normalizeU = TRUE,gridmult = 2.5,ruv = TRUE,Fuv = Fuv,
                 rho = prefit$rho,verbose = TRUE,
                 init = list(mu = prefit$mu,psi2 = prefit$psi2),
                 control = list(maxiter = 4))
cat("finish fitting poisson mash with ruv\n")
runtime <- proc.time() - start_time
res$runtime <- runtime
saveRDS(res,"pois_mash_ruv_fit.Rds")
print(sessionInfo())
