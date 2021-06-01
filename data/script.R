setwd("/scratch/midway2/yushaliu/multivariate_poisson/poisson.mash.alpha")
library(Matrix)
library(glmpca)
library(scran)
library(psych)
library(poilog)


################################ Read in functions and data #########################################
### load in functions to run poisson mash
source("code/util.R")
source("code/set_data.R")
source("code/pois_mash_ruv_prefit.R")
source("code/pois_cov_init.R")
source("code/pois_cov_canonical.R")
source("code/pois_cov_ed.R")
source("code/pois_mash_ruv.R")
source("code/pois_mash_posterior.R")


### read in a simulated single cell dataset
scdata <- readRDS("data/scdata.Rds")

### compute cell-specific size factors using scrna
clusters <- quickCluster(scdata$Y)
si <- computeSumFactors(scdata$Y, clusters=clusters)

### create a data object for poisson mash analysis
data <- pois_mash_set_data(scdata$Y, scdata$condition, si)



############################ Estimate the matrix of latent factors causing unwanted variation #####################################
### run glmpca to estimate matrix of latent factors while adjusting for cell-specific size factors and gene-specific, condition-specific means
start_time = proc.time()
print("##########################################")
print("start fitting GLMPCA to estimate matrix of latent factors")
design <- model.matrix(~scdata$condition)
fit.glmpca <- glmpca(Y=scdata$Y, X=design[,-1], L=10, fam="nb2", sz=si, ctl=list(verbose=TRUE, tol=1e-5))
Fuv <- as.matrix(fit.glmpca$loadings)
print("finish fitting GLMPCA to estimate matrix of latent factors")
runtime = proc.time() - start_time
fit.glmpca[["runtime"]] = runtime
saveRDS(fit.glmpca, file = "output/fit_glmpca.Rds")



############################ Prefit the model to initialize parameters #####################################
print("##########################################")
print("start prefit for parameter initialization")
prefit <- pois_mash_ruv_prefit(data, Fuv, verbose=TRUE)
print("finish prefit for parameter initialization")



############################ Estimate data-driven prior covariance matrices #####################################
### initialize data-driven prior covariance matrices
res.pca <- pois_cov_init(data, ruv=TRUE, Fuv=Fuv, rho=prefit$rho, npc=5)

### run the ED step
start_time = proc.time()
print("##########################################")
print("start fitting ED step")
fit.ed <- pois_cov_ed(data, subset=res.pca$subset, Ulist=res.pca$Ulist, ulist=res.pca$ulist, ruv=TRUE, Fuv=Fuv, verbose=TRUE, init=prefit)
print("finish fitting ED step")
runtime = proc.time() - start_time
fit.ed[["runtime"]] = runtime
saveRDS(fit.ed, file = "output/pois_mash_ruv_ed.Rds")



############################ Run Poisson mash ruv #####################################
### construct canonical prior covariance matrices
ulist.c <- pois_cov_canonical(data)

### combine all the rank-1 prior covariance matrices 
ulist <- c(fit.ed$ulist, ulist.c)

start_time = proc.time()
print("##########################################")
print("start fitting poisson mash with ruv")
res <- pois_mash(data=data, Ulist=fit.ed$Ulist, ulist=ulist, normalizeU=TRUE, gridmult=2.5, ruv=TRUE, Fuv=Fuv, rho=prefit$rho, verbose=TRUE,
                 init=list(mu=prefit$mu, psi2=prefit$psi2)) 
print("finish fitting poisson mash with ruv")
runtime = proc.time() - start_time
res[["runtime"]] = runtime
saveRDS(res, file = "output/pois_mash_ruv_fit.Rds")



print(sessionInfo())