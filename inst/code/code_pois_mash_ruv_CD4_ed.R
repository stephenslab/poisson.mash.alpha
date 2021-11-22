setwd("/scratch/midway2/yushaliu/poisson.mash.alpha")
devtools::load_all()

### specifiy the selected conditions measured on the same batch
trts <- c("Ctrl_2", "CCL20", "CXCL1", "CCL22", "CXCL5", "CCL11", "CCL4", "CCL17", "CCL5", "CXCL13", "CXCL10", "CXCL9",  
          "CXCL12", "GCSF", "MCSF", "GMCSF", "IFNg", "IL10", "IL12p70", "IL17a", "IL13", "IL15", "IL17f", "IL22",
          "IL18", "IL1a", "IL2", "IL3", "IL1b", "IL23", "IL21", "IL33", "IL25", "IL34", "IL27", "IL36a", 
          "IL4", "IL6", "IL5", "IL7", "IL9", "IL11", "TGFb", "CCL2", "CCL3", "TSLP") 
R <- length(trts)


### load in the aggregated count data and size factors
data.jr <- readRDS("inst/datafiles/data_jr_CD4.Rds")
data <- list(X=data.jr$data.jr, s=data.jr$s, subgroup=rep(1,R))
rm(data.jr)

### load in the GLM-PCA fit 
fit.glmpca <- readRDS("inst/datafiles/glmpca_CD4.Rds")

### get the J x D factor matrix related to unwanted variation
Fuv <- as.matrix(fit.glmpca$loadings)
sum(rownames(data$X)!=rownames(Fuv))



################################ Run prefit step by ignoring fixed effects #########################################
start_time = proc.time()
print("##########################################")
print("start prefit scdata from CD4 T cells")
prefit <- pois_mash_ruv_prefit(data, Fuv, verbose=TRUE, control=list(maxiter=500))
print("finish prefit scdata from CD4 T cells")
proc.time() - start_time
saveRDS(prefit, file = "inst/outputs/pois_mash_ruv_prefit_CD4.Rds")



################################ Run ED step #########################################
### initialize data-driven prior covariance matrices
res.pca <- pois_cov_init(data, ruv=FALSE, npc=5, cutoff=abs(qnorm(0.05/2/R)))

# rank-1 prior covariances
ulist.c <- pois_cov_canonical(data)
ulist <- c(res.pca$ulist, ulist.c)
ulist.dd <- c(rep(TRUE, length(res.pca$ulist)-1), rep(FALSE, R+1))

### run the ED step
start_time = proc.time()
print("##########################################")
print("start fitting ED of possion mash ruv for scdata from CD4 T cells")
fit.ed <- pois_cov_ed(data, subset=res.pca$subset, Ulist=res.pca$Ulist, ulist=ulist, ulist.dd=ulist.dd, ruv=TRUE, Fuv=Fuv, verbose=TRUE, 
                      control=list(maxiter=300, maxiter.q=25, maxpsi2=log(2), maxbias=1, tol.q=1e-2, tol.rho=1e-3, tol.stop=1e-1))
print("finish fitting ED of possion mash ruv for scdata from CD4 T cells")
runtime = proc.time() - start_time
fit.ed[["runtime"]] = runtime
saveRDS(fit.ed, file = "inst/outputs/pois_mash_ruv_ed_CD4.Rds")



print(sessionInfo())