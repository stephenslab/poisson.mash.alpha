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
prefit <- readRDS("inst/outputs/pois_mash_ruv_prefit_CD4.Rds")



################################ Run ED step #########################################
fit.ed <- readRDS("inst/outputs/pois_mash_ruv_ed_CD4.Rds")

# set ulist.epsilon2
H <- length(fit.ed$Ulist)
G <- length(fit.ed$ulist)
epsilon2.G <- rep(1e-8, G)
names(epsilon2.G) <- names(fit.ed$ulist)

for(g in 1:G){
  if(fit.ed$pi[H+g]>1e-2 & sum(fit.ed$ulist[[g]]!=0)>0){
    epsilon2.G[g] <- 1e-2
  }
}



################################ Run poisson mash #########################################
start_time = proc.time()
print("##########################################")
print("start fitting poisson mash with ruv for scdata from CD4 T cells")
res <- pois_mash(data=data, Ulist=fit.ed$Ulist, ulist=fit.ed$ulist, ulist.epsilon2=epsilon2.G, normalizeU=TRUE, gridmult=2, 
                 ruv=TRUE, Fuv=Fuv, rho=prefit$rho, verbose=TRUE, init=list(mu=prefit$mu, psi2=prefit$psi2), median_deviations=TRUE,
                 control=list(maxiter=300, maxiter.q=25, maxpsi2=log(2), maxbias=1, 
                              tol.q=1e-2, tol.rho=1e-3, tol.mu=1e-2, tol.psi2=2e-2, tol.bias=2e-2,
                              nc=1)) 
print("finish fitting poisson mash with ruv for scdata from CD4 T cells")
runtime = proc.time() - start_time
res[["runtime"]] = runtime
saveRDS(res, file = "inst/outputs/pois_mash_ruv_fit_CD4_nc1.Rds")



print(sessionInfo())