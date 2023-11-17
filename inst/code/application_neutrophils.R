library(Matrix)
library(glmpca)
library(scran)
library(poisson.mash.alpha)


################################ Select the raw count data and sample annotation information for neutrophils #########################################
### specifiy the 45 conditions measured on the same batch
trts <- c("Ctrl_2", "CCL20", "CXCL1", "CCL22", "CXCL5", "CCL11", "CCL4", "CCL17", "CCL5", "CXCL13", "CXCL10", "CXCL9",  
          "CXCL12", "GCSF", "MCSF", "GMCSF", "IFNg", "IL10", "IL12p70", "IL17a", "IL13", "IL15", "IL17f", "IL22",
          "IL18", "IL1a", "IL2", "IL3", "IL1b", "IL23", "IL21", "IL33", "IL25", "IL34", "IL36a", 
          "IL4", "IL6", "IL5", "IL7", "IL9", "IL11", "TGFb", "CCL2", "CCL3", "TSLP") 
R <- length(trts)

### load in the UMI count data and sample annotation information for all cells generated from this experiment
scdata <- readRDS("whole_cyto_scdata.rds")
sample.info <- read.csv("whole_cyto_annot.csv")
sum(colnames(scdata)!=sample.info$X0)

### select the neutrophils receiving one of the 45 chosen treatments
sample.info.sub <- sample.info[sample.info$cell_type=="Neutrophils" & sample.info$sample %in% trts,]
scdata.sub <- scdata[, sample.info$cell_type=="Neutrophils" & sample.info$sample %in% trts]
sum(colnames(scdata.sub)!=sample.info.sub$X0)

### remove genes with very few counts
idx.gene <- which(rowSums(scdata.sub)>=25)
scdata.sub <- scdata.sub[idx.gene,]

### save the neutrophils data
Y <- scdata.sub
conditions <- factor(sample.info.sub$sample, levels=trts)
save(Y, conditions, file = "data/neutrophils_whole.RData")



################################ Create a data object for Poisson mash analysis #########################################
### compute cell-specific size factors
clusters <- quickCluster(Y)
si <- calculateSumFactors(Y, clusters=clusters)
names(si) <- colnames(Y)

### create the poisson mash data object
data <- pois_mash_set_data(Y, conditions, si)



################################ Estimate latent factors capturing unwanted variation #########################################
### construct the design matrix for conditions
design <- model.matrix(~conditions)
design <- design[,-1]
colnames(design) <- trts[-1]

### run glm-pca while adjusting for cell-specific size factors and gene-specific, condition-specific means
start_time = proc.time()
print("##########################################")
print("start fitting GLM-PCA for scdata from Neutrophils")
fit.glmpca <- glmpca(Y=as.matrix(Y), X=design, L=4, fam="nb2", sz=si, ctl=list(verbose=TRUE, lr=0.01, tol=1e-6))
print("finish fitting GLM-PCA for scdata from Neutrophils")
runtime = proc.time() - start_time
print(runtime)

### get the J x D latent factor matrix related to unwanted variation
Fuv <- as.matrix(fit.glmpca$loadings)
sum(rownames(Fuv)!=rownames(Y))



################################ Prefit the Poisson mash model without DE effects to initialize parameters #########################################
start_time = proc.time()
print("##########################################")
print("start pre-fitting model for scdata from Neutrophils")
prefit <- pois_mash_prefit(data, ruv=TRUE, Fuv=Fuv, verbose=TRUE, control=list(maxiter=500))
print("finish pre-fitting model for scdata from Neutrophils")
runtime = proc.time() - start_time
print(runtime)



############################### Set up the prior covariance matrices, in particular estimate data-driven ones ######################################
### construct canonical rank-1 prior covariance matrices that model condition-specific effects for each condition
ulist.c <- pois_cov_canonical(data)

### initialize data-driven prior covariance matrices
res.pca <- pois_cov_init(data, npc=5, cutoff=abs(qnorm(0.05/2/R)))

### merge data-driven and canonical rank-1 prior covariance matrices, each of which is represented by an R x 1 vector
ulist <- c(res.pca$ulist, ulist.c)

### distinguish the data-driven rank-1 prior covariance matrices from the canonical ones
ulist.dd <- c(rep(TRUE, length(res.pca$ulist)-1), rep(FALSE, R+1))

### run the Extreme Eeconvolution (ED) algorithm to refine the initial estimates of data-driven covariance matrices
start_time = proc.time()
print("##########################################")
print("start fitting ED for scdata from Neutrophils")
fit.ed <- pois_cov_ed(data, subset=res.pca$subset, Ulist=res.pca$Ulist, ulist=ulist, ulist.dd=ulist.dd, ruv=TRUE, Fuv=Fuv, verbose=TRUE, 
                      init=prefit, control=list(maxiter=300, maxpsi2=log(2), maxbias=1, tol.stop=1e-6))
print("finish fitting ED for scdata from Neutrophils")
runtime = proc.time() - start_time
print(runtime)

### add epsilon2*I to each full-rank data-driven prior covariance matrix 
### all the full-rank covariance matrices are data-driven and need this modification
Ulist <- fit.ed$Ulist
H <- length(Ulist)
for(h in 1:H){
  Uh <- Ulist[[h]]
  Uh <- Uh/max(diag(Uh))
  Ulist[[h]] <- Uh + 1e-2*diag(R)
}

### add epsilon2*I to each rank-1 data-driven prior covariance matrix
### only a subset of the rank-1 covariance matrices are data-driven and need this modification
G <- length(fit.ed$ulist)
epsilon2.G <- rep(1e-8, G)
names(epsilon2.G) <- names(fit.ed$ulist)
epsilon2.G[ulist.dd] <- 1e-2



################################ Fit the Poisson mash model #########################################
start_time = proc.time()
print("##########################################")
print("start fitting poisson mash with ruv for scdata from Neutrophils")
res <- pois_mash(data=data, Ulist=Ulist, ulist=fit.ed$ulist, ulist.epsilon2=epsilon2.G,
                 ruv=TRUE, Fuv=Fuv, rho=prefit$rho, verbose=TRUE, median_deviations=TRUE, 
                 init=list(mu=prefit$mu, psi2=prefit$psi2),
                 control=list(maxiter=300, maxpsi2=log(2), maxbias=1, nc=3)) 
print("finish fitting poisson mash with ruv for scdata from Neutrophils")
runtime = proc.time() - start_time
res[["runtime"]] = runtime
save(res, file = "data/neutrophils_pois_mash_ruv_fit.RData")



print(sessionInfo())