#!/bin/bash

#SBATCH --job-name=cd4_pois_nc2
#SBATCH --output=cd4_pois_nc2.out
#SBATCH --error=cd4_pois_nc2.err
#SBATCH --partition=mstephens
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=2
#SBATCH --mem=72G
#SBATCH --time=24:00:00


module load R/4.1.0

echo "Running code_pois_mash_ruv_CD4_nc2.R"
export OPENBLAS_NUM_THREADS=1
Rscript /scratch/midway2/yushaliu/poisson.mash.alpha/inst/code/code_pois_mash_ruv_CD4_nc2.R  > pois_mash_ruv_CD4_nc2.Rout 2>&1


