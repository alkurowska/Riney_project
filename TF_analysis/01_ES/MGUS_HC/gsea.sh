#!/bin/bash -l 
## SLURM Resource requirement
#SBATCH --time 30:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH -c 8

module load R
Rscript GSEA_TF_abs.R
Rscript GSEA_TF_pos.R
Rscript GSEA_TF_neg.R