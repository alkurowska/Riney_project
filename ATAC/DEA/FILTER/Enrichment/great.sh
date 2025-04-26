#!/bin/bash -l 
## SLURM Resource requirement
#SBATCH --time 24:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH -c 8

module load R
Rscript 06_GREAT.R