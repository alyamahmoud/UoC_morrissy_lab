#!/bin/bash
# =======================================================================================
#SBATCH --job-name=pySCENIC  

#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40

#SBATCH --mem=900gb 
#SBATCH --time=05:00:00 

#SBATCH --partition=bigmem

# =======================================================================================
source activate.anaconda3
conda activate pySCENICgithub

./my-cluster2.py

# =======================================================================================
