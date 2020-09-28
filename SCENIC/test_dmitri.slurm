#!/bin/bash
# =======================================================================================
#SBATCH --job-name=pySCENIC  

#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40

#SBATCH --mem=0 
#SBATCH --time=05:00:00 

#SBATCH --partition=apophis-bf

# =======================================================================================
source activate.anaconda3
conda activate pySCENIC

./my-cluster2.py

# =======================================================================================
