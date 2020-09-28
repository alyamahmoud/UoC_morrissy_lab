#!/bin/bash
# =======================================================================================
#SBATCH --job-name=pySCENIC  
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail

#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40

#SBATCH --mem=1500000
#SBATCH --time=5:00:00 

#SBATCH --partition=bigmem

# =======================================================================================
module load python/anaconda3-2019.10-tensorflowgpu gcc/7.3.0 
#source activate.anaconda3
#conda activate pySCENIC
source activate scenic

#./my-cluster2.py
python 01_pyscenic.py
# =======================================================================================
