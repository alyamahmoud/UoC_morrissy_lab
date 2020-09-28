#!/bin/bash
#SBATCH --job-name=pySCENIC    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=85000                     # Job memory request
#SBATCH --time=5:00:00              # Time limit hrs:min:sec
#SBATCH --output=pySCENIC_%j.log   # Standard output and error log
#SBATCH --partition=cpu2019


module load python/anaconda3-2019.10-tensorflowgpu
module load R/3.5.3
module load lib/hdf5/1.10.0.1
cd /home/alyaa.mohamed/
source activate scenic 




cd /work/morrissylab/alyaa.mohamed/for_dmitri
python  pySCENIC.py
