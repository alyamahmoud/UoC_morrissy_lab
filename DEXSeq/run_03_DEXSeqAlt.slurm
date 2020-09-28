#!/bin/bash
#SBATCH --job-name=DEXSeqAlt   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=90000                    # Job memory request
#SBATCH --time=10:30:00              # Time limit hrs:min:sec
#SBATCH --output=DEXSeqAlt%j.log   # Standard output and error log
#SBATCH --partition=cpu2019

module load python/anaconda3-2019.10-tensorflowgpu
module load R
module load lib/hdf5/1.10.0.1
source activate splicing

cd /newwork/morrissy_lab/alyaa.mohamed/pedGBM/
xvfb-run R --no-save CMD BATCH 03_DEXSeqAlt.R
