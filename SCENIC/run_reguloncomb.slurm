#!/bin/bash
#SBATCH --job-name=reguloncomb  # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --ntasks=1 
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=185000                  # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --partition=cpu2019
#SBATCH --output=reguloncomb%j.log   # Standard output and error log


# direct temporary files to alyaa.mohamed/scratch
# keep large data files at /work
module load python/anaconda3-2019.10-tensorflowgpu 
source activate scenic
module load R
#module load lib/hdf5/1.10.0.1
xvfb-run R --no-save CMD BATCH Regulon_comb.R

