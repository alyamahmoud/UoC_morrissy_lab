#!/bin/bash
#SBATCH --job-name=SCENIC  # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --ntasks=1 
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --cpus-per-task=20
#SBATCH --mem=800000                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem
#SBATCH --output=SCENIC%j.log   # Standard output and error log


# direct temporary files to alyaa.mohamed/scratch
# keep large data files at /work
module load python/anaconda3-2019.10-tensorflowgpu 
source activate scenic
module load R
#module load lib/hdf5/1.10.0.1
xvfb-run R --no-save CMD BATCH 01_scenic.R

