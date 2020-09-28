#!/bin/bash
#SBATCH --job-name=xenofilteR    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=2000000        # Job memory request
#SBATCH --time=4:00:00              # Time limit hrs:min:sec
#SBATCH --output=xenofilteR_%j.log   # Standard output and error log
#SBATCH --partition=bigmem

module load R/3.5.3 gcc/6.1.0 lib/hdf5/1.10.0.1
xvfb-run R --no-save CMD BATCH xenofilter.R
