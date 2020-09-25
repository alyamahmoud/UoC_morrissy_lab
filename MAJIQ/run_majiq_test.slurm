#!/bin/bash
##SBATCH --job-name=majiq   # Job name
##SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
##SBATCH --nodes=1                   # Run on a single CPU
##SBATCH --ntasks=1
##SBATCH --ntasks-per-node=1
###SBATCH --cpus-per-task=1
##SBATCH --mem=8000                     # Job memory request
##SBATCH --time=48:00:00              # Time limit hrs:min:sec
##SBATCH --output=majiq%j.log   # Standard output and error log
##SBATCH --partition=cpu2019

module load python/anaconda3-2019.10-tensorflowgpu
source activate splicing

cd /work/morrissylab/alyaa.mohamed/splicing/workshop_example
majiq build DB.gff3 -c settings.ini -o build_groups --disable-denovo 
