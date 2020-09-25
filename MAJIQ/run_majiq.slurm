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

cd /newwork/morrissy_lab/alyaa.mohamed/splicing/
#majiq build gencode.v32.annotation.gff3 -c majiq_config_file -j 1 -o output
#majiq deltapsi -grp1 output/SMP18_P1_1.majiq output/SMP18_P1_2.majiq output/SMP18_P1_3.majiq output/SMP18_P1_4.majiq -grp2 output/SMP18_R1_1.majiq output/SMP18_R1_2.majiq output/SMP18_R1_3.majiq output/SMP18_R1_4.majiq -j 1  -n SMP18_P1 SMP18_R1  -o output/
voila view output/splicegraph.sql output/SMP18_P1_SMP18_R1.deltapsi.voila
