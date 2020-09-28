#!/bin/bash
#SBATCH --job-name=DEXSeq    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=85000                     # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=DEXSeq_%j.log   # Standard output and error log
#SBATCH --partition=cpu2019


module load python/anaconda3-2019.10-tensorflowgpu
module load R/3.5.3
module load lib/hdf5/1.10.0.1
#cd /home/alyaa.mohamed/
source activate splicing 




cd /work/morrissy_lab/alyaa.mohamed/splicing/
python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_P1_3.bam  exon_counts/SMP18_P1_3

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_P1_4.bam  exon_counts/SMP18_P1_4

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_R1_1.bam  exon_counts/SMP18_R1_1

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_R1_2.bam  exon_counts/SMP18_R1_2

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_R1_3.bam  exon_counts/SMP18_R1_3

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP18_R1_4.bam  exon_counts/SMP18_R1_4

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_P1_1.bam  exon_counts/SMP8_P1_1

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_P1_2.bam  exon_counts/SMP8_P1_2

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R1_1.bam  exon_counts/SMP8_R1_1

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R1_2.bam  exon_counts/SMP8_R1_2

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R1_3.bam  exon_counts/SMP8_R1_3

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R1_4.bam  exon_counts/SMP8_R1_4

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_1.bam  exon_counts/SMP8_R2_1

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_2.bam  exon_counts/SMP8_R2_2

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_3.bam  exon_counts/SMP8_R2_3

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_4.bam  exon_counts/SMP8_R2_4

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_5.bam  exon_counts/SMP8_R2_5

python /home/alyaa.mohamed/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos gencode.v32.annotation.DEXSeq.chr.gff data/SMP8_R2_6.bam  exon_counts/SMP8_R2_6


