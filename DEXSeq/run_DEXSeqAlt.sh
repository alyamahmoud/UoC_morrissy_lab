#! /usr/bin/env bash

#BSUB -J DEXSeqAlt
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -o DEXSeqAlt_%J.out
#BSUB -e DEXSeqAlt_%J.err

source activate splicing
cd /tiered/smorrissy/TFRI_2018_data/bamrnacopy/map-hg38-rel96/
Rscript 02_DEXSeqAlt_2.R   
