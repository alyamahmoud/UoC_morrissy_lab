#!/bin/bash
#SBATCH --job-name=star_index    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=3000000       # Job memory request
#SBATCH --time=5:00:00              # Time limit hrs:min:sec
#SBATCH --output=star_index_%j.log   # Standard output and error log
#SBATCH --partition=bigmem

export PATH=$PATH:/work/morrissy_lab/alyaa.mohamed/applications/STAR-2.7.5b/source/
STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir ./indices_spiked.am/ \
--genomeFastaFiles ./Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ./Homo_sapiens.GRCh38.99.gtf \
--sjdbOverhang 100 \
#--limitGenomeGenerateRAM 30000000000 \
#--genomeSAsparseD 3 \
#--genomeSAindexNbases 12 \
#--genomeChrBinNbits 16
#STAR --runThreadN 12 --runMode genomeGenerate --genomeSAsparseD 3 --genomeSAindexNbases 12 -- genomeChrBinNbits 14 --genomeDir /mnt/data/tg-hp/human_transcriptome --genomeFastaFiles /mnt/data/tg-hp/human_transcriptome/GCF_000001405.36_GRCh38.p10_rna.fna
