#!/bin/bash
#SBATCH --job-name=star_map    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=2000000        # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=star_map_%j.log   # Standard output and error log
#SBATCH --partition=bigmem

# this script runs all trios, soft clipping extend 5pR1 

export PATH=$PATH://work/morrissy_lab/alyaa.mohamed/applications/STAR-2.7.5b/source/

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__CellLine__RNA__A29397.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__CellLine__RNA__A29397.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CA2271__N__CellLine__RNA__A29397_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Tumor__RNA__A29398.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Tumor__RNA__A29398.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CA2271__N__Tumor__RNA__A29398_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Xenograft__RNA__A52018.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Xenograft__RNA__A52018.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CA2271__N__Xenograft__RNA__A52018_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A29399.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A29399.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__CellLine__RNA__A29399_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37792.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37792.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__CellLine__RNA__A37792_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37793.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37793.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__CellLine__RNA__A37793_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37794.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37794.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__Tumor__RNA__A37794_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37795.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37795.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__Tumor__RNA__A37795_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Xenograft__RNA__A52019.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Xenograft__RNA__A52019.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-CM1209__R__Xenograft__RNA__A52019_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__CellLine__RNA__A29404.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__CellLine__RNA__A29404.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-DS9789__R__CellLine__RNA__A29404_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Tumor__RNA__A29405.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Tumor__RNA__A29405.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-DS9789__R__Tumor__RNA__A29405_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Xenograft__RNA__A50724.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Xenograft__RNA__A50724.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-DS9789__R__Xenograft__RNA__A50724_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__CellLine__RNA__A29406.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__CellLine__RNA__A29406.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-EV3071__N__CellLine__RNA__A29406_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Tumor__RNA__A29407.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Tumor__RNA__A29407.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-EV3071__N__Tumor__RNA__A29407_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Xenograft__RNA__A50722.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Xenograft__RNA__A50722.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-EV3071__N__Xenograft__RNA__A50722_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__CellLine__RNA__A29419.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__CellLine__RNA__A29419.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-JE6375__N__CellLine__RNA__A29419_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Tumor__RNA__A29420.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Tumor__RNA__A29420.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-JE6375__N__Tumor__RNA__A29420_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Xenograft__RNA__A50727.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Xenograft__RNA__A50727.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-JE6375__N__Xenograft__RNA__A50727_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__CellLine__RNA__A29426.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__CellLine__RNA__A29426.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-LC3356__N__CellLine__RNA__A29426_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Xenograft__RNA__A52023.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Xenograft__RNA__A52023.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-LC3356__N__Xenograft__RNA__A52023_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37806.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37806.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__CellLine__RNA__A37806_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37807.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37807.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__CellLine__RNA__A37807_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A29433.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A29433.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__Tumor__RNA__A29433_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37808.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37808.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__Tumor__RNA__A37808_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37809.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37809.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__Tumor__RNA__A37809_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50723.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50723.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__Xenograft__RNA__A50723_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50732.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50732.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-MB9777__R__Xenograft__RNA__A50732_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A29438.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A29438.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-PV2594__N__CellLine__RNA__A29438_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52008.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52008.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-PV2594__N__CellLine__RNA__A52008_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Tumor__RNA__A29439.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Tumor__RNA__A29439.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-PV2594__N__Tumor__RNA__A29439_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Xenograft__RNA__A50726.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Xenograft__RNA__A50726.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-PV2594__N__Xenograft__RNA__A50726_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A2944.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A2944.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__CellLine__RNA__A2944_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1


STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52009.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52009.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__CellLine__RNA__A52009_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52010.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52010.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__CellLine__RNA__A52010_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A30194.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A30194.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__Tumor__RNA__A30194_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A37815.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A37815.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__Tumor__RNA__A37815_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Xenograft__RNA__A50730.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Xenograft__RNA__A50730.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-RW9277__NA__Xenograft__RNA__A50730_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__CellLine__RNA__A29446.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__CellLine__RNA__A29446.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-SJ5453__N__CellLine__RNA__A29446_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Tumor__RNA__A29447.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Tumor__RNA__A29447.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-SJ5453__N__Tumor__RNA__A29447_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/indices_spiked.am/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Xenograft__RNA__A50725.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Xenograft__RNA__A50725.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/bams/scratch_extend5pR1/PT-SJ5453__N__Xenograft__RNA__A50725_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000 \
--alignEndsType Extend5pOfRead1


