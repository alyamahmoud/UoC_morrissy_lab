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

export PATH=$PATH://work/morrissy_lab/alyaa.mohamed/applications/STAR-2.7.5b/source/

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AK7565__N__Tumor__RNA__A37788.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AK7565__N__Tumor__RNA__A3778.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AK7565__N__Tumor__RNA__A37788_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AL4257__N__CellLine__RNA__A29391.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AL4257__N__CellLine__RNA__A29391.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AL4257__N__CellLine__RNA__A29391_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AL4257__N__Tumor__RNA__A29392.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AL4257__N__Tumor__RNA__A29392.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AL4257__N__Tumor__RNA__A29392_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR3050__N__CellLine__RNA__A37789.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR3050__N__CellLine__RNA__A37789.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AR3050__N__CellLine__RNA__A37789_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR3050__N__Tumor__RNA__A37790.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR3050__N__Tumor__RNA__A37790.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AR3050__N__Tumor__RNA__A37790_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR5365__R__CellLine__RNA__A29393.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR5365__R__CellLine__RNA__A29393.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AR5365__R__CellLine__RNA__A29393_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR5365__R__Tumor__RNA__A29394.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-AR5365__R__Tumor__RNA__A29394.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-AR5365__R__Tumor__RNA__A29394_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-BK0248__N__CellLine__RNA__A37791.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-BK0248__N__CellLine__RNA__A37791.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-BK0248__N__CellLine__RNA__A37791_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-BM772__N__CellLine__RNA__A29395.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-BM772__N__CellLine__RNA__A29395.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-BM772__N__CellLine__RNA__A29395_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__CellLine__RNA__A29397.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__CellLine__RNA__A29397.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CA2271__N__CellLine__RNA__A29397_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Tumor__RNA__A29398.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Tumor__RNA__A29398.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CA2271__N__Tumor__RNA__A29398_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Xenograft__RNA__A52018.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CA2271__N__Xenograft__RNA__A52018.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CA2271__N__Xenograft__RNA__A52018_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A29399.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A29399.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__CellLine__RNA__A29399_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37792.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37792.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__CellLine__RNA__A37792_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37793.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__CellLine__RNA__A37793.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__CellLine__RNA__A37793_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37794.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37794.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__Tumor__RNA__A37794_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37795.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Tumor__RNA__A37795.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__Tumor__RNA__A37795_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Xenograft__RNA__A52019.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM1209__R__Xenograft__RNA__A52019.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM1209__R__Xenograft__RNA__A52019_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM3220__N__CellLine__RNA__A29400.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM3220__N__CellLine__RNA__A29400.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM3220__N__CellLine__RNA__A29400_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM3220__N__Xenograft__RNA__A52011.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-CM3220__N__Xenograft__RNA__A52011.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-CM3220__N__Xenograft__RNA__A52011_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DF5919__N__CellLine__RNA__A29401.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DF5919__N__CellLine__RNA__A29401.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DF5919__N__CellLine__RNA__A29401_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DF5919__N__Tumor__RNA__A29402.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DF5919__N__Tumor__RNA__A29402.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DF5919__N__Tumor__RNA__A29402_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DM9089__R__CellLine__RNA__A29403.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DM9089__R__CellLine__RNA__A29403.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DM9089__R__CellLine__RNA__A29403_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__CellLine__RNA__A29404.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__CellLine__RNA__A29404.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DS9789__R__CellLine__RNA__A29404_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Tumor__RNA__A29405.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Tumor__RNA__A29405.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DS9789__R__Tumor__RNA__A29405_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Xenograft__RNA__A50724.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-DS9789__R__Xenograft__RNA__A50724.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-DS9789__R__Xenograft__RNA__A50724_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__CellLine__RNA__A29406.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__CellLine__RNA__A29406.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-EV3071__N__CellLine__RNA__A29406_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Tumor__RNA__A29407.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Tumor__RNA__A29407.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-EV3071__N__Tumor__RNA__A29407_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Xenograft__RNA__A50722.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-EV3071__N__Xenograft__RNA__A50722.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-EV3071__N__Xenograft__RNA__A50722_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-FB6711__N__CellLine__RNA__A29408.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-FB6711__N__CellLine__RNA__A29408.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-FB6711__N__CellLine__RNA__A29408_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-FB6711__N__Tumor__RNA__A29409.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-FB6711__N__Tumor__RNA__A29409.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-FB6711__N__Tumor__RNA__A29409_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GB9186__N__CellLine__RNA__A29410.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GB9186__N__CellLine__RNA__A29410.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GB9186__N__CellLine__RNA__A29410_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000


STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GB9483__R__CellLine__RNA__A37796.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GB9483__R__CellLine__RNA__A37796.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GB9483__R__CellLine__RNA__A37796_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GC1519__N__CellLine__RNA__A37797.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GC1519__N__CellLine__RNA__A37797.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GC1519__N__CellLine__RNA__A37797_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GC1519__N__Tumor__RNA__A37798.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GC1519__N__Tumor__RNA__A37798.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GC1519__N__Tumor__RNA__A37798_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GE7528__R__Tumor__RNA__A29413.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GE7528__R__Tumor__RNA__A29413.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GE7528__R__Tumor__RNA__A29413_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GI2070__N__CellLine__RNA__A2941.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GI2070__N__CellLine__RNA__A2941.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GI2070__N__CellLine__RNA__A2941_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GJ3716__N__CellLine__RNA__A37799.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GJ3716__N__CellLine__RNA__A37799.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GJ3716__N__CellLine__RNA__A37799_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GJ3716__N__Tumor__RNA__A37800.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GJ3716__N__Tumor__RNA__A37800.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GJ3716__N__Tumor__RNA__A37800_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GR2309__N__CellLine__RNA__A37801.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-GR2309__N__CellLine__RNA__A37801.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-GR2309__N__CellLine__RNA__A37801_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HN6692__NA__CellLine__RNA__A37802.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HN6692__NA__CellLine__RNA__A37802.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-HN6692__NA__CellLine__RNA__A37802_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HN6692__NA__CellLine__RNA__A37802.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HN6692__NA__CellLine__RNA__A37802.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-HN6692__NA__CellLine__RNA__A37802_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HO0394__N__CellLine__RNA__A29415.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HO0394__N__CellLine__RNA__A29415.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-HO0394__N__CellLine__RNA__A29415_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HO0394__N__Xenograft__RNA__A52022.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-HO0394__N__Xenograft__RNA__A52022.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-HO0394__N__Xenograft__RNA__A52022_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JB1730__N__CellLine__RNA__A29417.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JB1730__N__CellLine__RNA__A29417.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JB1730__N__CellLine__RNA__A29417_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JB1730__N__Tumor__RNA__A29418.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JB1730__N__Tumor__RNA__A29418.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JB1730__N__Tumor__RNA__A29418_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__CellLine__RNA__A29419.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__CellLine__RNA__A29419.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JE6375__N__CellLine__RNA__A29419_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Tumor__RNA__A29420..R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Tumor__RNA__A29420..R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JE6375__N__Tumor__RNA__A29420._ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Xenograft__RNA__A50727.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JE6375__N__Xenograft__RNA__A50727.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JE6375__N__Xenograft__RNA__A50727_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JP2405__R__CellLine__RNA__A37803.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JP2405__R__CellLine__RNA__A37803.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JP2405__R__CellLine__RNA__A37803_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JP2405__R__Tumor__RNA__A37804.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JP2405__R__Tumor__RNA__A37804.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JP2405__R__Tumor__RNA__A37804_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JR9883__N__CellLine__RNA__A29421.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JR9883__N__CellLine__RNA__A29421.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JR9883__N__CellLine__RNA__A29421_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JW6420__N__CellLine__RNA__A29422.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JW6420__N__CellLine__RNA__A29422.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JW6420__N__CellLine__RNA__A29422_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JW6420__N__Tumor__RNA__A29423.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-JW6420__N__Tumor__RNA__A29423.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-JW6420__N__Tumor__RNA__A29423_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-KM5291__N__CellLine__RNA__A29424.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-KM5291__N__CellLine__RNA__A29424.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-KM5291__N__CellLine__RNA__A29424_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-KM5291__N__Tumor__RNA__A29425.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-KM5291__N__Tumor__RNA__A29425.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-KM5291__N__Tumor__RNA__A29425_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__CellLine__RNA__A29426.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__CellLine__RNA__A29426.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LC3356__N__CellLine__RNA__A29426_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Tumor__RNA__A29427.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Tumor__RNA__A29427.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LC3356__N__Tumor__RNA__A29427_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Xenograft__RNA__A52023.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC3356__N__Xenograft__RNA__A52023.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LC3356__N__Xenograft__RNA__A52023_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC6372__N__CellLine__RNA__A37805.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LC6372__N__CellLine__RNA__A37805.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LC6372__N__CellLine__RNA__A37805_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LR9369__R__CellLine__RNA__A29428.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LR9369__R__CellLine__RNA__A29428.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LR9369__R__CellLine__RNA__A29428_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LR9369__R__Tumor__RNA__A29429.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LR9369__R__Tumor__RNA__A29429.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LR9369__R__Tumor__RNA__A29429_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LS4891__N__CellLine__RNA__A29430.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LS4891__N__CellLine__RNA__A29430.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LS4891__N__CellLine__RNA__A29430_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LS4891__N__Tumor__RNA__A29431.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-LS4891__N__Tumor__RNA__A29431.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-LS4891__N__Tumor__RNA__A29431_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A29432.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A29432.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37806.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37806.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__CellLine__RNA__A37806_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37807.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__CellLine__RNA__A37807.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__CellLine__RNA__A37807_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A29433.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A29433.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__Tumor__RNA__A29433_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37808.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37808.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__Tumor__RNA__A37808_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37809.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Tumor__RNA__A37809.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__Tumor__RNA__A37809_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50723.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50723.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__Xenograft__RNA__A50723_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50732.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MB9777__R__Xenograft__RNA__A50732.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MB9777__R__Xenograft__RNA__A50732_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MD9088__R__CellLine__RNA__A37810.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MD9088__R__CellLine__RNA__A37810.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MD9088__R__CellLine__RNA__A37810_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MD9088__R__Tumor__RNA__A37811.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-MD9088__R__Tumor__RNA__A37811.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-MD9088__R__Tumor__RNA__A37811_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-ML9537__N__CellLine__RNA__A37812.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-ML9537__N__CellLine__RNA__A37812.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-ML9537__N__CellLine__RNA__A37812_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PD6881__R__CellLine__RNA__A29435.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PD6881__R__CellLine__RNA__A29435.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PD6881__R__CellLine__RNA__A29435_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PD6881__R__Tumor__RNA__A29436.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PD6881__R__Tumor__RNA__A29436.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PD6881__R__Tumor__RNA__A29436_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PR5617__N__CellLine__RNA__A29437.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PR5617__N__CellLine__RNA__A29437.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PR5617__N__CellLine__RNA__A29437_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A29438.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A29438.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PV2594__N__CellLine__RNA__A29438_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52007.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52007.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PV2594__N__CellLine__RNA__A52007_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52008.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__CellLine__RNA__A52008.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PV2594__N__CellLine__RNA__A52008_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Tumor__RNA__A29439.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Tumor__RNA__A29439.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PV2594__N__Tumor__RNA__A29439_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Xenograft__RNA__A50726.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-PV2594__N__Xenograft__RNA__A50726.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-PV2594__N__Xenograft__RNA__A50726_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RD1291__N__CellLine__RNA__A29416.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RD1291__N__CellLine__RNA__A29416.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RD1291__N__CellLine__RNA__A29416_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RD1291__N__Tumor__RNA__A37813.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RD1291__N__Tumor__RNA__A37813.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RD1291__N__Tumor__RNA__A37813_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL5404__N__CellLine__RNA__A29440.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL5404__N__CellLine__RNA__A29440.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RL5404__N__CellLine__RNA__A29440_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL5404__N__Tumor__RNA__A29441.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL5404__N__Tumor__RNA__A29441.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RL5404__N__Tumor__RNA__A29441_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL7940__N__CellLine__RNA__A29442.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL7940__N__CellLine__RNA__A29442.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RL7940__N__CellLine__RNA__A29442_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL7940__N__Tumor__RNA__A29443.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RL7940__N__Tumor__RNA__A29443.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RL7940__N__Tumor__RNA__A29443_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RV2286__R__CellLine__RNA__A37814.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RV2286__R__CellLine__RNA__A37814.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RV2286__R__CellLine__RNA__A37814_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A2944.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A2944.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__CellLine__RNA__A2944_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000


STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52009.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52009.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__CellLine__RNA__A52009_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52010.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__CellLine__RNA__A52010.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__CellLine__RNA__A52010_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A30194.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A30194.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__Tumor__RNA__A30194_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A37815.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Tumor__RNA__A37815.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__Tumor__RNA__A37815_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Xenograft__RNA__A50730..R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-RW9277__NA__Xenograft__RNA__A50730..R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-RW9277__NA__Xenograft__RNA__A50730._ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__CellLine__RNA__A29446.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__CellLine__RNA__A29446.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SJ5453__N__CellLine__RNA__A29446_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Tumor__RNA__A29447.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Tumor__RNA__A29447.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SJ5453__N__Tumor__RNA__A29447_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Xenograft__RNA__A50725.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SJ5453__N__Xenograft__RNA__A50725.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SJ5453__N__Xenograft__RNA__A50725_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SK0976__N__CellLine__RNA__A29448.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SK0976__N__CellLine__RNA__A29448.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SK0976__N__CellLine__RNA__A29448_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SK0976__N__Tumor__RNA__A37816.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SK0976__N__Tumor__RNA__A37816.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SK0976__N__Tumor__RNA__A37816_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SO0258__N__CellLine__RNA__A29449.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SO0258__N__CellLine__RNA__A29449.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SO0258__N__CellLine__RNA__A29449_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SO0258__N__Tumor__RNA__A29450.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SO0258__N__Tumor__RNA__A29450.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SO0258__N__Tumor__RNA__A29450_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SS3647__N__CellLine__RNA__A29451.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SS3647__N__CellLine__RNA__A29451.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SS3647__N__CellLine__RNA__A29451_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SS3647__N__Tumor__RNA__A29452.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-SS3647__N__Tumor__RNA__A29452.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-SS3647__N__Tumor__RNA__A29452_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-TM5196__R__CellLine__RNA__A37817.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-TM5196__R__CellLine__RNA__A37817.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-TM5196__R__CellLine__RNA__A37817_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-TM5196__R__Tumor__RNA__A37818.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-TM5196__R__Tumor__RNA__A37818.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-TM5196__R__Tumor__RNA__A37818_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT--VO7089__N__CellLine__RNA__A37819.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT--VO7089__N__CellLine__RNA__A37819.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT--VO7089__N__CellLine__RNA__A37819_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-VO7089__N__Tumor__RNA__A29454.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-VO7089__N__Tumor__RNA__A29454.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-VO7089__N__Tumor__RNA__A29454_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WP9124__R__CellLine__RNA__A37820.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WP9124__R__CellLine__RNA__A37820.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-WP9124__R__CellLine__RNA__A37820_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WP9124__R__Tumor__RNA__A37821.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WP9124__R__Tumor__RNA__A37821.R2.fastq.gz  \
--readFilesCommand zcat \ 
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-WP9124__R__Tumor__RNA__A37821_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WR7927__N__CellLine__RNA__A37822.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WR7927__N__CellLine__RNA__A37822.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000

STAR --genomeDir /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/indices_mus_spiked_new.97/ \
-- runThreadN 10 \
--readFilesIn /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WT4796__N__Tumor__RNA__A29456.R1.fastq.gz /work/morrissy_lab/TFRI_2018_data/bamrnacopy/newfastq/PT-WT4796__N__Tumor__RNA__A29456.R2.fastq.gz  \
--readFilesCommand zcat \
--outFileNamePrefix /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-mg38-rel97/bams/PT-WT4796__N__Tumor__RNA__A29456_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterType BySJout \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outSAMstrandField intronMotif \
--limitBAMsortRAM 70000000000



