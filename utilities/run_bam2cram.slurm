#!/bin/bash
#SBATCH --job-name=bam2cram  # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=alyaa.mohamed@ucalgary.ca     # Where to send mail  
#SBATCH --ntasks=1 
#SBATCH --nodes=1                   # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=185000                  # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --partition=cpu2019
#SBATCH --output=bam2cram%j.log   # Standard output and error log



export PATH=$PATH:/work/morrissy_lab/alyaa.mohamed/applications/samtools-1.10/
samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB6372__R__CellLine__RNA__A29389-Aligned.sortedByCoord.out.cram PT-AB6372__R__CellLine__RNA__A29389-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o 
PT-AB0029__N__Tumor__RNA__A29388-Aligned.sortedByCoord.out.cram PT-AB0029__N__Tumor__RNA__A29388-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB0029__N__Tumor__RNA__A29388-Aligned.toTranscriptome.out.cram PT-AB0029__N__Tumor__RNA__A29388-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB6372__R__CellLine__RNA__A29389-Aligned.sortedByCoord.out.cram PT-AB6372__R__CellLine__RNA__A29389-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB6372__R__CellLine__RNA__A29389-Aligned.toTranscriptome.out.cram PT-AB6372__R__CellLine__RNA__A29389-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB6372__R__Tumor__RNA__A37786-Aligned.sortedByCoord.out.cram PT-AB6372__R__Tumor__RNA__A37786-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AB6372__R__Tumor__RNA__A37786-Aligned.toTranscriptome.out.cram PT-AB6372__R__Tumor__RNA__A37786-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AH1410__N__CellLine__RNA__A29390-Aligned.sortedByCoord.out.cram PT-AH1410__N__CellLine__RNA__A29390-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AH1410__N__CellLine__RNA__A29390-Aligned.toTranscriptome.out.cram PT-AH1410__N__CellLine__RNA__A29390-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AK7565__N__CellLine__RNA__A37787-Aligned.sortedByCoord.out.cram PT-AK7565__N__CellLine__RNA__A37787-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AK7565__N__CellLine__RNA__A37787-Aligned.toTranscriptome.out.cram PT-AK7565__N__CellLine__RNA__A37787-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AK7565__N__Tumor__RNA__A37788-Aligned.sortedByCoord.out.cram PT-AK7565__N__Tumor__RNA__A37788-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AK7565__N__Tumor__RNA__A37788-Aligned.toTranscriptome.out.cram PT-AK7565__N__Tumor__RNA__A37788-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AL4257__N__CellLine__RNA__A29391-Aligned.sortedByCoord.out.cram PT-AL4257__N__CellLine__RNA__A29391-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AL4257__N__CellLine__RNA__A29391-Aligned.toTranscriptome.out.cram PT-AL4257__N__CellLine__RNA__A29391-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AL4257__N__Tumor__RNA__A29392-Aligned.sortedByCoord.out.cram PT-AL4257__N__Tumor__RNA__A29392-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AL4257__N__Tumor__RNA__A29392-Aligned.toTranscriptome.out.cram PT-AL4257__N__Tumor__RNA__A29392-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR3050__N__CellLine__RNA__A37789-Aligned.sortedByCoord.out.cram PT-AR3050__N__CellLine__RNA__A37789-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR3050__N__CellLine__RNA__A37789-Aligned.toTranscriptome.out.cram PT-AR3050__N__CellLine__RNA__A37789-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR3050__N__Tumor__RNA__A37790-Aligned.sortedByCoord.out.cram PT-AR3050__N__Tumor__RNA__A37790-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR3050__N__Tumor__RNA__A37790-Aligned.toTranscriptome.out.cram PT-AR3050__N__Tumor__RNA__A37790-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR5365__R__CellLine__RNA__A29393-Aligned.sortedByCoord.out.cram PT-AR5365__R__CellLine__RNA__A29393-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR5365__R__CellLine__RNA__A29393-Aligned.toTranscriptome.out.cram PT-AR5365__R__CellLine__RNA__A29393-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR5365__R__Tumor__RNA__A29394-Aligned.sortedByCoord.out.cram PT-AR5365__R__Tumor__RNA__A29394-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-AR5365__R__Tumor__RNA__A29394-Aligned.toTranscriptome.out.cram PT-AR5365__R__Tumor__RNA__A29394-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-BK0248__N__CellLine__RNA__A37791-Aligned.sortedByCoord.out.cram PT-BK0248__N__CellLine__RNA__A37791-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-BK0248__N__CellLine__RNA__A37791-Aligned.toTranscriptome.out.cram PT-BK0248__N__CellLine__RNA__A37791-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-BM772__N__CellLine__RNA__A29395-Aligned.sortedByCoord.out.cram PT-BM772__N__CellLine__RNA__A29395-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-BM772__N__CellLine__RNA__A29395-Aligned.toTranscriptome.out.cram PT-BM772__N__CellLine__RNA__A29395-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__CellLine__RNA__A29397-Aligned.sortedByCoord.out.cram PT-CA2271__N__CellLine__RNA__A29397-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__CellLine__RNA__A29397-Aligned.toTranscriptome.out.cram PT-CA2271__N__CellLine__RNA__A29397-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__Tumor__RNA__A29398-Aligned.sortedByCoord.out.cram PT-CA2271__N__Tumor__RNA__A29398-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__Tumor__RNA__A29398-Aligned.toTranscriptome.out.cram PT-CA2271__N__Tumor__RNA__A29398-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__Xenograft__RNA__A52018-Aligned.sortedByCoord.out.cram PT-CA2271__N__Xenograft__RNA__A52018-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CA2271__N__Xenograft__RNA__A52018-Aligned.toTranscriptome.out.cram PT-CA2271__N__Xenograft__RNA__A52018-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A29399-Aligned.sortedByCoord.out.cram PT-CM1209__R__CellLine__RNA__A29399-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A29399-Aligned.toTranscriptome.out.cram PT-CM1209__R__CellLine__RNA__A29399-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A37792-Aligned.sortedByCoord.out.cram PT-CM1209__R__CellLine__RNA__A37792-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A37792-Aligned.toTranscriptome.out.cram PT-CM1209__R__CellLine__RNA__A37792-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A37793-Aligned.sortedByCoord.out.cram PT-CM1209__R__CellLine__RNA__A37793-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__CellLine__RNA__A37793-Aligned.toTranscriptome.out.cram PT-CM1209__R__CellLine__RNA__A37793-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Tumor__RNA__A37794-Aligned.sortedByCoor.out.cram PT-CM1209__R__Tumor__RNA__A37794-Aligned.sortedByCoor.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Tumor__RNA__A37794-Aligned.toTranscriptome.out.cram PT-CM1209__R__Tumor__RNA__A37794-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Tumor__RNA__A37795-Aligned.sortedByCoord.out.cram PT-CM1209__R__Tumor__RNA__A37795-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Tumor__RNA__A37795-Aligned.toTranscriptome.out.cram PT-CM1209__R__Tumor__RNA__A37795-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Xenograft__RNA__A52019-Aligned.sortedByCoord.out.cram PT-CM1209__R__Xenograft__RNA__A52019-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM1209__R__Xenograft__RNA__A52019-Aligned.toTranscriptome.out.cram PT-CM1209__R__Xenograft__RNA__A52019-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM3220__N__CellLine__RNA__A29400-Aligned.sortedByCoord.out.cram PT-CM3220__N__CellLine__RNA__A29400-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM3220__N__CellLine__RNA__A29400-Aligned.toTranscriptome.out.cram PT-CM3220__N__CellLine__RNA__A29400-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM3220__N__Xenograft__RNA__A52011-Aligned.sortedByCoord.out.cram PT-CM3220__N__Xenograft__RNA__A52011-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-CM3220__N__Xenograft__RNA__A52011-Aligned.toTranscriptome.out.cram PT-CM3220__N__Xenograft__RNA__A52011-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DF5919__N__CellLine__RNA__A29401-Aligned.sortedByCoord.out.cram PT-DF5919__N__CellLine__RNA__A29401-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DF5919__N__CellLine__RNA__A29401-Aligned.toTranscriptome.out.cram PT-DF5919__N__CellLine__RNA__A29401-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DF5919__N__Tumor__RNA__A29402-Aligned.sortedByCoord.out.cram PT-DF5919__N__Tumor__RNA__A29402-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DF5919__N__Tumor__RNA__A29402-Aligned.toTranscriptome.out.cram PT-DF5919__N__Tumor__RNA__A29402-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DM9089__R__CellLine__RNA__A29403-Aligned.sortedByCoord.out.cram PT-DM9089__R__CellLine__RNA__A29403-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DM9089__R__CellLine__RNA__A29403-Aligned.toTranscriptome.out.cram PT-DM9089__R__CellLine__RNA__A29403-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__CellLine__RNA__A29404-Aligned.sortedByCoord.out.cram PT-DS9789__R__CellLine__RNA__A29404-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__CellLine__RNA__A29404-Aligned.toTranscriptome.out.cram PT-DS9789__R__CellLine__RNA__A29404-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__Tumor__RNA__A29405-Aligned.sortedByCoord.out.cram PT-DS9789__R__Tumor__RNA__A29405-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__Tumor__RNA__A29405-Aligned.toTranscriptome.out.cram PT-DS9789__R__Tumor__RNA__A29405-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__Xenograft__RNA__A50724-Aligned.sortedByCoord.out.cram PT-DS9789__R__Xenograft__RNA__A50724-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-DS9789__R__Xenograft__RNA__A50724-Aligned.toTranscriptome.out.cram PT-DS9789__R__Xenograft__RNA__A50724-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__CellLine__RNA__A29406-Aligned.sortedByCoord.out.cram PT-EV3071__N__CellLine__RNA__A29406-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__CellLine__RNA__A29406-Aligned.toTranscriptome.out.cram PT-EV3071__N__CellLine__RNA__A29406-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__Tumor__RNA__A29407-Aligned.sortedByCoord.out.cram PT-EV3071__N__Tumor__RNA__A29407-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__Tumor__RNA__A29407-Aligned.toTranscriptome.out.cram PT-EV3071__N__Tumor__RNA__A29407-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__Xenograft__RNA__A50722-Aligned.sortedByCoord.out.cram PT-EV3071__N__Xenograft__RNA__A50722-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-EV3071__N__Xenograft__RNA__A50722-Aligned.toTranscriptome.out.cram PT-EV3071__N__Xenograft__RNA__A50722-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-FB6711__N__CellLine__RNA__A29408-Aligned.sortedByCoord.out.cram PT-FB6711__N__CellLine__RNA__A29408-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-FB6711__N__CellLine__RNA__A29408-Aligned.toTranscriptome.out.cram PT-FB6711__N__CellLine__RNA__A29408-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-FB6711__N__Tumor__RNA__A29409-Aligned.sortedByCoord.out.cram PT-FB6711__N__Tumor__RNA__A29409-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-FB6711__N__Tumor__RNA__A29409-Aligned.toTranscriptome.out.cram PT-FB6711__N__Tumor__RNA__A29409-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9186__N__CellLine__RNA__A29410-Aligned.sortedByCoord.out.cram PT-GB9186__N__CellLine__RNA__A29410-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9186__N__CellLine__RNA__A29410-Aligned.toTranscriptome.out.cram PT-GB9186__N__CellLine__RNA__A29410-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9186__N__Tumor__RNA__A29411-Aligned.sortedByCoord.out.cram PT-GB9186__N__Tumor__RNA__A29411-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9186__N__Tumor__RNA__A29411-Aligned.toTranscriptome.out.cram PT-GB9186__N__Tumor__RNA__A29411-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9483__R__CellLine__RNA__A37796-Aligned.sortedByCoord.out.cram PT-GB9483__R__CellLine__RNA__A37796-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GB9483__R__CellLine__RNA__A37796-Aligned.toTranscriptome.out.cram PT-GB9483__R__CellLine__RNA__A37796-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GC1519__N__CellLine__RNA__A37797-Aligned.sortedByCoord.out.cram PT-GC1519__N__CellLine__RNA__A37797-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GC1519__N__CellLine__RNA__A37797-Aligned.toTranscriptome.out.cram PT-GC1519__N__CellLine__RNA__A37797-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GC1519__N__Tumor__RNA__A37798-Aligned.sortedByCoord.out.cram PT-GC1519__N__Tumor__RNA__A37798-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GC1519__N__Tumor__RNA__A37798-Aligned.toTranscriptome.out.cram PT-GC1519__N__Tumor__RNA__A37798-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GE7528__R__Tumor__RNA__A29413-Aligned.sortedByCoord.out.cram PT-GE7528__R__Tumor__RNA__A29413-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GE7528__R__Tumor__RNA__A29413-Aligned.toTranscriptome.out.cram PT-GE7528__R__Tumor__RNA__A29413-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GI2070__N__CellLine__RNA__A29414-Aligned.sortedByCoord.out.cram PT-GI2070__N__CellLine__RNA__A29414-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GI2070__N__CellLine__RNA__A29414-Aligned.toTranscriptome.out.cram PT-GI2070__N__CellLine__RNA__A29414-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GJ3716__N__CellLine__RNA__A37799-Aligned.sortedByCoord.out.cram PT-GJ3716__N__CellLine__RNA__A37799-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GJ3716__N__CellLine__RNA__A37799-Aligned.toTranscriptome.out.cram PT-GJ3716__N__CellLine__RNA__A37799-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GJ3716__N__Tumor__RNA__A37800-Aligned.sortedByCoord.out.cram PT-GJ3716__N__Tumor__RNA__A37800-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GJ3716__N__Tumor__RNA__A37800-Aligned.toTranscriptome.out.cram PT-GJ3716__N__Tumor__RNA__A37800-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GR2309__N__CellLine__RNA__A37801-Aligned.sortedByCoord.out.cram PT-GR2309__N__CellLine__RNA__A37801-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-GR2309__N__CellLine__RNA__A37801-Aligned.toTranscriptome.out.cram PT-GR2309__N__CellLine__RNA__A37801-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HN6692__NA__CellLine__RNA__A37802-Aligned.sortedByCoord.out.cram PT-HN6692__NA__CellLine__RNA__A37802-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HN6692__NA__CellLine__RNA__A37802-Aligned.toTranscriptome.out.cram PT-HN6692__NA__CellLine__RNA__A37802-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HO0394__N__CellLine__RNA__A29415-Aligned.sortedByCoord.out.cram PT-HO0394__N__CellLine__RNA__A29415-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HO0394__N__CellLine__RNA__A29415-Aligned.toTranscriptome.out.cram PT-HO0394__N__CellLine__RNA__A29415-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HO0394__N__Xenograft__RNA__A52022-Aligned.sortedByCoord.out.cram PT-HO0394__N__Xenograft__RNA__A52022-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-HO0394__N__Xenograft__RNA__A52022-Aligned.toTranscriptome.out.cram PT-HO0394__N__Xenograft__RNA__A52022-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JB1730__N__CellLine__RNA__A29417-Aligned.sortedByCoord.out.cram PT-JB1730__N__CellLine__RNA__A29417-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JB1730__N__CellLine__RNA__A29417-Aligned.toTranscriptome.out.cram PT-JB1730__N__CellLine__RNA__A29417-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JB1730__N__Tumor__RNA__A29418-Aligned.sortedByCoord.out.cram PT-JB1730__N__Tumor__RNA__A29418-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JB1730__N__Tumor__RNA__A29418-Aligned.toTranscriptome.out.cram PT-JB1730__N__Tumor__RNA__A29418-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__CellLine__RNA__A29419-Aligned.sortedByCoord.out.cram PT-JE6375__N__CellLine__RNA__A29419-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__CellLine__RNA__A29419-Aligned.toTranscriptome.out.cram PT-JE6375__N__CellLine__RNA__A29419-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__Tumor__RNA__A29420-Aligned.sortedByCoord.out.cram PT-JE6375__N__Tumor__RNA__A29420-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__Tumor__RNA__A29420-Aligned.toTranscriptome.out.cram PT-JE6375__N__Tumor__RNA__A29420-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__Xenograft__RNA__A50727-Aligned.sortedByCoord.out.cram PT-JE6375__N__Xenograft__RNA__A50727-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JE6375__N__Xenograft__RNA__A50727-Aligned.toTranscriptome.out.cram PT-JE6375__N__Xenograft__RNA__A50727-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JP2405__R__CellLine__RNA__A37803-Aligned.sortedByCoord.out.cram PT-JP2405__R__CellLine__RNA__A37803-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JP2405__R__CellLine__RNA__A37803-Aligned.toTranscriptome.out.cram PT-JP2405__R__CellLine__RNA__A37803-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JP2405__R__Tumor__RNA__A37804-Aligned.sortedByCoord.out.cram PT-JP2405__R__Tumor__RNA__A37804-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JP2405__R__Tumor__RNA__A37804-Aligned.toTranscriptome.out.cram PT-JP2405__R__Tumor__RNA__A37804-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JR9883__N__CellLine__RNA__A29421-Aligned.sortedByCoord.out.cram PT-JR9883__N__CellLine__RNA__A29421-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JR9883__N__CellLine__RNA__A29421-Aligned.toTranscriptome.out.cram PT-JR9883__N__CellLine__RNA__A29421-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JW6420__N__CellLine__RNA__A29422-Aligned.sortedByCoord.out.cram PT-JW6420__N__CellLine__RNA__A29422-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JW6420__N__CellLine__RNA__A29422-Aligned.toTranscriptome.out.cram PT-JW6420__N__CellLine__RNA__A29422-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JW6420__N__Tumor__RNA__A29423-Aligned.sortedByCoord.out.cram PT-JW6420__N__Tumor__RNA__A29423-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-JW6420__N__Tumor__RNA__A29423-Aligned.toTranscriptome.out.cram PT-JW6420__N__Tumor__RNA__A29423-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-KM5291__N__CellLine__RNA__A29424-Aligned.sortedByCoord.out.cram PT-KM5291__N__CellLine__RNA__A29424-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-KM5291__N__CellLine__RNA__A29424-Aligned.toTranscriptome.out.cram PT-KM5291__N__CellLine__RNA__A29424-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-KM5291__N__Tumor__RNA__A29425-Aligned.sortedByCoord.out.cram PT-KM5291__N__Tumor__RNA__A29425-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-KM5291__N__Tumor__RNA__A29425-Aligned.toTranscriptome.out.cram PT-KM5291__N__Tumor__RNA__A29425-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__CellLine__RNA__A29426-Aligned.sortedByCoord.out.cram PT-LC3356__N__CellLine__RNA__A29426-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__CellLine__RNA__A29426-Aligned.toTranscriptome.out.cram PT-LC3356__N__CellLine__RNA__A29426-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__Tumor__RNA__A29427-Aligned.sortedByCoord.out.cram PT-LC3356__N__Tumor__RNA__A29427-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__Tumor__RNA__A29427-Aligned.toTranscriptome.out.cram PT-LC3356__N__Tumor__RNA__A29427-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__Xenograft__RNA__A52023-Aligned.sortedByCoord.out.cram PT-LC3356__N__Xenograft__RNA__A52023-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC3356__N__Xenograft__RNA__A52023-Aligned.toTranscriptome.out.cram PT-LC3356__N__Xenograft__RNA__A52023-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC6372__N__CellLine__RNA__A37805-Aligned.sortedByCoord.out.cram PT-LC6372__N__CellLine__RNA__A37805-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LC6372__N__CellLine__RNA__A37805-Aligned.toTranscriptome.out.cram PT-LC6372__N__CellLine__RNA__A37805-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LR9369__R__CellLine__RNA__A29428-Aligned.sortedByCoord.out.cram PT-LR9369__R__CellLine__RNA__A29428-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LR9369__R__CellLine__RNA__A29428-Aligned.toTranscriptome.out.cram PT-LR9369__R__CellLine__RNA__A29428-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LR9369__R__Tumor__RNA__A29429-Aligned.sortedByCoord.out.cram PT-LR9369__R__Tumor__RNA__A29429-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LR9369__R__Tumor__RNA__A29429-Aligned.toTranscriptome.out.cram PT-LR9369__R__Tumor__RNA__A29429-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LS4891__N__CellLine__RNA__A29430-Aligned.sortedByCoord.out.cram PT-LS4891__N__CellLine__RNA__A29430-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LS4891__N__CellLine__RNA__A29430-Aligned.toTranscriptome.out.cram PT-LS4891__N__CellLine__RNA__A29430-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LS4891__N__Tumor__RNA__A29431-Aligned.sortedByCoord.out.out.cram PT-LS4891__N__Tumor__RNA__A29431-Aligned.sortedByCoord.out.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-LS4891__N__Tumor__RNA__A29431-Aligned.toTranscriptome.out.cram PT-LS4891__N__Tumor__RNA__A29431-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A29432-Aligned.sortedByCoord.out.cram PT-MB9777__R__CellLine__RNA__A29432-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A29432-Aligned.toTranscriptome..out.cram PT-MB9777__R__CellLine__RNA__A29432-Aligned.toTranscriptome..out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A37806-Aligned.sortedByCoord.out.cram PT-MB9777__R__CellLine__RNA__A37806-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A37806-Aligned.toTranscriptome.out.cram PT-MB9777__R__CellLine__RNA__A37806-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A37807-Aligned.sortedByCoord.out.cram PT-MB9777__R__CellLine__RNA__A37807-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__CellLine__RNA__A37807-Aligned.toTranscriptome.out.cram PT-MB9777__R__CellLine__RNA__A37807-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Tumor__RNA__A29433-Aligned.sortedByCoord.out.cram PT-MB9777__R__Tumor__RNA__A29433-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Tumor__RNA__A29433-Aligned.toTranscriptome.out.cram PT-MB9777__R__Tumor__RNA__A29433-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Tumor__RNA__A37808-Aligned.sortedByCoord.out.cram PT-MB9777__R__Tumor__RNA__A37808-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Tumor__RNA__A37808-Aligned.toTranscriptome.out.cram PT-MB9777__R__Tumor__RNA__A37808-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o T-MB9777__R__Tumor__RNA__A37809-Aligned.sortedByCoord..out.cram T-MB9777__R__Tumor__RNA__A37809-Aligned.sortedByCoord..out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Tumor__RNA__A37809-Aligned.toTranscriptome.out.cram PT-MB9777__R__Tumor__RNA__A37809-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Xenograft__RNA__A50723-Aligned.sortedByCoord.out.cram PT-MB9777__R__Xenograft__RNA__A50723-Aligned.sortedByCoord.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Xenograft__RNA__A50723-Aligned.toTranscriptome.out.cram PT-MB9777__R__Xenograft__RNA__A50723-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Xenograft__RNA__A50723-Aligned.toTranscriptome.out.cram PT-MB9777__R__Xenograft__RNA__A50723-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Xenograft__RNA__A50732-Aligned.sortedByCoord.out.cram PT-MB9777__R__Xenograft__RNA__A50732-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MB9777__R__Xenograft__RNA__A50732-Aligned.toTranscriptome.out.cram PT-MB9777__R__Xenograft__RNA__A50732-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MD9088__R__CellLine__RNA__A37810-Aligned.sortedByCoord.out.cram PT-MD9088__R__CellLine__RNA__A37810-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MD9088__R__CellLine__RNA__A37810-Aligned.toTranscriptome.out.cram PT-MD9088__R__CellLine__RNA__A37810-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MD9088__R__Tumor__RNA__A37811-Aligned.sortedByCoord.out.cram PT-MD9088__R__Tumor__RNA__A37811-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-MD9088__R__Tumor__RNA__A37811-Aligned.toTranscriptome.out.cram PT-MD9088__R__Tumor__RNA__A37811-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-ML9537__N__CellLine__RNA__A37812-Aligned.sortedByCoord.out.cram PT-ML9537__N__CellLine__RNA__A37812-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-ML9537__N__CellLine__RNA__A37812-Aligned.toTranscriptome.out.cram PT-ML9537__N__CellLine__RNA__A37812-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PD6881__R__CellLine__RNA__A29435-Aligned.sortedByCoord.out.cram PT-PD6881__R__CellLine__RNA__A29435-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PD6881__R__CellLine__RNA__A29435-Aligned.toTranscriptome.out.cram PT-PD6881__R__CellLine__RNA__A29435-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PD6881__R__Tumor__RNA__A29436-Aligned.sortedByCoord.out.cram PT-PD6881__R__Tumor__RNA__A29436-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PD6881__R__Tumor__RNA__A29436-Aligned.toTranscriptome.out.cram PT-PD6881__R__Tumor__RNA__A29436-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PR5617__N__CellLine__RNA__A29437-Aligned.sortedByCoord.out.cram PT-PR5617__N__CellLine__RNA__A29437-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PR5617__N__CellLine__RNA__A29437-Aligned.toTranscriptome.out.cram PT-PR5617__N__CellLine__RNA__A29437-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A29438-Aligned.sortedByCoord.out.cram PT-PV2594__N__CellLine__RNA__A29438-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A29438-Aligned.toTranscriptome.out.cram PT-PV2594__N__CellLine__RNA__A29438-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A52007-Aligned.sortedByCoord.out.cram PT-PV2594__N__CellLine__RNA__A52007-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A52007-Aligned.toTranscriptome.out.cram PT-PV2594__N__CellLine__RNA__A52007-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A52008-Aligned.sortedByCoord.out.cram PT-PV2594__N__CellLine__RNA__A52008-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__CellLine__RNA__A52008-Aligned.toTranscriptome.out.cram PT-PV2594__N__CellLine__RNA__A52008-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__Tumor__RNA__A29439-Aligned.sortedByCoord.out.cram PT-PV2594__N__Tumor__RNA__A29439-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__Tumor__RNA__A29439-Aligned.toTranscriptome.out.cram PT-PV2594__N__Tumor__RNA__A29439-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__Xenograft__RNA__A50726-Aligned.sortedByCoord.out.cram PT-PV2594__N__Xenograft__RNA__A50726-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-PV2594__N__Xenograft__RNA__A50726-Aligned.toTranscriptome.out.cram PT-PV2594__N__Xenograft__RNA__A50726-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RD1291__N__CellLine__RNA__A29416-Aligned.sortedByCoord.out.cram PT-RD1291__N__CellLine__RNA__A29416-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RD1291__N__CellLine__RNA__A29416-Aligned.toTranscriptome.out.cram PT-RD1291__N__CellLine__RNA__A29416-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RD1291__N__Tumor__RNA__A37813-Aligned.sortedByCoord.out.cram PT-RD1291__N__Tumor__RNA__A37813-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RD1291__N__Tumor__RNA__A37813-Aligned.toTranscriptome.out.cram PT-RD1291__N__Tumor__RNA__A37813-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL5404__N__CellLine__RNA__A29440-Aligned.sortedByCoord.out.cram PT-RL5404__N__CellLine__RNA__A29440-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL5404__N__CellLine__RNA__A29440-Aligned.toTranscriptome.out.cram PT-RL5404__N__CellLine__RNA__A29440-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL5404__N__Tumor__RNA__A29441-Aligned.sortedByCoord.out.cram PT-RL5404__N__Tumor__RNA__A29441-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL5404__N__Tumor__RNA__A29441-Aligned.toTranscriptome.out.cram PT-RL5404__N__Tumor__RNA__A29441-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL7940__N__CellLine__RNA__A29442-Aligned.sortedByCoord.out.cram PT-RL7940__N__CellLine__RNA__A29442-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL7940__N__CellLine__RNA__A29442-Aligned.toTranscriptome.out.cram PT-RL7940__N__CellLine__RNA__A29442-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL7940__N__Tumor__RNA__A29443-Aligned.sortedByCoord.out.cram PT-RL7940__N__Tumor__RNA__A29443-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RL7940__N__Tumor__RNA__A29443-Aligned.toTranscriptome.out.cram PT-RL7940__N__Tumor__RNA__A29443-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RV2286__R__CellLine__RNA__A37814-Aligned.sortedByCoord.out.cram PT-RV2286__R__CellLine__RNA__A37814-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RV2286__R__CellLine__RNA__A37814-Aligned.toTranscriptome.out.cram PT-RV2286__R__CellLine__RNA__A37814-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A29444-Aligned.sortedByCoord.out.cram PT-RW9277__NA__CellLine__RNA__A29444-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A29444-Aligned.toTranscriptome.out.cram PT-RW9277__NA__CellLine__RNA__A29444-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A52009-Aligned.sortedByCoord.out.cram PT-RW9277__NA__CellLine__RNA__A52009-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A52009-Aligned.toTranscriptome.out.cram PT-RW9277__NA__CellLine__RNA__A52009-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A52010-Aligned.sortedByCoord.out.cram PT-RW9277__NA__CellLine__RNA__A52010-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__CellLine__RNA__A52010-Aligned.toTranscriptome.out.cram PT-RW9277__NA__CellLine__RNA__A52010-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Tumor__RNA__A30194-Aligned.sortedByCoord.out.cram PT-RW9277__NA__Tumor__RNA__A30194-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Tumor__RNA__A30194-Aligned.toTranscriptome.out.cram PT-RW9277__NA__Tumor__RNA__A30194-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Tumor__RNA__A37815-Aligned.sortedByCoord.out.cram PT-RW9277__NA__Tumor__RNA__A37815-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Tumor__RNA__A37815-Aligned.toTranscriptome.out.cram PT-RW9277__NA__Tumor__RNA__A37815-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Xenograft__RNA__A50730-Aligned.sortedByCoord.out.cram PT-RW9277__NA__Xenograft__RNA__A50730-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-RW9277__NA__Xenograft__RNA__A50730-Aligned.toTranscriptome.out.cram PT-RW9277__NA__Xenograft__RNA__A50730-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__CellLine__RNA__A29446-Aligned.sortedByCoord.out.cram PT-SJ5453__N__CellLine__RNA__A29446-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__CellLine__RNA__A29446-Aligned.toTranscriptome.out.cram PT-SJ5453__N__CellLine__RNA__A29446-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__Tumor__RNA__A29447-Aligned.sortedByCoord.out.cram PT-SJ5453__N__Tumor__RNA__A29447-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__Tumor__RNA__A29447-Aligned.toTranscriptome.out.cram PT-SJ5453__N__Tumor__RNA__A29447-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__Xenograft__RNA__A50725-Aligned.sortedByCoord.out.cram PT-SJ5453__N__Xenograft__RNA__A50725-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SJ5453__N__Xenograft__RNA__A50725-Aligned.toTranscriptome.out.cram PT-SJ5453__N__Xenograft__RNA__A50725-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SK0976__N__CellLine__RNA__A29448-Aligned.sortedByCoord.out.cram PT-SK0976__N__CellLine__RNA__A29448-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SK0976__N__CellLine__RNA__A29448-Aligned.toTranscriptome.out.cram PT-SK0976__N__CellLine__RNA__A29448-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SK0976__N__Tumor__RNA__A37816-Aligned.sortedByCoord.out.cram PT-SK0976__N__Tumor__RNA__A37816-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SK0976__N__Tumor__RNA__A37816-Aligned.toTranscriptome.out.cram PT-SK0976__N__Tumor__RNA__A37816-Aligned.toTranscriptome.out.bam


samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SO0258__N__CellLine__RNA__A29449-Aligned.sortedByCoord..out.cram PT-SO0258__N__CellLine__RNA__A29449-Aligned.sortedByCoord..out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SO0258__N__CellLine__RNA__A29449-Aligned.toTranscriptome.out.cram PT-SO0258__N__CellLine__RNA__A29449-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SO0258__N__Tumor__RNA__A29450-Aligned.sortedByCoord.out.cram PT-SO0258__N__Tumor__RNA__A29450-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SO0258__N__Tumor__RNA__A29450-Aligned.toTranscriptome.out.cram PT-SO0258__N__Tumor__RNA__A29450-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SS3647__N__CellLine__RNA__A29451-Aligned.sortedByCoord.out.cram PT-SS3647__N__CellLine__RNA__A29451-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SS3647__N__CellLine__RNA__A29451-Aligned.toTranscriptome.out.cram PT-SS3647__N__CellLine__RNA__A29451-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SS3647__N__Tumor__RNA__A29452-Aligned.sortedByCoord.out.cram PT-SS3647__N__Tumor__RNA__A29452-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-SS3647__N__Tumor__RNA__A29452-Aligned.toTranscriptome.out.cram PT-SS3647__N__Tumor__RNA__A29452-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-TM5196__R__CellLine__RNA__A37817-Aligned.sortedByCoord.out.cram PT-TM5196__R__CellLine__RNA__A37817-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-TM5196__R__CellLine__RNA__A37817-Aligned.toTranscriptome.out.cram PT-TM5196__R__CellLine__RNA__A37817-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-TM5196__R__Tumor__RNA__A37818-Aligned.sortedByCoord.out.cram PT-TM5196__R__Tumor__RNA__A37818-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-TM5196__R__Tumor__RNA__A37818-Aligned.toTranscriptome.out.cram PT-TM5196__R__Tumor__RNA__A37818-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-VO7089__N__CellLine__RNA__A37819-Aligned.sortedByCoord.out.cram PT-VO7089__N__CellLine__RNA__A37819-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-VO7089__N__CellLine__RNA__A37819-Aligned.toTranscriptome.out.cram PT-VO7089__N__CellLine__RNA__A37819-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-VO7089__N__Tumor__RNA__A29454-Aligned.sortedByCoord.out.cram PT-VO7089__N__Tumor__RNA__A29454-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-VO7089__N__Tumor__RNA__A29454-Aligned.toTranscriptome.out.cram PT-VO7089__N__Tumor__RNA__A29454-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WP9124__R__CellLine__RNA__A37820-Aligned.sortedByCoord.out.cram PT-WP9124__R__CellLine__RNA__A37820-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WP9124__R__CellLine__RNA__A37820-Aligned.toTranscriptome.out.cram PT-WP9124__R__CellLine__RNA__A37820-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WP9124__R__Tumor__RNA__A37821-Aligned.sortedByCoord.out.cram PT-WP9124__R__Tumor__RNA__A37821-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WP9124__R__Tumor__RNA__A37821-Aligned.toTranscriptome.out.cram PT-WP9124__R__Tumor__RNA__A37821-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WR7927__N__CellLine__RNA__A37822-Aligned.sortedByCoord.out.cram PT-WR7927__N__CellLine__RNA__A37822-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WR7927__N__CellLine__RNA__A37822-Aligned.toTranscriptome.out.cram PT-WR7927__N__CellLine__RNA__A37822-Aligned.toTranscriptome.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WT4796__N__Tumor__RNA__A29456-Aligned.sortedByCoord.out.cram PT-WT4796__N__Tumor__RNA__A29456-Aligned.sortedByCoord.out.bam

samtools view -T /work/morrissy_lab/TFRI_2018_data/bamrnacopy/map-hg38-rel96/HS.GRCh38.rel96.fasta/Homo_sapiens.GRCh38.dna.chromosome.all.fa -C -o PT-WT4796__N__Tumor__RNA__A29456-Aligned.toTranscriptome.out.cram PT-WT4796__N__Tumor__RNA__A29456-Aligned.toTranscriptome.out.bam

