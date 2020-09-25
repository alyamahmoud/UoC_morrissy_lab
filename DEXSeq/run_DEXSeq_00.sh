#! /usr/bin/env bash

#BSUB -J htseq_00
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 48:00
#BSUB -o htseq_%J.out
#BSUB -e htseq_%J.err

cd /tiered/smorrissy/TFRI_2018_data/bamrnacopy/map-hg38-rel96/
python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-AB0029__N__Tumor__RNA__A29388-Aligned.sortedByCoord.out.bam exon_counts/A29388

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-KM5291__N__CellLine__RNA__A29424-Aligned.sortedByCoord.out.bam exon_counts/A29424

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-AB6372__R__CellLine__RNA__A29389-Aligned.sortedByCoord.out.bam exon_counts/A29389

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-KM5291__N__Tumor__RNA__A29425-Aligned.sortedByCoord.out.bam exon_counts/A29425

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-AB6372__R__Tumor__RNA__A37786-Aligned.sortedByCoord.out.bam exon_counts/A37786

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-LC3356__N__CellLine__RNA__A29426-Aligned.sortedByCoord.out.bam exon_counts/A29426

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-AH1410__N__CellLine__RNA__A29390-Aligned.sortedByCoord.out.bam exon_counts/A29390

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-LC3356__N__Tumor__RNA__A29427-Aligned.sortedByCoord.out.bam exon_counts/A29427

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-AK7565__N__CellLine__RNA__A37787-Aligned.sortedByCoord.out.bam exon_counts/A37787

python /tiered/smorrissy/home/amohamed/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos GRCh38_rel96_spiked_DEXSeq_chr.gff bams/PT-LC3356__N__Xenograft__RNA__A52023-Aligned.sortedByCoord.out.bam exon_counts/A52023

