library (tidyverse)
BiocManager::install('ComplexHeatmap')
library (ComplexHeatmap)#
library (DESeq2)

# brain span expression data
dat <- NULL; 
dat <- read_tsv('scratch/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt')

# TO DO use the raw counts and normalize them in DESeq rather than RPKMs
# mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count
# colnames refer temporal.spatial IDs
# convert rownames to gene symbols and summarise by mean for input to SCENIC
df <- dat; 
df$Geneid <- sapply (strsplit(df$Geneid, '\\|' ), function(x) paste(x[2]) )
df <- df %>%
  group_by(Geneid) %>%
  summarise_all(mean) %>%
  column_to_rownames(var = 'Geneid')
saveRDS(df, file = 'scratch/summarised_pyschENCODE.rds')
dat_psychENCODE <- df; 

# pyschENCODE annotation
annotation_psychENCODE <- read_tsv('scratch/mRNA-seq_Sample_metadata.txt')

# brain span regulons
auc_psychENCODE = readRDS('scratch/3.4_regulonAUC_psychENCODE.Rds')
regulons_psychENCODE = readRDS('scratch/2.6_regulons_asGeneSet_psychENCODE.Rds')

# heatmap
dat = NULL; 
dat = data.frame(assay(auc_psychENCODE))
dat = na.omit(t(scale(t(as.matrix(dat)), center = T, scale=T)))
colnames(dat) <- sapply(strsplit(colnames(dat), split='\\.'), function(x) x[[1]])
auc_psychENCODE <- dat; 


pdf ('scratch/plots/psychENCODE_regulons_genes_heatmap.pdf', width = 30, height = 18)
annotation_file <- annotation_psychENCODE %>%
  dplyr::select('Braincode', 'Age') %>%
  column_to_rownames(var= 'Braincode')
pheatmap::pheatmap(auc_psychENCODE, fontsize_row=3, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA, annotation = annotation_file, clustering_method = "ward.D2", fontsize_col = 3, show_rownames = FALSE)
dev.off()

# CBTTC expression data
expression_file <- 'scratch/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds'
expression_data <- readRDS(expression_file) %>%
	as.data.frame()
dat_cbttc <- expression_data

# CBTTC annotation
histologies <- read_tsv('scratch/pbta-histologies.tsv') 
histologies <- histologies[histologies$Kids_First_Biospecimen_ID %in% colnames(expression_data), ]

# CBTTC regulons
auc_cbttc <- readRDS('scratch/3.4_regulonAUC_cbttc.Rds')
regulons_cbttc <- readRDS('scratch/2.6_regulons_asGeneSet_cbttc.Rds')

# heatmap
dat = NULL; 
dat = data.frame(assay(auc_cbttc))
dat = na.omit(t(scale(t(as.matrix(dat)), center = T, scale=T)))
auc_cbttc <- dat; 

pdf ('scratch/plots/regulons_genes_full_dataset_heatmap.pdf', width = 12, height = 10)
histologies <- histologies %>%
  filter(short_histology == c('LGAT', 'HGAT', 'Medulloblastoma', 'Ependymoma')) %>%
  select (short_histology, Kids_First_Biospecimen_ID)
histologies <- histologies [!duplicated(histologies$Kids_First_Biospecimen_ID), ]
annotation_genes <- histologies %>%
  column_to_rownames(var = 'Kids_First_Biospecimen_ID')
dat <- dat[, colnames(dat) %in% histologies$Kids_First_Biospecimen_ID]
pheatmap::pheatmap(dat, fontsize_row=4, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA, annotation = annotation_genes, clustering_method = "ward.D2", fontsize_col = 3, show_rownames = FALSE)
dev.off()

# analysis of overlapping regulons
dat_tumor <- read_tsv('scratch/auc_cbttc_formatted.txt')
dat_normal <- read_tsv('scratch/auc_psychENCODE_formatted.txt')
# check overlapping regulons 
length (intersect(dat_tumor$regulon_name, dat_normal$regulon_name))
# remove _extended 
normal_regulons <- sapply(strsplit(dat_normal$regulon_name, split = '_'), function (x)x[[1]]); normal_regulons = unique(normal_regulons)
tumor_regulons <- sapply(strsplit(dat_tumor$regulon_name, split = '_'), function(x)x[[1]])

# check the targets for the overlapping and the private regulons
overlapping <- intersect(normal_regulons, tumor_regulons)
dat_tumor_overlapping <- dat_tumor[dat_tumor$regulon_name %in% overlapping, ]; colnames(dat_tumor_overlapping)[1:3] <- c('regulon_tumor', 'regulon_name', 'regulon_size_tumor')
dat_normal_overlapping <- dat_normal[dat_normal$regulon_name %in% overlapping, ]; 
colnames(dat_normal_overlapping)[1:3] <- c('regulon_normal', 'regulon_name', 'regulon_size_normal')
dat <- NULL; dat <- inner_join(dat_normal_overlapping, dat_tumor_overlapping, by = 'regulon_name')
df <- dat %>%
  select(-regulon_normal, -regulon_size_normal, -regulon_tumor, -regulon_size_tumor) %>%
  column_to_rownames(var = 'regulon_name')

# heatmap
pdf ('scratch/plots/regulons_CBTTC_vs_BrainSpan_overlapping_regulons_heatmap.pdf', width = 6, height = 3)
annotation <- data.frame(colnames(df))
annotation$dtype <- c(rep ('normal_developing_brain', '607'), rep ('pediatric_brain_tumor', 970))
annotation <- annotation %>%
  column_to_rownames(var = 'colnames.df.')
pheatmap::pheatmap(as.matrix(df), fontsize_row=2, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA, annotation = annotation, clustering_method = "ward.D2", fontsize_col = 2, show_rownames = FALSE, show_colnames = FALSE)
dev.off()

# which tumor types cluster with which brain development stage
# add another annotation level to indicate tumor type and brain development stage
pdf ('scratch/plots/regulons_CBTTC_vs_BrainSpan_overlapping_regulons_heatmap_annotated.pdf', width = 10, height = 14)
normal <- sapply(strsplit(rownames(annotation)[1:607], split = '_'), function (x) x[[1]])
normal <- data.frame(normal); colnames(normal) <- c('Braincode')
normal_annotation <- left_join(normal, data.frame(annotation_psychENCODE), by = 'Braincode') %>%
  select('Braincode', 'Age') 

tumor <- data.frame(rownames(annotation)[608:1577])
colnames(tumor) <- c('Kids_First_Biospecimen_ID')
tumor_annotation <- left_join(tumor, histologies, by = 'Kids_First_Biospecimen_ID') %>%
  select('Kids_First_Biospecimen_ID', 'short_histology') 
  
pheatmap::pheatmap(as.matrix(df), fontsize_row=4, color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),treeheight_row=10, treeheight_col=10, border_color=NA, annotation = annotation, clustering_method = "ward.D2", fontsize_col = 2, show_rownames = TRUE, show_colnames = FALSE)
dev.off()

# complexHeatmap
# k means clustering
library (ComplexHeatmap)
pdf ('scratch/plots/cbttc_brain_span_clustering_full_annotation.pdf', width = 12, height = 10)
Heatmap(df, row_km = 3, column_km = 2, row_km_repeats = 100, column_km_repeats = 100, show_parent_dend_line = FALSE)
# split heatmap by the original dendrogram
dend = hclust(dist(df), method = "ward.D2")
#dend = color_branches(dend, k = 3)
dend_r <- hclust(dist(t(df)), method = "ward.D2")
col <- list(dtype = c("normal_developing_brain" =  "red", "pediatric_brain_tumor" = "green"))
ha <- HeatmapAnnotation(df = annotation, col = col)
Heatmap(df, name = "mat",  cluster_rows = dend, row_split = 3, column_split = 6, cluster_columns = dend_r, border = TRUE, top_annotation = ha)
# append purity to the heatmap
# read tumor and nomral samples purity
out_brain_span <- readRDS('scratch/out_brain_span.Rds')
out_cbttc <- readRDS('scratch/out_cbttc.Rds')
out_brain_span_tref <- readRDS('scratch/out_brain_span_TRef.Rds')
out_cbttc_tref <- readRDS('scratch/out_cbttc_TRef.Rds')
tumor_purity <- data.frame(out_cbttc$cellFractions) %>%
 dplyr:: select('otherCells') %>%
  rownames_to_column(var = 'sample_ID')
colnames(tumor_purity)[2] <- c('fraction_cancer_cells')
normal_purity <- data.frame(out_brain_span$cellFractions) %>%
  dplyr::select('otherCells') %>%
  rownames_to_column(var = 'sample_ID') 
colnames(normal_purity)[2] <- c('fraction_cancer_cells')  
purity_norm_brain <- rbind(normal_purity, tumor_purity)
tumor_purity <- data.frame(out_cbttc_tref$cellFractions) %>%
 dplyr:: select('otherCells') %>%
  rownames_to_column(var = 'sample_ID')
colnames(tumor_purity)[2] <- c('fraction_cancer_cells')
normal_purity <- data.frame(out_brain_span_tref$cellFractions) %>%
  dplyr::select('otherCells') %>%
  rownames_to_column(var = 'sample_ID') 
colnames(normal_purity)[2] <- c('fraction_cancer_cells')  
purity_tref <- rbind (normal_purity, tumor_purity)
purity_tref$fraction_cancer_cells <- 1-(purity_tref$fraction_cancer_cells)
purity_norm_brain$fraction_cancer_cells <- 1-(purity_norm_brain$fraction_cancer_cells)
colnames(purity_norm_brain)[2] <- c('purity_norm_brain')
colnames(purity_tref)[2] <- c('purity_tref')
purity <- inner_join(purity_norm_brain, purity_tref, by = 'sample_ID')
purity$sample_ID <- sapply(strsplit(purity$sample_ID, split = '\\.'), function (x)x[[1]])
purity$non_immune_non_brain <- 1-(purity$purity_norm_brain - purity$purity_tref)
# add age group, region and tumor type to the annotation
annotation <- read_tsv('scratch/annotation_for_complex_heatmap.txt')
library(circlize)
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
column_ha <- HeatmapAnnotation(sample_type = annotation$dtype, age   = annotation$age_grouped, tumor = annotation$tumor_type_grouped, purity = anno_barplot(purity$non_immune_non_brain), col = list(sample_type = c('normal_developing_brain' = 'red', 'pediatric_brain_tumor' = 'blue'), age = c('adult' = 'aquamarine', 'infant' = 'bisque', 'kid' = 'blueviolet', 'pediatric_brain_tumor' = 'blue', 'prenatal_12_19' = 'darkgreen', 'youth' = 'darkorange', 'prenatal_35_37' = 'springgreen', 'prenatal_8_9' = 'violetred', 'prenatal_21_22' = 'yellow'), tumor = c('ATRT' = 'lightyellow', 'Craniopharyngioma' = 'darkslategray1', 'Ependymoma' = 'goldenrod2', 'ETMR' = 'deeppink', 'HGAT' = 'lawngreen', 'LGAT' = 'yellow', 'Medulloblastoma' = 'seagreen', 'Meningioma' = 'orchid', 'normal_developing_brain' = 'red', 'Oligodendroglioma' = 'rosybrown1', 'Other' = 'black')))            
Heatmap(df, top_annotation = column_ha, cluster_rows = dend, row_split = 6, column_split = 6, cluster_columns = dend_r, border = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 5))
dev.off()
