library(NETI2)
expression_data <- readRDS('data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds')
scenic_selected_genes <- readRDS ('/work/morrissy_lab/alyaa.mohamed/OpenPBTA/OpenPBTA-analysis/selected_genes_NETI2_input.Rds')
expression_data <- expression_data[rownames(expression_data) %in% scenic_selected_genes, ]

histologies <- read.delim ('data/pbta-histologies.tsv')
keep <- c('Kids_First_Biospecimen_ID', 'short_histology', 'broad_histology', 'broad_composition', 'Kids_First_Participant_ID')
histologies <- histologies[, keep]
keep <- c('HGAT', 'LGAT', 'Ependymoma', 'Medulloblastoma')
histologies <- histologies[histologies$short_histology %in% keep, ]

# split data by tumor types
HGAT <- histologies$Kids_First_Biospecimen_ID[histologies$short_histology == 'HGAT' ]
LGAT <- histologies$Kids_First_Biospecimen_ID[histologies$short_histology == 'LGAT' ]
Medulloblastoma <- histologies$Kids_First_Biospecimen_ID[histologies$short_histology == 'Medulloblastoma' ]
Ependymoma <- histologies$Kids_First_Biospecimen_ID[histologies$short_histology == 'Ependymoma' ]

HGAT_dat <- expression_data[, colnames(expression_data) %in% HGAT]
LGAT_dat <- expression_data[, colnames(expression_data) %in% LGAT]
Medulloblastoma_dat <- expression_data[, colnames(expression_data) %in% Medulloblastoma]
Ependymoma_dat <- expression_data[, colnames(expression_data) %in% Ependymoma]

HGAT_dat <- t(HGAT_dat)
LGAT_dat <- t(LGAT_dat)
Medulloblastoma_dat <- t(Medulloblastoma_dat)
Ependymoma_dat <- t(Ependymoma_dat)

cbttc_dat <- as.matrix(list (HGAT_dat, LGAT_dat, Medulloblastoma_dat, Ependymoma_dat))

purity <- readRDS('out_cbttc.Rds')
purity <- data.frame(purity$cellFractions)
purity_HGAT <- purity$otherCells[rownames(purity) %in% rownames(HGAT_dat)]
purity_LGAT <- purity$otherCells[rownames(purity) %in% rownames(LGAT_dat)]
purity_Medulloblastoma <- purity$otherCells[rownames(purity) %in% rownames(Medulloblastoma_dat)]
purity_Ependymoma <- purity$otherCells[rownames(purity) %in% rownames(Ependymoma_dat)]

cbttc_purity <- as.matrix(list(purity_HGAT, purity_LGAT, purity_Medulloblastoma, purity_Ependymoma))


# run NETI2
out <- NETI2(cbttc_dat,cbttc_purity, lambda = 0.6, tau = 0.4,delta = 0.2)
saveRDS(out, file = 'data/NETI2_out_cbttc.Rds')
