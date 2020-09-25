suppressPackageStartupMessages({
  library (data.table)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
library (Matrix)
})
BiocManager::install('tidyverse')
library (tidyverse)

# read in RSEM FPKM collapsed stranded 
expression_file <- "/work/morrissy_lab/alyaa.mohamed/OpenPBTA/OpenPBTA-analysis/data/expression_data_matched_proteomics.Rds" 

# Read in expression data
expression_data <- readr::read_rds(expression_file) %>%
  as.data.frame()	
dim(expression_data)
#expression_data <- column_to_rownames(expression_data, var = 'gene_id')

genes_to_keep <- rowSums(expression_data) >= 100
expression_data <- expression_data[genes_to_keep, ]

# metadata
histologies_file <- "/work/morrissy_lab/alyaa.mohamed/OpenPBTA/OpenPBTA-analysis/data/pbta-histologies.tsv"
histologies <- read_tsv(histologies_file) 
histologies <- histologies[histologies$Kids_First_Biospecimen_ID %in% colnames(expression_data), ]
histologies <- histologies %>%
  dplyr::select(Kids_First_Biospecimen_ID, short_histology, broad_histology, broad_composition) %>%
  column_to_rownames(var = 'Kids_First_Biospecimen_ID')

# SCENIC
scenic_in <- expression_data
dir.create("/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int")
saveRDS(scenic_in, file = "/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/scenic_in.Rds")
# set scenicoptions
org <- "hgnc" # or mgi, or dmel
dbDir <- "/work/morrissy_lab/alyaa.mohamed/OpenPBTA/cisTarget_databases"
myDatasetTitle <- "SCENIC on CBTTC" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20)

cellInfo <- histologies
saveRDS(cellInfo, file="/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/cellInfo.Rds")
saveRDS(scenicOptions, file ="/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/scenicOptions.Rds")

scenicOptions@inputDatasetInfo$cellInfo <- "/work/morrissylab/alyaa.mohamed/OpenPBTA/int/cellInfo.Rds"
saveRDS(scenicOptions, file ="/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/scenicOptions.Rds")


### Co-expression network
#genesKept <- geneFiltering(exprMat, scenicOptions)
#exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- scenic_in
exprMat <- scenic_in
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
exprMat_filtered_log <- as.matrix(exprMat_filtered_log)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
logMat <- log2(exprMat+1)
logMat <- as.matrix(logMat)
scenicOptions <- readRDS("/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 20
scenicOptions@settings$seed <- 123
scenicOptions@inputDatasetInfo$cellInfo <- "/work/morrissy_lab/alyaa.mohamed/OpenPBTA/int/cellInfo.Rds"

# For toy run
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
#runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, logMat)

