Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
library (DEXSeq)
library (tidyverse)
library (dplyr)
library('DEXSeqAlt')

# set params

# read data
filenames=list.files(path = 'exon_counts/', pattern="SM*")
countFiles = filenames
countFiles <- countFiles[-1]
# sampleTable
# sampleTable = read_tsv('bam2ptID2metadata.txt')
sampleTable <- data.frame (rownames=countFiles, condition = c(rep ('SMP18', 6), rep ('SMP8', 12)), clinical_status = c(rep ('P1', 2), rep ('R1', 4), rep ('P1', 2), rep ('R1', 4), rep ('R2', 6))) 


# flattenedFile 
flattenedFile = 'gencode.v32.annotation.DEXSeq.chr.gff'

# design
#formulaFullModel <- ~ sample + exon + progress:exon + condition:exon
#formulaReducedModel <- ~ sample + exon + condition:exon
# formulaReducedModel <- ~ sample + exon + type:exon
formulaReducedModel <- ~ sample + exon 
formulaFullModel <- ~ sample + exon + condition:exon


# create DEXSeqDataSet object
setwd('exon_counts/')
dxd <- DEXSeqDataSetFromHTSeq(countFiles,  sampleData=sampleTable, design=  ~ sample + exon + condition:exon, flattenedfile = flattenedFile )

# remove _ambiguous, _ambiguous_readpair_position, _empty, _lowaqual , _notaligned 
remove <- c('_ambiguous', '_ambiguous_readpair_position', '_empty', '_lowaqual', '_notaligned')
rownames(dxd) <- rownames(dxd)[!rownames(dxd) %in% remove]


# normalisation
dxd <- estimateSizeFactors(dxd)


# dxd <- estimateDispersions(dxd, formula=formulaFullModel)
dxd <- estimateDispersionsAlt(dxd, verbose=TRUE)
saveRDS(dxd, file = '../output/dxd_BTIC_DEXSeqAlt.Rds')

# testing for differential exon usage
dxd <- testForDEUAlt(dxd, reducedModel=formulaReducedModel, fullModel=formulaFullModel)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition")
saveRDS(dxd, file = '../output/dxdtestForDEUAlt.Rds')
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
dxr1 <- DEXSeqResultsAlt(dxd)
# save R object
saveRDS(dxr1, file= '../output/dxr1_BTIC_DEXSeqAlt.Rds')
dxr1_tab <- cbind( dxr, mcols(dxd)[,which(elementMetadata(mcols(dxd))$type == "DEXSeq results")] )
saveRDS(dxr1_tab, file = '../output/dxr1_tab_BTIC_DEXSeqAlt.Rds')
# how many exonic regions are signficant with false discovery rate < 10%
table(dxr1$padj < 0.1)

# how many genes are affected
table(tapply(dxr1$padj < 0.1, dxr1$groupID, any))

pdf ('../plots/DEXSeq_BTIC_DEXSeqAlt.pdf')
plotDispEsts(dxd)
# MA plot
plotMA(dxr1, cex = 0.8)
dev.off()

# save R object
saveRDS(dxr1, file= '../output/dxr1_BTIC_DEXSeqAlt.Rds')
