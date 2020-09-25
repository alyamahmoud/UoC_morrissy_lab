#BiocManager::install ('tidyverse')
BiocManager::install('NMF')
#library(tidyverse)
library(NMF)
df <- read.delim('/work/morrissy_lab/alyaa.mohamed/OpenPBTA/OpenPBTA-analysis/overlapping_regulons_activity.txt', stringsAsFactors = FALSE, row.names = 1)
# convert matrix to non-negative values
dat <- nneg(as.matrix(df), method = 'min')
# NMF
# run at multiple ranks and choose the one at which the cophenetic value starts to decrease
# Estimation of the rank: Consensus matrices computed from 10 runs for each value of r
seed= 123456
nr=2000
sample <- 'overlapping_regulons'
res <- nmf (dat, seq(2,10), nrun = nr, seed = seed, method='brunet', .options='tvp50')
save(res, file=paste(sample, genes, nr, "NMFresult.Robj", sep="_"))

# permute data and recalculate (to estimate if overfitting is happening):
V.random <- randomize(dat)
estim.r.random <- nmf(V.random, seq(2,20), method='brunet', nrun = nr, seed = seed, .options='tvp50')
save(estim.r.random, file=paste(sample, genes, nr, "NMFresult_randomized.Robj", sep="_"))
# plot summary information
pdf(paste(sample, genes, nr, "NMFstats_randomized.pdf", sep="_"))
plot(res, estim.r.random)
dev.off()
pdf(paste(sample, genes, nr, "NMFstats.pdf", sep="_"), width=10, height=7)
plot(res)
dev.off()
pdf(paste(sample, genes, nr, "NMFconsensusmaps.pdf", sep="_"), width=28, height=24, onefile=TRUE)
consensusmap(res, Rowv=TRUE, Colv=TRUE)
dev.off()

