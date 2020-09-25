# compute incidence and percentage relative to AUC activity scores
incidMat <- readr::read_rds('/work/morrissy_lab/alyaa.mohamed/psychENCODE/int/2.6_regulons_asIncidMat.Rds')
Reg <- incidMat   #TF x genes #incidMat
R_comb <- combn(rownames(Reg), m=2)  #Create pairwise combinations

# Create required dataframe 
Reg1 = Reg2  = nR1_genes = nR2_genes = nR1_R2_genes = R1_R2_genes = NULL;
 
for (comb in 1:ncol(R_comb)) {
  #two_genes <- print(as.character(gene_comb2[,combo]))
 
  two_R <- as.character(R_comb[,comb])
  R1 <- data.frame(Reg[which(rownames(Reg) == two_R[1]),])
  R2 <- data.frame(Reg[which(rownames(Reg) == two_R[2]),])
  R1_2 <- cbind(R1, R2)
  colnames(R1_2) <- c("R1", "R2")
  nR1_2 <- length(which(R1_2[,1] > 0 & R1_2[,2] > 0))
  R1R2_genes <- rownames(R1_2[which(R1_2[,1] > 0 & R1_2[,2] > 0),])
  R1R2_genes <- paste(R1R2_genes,collapse=",")
  #common_Regs <- paste(cat(rownames(samples1_2[which(samples1_2[,1] > 0 & samples1_2[,2] > 0),]), sep = “,”))
  Reg1 = c(Reg1, two_R[1])
  Reg2 = c(Reg2, two_R[2])
  nR1_R2_genes = c(nR1_R2_genes, nR1_2)
  nR1_genes = c(nR1_genes, sum(R1_2[,1]))
  nR2_genes = c(nR2_genes, sum(R1_2[,2]))
 R1_R2_genes = c(R1_R2_genes, R1R2_genes[1])
 }
 Regulon_comb <- data.frame(cbind(Reg1 ,Reg2 ,nR1_R2_genes, nR1_genes , nR2_genes ,R1_R2_genes))
 Regulon_comb$R1_R2_genes <- as.character(Regulon_comb$R1_R2_genes)
 Regulon_comb$R1_R2_genes[Regulon_comb$R1_R2_genes == ""] <- "NA"
 saveRDS(Regulon_comb, file = '/work/morrissy_lab/alyaa.mohamed/psychENCODE/int/Regulon_comb_brainspan.Rds')
