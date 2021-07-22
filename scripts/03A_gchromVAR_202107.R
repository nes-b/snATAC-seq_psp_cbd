###########################.
## TITLE: Identification of brain cell populations with GWAS SNP enrichment with gchromVAR
## Author: Nils Briel, Feodor-Lynen-Strasse 17, 81377 Munich, Bavaria
## Date: "22/07/2021"
#@ In this part, we determine population-specific GWAS SNP enrichment of brain cells in PSP and CBD frontal cortex samples.
#@ Therefore, we include data from two previously published and reproduced GWAS studies involving either PSP (https://www.niagads.org/summary-statistics-psp) or CBD individuals. 
  # download reference GWAS data from:
  # https://www.ebi.ac.uk/gwas/efotraits/Orphanet_683 --> PSP
  # https://www.ebi.ac.uk/gwas/efotraits/Orphanet_278 --> CBD
  # https://www.ebi.ac.uk/gwas/efotraits/Orphanet_282 --> FTD
  # https://www.ebi.ac.uk/gwas/efotraits/EFO_0000249 --> AD
  # https://www.ebi.ac.uk/gwas/efotraits/EFO_0002508 --> PD
  # https://www.ebi.ac.uk/gwas/efotraits/EFO_1001050 --> MSA
  # https://www.ebi.ac.uk/gwas/efotraits/EFO_0006792 --> LBD
  # https://www.ebi.ac.uk/gwas/efotraits/EFO_0000253 --> ALS
###########################.
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(SnapATAC)
library(data.table)
set.seed(123)
setwd('..')

#### 1 Prepare input ####
# 1.1 Prepare peaks
# First, we need to construct a SumamrizedExperiment object, fed with peak counts from previous peak calling in SnapATAC:
    
# take non-binarized pmat from SnapATAC and prepare input with cluster assignment
x.sp <- readRDS("psp_cbd_ctrl_peaks_rev.Rds") 
# exclude unmapped regions

newdata <- ! seqnames(x.sp@peak) %in% c('s37d5', unique(as.vector(seqnames(x.sp@peak[grepl('GL', seqnames(x.sp@peak))]))))

# preaparation function
preparePeakMatFromXSP <- function(SnapATACFile, newdata){
  ct <- as.factor(SnapATACFile@metaData$celltype)
  # construct per cell type sums of count matrix
  x <- SnapATACFile@pmat[,as.vector(newdata)]
  l = list()
  for(i in levels(ct)){
    l[[i]] <- DelayedMatrixStats::colSums2(x, rows = ct==i)
    }
  x <- t(Reduce(cbind, l))
  rownames(x) <- levels(ct)

  compl_peaks <- which(Matrix::colSums(x) > 0)
  
  x <- x[,compl_peaks] %>% t()
  
  rowD <- x.sp@peak[newdata] %>% .[compl_peaks]
  
  # call functions from SummarizedExperiment to construct SE 
  SE <- SummarizedExperiment(assays = list(counts = x),
                             rowData = rowD, 
                             colData = DataFrame(names = colnames(x)))
  
  # ... and chromVAR to correct for hg19 GCbias
  SE <- suppressWarnings(addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19))
  return(SE)
  print('SE object construciton successfull')
}

SE_tau <- preparePeakMatFromXSP(x.sp[which(x.sp@metaData$disease %in% c('PSP', 'CBD')),], newdata = newdata)
SE_tau@metadata$dis = 'tauopathy'

saveRDS(SE_tau, 'gchromVAR/SE_tau.Rds')

# clean R.envir
rm(list=setdiff(ls(), c("SE_tau")))
 
# 1.2 Prepare GWAS summary
setwd('gchromVAR')

# define function for input formatting - produce gchromVAR compatible GWAS summary input
 rearrGWASsum <- function(path, format, pvalue, chr=NULL, bp=NULL, psig=NULL){
  
  if(format == 'xlsx'){
    table <- readxl::read_xlsx(path) %>% 
      dplyr::select(chr, bp, psig)
  }else if(format == 'xls'){
    table <- readxl::read_xls(path) %>% 
      dplyr::select(chr, bp, psig)
  }else{
    table <- read_delim(file = path, 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      dplyr::select('CHR_ID', 'CHR_POS', 'P-VALUE')
  }
  colnames(table) = c('chr.1', 'location.bp.1', 'pval')
  table$chr.1 = as.factor(paste0('chr', table$chr.1))
  table$loc2 = as.integer(as.integer(table$location.bp.1)+1)
  table$location.bp.1 = as.integer(table$location.bp.1)
  colnames(table) = c('V1', 'V2', 'V3', 'V4') 
  table$V5 = as.factor('region1')
  table <- table %>% dplyr::select(1,2,4,5,3)     
  colnames(table) = c('V1', 'V2', 'V3', 'V4', 'V5')
  subset(table, V5 <= pvalue)
  table <- table[complete.cases(table),]
  return(table)
}

    # apply function on PSP GWAS
    GWASsumPSP <-  rearrGWASsum('gwas-association-downloaded_2021-01-07-Orphanet_683.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumPSP, 'GWASsumPSP.bed')
    
    # apply function on CBD GWAS
    GWASsumCBD <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-Orphanet_278.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumCBD, 'GWASsumCBD.bed')
    
    # apply function on FTD GWAS
    GWASsumFTD <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-Orphanet_683.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumFTD, 'GWASsumFTD.bed')
    
    # apply function on AD GWAS
    GWASsumAD <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-EFO_1001870-withChildTraits.tsv', 
                               format='tsv', 
                               pvalue = 0.05)
    write.table(GWASsumAD, 'GWASsumAD.bed')
    
    # apply function on PD GWAS
    GWASsumPD <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-Orphanet_2828-withChildTraits.tsv', 
                               format='tsv', 
                               pvalue = 0.05)
    write.table(GWASsumPD, 'GWASsumPD.bed')
    
    # apply function on MSA GWAS
    GWASsumMSA <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-EFO_1001050.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumMSA, 'GWASsumMSA.bed')
    
    # apply function on LBD GWAS
    GWASsumLBD <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-EFO_0006792.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumLBD, 'GWASsumLBD.bed')
    
    # apply function on ALS GWAS
    GWASsumALS <-  rearrGWASsum(path='gwas-association-downloaded_2021-01-07-EFO_0001357-withChildTraits.tsv', 
                                format='tsv', 
                                pvalue = 0.05)
    write.table(GWASsumALS, 'GWASsumALS.bed')
    
    # clean R.envir
    rm(list=setdiff(ls(), c("SE_tau")))
    files <- list.files('../gchromVAR',
                        full.names = TRUE, pattern = ".bed$")
    head(read.table(files[6]))
     
  # loop over PSP and CBD 
   
#### 2 Construct weight scores ####
bedscores <- importBedScore(rowRanges(SE_tau), files, colidx = 5)

#Compute weights deviation:
        
    wDEV <- computeWeightedDeviations(SE_tau, bedscores)
    zdf <- reshape2::melt(t(assays(wDEV)[["z"]]))
    zdf[,2] <- gsub("_PP001", "", zdf[,2])
    colnames(zdf) <- c("ct", "tr", "Zscore")
    head(zdf)
    zdf <- zdf[complete.cases(zdf),]
 
#### 3 Lineage-specific enrichment test ####
## Annotate your single cell lineages..
    Ast <- c('Ast')
    Mic <- c('Mic')
    Neu <- c('Exc. ULN', 'Exc. DLN', 'Exc. N. NEFM|BCL11B', 'Inh. N. SST|PVALB', 'Inh. N. VIP')
    Oli <- c('Oli #1 MOBP',  ' Oli #2',   'Oli #3' , 'OPC')
 
# annotate df lineage-specific (LS) enrichments:
zdf$LS  <-
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumMSA" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumMSA" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumMSA" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumMSA" & zdf$ct %in% Oli) +
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Oli) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Oli)
 
# Mann-Whitney Rank-sum statistic for relative enrichment
  set.seed(123)
  zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
  gchromVAR_ranksum <- sum(1:dim(zdf)[1]*zdf[order(zdf$gchromVAR_pvalue, decreasing = FALSE), "LS"])
  permuted <- sapply(1:10000, function(i) sum(1:dim(zdf)[1] * sample(zdf$LS, length(zdf$LS))))
  pnorm((mean(permuted) - gchromVAR_ranksum)/sd(permuted), lower.tail = FALSE)
   
#### 4 Extracting individual enrichment: ####
  zdf$adj.pvalue <- p.adjust(zdf$gchromVAR_pvalue, method = 'bonferroni')
  zdf$sign <- ifelse(zdf$adj.pvalue < .05, 'sign', 'n.s.') %>% as.factor()
  saveRDS(zdf, paste0('../output/gchromVAR_z_df_',SE_tau@metadata$dis,'.Rds'))


print("Part 07 is done: Now continue with '03B_ast_clus_analysis_202107.R'")
   
   