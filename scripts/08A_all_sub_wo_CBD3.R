###########################.
## TITLE: SnapATAC analysis - Part 2
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "26/10/2021"
## Description: Here we assess differentially acessible regions (DARs) of each cluster and subject these findings to a Gene Ontology database fetching package.
##              We do also apply chromVAR-motif to extract TF motif enrichments in the peaks, and finally calculate z-scores of protein degradation involvement. 
###########################.

Sys.setenv(RETICULATE_PYTHON = "/home/nes/miniconda3/bin/python") ## modify here to navigate to the specific directory
library(reticulate)
path_to_python <- "/home/nes/miniconda3/bin/python" ## modify here to navigate to the specific directory
use_python(path_to_python, required = T)
py_config()
library(SnapATAC)
library(GenomicRanges)
library(viridisLite)
library(tidyverse)
library(leiden)
library(SummarizedExperiment)
library(ggpubr)
library(data.table)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
setwd("..")
source('scripts/plotViz2.R')

### ViZ w/o CBD3
    x.sp <- readRDS("psp_cbd_ctrl_motifs_rev.Rds")
    
    x.sp@metaData$disease <- factor(x.sp@metaData$disease, levels = c("Ctrl", "PSP", "CBD"))
    dis_cols <- c("#377EB8","#4DAF4A","#E41A1C")
    
    x.sp <- x.sp[which(x.sp@metaData$case != 'num112'),]
    
    
    ##### 2.3 Batch effect normalization ####
    library(harmony);
    x.after.sp = runHarmony(obj=x.sp,eigs.dim=1:25, 
                            meta_data=x.sp@sample # sample index
    );
    
    ##### 2.4 Graph-based clustering ####
    x.after.sp = runKNN(obj=x.after.sp,
                        eigs.dims=1:25,k=15
    );
    
    x.after.sp=runCluster(obj=x.after.sp,tmp.folder=tempdir(),
                          louvain.lib="leiden", # sophisticated community detection algorithm
                          seed.use=10,resolution=1
    );
    x.after.sp@metaData$cluster = x.after.sp@cluster
   # clus_cols <- c("#9970ff", "#458cff", "#85db8b",  "#1a8239" ,"#fa7878", "#afdb85", "#00b8b1", "#fcc44c",  "#3efae7", "#e695dc", 
   #                "#ccb8ff", "#deb8ff", "#a491c4", "#ecb8ff","#ffa6ff",  "#ff96ff")
    ct_cols <- c('#fa7878',  '#afdb85',  '#85db8b',  '#1a8239',  '#00b8b1',  '#3efae7',  '#fcc44c',  '#458cff',  '#7370ff',  '#9970ff',  '#e695dc')
    
    cairo_pdf(file = "output/Ast_sub_wo_CBD3/Fig2ct_UMAP_20210828.pdf", width = 8, height = 8)
    plotViz2(obj=x.sp, method="umap", main="",
            point.color='celltype', text.add=F, legend.add=TRUE
    ) + scale_fill_manual(values = ct_cols)
    dev.off()
    
    dir.create('output/Ast_sub_wo_CBD3')
    saveRDS(x.after.sp, 'output/Ast_sub_wo_CBD3/psp_cbd_ctrl_motifs_wo_CBD3.Rds')  

    gc()
    
##### 2.9 Identify differentially accessible regions ####
##### 2.10 chromVAR-motif ####
   set.seed(42)
   x.cbd <- subset(x.after.sp@metaData, disease == 'CBD')
   x.psp <- subset(x.after.sp@metaData, disease == 'PSP')
   x.ctrl <- subset(x.after.sp@metaData, disease == 'Ctrl')
   
   x.ctrl <- x.ctrl %>% dplyr::count(as.character(.$celltype)) %>%     
     mutate(prop = prop.table(n))
   
   x.psp_list <- list()
   x.psp_list <- mclapply(levels(as.factor(x.psp$case)), mc.cores = 4, function(i){
     x <- subset(x.psp, case == i)
     x <- x %>% dplyr::count(as.character(.$celltype)) %>%     
       mutate(prop = prop.table(n))
     if(nrow(x) == nrow(x.ctrl)){
       x$prop_to_ctrl <- x$prop/x.ctrl$prop
       x$prop_ctrl <- x.ctrl$prop
       
     }else{
       x <- left_join(select(x.ctrl, c(-prop, -n)), x, by = 'as.character(.$celltype)')
       x[is.na(x)] <- 0
       x$prop_to_ctrl <- x$prop/x.ctrl$prop
       x$prop_ctrl <- x.ctrl$prop
     }
     colnames(x)[1] <- 'celltype'
     x$case <- i
     x.psp_list[[i]] <- x
   }
   )
   x.psp <- data.table::rbindlist(x.psp_list)
   
   x.cbd_list <- list()
   x.cbd_list <- mclapply(levels(as.factor(x.cbd$case)), 
                          mc.cores = 4, 
                          function(i){
                            x <- base::subset(x.cbd, case == i)
                            x <- x %>% dplyr::count(as.character(.$celltype)) %>%     
                              dplyr::mutate(prop = prop.table(n))
                            if(nrow(x) == nrow(x.ctrl)){
                              x$prop_to_ctrl <- x$prop/x.ctrl$prop
                              x$prop_ctrl <- x.ctrl$prop
                              
                            }else{
                              x <- left_join(dplyr::select(x.ctrl, c(-prop, -n)), x, by = 'as.character(.$celltype)')
                              x[is.na(x)] <- 0
                              x$prop_to_ctrl <- x$prop/x.ctrl$prop
                              x$prop_ctrl <- x.ctrl$prop
                            }
                            colnames(x)[1] <- 'celltype'
                            x$case <- i
                            x.cbd_list[[i]] <- x
                          }
   )
   
   x.cbd <- data.table::rbindlist(x.cbd_list)
   
   x <- data.frame(table(x.after.sp@metaData$celltype))
   colnames(x)[1] <- 'Celltype'
   
   x.psp_p <- ggpubr::ggboxplot(x.psp, x = 'celltype', y = 'prop_to_ctrl', fill = 'celltype', 
                                xlab = '', ylab = 'Relative frequency') +
     geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) + ylim(0,4.5) + #scale_x_discrete(limits=rev) +
     coord_flip() + 
     stat_compare_means(label = "p.signif", method = "t.test",
                        ref.group = ".all.", hide.ns = T)+
     ggthemes::theme_few(base_family="Arial", base_size = 14)  + scale_fill_manual(values = ct_cols) 
   
   x.cbd_p <- ggpubr::ggboxplot(x.cbd, x = 'celltype', y = 'prop_to_ctrl', fill = 'celltype', 
                                xlab = '', ylab = 'Relative frequency') +
     geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) + ylim(0,4.5) + #scale_x_discrete(limits=rev) + 
     coord_flip() + stat_compare_means(label = "p.signif", method = "t.test",
                                       ref.group = ".all.", hide.ns = T) + scale_fill_manual(values = ct_cols) + 
     ggthemes::theme_few(base_family="Arial", base_size = 14) + theme(axis.text.y=element_blank())
   
   tot_ct <- ggplot(x, aes(x = Freq, y = Celltype, fill = Celltype)) + geom_col(color = 'black') + labs(x = 'Total frequency', y = '', fill = 'Cell type') +
     scale_fill_manual(values = ct_cols) + geom_text(aes(label=Freq),color = 'black',size = 4,hjust = -0.15) + scale_x_continuous(limits = c(0,1.1e4)) + scale_y_discrete(limits=rev) + 
     ggthemes::theme_tufte(base_family="Arial", base_size = 14) + theme(axis.text.y=element_blank())
   
   ct_frqs <- ggpubr::ggarrange(plotlist = list(x.psp_p, x.cbd_p, tot_ct), common.legend = T, ncol = 3, 
                                legend = 'right', labels = c('PSP', 'CBD', ''),label.y = 1.03,label.x=c(0.37,0,0),widths = c(1.4,0.8, 0.6))
   cairo_pdf(file = "output/Ast_sub_wo_CBD3/Fig2ct_frqs_20210828.pdf", width = 12, height = 8)
   plot(ct_frqs)
   dev.off()

# 2.10.1 chromVAR-motif visualization

# 2.10.2 cluster TF activity - pheatmap
   ##### 2.12 Assess protein degradation on system-level #####
    library(pheatmap)
    library(ggpubr)
    library(matrixTests)
    set.seed(42)
    x.sp <- readRDS("~/projs/scATACseq_MLL/psp_cbd_ctrl_protein_degradation_genes.Rds")
    rm(x.after.sp)
    
    cma_amiGO <- read.delim("refs/cma_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    upr_amiGO <- read.delim("refs/upr_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    ups_amiGO <- read.delim("refs/ups_genes_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    cma_amiGO$data <- c('CMA')
    upr_amiGO$data <- c('UPR')
    ups_amiGO$data <- c('UPS')

    set.seed(42)
    plot_list <- list(ups_amiGO, cma_amiGO,upr_amiGO)
    data_list <- readRDS('degra_data_list.Rds')
    
    # identify and remove num112/CBD3 cells from snap-obj
    idy <- rownames(x.sp@metaData[which(x.sp@metaData$case == 'num112' | x.sp@metaData$celltype == 'Oli #3'),]) %>% as.numeric()
    idy <- !(as.numeric(rownames(x.sp@metaData) %in% idy))
    idy <- rep(idy,3)
    
    Degra_table <- data.table::rbindlist(data_list) %>% as.data.frame() #.[which(complete.cases(.)),] 
    Degra_table <- Degra_table[idy,] %>% as.data.table()
    
    for(pathway in levels(as.factor(Degra_table$path))){
      degra_ref <- subset(Degra_table, path %in% pathway)
      df <- data.frame()
      for(ctype in levels(as.factor(degra_ref$ct))){
        p <- subset(degra_ref, ct %in% ctype) %>% group_by(dis) %>%
          summarise(mean = mean(.), sd = sd(.))
        p$ct <- ctype
        p$mean <- p$mean - p$mean[p$dis == 'Ctrl']
        v2 <- t.test(subset(degra_ref, ct %in% ctype & dis == 'PSP')$., subset(degra_ref, ct %in% ctype & dis == 'Ctrl')$.) 
        v3 <- t.test(subset(degra_ref, ct %in% ctype & dis == 'CBD')$., subset(degra_ref, ct %in% ctype & dis == 'Ctrl')$.) 
        p$pvalue = c(1, v2$p.value, v3$p.value)
        df <- rbind(df, p)
      }
      df$adj.pvalue <- p.adjust(df$pvalue, method = 'fdr', n = length(df$pvalue)) 
      df$updown <- ifelse(df$mean > 0, 'up', 'down')
      
      plot_list[[paste0(pathway,'_heat')]] <- plot(ggplot(df, aes(
        x = ct,
        y = dis,
        fill = mean)) + labs(fill = 'Z-score', y = '', x = '') + 
          colorspace::scale_fill_continuous_diverging("Blue-Red 3") + geom_tile(color = 'black',size = 0.5) +
          geom_text(aes(label=ifelse(pvalue <= 5e-2, signif((pvalue), 2),'')), nudge_y = 0.0, size = 3, fontface=2)+
          ggthemes::theme_tufte(base_family = 'arial', base_size = 12) + 
          theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))
      )
    }

    cairo_pdf(paste0('output/Ast_sub_wo_CBD3/Fig3_heatmap_degra_20210828.pdf'), width = 8, height = 12)
    plot(ggarrange(plotlist = plot_list[4:6], ncol = 1, labels = c('CMA', 'UPR', 'UPS'), common.legend = F))
    dev.off()
    
  
    
#### 2.11 chromVAR-motif visualization ####
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
 
  ## 2.11.1 Prepare input 
  # 2.11.1.1 Prepare peaks
  # First, we need to construct a SumamrizedExperiment object, fed with peak counts from previous peak calling in SnapATAC:
  
  # take non-binarized pmat from SnapATAC and prepare input with cluster assignment
  x.sp <- readRDS("output/Ast_sub_wo_CBD3/psp_cbd_ctrl_motifs_wo_CBD3.Rds") 
  
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
  
  dir.create('output/Ast_sub_wo_CBD3/gchromVAR')
  saveRDS(SE_tau, 'output/Ast_sub_wo_CBD3/gchromVAR/SE_tau_sub_wo_CBD3.Rds')
  
  # clean R.envir
  rm(list=setdiff(ls(), c("SE_tau")))
  gc()
  
  # 2.11.1.2 Prepare GWAS summary
  #setwd('output/Ast_sub_wo_CBD3/gchromVAR')
  
  # use same GWAS sum tables as in preceding gchromVAR analysis
  # clean R.envir
  files <- list.files('gchromVAR',
                      full.names = TRUE, pattern = ".bed$")
  head(read.table(files[6]))
  
  # loop over PSP and CBD 
  
 ## 2.11.2 Construct weight scores 
  bedscores <- importBedScore(rowRanges(SE_tau), files, colidx = 5)
  
  #Compute weights deviation:
  wDEV <- computeWeightedDeviations(SE_tau, bedscores)
  zdf <- reshape2::melt(t(assays(wDEV)[["z"]]))
  zdf[,2] <- gsub("_PP001", "", zdf[,2])
  colnames(zdf) <- c("ct", "tr", "Zscore")
  head(zdf)
  zdf <- zdf[complete.cases(zdf),]
  
 ## 2.11.3 Lineage-specific enrichment test 
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
  
 ## 2.11.4 Extracting individual enrichment
  zdf$adj.pvalue <- p.adjust(zdf$gchromVAR_pvalue, method = 'bonferroni')
  zdf$sign <- ifelse(zdf$adj.pvalue < .05, 'sign', 'n.s.') %>% as.factor()
  saveRDS(zdf, paste0('output/Ast_sub_wo_CBD3/gchromVAR/gchromVAR_z_df_',SE_tau@metadata$dis,'.Rds'))
  
  gchrom <- ggplot(zdf, aes(
    x = ct,
    y = tr,
    fill = Zscore)) + labs(fill = 'Z-score', y = '', x = '') + 
    colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red') + coord_equal() + geom_tile() + 
    geom_text(aes(label=signif(adj.pvalue, 3)), nudge_y = 0.2, size = 3, fontface=2) + 
    geom_text(aes(label=signif(gchromVAR_pvalue, 3)), nudge_y = -0.2, size = 2.5, fontface=3)+
    ggthemes::theme_tufte(base_family = 'Helvetica', base_size = 12) + theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))
  
  cairo_pdf(paste0('output/Ast_sub_wo_CBD3/gchromVAR/Fig3_heatmap_gchrom_sub_wo_CBD3.pdf'), width = 8, height = 6)
  plot(gchrom)
  dev.off()
  

################## finished #####################.
print("Part 8A is done: now continue with '08B_ast_sub_wo_CBD3.R'")
sessionInfo()

