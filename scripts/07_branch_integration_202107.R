###########################.
## TITLE: branch integration
## Author: Nils Briel, Feodor-Lynen-Strasse 17, 81377 Munich, Bavaria
## Date: "22/07/2021"
## Descriptor: In this R-script we'll load TF hits from:
#  1 RTN-TA-gsea1 results (obtained from Allen et al. 2018, Acta Neuropathologica, and analysed with Bioconductor::RTN-package on PSP TCX bulkRNA-seq data)
#  2 tradeSeq's implemented tests' results (for lineage association, pattern-, startVsend-, endVsend-differences and visual clustering-results)
#  3 SnapATAC's/chromVAR-motif pair-wise two_sample_wilcoxon-test results (significantly scoring between triangular comparison of Ctrl-PSP-CBD[list#1], 
#                                                                          and exclusively significant candidates not shared with other comparisons[list#2])
#  4 XGB-model disease-specific weights (among the highest disease-specific feature_values & highest enrichment compared to the other groups)
# 
# and describe the intersections of TF candidates, which resulted from the particular analysis branches. 
###########################.
library(SnapATAC)
library(monocle3)
library(cicero)
library(Biobase)
library(tidyverse)
library(lime)
library(UpSetR)
library(VennDiagram)
library(data.table)
library(ggcorrplot)
library(finalfit)
library(parallel)
library(matrixTests)
set.seed(42)
setwd("..")
###########################.

#### 1 Combine single analysis branches ####
  ##### 1.1 PSP-Allen2018-dataset: RTN associated with TA-scores ####
  # keep only significant ones 
  gsea1_TA <- readRDS("~/projs/03_Knit_RTN/RTN/final/gsea1_final_TA.Rds") 
  gsea1_TA_sig <- subset(gsea1_TA, Adjusted.Pvalue < 0.05) %>% rownames()
  mra <- read.csv("~/projs/03_Knit_RTN/mra_regulons.csv", row.names = 1)
  mra_sig <- subset(mra, Adjusted.Pvalue < 0.05) %>% rownames()
  
  ##### 1.2 tradeSeq data ####
  tf_int <- readRDS("output/Ast_rev_cbd/5.2.3_TF_int_tradeseq.Rds")
  TF_int_1 <- tf_int[[1]] #Tradeseq's waldtest-results
  TF_int_2 <- tf_int[[2]] #Tradeseq's visual clustering results  
  TF_int_1_sig <- as.data.frame(TF_int_1) %>% separate(col = 'TF_int_1', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF
  TF_int_2_sig <- as.data.frame(TF_int_2) %>% separate(col = 'TF_int_2', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF
  
  ##### 1.3 Group-wise comparisons (triangular) SnapATAC / chromVAR-motifs results ####
  tf_mat_snap <- readRDS("TF_mat_comparison_snap.Rds") %>% .[[1]] %>% .$TF
  
  ##### 1.4 XGB-model lime explainer ####
  lime_explanation <- readRDS("output/modelling_rev/lime_explanation.Rds")
  
  # Function to compute rank stats in lime output:
  ## that score highest in terms of enrichment values among all three groups,
  comp_rank <- function(tab){
    temp <- tab
    #aggregate data by feature(TF) and label(disease) of interest and compute medians
    ref <- lime_explanation[which(lime_explanation$label != tab[1,'label'][[1]] & lime_explanation$feature == tab[1,'feature'][[1]] ),c('label', 'feature_value')] %>% 
      group_by(label) %>% 
      dplyr::summarise(median = median(feature_value))
    ref$feature <- tab[1,'feature'][[1]]
    ref <- ref %>% pivot_wider(names_from = 'label', values_from = 'median')
    # compute the same aggregation stats for the remaining labels
    for(f in tab[,'feature'][[1]][2:length(tab[,'feature'][[1]])]){
      ref_ <- lime_explanation[which(lime_explanation$label != tab[1,'label'][[1]] & lime_explanation$feature == f),c('label', 'feature_value')] %>% 
        group_by(label) %>% 
        dplyr::summarise(median = median(feature_value))
      ref_$feature <- f
      ref_ <- ref_ %>% pivot_wider(names_from = 'label', values_from = 'median')
      ref <- rbind(ref, ref_)
    }
    # join these tables, calculate feature_value ranks for each feature(TF) and report the label with the highest score
    tab <- left_join(temp, ref, by = 'feature',keep = F) %>% unique()
    x <- t(apply(-tab[,c(9,14,15)],1,rank))
    # return table which only contains the highest scoring ones, ordered by feature_weight
    tab <- tab[which(x[,1] == 1),] %>% .[order(.[,"feature_weight"]),]
    return(tab)
  }
  
  # We select only those TFs that are among the top 30 by their feature_weight for a given disease/label in lime's output and subject them to the comp_rank function
  expl_psp <- lime_explanation[which(lime_explanation$label == 'PSP'),] %>% .[order(.[,"feature_weight"], decreasing = T),] %>% .[1:30,] %>% comp_rank()
  expl_cbd <- lime_explanation[which(lime_explanation$label == 'CBD'),] %>% .[order(.[,"feature_weight"], decreasing = T),] %>% .[1:30,] %>% comp_rank()
  expl_ctrl <- lime_explanation[which(lime_explanation$label == 'Ctrl'),] %>% .[order(.[,"feature_weight", decreasing = T]),] %>% .[1:30,] %>% comp_rank()
  # retrieve the actual TFs without motif-identifier
  xgb_psp <- as.data.frame(expl_psp) %>% separate(col = 'feature', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF %>% .[which(!is.na(.))] %>% unique()
  xgb_cbd <- as.data.frame(expl_cbd) %>% separate(col = 'feature', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF %>% .[which(!is.na(.))] %>% unique()
  xgb_ctrl <- as.data.frame(expl_ctrl) %>% separate(col = 'feature', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF %>% .[which(!is.na(.))] %>% unique()

  
#### 2 1st stage integration ####
  # to show obtain overlap between
  #  1 RTN-TA-gsea1 results 
  #  2 tradeSeq's wald tests and visual clustering-results 
  #  3 SnapATAC's/chromVAR-motif disease-wise wilcoxon-test results 
  #  4 XGB-model disease-specific weights   
  # which suggest tau-associatied changes
  
  x = list(
    RTN = c(gsea1_TA_sig, mra_sig),
    trajectory = c(TF_int_1_sig, TF_int_2_sig),
    Pw.TFME = tf_mat_snap,
    XGB_model = c(xgb_psp, xgb_cbd, xgb_ctrl))
  
  saveRDS(x, 'output/upset_list_all.Rds')
  
  tiff(paste0('output/Ast_rev/4.6_upset_tau.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2000, height = 1000)
  upset(fromList(x), order.by = 'degree', sets.bar.color = 'lightblue')
  dev.off()
  
  int_1 <- Reduce(intersect, x)
  int_2 <- Reduce(intersect, list(x[['XGB_model']], x[['trajectory']], x[['Pw.TFME']])) %>% .[!(. %in% int_1)]
  int_3 <- Reduce(intersect, list(x[['XGB_model']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2)]
  int_4 <- Reduce(intersect, list(x[['trajectory']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2 | . %in% int_3)]
  xplan <-  as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(c(int_1, int_2, int_3, int_4))) %>% colnames() 

  
#### 3 2nd stage integration ####
  # to show divergence through
  #  1 SnapATAC's/chromVAR-motif disease-wise Wilcoxon-test results 
  #  2 XGB-model disease-specific weights   
  
  # SnapATAC / chromVAR-motifs results
  tf_mat_snap <- readRDS("TF_mat_comparison_snap.Rds") %>% .[[1]] 
  tf_mat_snap <- subset(tf_mat_snap, comparison %in% "PSP vs. CBD") %>% .$TF
  
  y = list(Pw.TFME_psp.vs.cbd = tf_mat_snap, XGB_PSP = xgb_psp, XGB_CBD = xgb_cbd)
  saveRDS(y, 'output/upset_list_psp_cbd.Rds')
  tiff(paste0('output/Ast_rev/4.6_upset_tau_diff.tiff'), units = 'px',compression = 'lzw', res = 300, width = 1800, height = 1000)
  upset(fromList(y), order.by = 'degree', sets.bar.color = 'lightblue')
  dev.off()
  
  #int_1 <- Reduce(intersect, list(xgb_cbd, tf_mat_snap))
  #int_2 <- Reduce(intersect, list(xgb_psp, tf_mat_snap))
  int_5 <- Reduce(intersect, y)
  int_6 <- Reduce(intersect, list(y[['XGB_CBD']], y[['Pw.TFME_psp.vs.cbd']]))
  int_7 <- Reduce(intersect, list(y[['XGB_PSP']], y[['Pw.TFME_psp.vs.cbd']]))
  explanatory_c = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_6)) %>% colnames()
  explanatory_p = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_7)) %>% colnames()
  
  rm(list = ls())


#### 4 Biological pathways related to resulting candidates from integration ####
  ##### 4.1 TF pathway analysis ####
  library("pathfindR")
  library(colorspace)
  library(ggpubr)
  source('scripts/04_cicero_analysis_functions.R')
  x.sp <- readRDS("output/Ast_rev/psp_cbd_ast.Rds")  
  annot <- getTFannotations(unique(colnames(x.sp@mmat)))
  
  levs <- levels(x.sp@metaData[,'disease'])[2:3]
  temp_list <- list()
  temp_list <- mclapply(levs, mc.cores = 4, mc.set.seed = 42, function(dis){
    mat <- x.sp@mmat[which(x.sp@metaData[,'disease'] %in% dis),] 
    mat <- t(mat[,annot$motif])
    mat_bg <- t(x.sp.ctrl@mmat[,annot$motif])
    x <- matrixTests::row_wilcoxon_twosample(as.matrix(mat), as.matrix(mat_bg)) # Wilcoxon rank-sum test
    x$FDR <- p.adjust(x$pvalue, method = 'BH', n = length(x$pvalue)) # Benjamini-Hochberg False discovery correction applied
    x$median.diff <- matrixStats::rowMedians(as.matrix(mat)) - matrixStats::rowMedians(as.matrix(mat_bg))
    x$comp <- paste0(dis, ' vs. Ctrl - adjust. method: Benjamini-Hochmerg')
    x$motif <- rownames(x)
    temp_list[[dis]] <- x
  } )
  
  pathPSP <- left_join(annot,temp_list[[1]], by = 'motif') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
  pathCBD <- left_join(annot,temp_list[[2]], by = 'motif') %>% dplyr::select('gene_name', 'median.diff', 'FDR')
  
  # Fig.6f-g
  set.seed(42)
  annot <- getTFannotations(unique(explanatory_p))
  output_psp_ <- run_pathfindR(pathPSP[which(pathPSP$gene_name %in% annot$gene_name),], output_dir = 'fig_06_PSP_explanatory')
  heatp <- term_gene_heatmap(output_psp_[order(output_psp_$Fold_Enrichment),][1:25,], use_description = T, genes_df = pathPSP, num_terms=15,
                             low = diverge_hcl(n = 3, palette = 'Blue-Red 3')[1],
                             mid = diverge_hcl(n = 3, palette = 'Blue-Red 3')[2],
                             high = diverge_hcl(n = 3, palette = 'Blue-Red 3')[3]) + coord_fixed()
  annot <- getTFannotations(unique(explanatory_c))
  output_cbd_ <- run_pathfindR(pathCBD[which(pathCBD$gene_name %in% annot$gene_name),], output_dir = 'fig_06_CBD_explanatory')
  heatc <- term_gene_heatmap(output_cbd_, use_description = T, genes_df = pathCBD, num_terms=15,
                             low = diverge_hcl(n = 3, palette = 'Blue-Red 3')[1],
                             mid = diverge_hcl(n = 3, palette = 'Blue-Red 3')[2],
                             high = diverge_hcl(n = 3, palette = 'Blue-Red 3')[3]) + coord_fixed()
  cairo_pdf(file = "output/Fig06_pathways.pdf", width = 12, height = 10)
  plot(plot_grid(heatp,heatc,ncol = 2, labels = c('PSP','CBD')))
  dev.off()
  
  ##
  mat <- x.sp@mmat[,annot$motif] %>% scales::rescale(to = c(0,1), from = range(., na.rm = TRUE, finite = TRUE)) %>% log1p()
  colnames(mat) <- annot$gene_name
  rownames(mat) <- paste('id',1:nrow(mat))
  score_matrix <- score_terms(enrichment_table = clustered_psp[clustered_psp$Status == "Representative", ],
                              exp_mat = t(mat[which(x.sp@metaData$disease %in% c('Ctrl', 'PSP')),]),
                              cases = rownames(mat)[which(x.sp@metaData$disease %in% 'PSP')],
                              use_description = TRUE,
                              label_samples = FALSE, 
                              case_title = "PSP",  
                              control_title = "Control", 
                              low = 'black',
                              mid = 'white',
                              high = 'blue')
  hmap <- plot_scores(score_matrix, label_samples = FALSE) + geom_tile(colour=NA,width=0, alpha = 0)
  combined_df <- combine_pathfindR_results(result_A = output_psp, 
                                           result_B = output_cbd, 
                                           plot_common = FALSE)
  
  library(ggraph)
  sel1 <- combined_df[which(combined_df$status %in% 'common'),]
  sel2 <- sel1[order(sel1$Fold_Enrichment_A),][1:10,'Term_Description']
  graph <- combined_results_graph(combined_df, use_description = T,selected_terms = c(sel2)) 

  ##### 4.2 Differentiate protein homeostasis pathways #####
    library(pheatmap)
    library(ggpubr)
    x.sp <- readRDS("psp_cbd_ast.Rds")
    
    cma_amiGO <- read.delim("refs/cma_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    upr_amiGO <- read.delim("refs/upr_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    ups_amiGO <- read.delim("refs/ups_genes_amiGO.txt", header=FALSE, col.names = c('uniprotkb', 'gene')) %>% as.data.table() %>% dplyr::filter(gene %in% colnames(x.sp@gmat)) 
    cma_amiGO$data <- c('CMA')
    upr_amiGO$data <- c('UPR')
    ups_amiGO$data <- c('UPS')
    int <- data.table::rbindlist(list(cma_amiGO,upr_amiGO, ups_amiGO)) %>% dplyr::filter(gene %in% colnames(x.sp@gmat))
    
    ## heatmaps ##
    plot_list_heat <- list(ups_amiGO, cma_amiGO,upr_amiGO)
    for(i in 1:3){
      # assign degradation pathway-specific gene set to 'p'
      p <- plot_list_heat[[i]] 
      # prepare pheatmap input for every degradation gene set
      ref1 <- col_t_welch(x = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'PSP'),p$gene]),
                          y = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'Ctrl'),p$gene]))
      ref1$bonf <- p.adjust(ref1$pvalue, method = 'bonferroni', n = length(ref1$pvalue)) 
      ref2 <- col_t_welch(x = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'CBD'),p$gene]),
                          y = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'Ctrl'),p$gene]))
      ref2$bonf <- p.adjust(ref2$pvalue, method = 'bonferroni', n = length(ref2$pvalue)) 
      idy <- which(ref1$bonf < 5e-2 | ref2$bonf < 5e-2 )
      
      mat <- t(x.sp@gmat[,p$gene[idy]]) %>% .[rowSums(.)>0,]
      
      mat_row <- data.frame(row.names = paste0('gene_',rownames(mat)), 
                            Category = as.factor(subset(p, gene %in% rownames(mat))$data))
      colnames(mat) <- paste0('cell_',1:nrow(x.sp@gmat))
      rownames(mat) <- paste0('gene_',rownames(mat))
      mat_col <- data.frame(row.names = colnames(mat),
                            Diagnosis = x.sp@metaData$disease, 
                            Celltype = x.sp@metaData$celltype)
      mat_col<- mat_col[order(as.character(mat_col$Diagnosis)),]
      mat <- mat[,rownames(mat_col)]
      mat_colors <- list(Celltype = 'red',
                         Diagnosis = c("#E41A1C", "#377EB8","#4DAF4A"),
                         Category =  rep(c('#a81675'), nrow(mat))
      )
      names(mat_colors$Celltype) <- unique(mat_col$Celltype)
      names(mat_colors$Diagnosis) <- unique(mat_col$Diagnosis)
      names(mat_colors$Category) <- unique(mat_row$Category)
      plot_list_heat[[paste(p$data[1])]] <- list(mat = mat,
                                                 mat_row = mat_row,
                                                 mat_col = mat_col,
                                                 mat_colors = mat_colors)
    }
    
    ups <- plot(ggplotify::as.ggplot(pheatmap(log1p(plot_list_heat[[4]]$mat), 
                                              main="UPS",
                                              annotation_col = plot_list_heat[[4]]$mat_col, 
                                              annotation_row = plot_list_heat[[4]]$mat_row,
                                              annotation_colors = plot_list_heat[[4]]$mat_colors,
                                              show_colnames = F, show_rownames = T, annotation_legend = TRUE,
                                              treeheight_row = 30, 
                                              color = colorspace::diverge_hcl(palette = 'Berlin', 7),
                                              border_color = 'black',
                                              clustering_method = "ward.D2", fontsize = 7,
                                              cluster_rows = T,
                                              cluster_cols = F,
                                              scale = 'row',
                                              clustering_distance_row = 'manhattan')))
    cma <- plot(ggplotify::as.ggplot(pheatmap(log1p(plot_list_heat[[5]]$mat), 
                                              main="CMA",
                                              annotation_col = plot_list_heat[[5]]$mat_col, 
                                              annotation_row = plot_list_heat[[5]]$mat_row,
                                              annotation_colors = plot_list_heat[[5]]$mat_colors,
                                              show_colnames = F, show_rownames = T, annotation_legend = TRUE,
                                              treeheight_row = 30, 
                                              color = colorspace::diverge_hcl(palette = 'Berlin', 7),
                                              border_color = 'black',
                                              clustering_method = "ward.D2", fontsize = 7,
                                              cluster_rows = T,
                                              cluster_cols = F,
                                              scale = 'row',
                                              clustering_distance_row = 'manhattan')))
    upr <- plot(ggplotify::as.ggplot(pheatmap(log1p(plot_list_heat[[6]]$mat), 
                                              main="UPR",
                                              annotation_col = plot_list_heat[[6]]$mat_col, 
                                              annotation_row = plot_list_heat[[6]]$mat_row,
                                              annotation_colors = plot_list_heat[[6]]$mat_colors,
                                              show_colnames = F, show_rownames = T, annotation_legend = TRUE,
                                              treeheight_row = 30, 
                                              color = colorspace::diverge_hcl(palette = 'Berlin', 7),
                                              border_color = 'black',
                                              clustering_method = "ward.D2", fontsize = 7,
                                              cluster_rows = T,
                                              cluster_cols = F,
                                              scale = 'row',
                                              clustering_distance_row = 'manhattan')))
  tiff(paste0('output/supplFig14_heats.tiff'), units = 'px',compression = 'lzw', res = 450, width = 5000, height = 5500)
  cowplot::plot_grid(ups, cma, upr,ncol = 1)
  dev.off()

  ##### 4.3 GA changes protein-homeostasis candidates along pseudotime ####
  library(tradeSeq)
  sce <- readRDS('output/Ast_rev_cbd/fitGAM_ga.Rds')
  ts_ga <- readRDS('output/Ast_rev_cbd/5.1.3_GA_int_tradeseq.Rds')
  plot_list_traj <- list(ups_amiGO, cma_amiGO,upr_amiGO)
  df <- data.frame()
  for(i in 1:3){
    # assign degradation pathway-specific gene set to 'p'
    p <- plot_list_traj[[i]] 
  ref2 <- col_t_welch(x = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'CBD'),p$gene]),
                      y = as.matrix(x.sp@gmat[which(x.sp@metaData$disease == 'Ctrl'),p$gene]))
  ref2$bonf <- p.adjust(ref2$pvalue, method = 'bonferroni', n = length(ref2$pvalue)) 
  ref2$data <- as.factor(p$data[1])
  ref2 <- subset(ref2, bonf < 5e-2)
   df <- rbind(df, ref2)
  }
  
  idy <- df[intersect(rownames(df), ts_ga[[1]]),]
  mat <- as.matrix(t(x.sp@gmat[,rownames(idy)]))
  yhatSmooth <- tradeSeq::predictSmooth(sce, gene =rownames(idy), nPoints = 256, tidy = F)
  heat2 <- plot(ggplotify::as.ggplot(pheatmap::pheatmap(t(scale(t(yhatSmooth))),
                                                        color = viridis::viridis(256),cluster_rows = F,
                                                        cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize= 8)))
  
  cairo_pdf(paste0('output/supplFig14_trajs_cbd_asts.pdf'), 
            width = 7, height = 4)
  cowplot::plot_grid(heat2,ncol = 1)
  dev.off()


###################### finished script ############################.
print("Part 08 is done.")
sessionInfo()



