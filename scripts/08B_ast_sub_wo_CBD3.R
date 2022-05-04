###########################.
## TITLE: Cicero analysis + trajectories
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "26/10/2021"
## Description: This R script is written to convert the astrocytes snap-file to cicero (1) and to compute cicero-connections on it(2).
##              In the next step (3) Seurat's/Signac's marker genes plot is used to visualize DARs/closest genes in the underlying dataset.
##              The complete steps summarized under (4) compute cicero/monocle3 dimension reduction, UMAP embeddings, clustering and pseudotime trajectories 
##              only for CBD and Ctrls astrocytes. tradeSeq analysis is conducted in (5).
##              Then the same workflow is also computed for PSP astrocytes (6).
###########################.
  
  library(SnapATAC)
  library(monocle3)
  library(cicero)
  library(Biobase)
  library(tidyverse)
  library(BiocGenerics)
  library(chromVAR)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(data.table)
  library(plotly)
  setwd("..")
  source('scripts/plotViz2_v3.R')
  source('scripts/04_cicero_analysis_functions.R')
############################.
  
  x.sp <- readRDS('psp_cbd_ctrl_motifs_rev.Rds')
  
  x.sp@metaData$disease <- factor(x.sp@metaData$disease, levels = c("Ctrl", "PSP", "CBD"))
  dis_cols <- c("#377EB8","#4DAF4A","#E41A1C")
  clus_cols <- ggpubr::get_palette(palette = "Dark2",length(levels(x.sp@cluster)))
  
  x.sp <- x.sp[which(x.sp@cluster==4 & x.sp@metaData$case != 'num112'),]
  
  x.sp@metaData$sample <- x.sp@sample
  set.seed(42)
  
  #### Astrocytes: Preprocessing and UMAP embedding  ####
  set.seed(100)
  cicero_cds <- snapToCicero(x.sp, preprocces_cic = F, preprocess_traj = F)
  cicero_cds@metadata$line <- 'pspcbdctrl'
  cicero_cds <- preprocess_cds(cicero_cds, method = "PCA")
  cicero_cds <- align_cds(cicero_cds, alignment_group = c('case'), preprocess_method = 'PCA')
  cicero_cds <- reduce_dimension(cicero_cds, preprocess_method = "Aligned",
                                 umap.metric = "cosine",umap.min_dist = 0.33, umap.n_neighbors = 18,
                                 cores = 8) 
  monocle3::plot_cells(cicero_cds, color_cells_by = 'disease', cell_size = 1) + scale_color_manual(values = dis_cols)
  saveRDS(cicero_cds, 'output/Ast_sub_wo_CBD3/cicero_cds_ast_rev3.Rds')
  
  # transferring Ast data from cicero_cds to snap-object, and perform k-means clustering
  rown <- rownames(cicero_cds@int_colData@listData$reducedDims@listData$UMAP)
  cicero_cds@int_colData@listData$reducedDims@listData$UMAP -> x.sp@umap[,]
  
  x.sp@metaData$clusterkmeans <- as.factor(kmeans(x.sp@umap, centers = 7, nstart = 40, iter.max = 1000)$cluster)
  x.sp@cluster <- x.sp@metaData$clusterkmeans
  
  # Plots differentiated by group entity 
  library(colorspace)
  p1 <- plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'PSP'),], method="umap", main="PSP", 
                 point.color='disease',point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
  ) + scale_fill_manual(values = '#4DAF4A')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
  p2 <-  plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'CBD'),],method="umap", main="CBD",
                  point.color='disease',point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
  ) + scale_fill_manual(values = '#E41A1C')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
  p3 <- plotViz2(obj=x.sp[which(x.sp@metaData$disease %in% 'Ctrl'),], method="umap", main="Ctrl",
                 point.color='disease', point.size=1.2,text.add=FALSE,legend.add=TRUE, stroke = 0.09, point.alpha = 1
  ) + scale_fill_manual(values = '#377EB8')+coord_equal()+theme(title = element_text(size = 14), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
  p4 <- plotViz2(obj=x.sp, method="umap", main="Clusters",
                 point.color='clusterkmeans', point.size=1.2,point.alpha=1,
                 text.add=TRUE, text.size=4.5,legend.add=FALSE, stroke = 0.09)+coord_equal() + scale_fill_manual(values = sequential_hcl(n=length(levels(x.sp@metaData$clusterkmeans)), palette = 'Reds')) +
    theme(title = element_text(size = 14), legend.text = element_text(size=10), legend.position = 'none')+  coord_fixed(xlim = c(-10,10), ylim = c(-10,10))
  
  tiff(paste0('output/Ast_sub_wo_CBD3/Fig4_ast_metadata_cluster_viz_rev20210518.tiff'), units = 'px',compression = 'lzw', res = 500, width = 7000, height = 7000)
  plot(cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), 
                          ncol = 4, nrow = 1, axis = 'tblr', align = 'hv') )
  dev.off()
  cairo_pdf(paste0('output/Ast_sub_wo_CBD3/Fig4_ast_metadata_cluster_viz_rev20210518.pdf'), width = 16, height = 4)
  plot(cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), 
                          ncol = 4, nrow = 1, axis = 'tblr', align = 'hv') )
  dev.off()
  
  saveRDS(x.sp, 'output/Ast_sub_wo_CBD3/psp_cbd_ast.Rds')
  
  #### 3B.1 Astclust: chromVAR-motif Group-wise Comparisons ####
  library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(gridExtra)
  
  library(matrixTests)
  set.seed(42)
  cpsp <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'Ctrl'),], x.sp@mmat[which(x.sp@metaData$disease == 'PSP'),])
  ccbd <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'Ctrl'),], x.sp@mmat[which(x.sp@metaData$disease == 'CBD'),])
  pspcbd <- col_wilcoxon_twosample(x.sp@mmat[which(x.sp@metaData$disease == 'CBD'),], x.sp@mmat[which(x.sp@metaData$disease == 'PSP'),])
  
  cpsp$bh <- p.adjust(cpsp$pvalue, 'BH', n = nrow(cpsp))
  ccbd$bh <- p.adjust(ccbd$pvalue, 'BH', n = nrow(ccbd))
  pspcbd$bh <- p.adjust(pspcbd$pvalue, 'BH', n = nrow(pspcbd))
  cpsp$feature <- rownames(cpsp)
  ccbd$feature <- rownames(ccbd)
  pspcbd$feature <- rownames(pspcbd)
  cpsp_s <- cpsp[which(cpsp$bh < 0.05),]
  ccbd_s <- ccbd[which(ccbd$bh < 0.05),]
  pspcbd_s <- pspcbd[which(pspcbd$bh < 0.05),] 
  cpsp_s$comparison <- 'PSP vs. Ctrl'
  ccbd_s$comparison <- 'CBD vs. Ctrl'
  pspcbd_s$comparison <- 'PSP vs. CBD'
  cpsp_s <- separate(cpsp_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  ccbd_s <- separate(ccbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  pspcbd_s <- separate(pspcbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  TF_mat_3comp <- data.table::rbindlist(list(cpsp_s, ccbd_s, pspcbd_s))
  TF_mat_3comp <- TF_mat_3comp[order(TF_mat_3comp$bh),]
  write.csv(TF_mat_3comp, 'output/Ast_sub_wo_CBD3/scatac_tfs_psp_cbd_allcomparisons.csv')
  
  cpsp_s <- cpsp[which(cpsp$bh < 0.05 & ccbd$bh > 0.05 & pspcbd$bh > 0.05),] 
  ccbd_s <- ccbd[which(ccbd$bh < 0.05 & cpsp$bh > 0.05 & pspcbd$bh > 0.05),]
  pspcbd_s <- pspcbd[which(pspcbd$bh < 0.05 & cpsp$bh > 0.05 & ccbd$bh > 0.05),]
  cpsp_s$comparison <- 'PSP vs. Ctrl'
  ccbd_s$comparison <- 'CBD vs. Ctrl'
  pspcbd_s$comparison <- 'PSP vs. CBD'
  cpsp_s <- separate(cpsp_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  ccbd_s <- separate(ccbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  pspcbd_s <- separate(pspcbd_s, col = 'feature', into =c('prefix', 'TF'), sep = c("_"), remove = F)
  TF_mat_3exclusive<- data.table::rbindlist(list(cpsp_s, ccbd_s, pspcbd_s))
  write.csv(TF_mat_3exclusive, 'output/Ast_sub_wo_CBD3/scatac_tfs_psp_cbd_exclcomparisons.csv')
  

  ##### Modelling disease ####
  library(caret)
  library(GenomicRanges)
  library(viridisLite)
  library(tidyverse)
  library(ggpubr)
  library(cowplot)
  library(grid)
  library(parallel)
  library(fastglm)
  library(ggrepel) 
  library(lime)
  # define label data
  outcome <- as.numeric(as.factor(x.sp@metaData[,c('disease')])) %>% scales::rescale(c(0,1))
  # define predictors
  mat_t <- cbind(as.matrix(x.sp@mmat),outcome)
  # train-test split
  idy <- sample(nrow(mat_t), nrow(mat_t)*0.8)
  mat_tr <- mat_t[idy,]
  mat_te <- mat_t[-idy,]
  
  # Train Xtreme Gradient Boosting Machine
  fitControl <- trainControl(method = "repeatedcv", number=10, repeats=3, classProbs = T)
  outcome <- as.factor(x.sp@metaData[,c('disease')])
  
  ## xgb Tree
  xgbTree_mod <- list()
  set.seed(825)
  mod <- train(y = outcome[idy], x = mat_tr[,-ncol(mat_tr)], 
               method = "xgbTree", 
               trControl = fitControl, 
               verbose = T, 
               metric = "ROC")
  
  xgbTree_mod[['summa']] <- summary(mod)
  xgbTree_mod[['predy']] <- stats::predict(mod, mat_te[,-ncol(mat_te)], type = 'raw')
  summary(xgbTree_mod[['predy']])
  xgbTree_mod[['confmat']] <- caret::confusionMatrix(xgbTree_mod[['predy']], as.factor(x.sp@metaData[-idy,c('disease')]))
  xgbTree_mod[['ROC_curves']] <- caTools::colAUC(as.numeric(as.factor(xgbTree_mod[['predy']])), as.factor(x.sp@metaData[-idy,c('disease')]), plotROC = TRUE)
  xgbTree_mod[['model']] <- mod
  
  
  ## Assess performance and confusion matrix
  tiff(paste0('output/Ast_sub_wo_CBD3/5_xgbTreeaccroc_dis_models.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2300, height = 1600) 
  caTools::colAUC(as.numeric(as.factor(xgbTree_mod[['predy']])), as.factor(x.sp@metaData[-idy,c('disease')]), plotROC = TRUE)
  dev.off()
  
  library(cvms)
  data <- data.frame(
    "target" = as.character(x.sp@metaData[-idy,c('disease')]),
    "prediction" = as.character(xgbTree_mod[['predy']]),
    stringsAsFactors = FALSE
  )
  eval <- evaluate(
    data = data,
    target_col = "target",
    prediction_cols = "prediction",
    type = 'multinomial'
  )
  
  tiff(paste0('output/Ast_sub_wo_CBD3/5_xgb_confmat.tiff'), units = 'px',compression = 'lzw', res = 300, width = 1400, height = 1400) 
  plot_confusion_matrix(eval[["Confusion Matrix"]][[1]], add_sums = TRUE, add_row_percentages = F, add_col_percentages = F,
                        sums_settings = sum_tile_settings(
                          label = "Total",
                          tc_tile_border_color = "black"
                        ))
  dev.off()
  pdf("output/Ast_sub_wo_CBD3/5_xgb_confmat.pdf")       # Export PDF
  x <- as.data.frame(eval[1:8]) %>% t()
  x[,1] <- round(x[,1], 3)
  gridExtra::grid.table(x)
  dev.off()
  
  # Use LIME to explain the model
  # Create explainer object
  components_lime <- lime(
    x = as.data.frame(mat_tr),
    model = xgbTree_mod[['model']], 
    n_bins = 10
  )
  summary(components_lime)
  
  psp_ob <- subset(mat_te, outcome = 1)[1,] %>% as.data.frame()
  ctrl_ob <- subset(mat_te, outcome = 0.5)[1,] %>% as.data.frame()
  cbd_ob <- subset(mat_te, outcome = 0)[1,] %>% as.data.frame()
  
  lime_explanation <- lime::explain(
    x = as.data.frame(mat_te), 
    explainer = components_lime, 
    n_permutations = 100,
    dist_fun = "gower",
    bin_continuous = T,
    kernel_width = 0.5,
    n_features = 10, 
    n_labels = 3,
    feature_select = "highest_weights"
  )
  
  best_pred_PSP <- subset(lime_explanation, label %in% 'PSP') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
  best_pred_CBD <- subset(lime_explanation, label %in% 'CBD') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
  best_pred_Ctrl <- subset(lime_explanation, label %in% 'Ctrl') %>% .[.$model_prediction == max(.$model_prediction),] %>% .$case %>% .[1]
  
  cairo_pdf(paste0('output/Ast_sub_wo_CBD3/5_lime_bestmodelxgb_featimportance.pdf'), width = 20, height = 10) 
  plot_features(lime_explanation, ncol = 3, cases = c(best_pred_CBD, best_pred_Ctrl, best_pred_PSP))
  dev.off()
  
  top_TFs_ast <- list()
  top_TFs_ast[['top_pred_PSP']] <- subset(lime_explanation, label %in% 'PSP') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']
  top_TFs_ast[['top_pred_CBD']] <- subset(lime_explanation, label %in% 'CBD') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),]  %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']
  top_TFs_ast[['top_pred_Ctrl']] <- subset(lime_explanation, label %in% 'Ctrl') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:5, 'feature']
  
  tiff(paste0('output/Ast_sub_wo_CBD3/5_featureplot.tiff'), units = 'px',compression = 'lzw', res = 450, width = 6000, height = 3600) 
  plot(caret::featurePlot(mat_t[,as.vector(unlist(top_TFs_ast))], outcome, "box", jitter = F))
  dev.off()
  
  top_TFs_ast <- list()
  top_TFs_ast[['top_pred_PSP']] <- subset(lime_explanation, label %in% 'PSP') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]
  top_TFs_ast[['top_pred_CBD']] <- subset(lime_explanation, label %in% 'CBD') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),]  %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]
  top_TFs_ast[['top_pred_Ctrl']] <- subset(lime_explanation, label %in% 'Ctrl') %>% .[(.$model_prediction > 0.6 & .$feature_weight > 0.1),] %>% .[order(desc(.$feature_weight)),] %>% .[1:10, c('feature','feature_desc')]
  
  saveRDS(top_TFs_ast, 'output/Ast_sub_wo_CBD3/top_TF_ast.Rds')
  saveRDS(components_lime, 'output/Ast_sub_wo_CBD3/components_lime.Rds')
  saveRDS(lime_explanation, 'output/Ast_sub_wo_CBD3/lime_explanation.Rds')
  saveRDS(xgbTree_mod, 'output/Ast_sub_wo_CBD3/xgbTree_mod.Rds')

  
    #### 1 Convert SnapATAC input ####
  #### 2 Co-Accessibility ####
  #### 3 Pseudotime Trajectories ####
    
    # Use CICERO for  prepare data set only for Ctrl and CBD!
    ## 3.1 Preprocessing and Traj Plot
    # Run PCA then UMAP on the data
    x.sp <- x.sp[which(x.sp@metaData$disease %in% c('Ctrl', 'CBD')),]
    cicero_cds <- snapToCicero(x.sp,preprocces_cic = F, preprocess_traj = T, umap.min_dist = 0.25, umap.n_neighbors = 15, cores = 8, cell_bin = NULL)
    cicero_cds@metadata$line <- 'cbdctrl'
    cicero_cds@colData$TF <- x.sp@mmat[,'MA0886.1_EMX2']
    
    set.seed(22)
    learn_graph_control_params <- list(
      euclidean_distance_ratio = 20, 
      geodesic_distance_ratio = 10, 
      minimal_branch_len = 20,
      orthogonal_proj_tip = F,prune_graph = T, rann.k = 10)
    cicero_cds <- learn_graph(cicero_cds, learn_graph_control = learn_graph_control_params, use_partition = F)
    
    # We find all the cells that are close to the starting point
    cicero_cds <- order_cells(cicero_cds, root_pr_nodes=get_earliest_principal_node(cicero_cds,mmat = x.sp@mmat,TF = 'EMX2', TF_cut_off = 0.01))
    cicero_cds@colData$Pseudotime <- as.vector(pseudotime(cicero_cds))
    
     clus_cols <- ggpubr::get_palette(palette = "Dark2",length(levels(cicero_cds@clusters$UMAP$clusters)))
    dis_cols <- c("#377EB8", # blue
                  "#E41A1C" # red
                  )
    
    # Plot the graph
    plot_cells(cicero_cds, 
               reduction_method = 'UMAP',
               color_cells_by = "Pseudotime",
               group_cells_by = "cluster",
               graph_label_size = 3,
               label_branch_points = F,
               label_cell_groups = F,
               cell_size = 0.75) + theme_void() + theme(legend.position="bottom") #+ scale_color_manual(values = dis_cols)
  saveRDS(cicero_cds, 'output/Ast_sub_wo_CBD3/cicero_cds_ast_cbd.Rds')
  saveRDS(x.sp, 'output/Ast_sub_wo_CBD3/snap_cbd_ast.Rds')
  
  p1 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "Pseudotime",
                   group_cells_by = "cluster",graph_label_size = 5,label_branch_points = F,label_cell_groups = F,
                   cell_size = 1.1) + #geom_point(shape=21, color ='grey21',alpha =1, stroke = 0.01) +
    theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom") + coord_equal() + 
    scale_color_viridis_c(option = 'cividis') + labs(color = 'Pseudotime')
  cicero_cds@colData$EMX2 <- log2(x.sp@mmat[,colnames(dplyr::select(as.data.frame(x.sp@mmat),contains('EMX2')))]+1)
  p2 <- plot_cells(cicero_cds,reduction_method = 'UMAP',color_cells_by = "EMX2",
                   group_cells_by = "cluster",group_label_size = 1.5,label_branch_points = F,label_cell_groups = F,
                   cell_size = 1.1) + theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom") + coord_equal()+
    colorspace::scale_color_continuous_diverging(palette = 'Purp', p2 = .3,rev = T) + labs(color = "log2(TFME:EMX2)")
  
  p3 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "disease",
                   group_cells_by = "cluster",graph_label_size = 5,label_branch_points = F, label_cell_groups = F,
                   cell_size = 1.1) + scale_color_manual(values =   c("#377EB8","#E41A1C")) + #geom_point(shape=21, color ='grey21',alpha =1, stroke = 0.01) +
    theme_void(base_family="Arial", base_size = 10) + theme(legend.position="bottom")+ coord_equal() + labs(fill = "Entity")
  
  cairo_pdf(file = "output/Ast_sub_wo_CBD3/asts_umap_cbd_wo_CBD3.pdf", width = 18, height = 18)
  plot(cowplot::plot_grid(p2, p1, p3,ncol = 4, nrow = 1))
  dev.off()
  
  # re-organize data environment 
   rm(list=setdiff(ls(), c("x.sp", 'cicero_cds', 'conns', 'clus_cols', 'dis_cols', 'get_earliest_principal_node', 'learn_graph_control_params')))
   source('scripts/04_cicero_analysis_functions.R')
   
  #### 4 Transcription Factors ####
   # 4.1 Trajectory TF motif enrichment tile gradients
    myTFs<- c("HESX1","FOS","JUN", "RFX4","EMX2")
    ifelse(dir.exists(paste0('output/Ast_sub_wo_CBD3/',cicero_cds@metadata$line)),
             stop(paste0("A storage directory for this project ",cicero_cds@metadata$line,"already exists. Please edit first.")),
             dir.create(paste0('output/Ast_sub_wo_CBD3/',cicero_cds@metadata$line)))
    plotTFenrichment(cds = cicero_cds,
                     TF = myTFs, pt_cut_off = 2, breaks = 5, cut_off = 10, explore = 0, legend_position = 'bottom', save = T, path = paste0('output/Ast_sub_wo_CBD3/',cicero_cds@metadata$line,'/4.4.1_traj_'))
   
    rm(list=setdiff(ls(), c("x.sp", 'cicero_cds', 'cic_ls', 'conns', 'clus_cols', 'dis_cols', 'c3d')))
    
  #### 5 tradeSeq: CBD astrocytes ####
    library(tradeSeq)
    library(SingleCellExperiment)
    library(slingshot)
    library(magrittr)
    library(monocle3)
    library(clusterExperiment)
    source('scripts/04_cicero_analysis_functions.R')
    set.seed(42)
    cds <- cicero_cds
    rm(cicero_cds)
    
    meta <- extractPtRootWeigths(cds)
    pseudotime <- meta[['pseudotime']]
    cellWeights <- meta[['cellWeights']]
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- 4 

    ##### 5.1 Gene accessibility in tradeSeq: CBD #### 
    
    ##### 5.2 TFME in tradeSeq: CBD ####
    # prep matrix, with nonnegative numbers
    mat <- x.sp@mmat %>% scales::rescale(to = c(0,1), from = range(., na.rm = TRUE, finite = TRUE)) %>% t()
    
    # search for best K in 
    set.seed(42)
    icMat <- evaluateK(counts = as.matrix(mat), pseudotime = pseudotime, 
                       cellWeights = cellWeights, 
                       k = c(5:13), nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
    
   ## 5.1.1 fit the binomial GAM
    sce_tf <- fitGAM(counts = as.matrix(mat),
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 7,#since AIC in icMat is lowest at 5 or 7
                  parallel = T,
                  BPPARAM = BiocParallel::bpparam())
    
    saveRDS(sce_tf, 'output/Ast_sub_wo_CBD3/fitGAM_tf.Rds')
    
  ## 5.2.2 Cluster all TF-candidates based on their TFME/pt curves 
    set.seed(42)
    nPointsClus <- 20
    clusPat <- clusterExpressionPatterns(sce_tf, nPoints = nPointsClus, ncores = 6, genes = colnames(x.sp@mmat))
    ref_tbl <- data.frame(tf = names(sce_tf), tf_clus = primaryCluster(clusPat$rsec))
    saveRDS(clusPat, 'output/Ast_sub_wo_CBD3/clusPat_tf.Rds')
    cUniq <- unique(ref_tbl$tf_clus)
    cUniq <- cUniq[!cUniq == -1] # remove unclustered TFs
    
    cluster_tf_plot <- list()
    for (xx in cUniq) {
      cId <- which(ref_tbl$tf_clus == xx)
      p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                    y = rep(range(clusPat$yhatScaled[cId, ]),
                                            nPointsClus / 2)),
                  aes(x = x, y = y)) +
        geom_point(alpha = 0) +
        labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "TFME") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
      for (ii in 1:length(cId)) {
        geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
        p <- p +
          geom_line(data = data.frame(x = rep(1:nPointsClus, 1),#change integer here, depending on n lineages
                                      y = clusPat$yhatScaled[geneId, ],
                                      lineage = rep(1:1, each = nPointsClus)), #change integer here, depending on n lineages
                    aes(col = as.character(lineage), group = lineage), lwd = 1.5)
      }
      p <- p + guides(color = FALSE) +
        scale_color_manual(values = dis_cols
        )  
      #print(p)
      cluster_tf_plot[[xx]] <- p 
    }
    
    tiff(paste0('output/Ast_sub_wo_CBD3/5.2.2_tradeseq_cluster_tf.tiff'), units = 'px',compression = 'lzw', res = 450, width = 6000, height = 6000)
    ggpubr::ggarrange(plotlist = cluster_tf_plot, nrow = round(sqrt(length(cluster_tf_plot)))+1,ncol=round(sqrt(length(cluster_tf_plot))), common.legend = T)
    dev.off()
    tf_int <- subset(ref_tbl, tf_clus %in% c(2,3))$tf %>% as.character()
    tf_down <- subset(ref_tbl, tf_clus %in% c(1))$tf %>% as.character()
    
  ## 5.2.3 Test TFME/Pt curves and changes
    # Within-lineage comparison
    # function to calculate the false discovery rate 
    extent_tradeseq_output <- function(test_res, pval){
      test_res$pvalue_fdr <- p.adjust(test_res[,pval[1]], 'fdr', n = nrow(test_res))
      return(test_res)
    }
    
    assoTF <- associationTest(sce_tf) 
    sigAsso <- subset(assoTF, pvalue < 0.05) %>% extent_tradeseq_output(pval=c('pvalue'))
    
    # Start Res: progenitor and differentiating marker genes
    startRes <- startVsEndTest(sce_tf)
    sigGeneStart <- subset(startRes, pvalue < 0.1) %>% extent_tradeseq_output(pval=c('pvalue'))
    
    # retrieve all top scoring genes
      TF_int_1<- c(rownames(sigGeneStart), rownames(sigAsso)) %>% unique()
      TF_int_2 <- c(tf_int, tf_down)
      
      saveRDS(list(TF_int_1, TF_int_2), 'output/Ast_sub_wo_CBD3/5.2.3_TF_int_tradeseq.Rds')
      
      VennDiagram::venn.diagram(x = list(A = rownames(sigGeneStart),
                                        B = rownames(sigAsso),
                                        C = tf_int,
                                        D = tf_down), col = "transparent",
                                fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                                alpha = 0.50, filename = 'output/Ast_sub_wo_CBD3/5.2.3_venny_tf_int')
      
    # plot 'em all
    for(i in intersect(TF_int_1, TF_int_2)){
      cds@colData$sigGene <- x.sp@mmat[,i]
      a <- plot_cells(cds, color_cells_by = 'disease', cell_size = 1, graph_label_size = 2.5, label_cell_groups = F) + scale_color_manual(values = dis_cols) + coord_equal() + labs(title = paste0(i))
      b <- plotSmoothers(sce_tf[[i]])+labs(y='TFME score')
      tiff(paste0('output/Ast_sub_wo_CBD3/5.2.3_tradeseq_tf_',i,'.tiff'), units = 'px',compression = 'lzw', res = 400, width = 4000, height = 1800)
      plot(ggpubr::ggarrange(a, b, nrow = 1,ncol=2))
    dev.off()
      }
    
    # TFs of interest: plot smoothers differentiated by lineage
    tf_plot <- list()
    tf_plot <- mclapply(intersect(TF_int_1, TF_int_2), mc.cores = 8, function(i){ 
      motif_i <- x.sp@mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(i)) %>% colnames()
      dat <- data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i])
      
      if(!is.null(dat$y[!is.na(dat$y)])){
        p <- plotSmoothers(sce_tf[[motif_i]], alpha = 1, border = TRUE) + ggtitle(paste0(motif_i)) + labs(y = 'TFME')
        tf_plot[[motif_i]] <- p
      }else{
        paste0(i,'is not found in the motif enrichment table.')
      }
    })
    
    tiff(paste0('output/Ast_sub_wo_CBD3/5.2.3_tradeseq_tf_interest.tiff'), units = 'px',compression = 'lzw', res = 300, width = 3500, height = 4000)
    ggpubr::ggarrange(plotlist = tf_plot, nrow = round(sqrt(length(tf_plot)))+1,ncol=round(sqrt(length(tf_plot))), common.legend = T)
    dev.off()
    
  # 5.2.4 Visualize in heatmap
    # Pseudotime gradient-heatmaps
    library(ggthemes)
    sce_tf@int_colData$disease <- x.sp@metaData$disease
    set.seed(42)

    add <- x.sp@mmat %>% as.data.frame() %>% dplyr::select(contains(c('LHX9','SHOX','EMX','HESX','RFX4','IRF','TFEB','CREB'))) %>% colnames()
    i <- c(TF_int_1, add)
      yhatSmooth <- matrix(nrow = length(i), ncol = 200)
      rownames(yhatSmooth) <- i
      for(k in i){yhatSmooth[k,] <- predictSmooth(sce_tf[[k]], nPoints = 200)}
    hc <- hclust(dist(yhatSmooth), method = 'ward.D2')
    ord <- rownames(yhatSmooth)[hc$order] %>% unique()
    tiff(paste0('output/Ast_sub_wo_CBD3/5.2.4_tradeseq_tf_lin1_hctree',i,'.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2200, height = 4000)
    plot(ape::as.phylo(hc), cex=0.9, label.offset=1)
    dev.off()
    
    yhatSmooth <- scales::rescale(yhatSmooth, to = c(0,1), from=range(yhatSmooth)) %>% as.data.frame()
    yhatSmooth$gene <- factor(rownames(yhatSmooth), levels = ord)
    yhatSmooth <- yhatSmooth[order(yhatSmooth$gene),]
    Pseudotime <- as.numeric(1:200)
    yhatSmooth <- rbind(yhatSmooth[,1:200],Pseudotime)
    yhatSmooth <- pivot_longer(as.data.frame(t(yhatSmooth)),  
                 cols = starts_with("MA"),
                 names_to = "TF",
                 values_to = "yhat",
                 values_drop_na = TRUE)
    colnames(yhatSmooth)[1] <- 'time'
    
    # smoothers
    smoo1 <- plotSmoothers(sce_tf[['MA0099.2_FOS::JUN']], alpha = 1, border = TRUE) + ggtitle(paste0('MA0099.2_FOS::JUN')) + labs(y = 'TFME')+
      theme_base(base_family="Arial", base_size = 12) 
    smoo2 <- plotSmoothers(sce_tf[['MA0841.1_NFE2']], alpha = 1, border = TRUE) + ggtitle(paste0('MA0841.1_NFE2')) + labs(y = 'TFME')+
      theme_base(base_family="Arial", base_size = 12)
    # heatmap
    
    heat1 <- plot(ggplotify::as.ggplot(pheatmap::pheatmap(t(scale(t(yhatSmooth))),
                                                          color = colorspace::diverging_hcl(n = 200, palette = 'Blue-Red 3',power = 0.75),
                                                          fontsize_row = 7,
                                                          cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize = 7)))
    cairo_pdf(file = "output/Ast_sub_wo_CBD3/supplFig_fig4_wo_CBD3_tfs.pdf", width = 10, height = 10)
    plot(cowplot::plot_grid(smoo1, smoo2, heat1, nrow = 2, axis = 'tblr'))
    dev.off()

    
    
  #### 6 Integration and Ast Signature ##### 
    x.sp <- readRDS('output/Ast_sub_wo_CBD3/psp_cbd_ast.Rds')
    ##### 6.1 Combine single analysis branches #####
    ##### 6.1.1 PSP-Allen2018-dataset: RTN associated with TA-scores ####
    # keep only significant ones 
    gsea1_TA <- readRDS("~/projs/03_Knit_RTN/RTN/final/gsea1_final_TA.Rds") 
    gsea1_TA_sig <- subset(gsea1_TA, Adjusted.Pvalue < 0.05) %>% rownames()
    mra <- read.csv("~/projs/03_Knit_RTN/mra_regulons.csv", row.names = 1)
    mra_sig <- subset(mra, Adjusted.Pvalue < 0.05) %>% rownames()
    
    ##### 6.1.2 tradeSeq data ####
    tf_int <- readRDS("output/Ast_sub_wo_CBD3/5.2.3_TF_int_tradeseq.Rds")
    TF_int_1 <- tf_int[[1]] #Tradeseq's waldtest-results
    TF_int_2 <- tf_int[[2]] #Tradeseq's visual clustering results  
    TF_int_1_sig <- as.data.frame(TF_int_1) %>% separate(col = 'TF_int_1', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF
    TF_int_2_sig <- as.data.frame(TF_int_2) %>% separate(col = 'TF_int_2', into =c('prefix', 'TF'), sep = c("_")) %>% .$TF
    
    ##### 6.1.3 Group-wise comparisons (triangular) SnapATAC / chromVAR-motifs results ####
    tf_mat_snap <- readRDS("output/Ast_sub_wo_CBD3/TF_mat_comparison_snap.Rds") %>% .[[1]] %>% .$TF
    
    ##### 6.1.4 XGB-model lime explainer ####
    lime_explanation <- readRDS("output/Ast_sub_wo_CBD3/lime_explanation.Rds")
    
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
    
    
    ##### 6.2 1st stage integration ####
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
    
    saveRDS(x, 'output/Ast_sub_wo_CBD3/upset_list_all.Rds')

    int_1 <- Reduce(intersect, x)
    int_2 <- Reduce(intersect, list(x[['XGB_model']], x[['trajectory']], x[['Pw.TFME']])) %>% .[!(. %in% int_1)]
    int_3 <- Reduce(intersect, list(x[['XGB_model']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2)]
    int_4 <- Reduce(intersect, list(x[['trajectory']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2 | . %in% int_3)]
    xplan <-  as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(c(int_1, int_2, int_3, int_4))) %>% colnames() 
    
    ##### 6.3 2nd stage integration ####
    # to show divergence through
    #  1 SnapATAC's/chromVAR-motif disease-wise Wilcoxon-test results 
    #  2 XGB-model disease-specific weights   
    
    # SnapATAC / chromVAR-motifs results
    tf_mat_snap <- readRDS("output/Ast_sub_wo_CBD3/TF_mat_comparison_snap.Rds") %>% .[[1]] 
    tf_mat_snap <- subset(tf_mat_snap, comparison %in% "PSP vs. CBD") %>% .$TF
    
    y = list(Pw.TFME_psp.vs.cbd = tf_mat_snap, XGB_PSP = xgb_psp, XGB_CBD = xgb_cbd)
    saveRDS(y, 'output/upset_list_psp_cbd.Rds')

    int_5 <- Reduce(intersect, y)
    int_6 <- Reduce(intersect, list(y[['XGB_CBD']], y[['Pw.TFME_psp.vs.cbd']]))
    int_7 <- Reduce(intersect, list(y[['XGB_PSP']], y[['Pw.TFME_psp.vs.cbd']]))
    explanatory_c = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_6)) %>% colnames()
    explanatory_p = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_7)) %>% colnames()
    
    rm(list = ls())
    library(UpSetR)
    library(ggplotify)
    x <- readRDS('output/Ast_sub_wo_CBD3/upset_list_all.Rds')
    y <- readRDS('output/upset_list_psp_cbd.Rds')
    x.sp <- readRDS("output/Ast_sub_wo_CBD3/psp_cbd_ast.Rds")  
    int_1 <- Reduce(intersect, x)
    int_2 <- Reduce(intersect, list(x[['XGB_model']], x[['trajectory']], x[['Pw.TFME']])) %>% .[!(. %in% int_1)]
    int_3 <- Reduce(intersect, list(x[['XGB_model']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2)]
    int_4 <- Reduce(intersect, list(x[['trajectory']], x[['RTN']], x[['Pw.TFME']])) %>% .[!(. %in% int_1 | . %in% int_2 | . %in% int_3)]
    xplan <-  as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(c(int_1, int_2, int_3, int_4))) %>% colnames() 
    
    tf_mat_snap <- readRDS("TF_mat_comparison_snap.Rds") %>% .[[1]] 
    tf_mat_snap <- subset(tf_mat_snap, comparison %in% "PSP vs. CBD") %>% .$TF
    
    int_5 <- Reduce(intersect, y)
    int_6 <- Reduce(intersect, list(y[['XGB_CBD']], y[['Pw.TFME_psp.vs.cbd']]))
    int_7 <- Reduce(intersect, list(y[['XGB_PSP']], y[['Pw.TFME_psp.vs.cbd']]))
    explanatory_c = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_6)) %>% colnames()
    explanatory_p = as.data.frame(x.sp@mmat[1:2,]) %>% dplyr::select(contains(int_7)) %>% colnames()
    
    # Fig6a
    cairo_pdf(file = "output/Fig6_up1.pdf", width = 10, height = 10)
    upset(fromList(x), order.by = 'degree', sets.bar.color = 'lightblue')
    dev.off()
    
    # Fig6b
    source('scripts/plotViz2.R')
    x.sp.ctrl <- x.sp[which(x.sp@metaData$disease %in% 'Ctrl'),]
    triang1 <- plotTriang(mat = x.sp@mmat, 
                          mat_bg = x.sp.ctrl@mmat, 
                          filter_obj = unique(xplan),
                          facet_x = 'disease', 
                          xlab = '',
                          ylab = 'TFME',
                          title = '', 
                          aggregate_fun = median,
                          order = unique(xplan),
                          save = F,
                          write_table = T)
    # Fig6c
    cairo_pdf(file = "output/Fig6_up2.pdf", width = 10, height = 10)
    upset(fromList(y), order.by = 'degree', sets.bar.color = 'darkorange')
    dev.off()
    
    # Fig6d-e
    triang2 <- plotTriang(mat = x.sp@mmat, 
                          mat_bg = x.sp.ctrl@mmat, 
                          filter_obj = unique(explanatory_c),
                          facet_x = 'disease', 
                          xlab = '',
                          ylab = 'TFME',
                          title = '', 
                          aggregate_fun = median,
                          order = unique(explanatory_c),
                          save = F,
                          write_table = F)
    triang3 <- plotTriang(mat = x.sp@mmat, 
                          mat_bg = x.sp.ctrl@mmat, 
                          filter_obj = unique(explanatory_p),
                          facet_x = 'disease', 
                          xlab = '',
                          ylab = 'TFME',
                          title = '', 
                          aggregate_fun = median,
                          order = unique(explanatory_p),
                          save = F,
                          write_table = F)
    
    tri <- ggarrange(triang3, triang2, common.legend = T, ncol = 1, legend = 'none')
    cairo_pdf(file = "output/Fig6_triang.pdf", width = 8, height = 10)
    plot(ggarrange(triang1, tri, common.legend = T, ncol = 2, legend = 'right'))
    dev.off()
    
################## finished #####################.    
print("Part 08B is done")
sessionInfo()     
     