###########################.
## TITLE: Cicero analysis + trajectories
## Author: Nils Briel, Feodor-Lynen-Strasse 23, 81377 Munich, Bavaria
## Date: "23/08/2021"
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
  source('scripts/04_cicero_analysis_functions.R')
############################.
  
  # load snap object
  x.sp <- readRDS("output/Ast_rev/psp_cbd_ast.Rds") 
  
  #### 1 Convert SnapATAC input ####
      cds_PSP <- snapToCicero(x.sp[which(x.sp@metaData$disease =='PSP'),], preprocces_cic = T, cell_bin = 2)
      cds_CBD <- snapToCicero(x.sp[which(x.sp@metaData$disease =='CBD'),], preprocces_cic = T, cell_bin = 2)
      cds_Ctrl <- snapToCicero(x.sp[which(x.sp@metaData$disease =='Ctrl'),], preprocces_cic = T, cell_bin = 2)
      
      cds_PSP@metadata$dis <- 'PSP'
      cds_CBD@metadata$dis <- 'CBD'
      cds_Ctrl@metadata$dis <- 'Ctrl'
  
  #### 2 Co-Accessibility ####
    register(MulticoreParam(8))
    data("human.hg19.genome")
 
    ## 2.2 Visualization
    temp <- tempfile()
    download.file("ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz", temp)
    gene_anno <- rtracklayer::readGFF(temp)
    unlink(temp)
    
    # rename some columns to match requirements
    gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
    gene_anno$gene <- gene_anno$gene_id
    gene_anno$transcript <- gene_anno$transcript_id
    gene_anno$symbol <- gene_anno$gene_name
    
    register(MulticoreParam(8))
    data(gene_annotation_sample)
    
    conns_plot(cds_PSP, conns = NULL, compare = NULL, genes = c('MAPT', 'STX6', 'MOBP'), gene_anno = gene_anno, CCANs = F)
    conns_plot(cds_CBD, conns = NULL, compare = NULL, genes = c('MAPT', 'STX6', 'MOBP'), gene_anno = gene_anno, CCANs = F)
    conns_plot(cds_Ctrl, conns = NULL, compare = NULL, genes = c('MAPT', 'STX6', 'MOBP'), gene_anno = gene_anno, CCANs = F)
    rm(cds_PSP, cds_CBD, cds_Ctrl)
  
    cicero_cds <- readRDS('output/Ast_rev/cicero_cds_ast_rev.Rds')
    
    conns <- run_cicero(cicero_cds, 
                        human.hg19.genome, 
                        sample_num = 100)
    head(dplyr::arrange(conns, desc(coaccess)))  
    
    saveRDS(conns, paste0('conns_',cicero_cds@metadata$dis,'.Rds'))
    
    
  #### 3 Pseudotime Trajectories ####
    
    # Use CICERO for  prepare data set only for Ctrl and CBD!
    ## 3.1 Preprocessing and Traj Plot
    # Run PCA then UMAP on the data
    x.sp <- x.sp[which(x.sp@metaData$disease %in% c('Ctrl', 'CBD')),]
    cicero_cds <- snapToCicero(x.sp, preprocess_traj = T, umap.min_dist = 0.25, umap.n_neighbors = 15, cores = 8, cell_bin = NULL)
    cicero_cds@metadata$line <- 'cbdctrl'
    
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
               cell_size = 0.75) + theme_void() + theme(legend.position="bottom")
  saveRDS(cicero_cds, 'output/Ast_rev_cbd/cicero_cds_ast_cbd.Rds')
  saveRDS(x.sp, 'output/Ast_rev_cbd/snap_cbd_ast_rev.Rds')
  
  # re-organize data environment 
   rm(list=setdiff(ls(), c("x.sp", 'cicero_cds', 'conns', 'clus_cols', 'dis_cols', 'gene_anno', 'get_earliest_principal_node', 'learn_graph_control_params')))
   source('scripts/04_cicero_analysis_functions.R')
   
  #### 4 Transcription Factors ####
   # 4.1 Trajectory TF motif enrichment tile gradients
    myTFs<- c("HESX1","FOS","JUN", "RFX4","EMX2")
    ifelse(dir.exists(paste0('output/Ast_rev_cbd/',cicero_cds@metadata$line)),
             stop(paste0("A storage directory for this project ",cicero_cds@metadata$line,"already exists. Please edit first.")),
             dir.create(paste0('output/Ast_rev_cbd/',cicero_cds@metadata$line)))
    plotTFenrichment(cds = cicero_cds,
                     TF = myTFs, pt_cut_off = 2, breaks = 5, cut_off = 10, explore = 0, legend_position = 'bottom', save = T, path = paste0('output/Ast_rev_cbd/',cicero_cds@metadata$line,'/4.4.1_traj_'))
   
    rm(list=setdiff(ls(), c("x.sp", 'cicero_cds', 'cic_ls', 'conns', 'clus_cols', 'dis_cols', 'c3d')))
    
  #### 5 tradeSeq: CBD astrocytes ####
    library(tradeSeq)
    library(SingleCellExperiment)
    library(slingshot)
    library(magrittr)
    library(monocle3)
    library(clusterExperiment)
    set.seed(42)
    cds <- cicero_cds
    rm(cicero_cds)
    
    meta <- extractPtRootWeigths(cds)
    pseudotime <- meta[['pseudotime']]
    cellWeights <- meta[['cellWeights']]
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- 4 

    ##### 5.1 Gene accessibility in tradeSeq: CBD #### 
    set.seed(42)
    icMat <- evaluateK(counts = as.matrix(t(x.sp@gmat)), pseudotime = pseudotime, 
                       cellWeights = cellWeights, 
                       k = 3:15, nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
    
    ## 5.1.1 fit the binomial GAM
    idy <- unique(colnames(x.sp@gmat))
    mat <- as.matrix(t(x.sp@gmat[,idy]))
    sce <- fitGAM(counts = mat,
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 13,
                  parallel = T,
                  BPPARAM = BiocParallel::bpparam())
    
    saveRDS(sce, 'output/Ast_rev_cbd/fitGAM_ga.Rds')
    
    cds@colData <- cbind(cds@colData, cellWeights)
    colnames(cellWeights)
    
  # 5.1.2 Cluster plot 
    nPointsClus <- 20
    clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus, ncores = 6,genes = colnames(x.sp@gmat))
    clusterLabels <- primaryCluster(clusPat$rsec)
    ref_tbl <- data.frame(gene = rownames(sce)[1:100], tf_clus = primaryCluster(clusPat$rsec))
    cUniq <- unique(clusterLabels)
    cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
    
    cluster_plot <- list()
    for (xx in cUniq) {
      cId <- which(clusterLabels == xx)
      p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                    y = rep(range(clusPat$yhatScaled[cId, ]),
                                            nPointsClus / 2)),
                  aes(x = x, y = y)) +
        geom_point(alpha = 0) +
        labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Gene activity") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
      for (ii in 1:length(cId)) {
        geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
        p <- p +
          geom_line(data = data.frame(x = rep(1:nPointsClus, 1),
                                      y = clusPat$yhatScaled[geneId, ],
                                      lineage = rep(1:1, each = nPointsClus)),
                    aes(col = as.character(lineage), group = lineage), lwd = 1.5)
      }
      p <- p + guides(color = FALSE) +
        scale_color_manual(values = dis_cols
        )  
      #print(p)
      cluster_plot[[xx]] <- p 
    }
    ga_int <- subset(ref_tbl, tf_clus %in% c(3))$gene %>% as.character()
    ga_down <- subset(ref_tbl, tf_clus %in% c(8))$gene %>% as.character()
    
    tiff(paste0('output/Ast_rev_cbd/5.1.2_tradeseq_cluster_top100.tiff'), units = 'px',compression = 'lzw', res = 450, width = 7000, height = 6000)
    ggpubr::ggarrange(plotlist = cluster_plot, nrow = 4,ncol=3, common.legend = T)
    dev.off()
    
  # 5.1.3 Test GA/Pt curves and changes
    # function to calculate the false discovery rate 
    extent_tradeseq_output <- function(test_res, pval){
      test_res$pvalue_fdr <- p.adjust(test_res[,pval[1]], 'fdr', n = nrow(test_res))
      return(test_res)
    }
    assoGA <- associationTest(sce, ) %>% extent_tradeseq_output(pval=c('pvalue'))
    sigAsso <- subset(assoGA, pvalue_fdr < 0.05 & waldStat > mean(waldStat, na.rm = T) & meanLogFC > mean(meanLogFC, na.rm = T)) %>% .[order(.$waldStat),] %>% .[1:250,]
    
    # Start Res: progenitor and differentiating marker genes
    startResGA <- startVsEndTest(sce) %>% extent_tradeseq_output(pval=c('pvalue'))
    sigGeneStart <- subset(startResGA, pvalue < 0.05 & waldStat > mean(waldStat, na.rm = T) & logFClineage1 > mean(logFClineage1, na.rm = T)) %>% .[order(.$waldStat),] %>% .[1:250,]
    
    # retrieve all top scoring genes
    ga_int_1<- c(rownames(sigGeneStart), rownames(sigAsso)) %>% unique()# rownames(assoTF[which(assoTF$pvalue_1_fdr<.05|assoTF$pvalue_2_fdr<.05|assoTF$pvalue_3_fdr<.05),]))
    ga_int_2 <- c(ga_int, ga_down)
    
    genes_interest <- c(intersect(rownames(sigGeneStart),rownames(sigAsso))) 
    genes_interest2 <- unique(c(genes_interest, c('MAPT', 'APP', 'EIF2AK3', 'APOE','LRRK2', 'SNCA', 'TARDBP', 'RBMS3', 'TREM2', 'BDNF', 'NRXN1', 'NRXN3', 'SMUG1')))
  
    saveRDS(list(ga_int_1, ga_int_2, genes_interest, genes_interest2), 'output/Ast_rev_cbd/5.1.3_GA_int_tradeseq.Rds')
    

    ##### 5.2 TFME in tradeSeq: CBD ####
    # prep matrix, with nonnegative numbers
    mat <- x.sp@mmat %>% scales::rescale(to = c(0,1), from = range(., na.rm = TRUE, finite = TRUE)) %>% t()
    
    # search for best K in 
    set.seed(42)
    icMat <- evaluateK(counts = as.matrix(mat), pseudotime = pseudotime, 
                       cellWeights = cellWeights, 
                       k = c(5:13), nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
    
   ## 5.1.1 fit the binomial GAM
    sce_tf <- fitGAM(counts = Seurat::as.sparse(mat),
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 7,#since AIC in icMat is lowest at 5 or 7
                  parallel = T,sce = TRUE,
                  BPPARAM = BiocParallel::bpparam())
    
    saveRDS(sce_tf, 'output/Ast_rev_cbd/fitGAM_tf.Rds')
    
  ## 5.2.2 Cluster all TF-candidates based on their TFME/pt curves 
    set.seed(42)
    nPointsClus <- 20
    clusPat <- clusterExpressionPatterns(sce_tf, nPoints = nPointsClus, ncores = 6, genes = colnames(x.sp@mmat))
    ref_tbl <- data.frame(tf = rownames(sce_tf), tf_clus = primaryCluster(clusPat$rsec))
    saveRDS(clusPat, 'output/Ast_rev_cbd/clusPat_tf.Rds')
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
    
    tiff(paste0('output/Ast_rev_cbd/5.2.2_tradeseq_cluster_tf.tiff'), units = 'px',compression = 'lzw', res = 450, width = 6000, height = 6000)
    ggpubr::ggarrange(plotlist = cluster_tf_plot, nrow = round(sqrt(length(cluster_tf_plot)))+1,ncol=round(sqrt(length(cluster_tf_plot))), common.legend = T)
    dev.off()
    tf_int <- subset(ref_tbl, tf_clus %in% c(3,4))$tf %>% as.character()
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
      
      saveRDS(list(TF_int_1, TF_int_2), 'output/Ast_rev_cbd/5.2.3_TF_int_tradeseq.Rds')
      
      VennDiagram::venn.diagram(x = list(A = rownames(sigGeneStart),
                                        B = rownames(sigAsso),
                                        C = tf_int,
                                        D = tf_down), col = "transparent",
                                fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                                alpha = 0.50, filename = 'output/Ast_rev_cbd/5.2.3_venny_tf_int')
      
    # plot 'em all
    for(i in intersect(TF_int_1, TF_int_2)){
      cds@colData$sigGene <- x.sp@mmat[,i]
      a <- plot_cells(cds, color_cells_by = 'disease', cell_size = 1, graph_label_size = 2.5, label_cell_groups = F) + scale_color_manual(values = dis_cols) + coord_equal() + labs(title = paste0(i))
      b <- plotSmoothers(sce_tf, counts=mat, i)+labs(y='TFME score')
      tiff(paste0('output/Ast_rev_cbd/5.2.3_tradeseq_tf_',i,'.tiff'), units = 'px',compression = 'lzw', res = 400, width = 4000, height = 1800)
      plot(ggpubr::ggarrange(a, b, nrow = 1,ncol=2))
    dev.off()
      }
    
    # TFs of interest: plot smoothers differentiated by lineage
    tf_plot <- list()
    tf_plot <- mclapply(intersect(TF_int_1, TF_int_2), mc.cores = 8, function(i){ 
      motif_i <- x.sp@mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(i)) %>% colnames()
      dat <- data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i])
      
      if(!is.null(dat$y[!is.na(dat$y)])){
        p <- plotSmoothers(sce_tf, gene = motif_i, counts = mat, alpha = 1, border = TRUE) + ggtitle(paste0(motif_i)) + labs(y = 'TFME')
        tf_plot[[motif_i]] <- p
      }else{
        paste0(i,'is not found in the motif enrichment table.')
      }
    })
    
    tiff(paste0('output/Ast_rev_cbd/5.2.3_tradeseq_tf_interest.tiff'), units = 'px',compression = 'lzw', res = 300, width = 3500, height = 4000)
    ggpubr::ggarrange(plotlist = tf_plot, nrow = round(sqrt(length(tf_plot)))+1,ncol=round(sqrt(length(tf_plot))), common.legend = T)
    dev.off()
    
  # 5.2.4 Visualize in heatmap
    # Pseudotime gradient-heatmaps
    sce_tf@int_colData$disease <- x.sp@metaData$disease
    set.seed(42)
    
    for(i in list(TF_int_1, tf_int, tf_down)){
    yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =i, nPoints = 20, tidy = F)
    hc <- hclust(dist(yhatSmooth), method = 'ward.D2')
    ord <- rownames(yhatSmooth)[hc$order] %>% unique()
    tiff(paste0('output/Ast_rev_cbd/5.2.4_tradeseq_tf_lin1_hctree',i,'.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2200, height = 4000)
    plot(ape::as.phylo(hc), cex=0.9, label.offset=1)
    dev.off()
    yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =i, nPoints = 200, tidy = T)
    yhatSmooth$yhat <- scales::rescale(yhatSmooth$yhat, to = c(0,1), from=range(yhatSmooth$yhat))
    yhatSmooth$gene <- factor(yhatSmooth$gene, levels = ord)
    yhatSmooth <- yhatSmooth[order(yhatSmooth$gene),]
    
    tiff(paste0('output/Ast_rev_cbd/5.2.4_tradeseq_tf_lin1tau_cbd',i,'.tiff'), units = 'px',compression = 'lzw', res = 400, width = 2200, height = 3500)
    p1 <- plot(ggplot(yhatSmooth[complete.cases(yhatSmooth),], aes(x=time, y=gene, fill = yhat)) + 
      scale_fill_viridis_c(option = 'inferno', begin = min(yhatSmooth$yhat), end = max(yhatSmooth$yhat)) +
      geom_tile(color="white", size=0.1, linetype='blank') + 
      coord_equal() + labs(x = 'pseudotime', y = '', fill = 'TFME') +
      ggthemes::theme_tufte(base_family="Helvetica") +
      theme(legend.position="bottom") + scale_y_discrete(limit = unique(yhatSmooth$gene)) + 
      theme(axis.text.x = element_text(angle = 0, hjust = 1)))
    dev.off()
    }
    
    ###########################.
    #### 6 tradeSeq: PSP astrocytes ####
    source('scripts/04_cicero_analysis_functions.R')
    # create new output sub-folder
    ifelse(dir.exists(paste0('output/Ast_rev_psp/')),
           stop(paste0("A storage directory for this project 'Ast_rev_psp' already exists. Please edit first.")),
           dir.create('output/Ast_rev_psp/'))
    ##########.
    # load snap object again
    x.sp <- readRDS("output/Ast_rev/psp_cbd_ast.Rds") 
    x.sp <- x.sp[which(x.sp@metaData$disease %in% c('Ctrl', 'PSP')),]
    cicero_cds <- snapToCicero(x.sp, preprocess_traj = T, prep = 'LSI', umap.min_dist = 0.5, umap.n_neighbors = 30, cores = 8, cell_bin = NULL)
    cicero_cds@metadata$line <- 'pspctrl'

    learn_graph_control_params <- list(
      euclidean_distance_ratio = 40, 
      geodesic_distance_ratio = 5, 
      minimal_branch_len = 30,
      orthogonal_proj_tip = F,prune_graph = T, rann.k = 30)
      cicero_cds <- learn_graph(cicero_cds,learn_graph_control = learn_graph_control_params,use_partition = F, close_loop = F)
      cicero_cds <- order_cells(cicero_cds, root_pr_nodes=get_earliest_principal_node(cicero_cds,mmat = x.sp@mmat,TF = 'EMX2', TF_cut_off = 0.1))
      
    clus_cols <- ggpubr::get_palette(palette = "Dark2",length(levels(cicero_cds@clusters$UMAP$clusters)))
    dis_cols <- c("#377EB8", # blue
                  "#4DAF4A" # green
    )
    
    # Plot the graph
    plot_cells(cicero_cds, 
               reduction_method = 'UMAP',
               color_cells_by = "Pseudotime",
               group_cells_by = "cluster",
               graph_label_size = 3,
               label_branch_points = F,
               label_cell_groups = F,
               cell_size = 0.75) + theme_void() + theme(legend.position="bottom")
    saveRDS(cicero_cds, 'output/Ast_rev_psp/cicero_cds_ast_psp.Rds')
    saveRDS(x.sp, 'output/Ast_rev_psp/snap_ast_psp.Rds')
    
    # Plotting
    tiff(paste0('output/Ast_rev_psp/6.1_pt_sub_overview_part',cicero_cds@metadata$line,'.tiff'), units = 'px',compression = 'lzw', res = 300, width = 4000, height = 1200) 
    p1 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "Pseudotime",
                     group_cells_by = "cluster",graph_label_size = 3,label_branch_points = F,label_cell_groups = F,
                     cell_size = 0.75) + 
      theme_void() + theme(legend.position="bottom") + coord_equal() + scale_color_viridis_c(option = 'cividis') + labs(color = 'Pseudotime')
    cicero_cds@colData$EMX2 <- x.sp@mmat[,colnames(dplyr::select(as.data.frame(x.sp@mmat),contains('EMX2')))]
    p2 <- plot_cells(cicero_cds,reduction_method = 'UMAP',color_cells_by = "EMX2",
                     group_cells_by = "cluster",group_label_size = 3,label_branch_points = F,label_cell_groups = F,
                     cell_size = 0.75) + 
      theme_void() + theme(legend.position="bottom") + coord_equal() + 
      scale_colour_gradient2(low = "white",mid = 'grey90', high = "blue", limits = c(-0.5,0.7)) + labs(color = "TFME: EMX2")
    p3 <- plot_cells(cicero_cds, reduction_method = 'UMAP',color_cells_by = "disease",
                     group_cells_by = "cluster",graph_label_size = 3,label_branch_points = F, label_cell_groups = F,
                     cell_size = 0.75) + 
      scale_color_manual(values = dis_cols) + theme_void() + theme(legend.position="bottom")+ coord_equal()
    
    print(cowplot::plot_grid(p2, p1, p3, p4, ncol = 4, nrow = 1, greedy = T, rel_widths = c(1, 1, 1, 1)))
    dev.off()
    
    ##### 6.2 TFME in tradeSeq: PSP #####
    set.seed(42)
    cds <- cicero_cds
    rm(cicero_cds)
    meta <- extractPtRootWeigths(cds)
    pseudotime <- meta[['pseudotime']]
    cellWeights <- meta[['cellWeights']]
    
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- 4 
    
    # prep matrix, with nonnegative numbers
    mat <- x.sp@mmat %>% scales::rescale(to = c(0,1), from = range(., na.rm = TRUE, finite = TRUE)) %>% t()
    
    # search for best K in 
    set.seed(42)
    icMat <- evaluateK(counts = as.matrix(mat), pseudotime = pseudotime, 
                       cellWeights = cellWeights, 
                       k = c(5:13), nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
    
    ## 6.2.1 fit the binomial GAM
    sce_tf <- fitGAM(counts = Seurat::as.sparse(mat),
                     pseudotime = pseudotime,
                     cellWeights = cellWeights,
                     nknots = 6,
                     parallel = T,sce = TRUE,
                     BPPARAM = BiocParallel::bpparam())
    
    saveRDS(sce_tf, 'output/Ast_rev_psp/fitGAM_tf.Rds')
    
    ## 6.2.2 Cluster all TF-candidates based on their TFME/pt curves 
    nPointsClus <- 20
    clusPat <- clusterExpressionPatterns(sce_tf, nPoints = nPointsClus, ncores = 6, genes = colnames(x.sp@mmat))
    ref_tbl <- data.frame(tf = rownames(sce_tf), tf_clus = primaryCluster(clusPat$rsec))
    saveRDS(clusPat, 'output/Ast_rev_psp/clusPat_tf.Rds')
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
    
    tiff(paste0('output/Ast_rev_psp/5.2.2_tradeseq_cluster_tf.tiff'), units = 'px',compression = 'lzw', res = 450, width = 6000, height = 6000)
    ggpubr::ggarrange(plotlist = cluster_tf_plot, nrow = round(sqrt(length(cluster_tf_plot))),ncol=round(sqrt(length(cluster_tf_plot))), common.legend = T)
    dev.off()
    tf_int <- subset(ref_tbl, tf_clus %in% c(2))$tf %>% as.character()
    tf_down <- subset(ref_tbl, tf_clus %in% c(3))$tf %>% as.character()
    
    ## 6.2.3 Test TFME/Pt curves and changes
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
    sigGeneStart <- subset(startRes, pvalue < 0.05) %>% extent_tradeseq_output(pval=c('pvalue'))
    
    # retrieve all top scoring genes
    TF_int_1<- c(rownames(sigGeneStart), rownames(sigAsso)) %>% unique()# rownames(assoTF[which(assoTF$pvalue_1_fdr<.05|assoTF$pvalue_2_fdr<.05|assoTF$pvalue_3_fdr<.05),]))
    TF_int_2 <- c(tf_int, tf_down)
    
    saveRDS(list(TF_int_1, TF_int_2), 'output/Ast_rev_psp/5.2.3_TF_int_tradeseq.Rds')
    
    
    # plot 'em all
    for(i in TF_int_2){
      cds@colData$sigGene <- x.sp@mmat[,i]
      a <- plot_cells(cds, color_cells_by = 'disease', cell_size = 1, graph_label_size = 2.5, label_cell_groups = F) + scale_color_manual(values = dis_cols) + coord_equal() + labs(title = paste0(i))#scale_color_manual(values = dis_cols)
      b <- plotSmoothers(sce_tf, counts=mat, i)+labs(y='TFME score')
      tiff(paste0('output/Ast_rev_psp/5.2.3_tradeseq_tf_',i,'.tiff'), units = 'px',compression = 'lzw', res = 400, width = 4000, height = 1800)
      plot(ggpubr::ggarrange(a, b, nrow = 1,ncol=2))
      dev.off()
    }
    
    # TFs of interest: plot smoothers differentiated by lineage
    tf_plot <- list()
    tf_plot <- mclapply(TF_int_2, mc.cores = 8, function(i){ 
      motif_i <- x.sp@mmat[1:2,] %>% as.data.frame() %>% dplyr::select(contains(i)) %>% colnames()
      dat <- data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i])
      
      if(!is.null(dat$y[!is.na(dat$y)])){
        p <- plotSmoothers(sce_tf, gene = motif_i, counts = mat, alpha = 1, border = TRUE) + ggtitle(paste0(motif_i)) + labs(y = 'TFME')# + scale_color_manual(values = dis_cols)
        tf_plot[[motif_i]] <- p
      }else{
        paste0(i,'is not found in the motif enrichment table.')
      }
    })
    
    tiff(paste0('output/Ast_rev_psp/5.2.3_tradeseq_tf_interest.tiff'), units = 'px',compression = 'lzw', res = 300, width = 3500, height = 4000)
    ggpubr::ggarrange(plotlist = tf_plot, nrow = round(sqrt(length(tf_plot)))+1,ncol=round(sqrt(length(tf_plot))), common.legend = T)
    dev.off()
    
    # 6.2.4 Vizualize in heatmap
    # Pseudotime gradient-heatmaps
    sce_tf@int_colData$disease <- x.sp@metaData$disease
    set.seed(42)
    yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =c(tf_int, tf_down), nPoints = 20, tidy = F)
    hc <- hclust(dist(yhatSmooth), method = 'centroid')
    ord <- rownames(yhatSmooth)[hc$order] %>% unique()
    tiff(paste0('output/Ast_rev_psp/5.2.4_tradeseq_tf_lin1_hctree.tiff'), units = 'px',compression = 'lzw', res = 300, width = 2200, height = 4000)
    plot(ape::as.phylo(hc), cex=0.9, label.offset=1)
    dev.off()
    yhatSmooth <- tradeSeq::predictSmooth(sce_tf, gene =c(TF_int_1, TF_int_2), nPoints = 200, tidy = T)
    yhatSmooth$yhat <- scales::rescale(yhatSmooth$yhat, to = c(0,1), from=range(yhatSmooth$yhat))
    yhatSmooth$gene <- factor(yhatSmooth$gene, levels = ord)
    yhatSmooth <- yhatSmooth[order(yhatSmooth$gene),]
    
    tiff(paste0('output/Ast_rev_psp/5.2.4_tradeseq_tf_lin1tau_psp.tiff'), units = 'px',compression = 'lzw', res = 400, width = 2200, height = 2000)
    p1 <- plot(ggplot(yhatSmooth[complete.cases(yhatSmooth),], aes(x=time, y=gene, fill = yhat)) + 
                 scale_fill_viridis_c(option = 'inferno', begin = min(yhatSmooth$yhat), end = max(yhatSmooth$yhat)) +
                 geom_tile(color="white", size=0.1, linetype='blank') + 
                 coord_equal() + labs(x = 'pseudotime', y = '', fill = 'TFME') +
                 ggthemes::theme_tufte(base_family="Helvetica") +
                 theme(legend.position="bottom") + scale_y_discrete(limit = unique(yhatSmooth$gene)) + 
                 theme(axis.text.x = element_text(angle = 0, hjust = 1)))
    dev.off()
    
    
    
    ##### 6.3 Gene accessibility in tradeSeq: PSP ####
    set.seed(42)
    icMat <- evaluateK(counts = as.matrix(t(x.sp@gmat)), pseudotime = pseudotime, 
                       cellWeights = cellWeights, 
                       k = 3:15, nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
    
    ## 6.3.1 fit the binomial GAM
    idy <- unique(colnames(x.sp@gmat))
    mat <- as.matrix(t(x.sp@gmat[,idy]))
    sce <- fitGAM(counts = mat,
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 10,
                  parallel = T,
                  BPPARAM = BiocParallel::bpparam())
    
    saveRDS(sce, 'output/Ast_rev_psp/fitGAM_ga.Rds')
    
    cds@colData <- cbind(cds@colData, cellWeights)
    colnames(cellWeights)
    
    # 6.3.2 Cluster plot 
    nPointsClus <- 20
    clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus, ncores = 6,genes = colnames(x.sp@gmat)[1:100])
    clusterLabels <- primaryCluster(clusPat$rsec)
    cUniq <- unique(clusterLabels)
    cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
    
    cluster_plot <- list()
    for (xx in cUniq) {
      cId <- which(clusterLabels == xx)
      p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                    y = rep(range(clusPat$yhatScaled[cId, ]),
                                            nPointsClus / 2)),
                  aes(x = x, y = y)) +
        geom_point(alpha = 0) +
        labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Gene activity") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
      for (ii in 1:length(cId)) {
        geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
        p <- p +
          geom_line(data = data.frame(x = rep(1:nPointsClus, 1),
                                      y = clusPat$yhatScaled[geneId, ],
                                      lineage = rep(1:1, each = nPointsClus)),
                    aes(col = as.character(lineage), group = lineage), lwd = 1.5)
      }
      p <- p + guides(color = FALSE) +
        scale_color_manual(values = dis_cols
        )  
      #print(p)
      cluster_plot[[xx]] <- p 
    }
    f_int <- subset(ref_tbl, tf_clus %in% c(1))$tf %>% as.character()
    tf_down <- subset(ref_tbl, tf_clus %in% c(4))$tf %>% as.character()
    
    tiff(paste0('output/Ast_rev_psp/5.1.2_tradeseq_cluster_top100.tiff'), units = 'px',compression = 'lzw', res = 450, width = 7000, height = 6000)
    ggpubr::ggarrange(plotlist = cluster_plot, nrow = 4,ncol=3, common.legend = T)
    dev.off()
    
    ## 6.3.3 Test GA/Pt curves nd changes
    # function to calculate the false discovery rate 
    extent_tradeseq_output <- function(test_res, pval){
      test_res$pvalue_fdr <- p.adjust(test_res[,pval[1]], 'fdr', n = nrow(test_res))
      return(test_res)
    }
    assoGA <- associationTest(sce, ) %>% extent_tradeseq_output(pval=c('pvalue'))
    sigAsso <- subset(assoGA, pvalue_fdr < 0.05 & waldStat > mean(waldStat, na.rm = T) & meanLogFC > mean(meanLogFC, na.rm = T)) %>% .[order(.$waldStat),] %>% .[1:250,]
    
    # Start Res: progenitor and differentiating marker genes
    startResGA <- startVsEndTest(sce) %>% extent_tradeseq_output(pval=c('pvalue'))
    sigGeneStart <- subset(startResGA, pvalue < 0.05 & waldStat > mean(waldStat, na.rm = T) & logFClineage1 > mean(logFClineage1, na.rm = T)) %>% .[order(.$waldStat),] %>% .[1:250,]
    
    # retrieve all top scoring genes
    TF_int_1<- c(rownames(sigGeneStart), rownames(sigAsso)) %>% unique()# rownames(assoTF[which(assoTF$pvalue_1_fdr<.05|assoTF$pvalue_2_fdr<.05|assoTF$pvalue_3_fdr<.05),]))
    TF_int_2 <- c(tf_int, tf_down)
    
    saveRDS(list(TF_int_1, TF_int_2), 'output/Ast_rev_psp/5.1.3_GA_int_tradeseq.Rds')
    
    genes_interest <- c(intersect(rownames(sigGeneStart),rownames(sigAsso))) 
    genes_interest <- unique(c(genes_interest, c('MAPT', 'APP', 'EIF2AK3', 'APOE','LRRK2', 'SNCA', 'TARDBP', 'RBMS3', 'TREM2', 'BDNF', 'NRXN1', 'NRXN3', 'SMUG1')))
    
    genes_plot <- list()
    genes_plot <- mclapply(genes_interest, mc.cores = 4, function(i){ 
      p <- plotSmoothers(sce, counts = mat, gene = i, alpha = 1, border = TRUE) + ggtitle(paste0(i)) + labs(y = 'Gene activity') #+ scale_color_manual(values = dis_cols)
      genes_plot[[i]] <- p
    })
    
    tiff(paste0('output/Ast_rev_psp/5.1.3_tradeseq_genes_interest.tiff'), units = 'px',compression = 'lzw', res = 400, width = 8000, height = 8000)
    ggpubr::ggarrange(plotlist = genes_plot, nrow = 5,ncol=3, common.legend = T)
    dev.off()
    
    # 6.3.4 Vizualize in heatmap
    library(pheatmap)
    yhatSmooth <- tradeSeq::predictSmooth(sce, gene = genes_interest, nPoints = 200, tidy = F)
    
    tiff(paste0('output/Ast_rev_psp/5.1.4_tradeseq_ga_lin1tau_psp.tiff'), units = 'px',compression = 'lzw', res = 400, width = 2200, height = 2000)
    heatSmooth <- pheatmap(t(scale(t(yhatSmooth))),
                           color = viridis::viridis(256),
                           cluster_cols = FALSE,show_rownames = T,show_colnames = FALSE, scale = 'row', fontsize = 7, )
    dev.off()
    
  
################## finished #####################.    
print("Part 04 is done: now continue with '05_modelling_disease_202107.R'")
sessionInfo()     
     